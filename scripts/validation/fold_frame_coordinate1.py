"""Second fold frame coordinate on a known folded geometry.

The first coordinate of a fold frame is a folded scalar field c0. The second
coordinate c1 must increase monotonically along the isosurfaces of c0, so it
measures distance along the folded surfaces. The builder asks for this with
``grad(c1) . grad(c0) = 0`` and a few values of c1. That equation is also
satisfied by a field with a vanishing gradient, which is the gradient
collapse this module measures.

The ground truth is a similar fold of a superposed profile curve ``z = f(x)``,
extruded along y:

* ``c0 = z - f(x)``, so every isosurface is the profile shifted in z,
* ``c1`` is the distance along the reference isosurface ``c0 = 0`` of the point
  reached by following the orthogonal trajectories of c0.

``grad(c1)`` is exactly orthogonal to ``grad(c0)``. The trajectories obey
``dx/dc0 = -f' / (1 + f'^2)`` and never cross, so the truth is defined everywhere in the
mesh and ``c1`` increases along every isosurface of ``c0``. Unlike a parallel
fold this has no focal curves, where the nearest point on the profile jumps.

Metrics are compared with a ramp ``c1 = x`` that ignores the folds, because a
field can be monotonic without following the folded surfaces at all.

Run from the repository root to print the comparison table and write the
figures to ``validation/fold_frame_coordinate1``::

    python -m scripts.validation.fold_frame_coordinate1
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import logging
from pathlib import Path
import time

import numpy as np
from loop_interpolation import P1Interpolator, P2Interpolator, P2UnstructuredTetMesh, TetMesh
from loop_interpolation._quadrature import tet_quadrature_points

BOX = (np.array([0.0, 0.0, -3.0]), np.array([12.0, 2.0, 3.0]))
SOLVER = {"maxiter": 20000, "rtol": 1e-10}
# interior of the box where metrics are evaluated, away from the mesh boundary
EVAL_LO = np.array([1.5, 0.5, -1.5])
EVAL_HI = np.array([10.5, 1.5, 1.5])


class SimilarFold:
    """Similar fold of a profile curve ``z = f(x)``, extruded along y."""

    def __init__(self, a1=1.2, l1=7.0, a2=0.2, l2=3.5, phase=0.7, n_steps=40):
        self.k1, self.k2 = 2 * np.pi / l1, 2 * np.pi / l2
        self.a1, self.a2, self.phase, self.n_steps = a1, a2, phase, n_steps
        x = np.linspace(-10.0, 22.0, 64001)
        dx = np.sqrt(1 + self.slope(x) ** 2)
        self._table = (x, np.r_[0.0, np.cumsum((dx[1:] + dx[:-1]) / 2 * np.diff(x))])

    def profile(self, x):
        return self.a1 * np.sin(self.k1 * x) + self.a2 * np.sin(self.k2 * x + self.phase)

    def slope(self, x):
        return self.a1 * self.k1 * np.cos(self.k1 * x) + self.a2 * self.k2 * np.cos(
            self.k2 * x + self.phase
        )

    def curvature(self, x):
        return -self.a1 * self.k1**2 * np.sin(self.k1 * x) - self.a2 * self.k2**2 * np.sin(
            self.k2 * x + self.phase
        )

    def _c1(self, x, z):
        """Distance along ``c0 = 0`` of the point where the orthogonal trajectory
        through ``(x, z)`` meets it.

        With ``u = c0`` as the variable along a trajectory, ``dx/du = -f' / (1 + f'^2)``.
        This is integrated from ``u = c0`` to 0 with Runge-Kutta 4.
        """
        x = np.array(x, dtype=float)
        h = -(z - self.profile(x)) / self.n_steps

        def rate(xi):
            slope = self.slope(xi)
            return -slope / (1 + slope**2)

        for _ in range(self.n_steps):
            k1 = rate(x)
            k2 = rate(x + h / 2 * k1)
            k3 = rate(x + h / 2 * k2)
            k4 = rate(x + h * k3)
            x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return np.interp(x, *self._table), x

    def truth(self, xyz):
        """Return ``(c0, c1, rate, grad_c0, tangent)`` at ``xyz``.

        ``rate`` is the derivative of c1 along the unit ``tangent`` of the folded
        surface. ``grad_c0`` is the unit normal, the unit gradient of c0. Both
        vectors have zero y component.
        """
        xyz = np.atleast_2d(xyz)
        x, z = xyz[:, 0], xyz[:, 2]
        slope = self.slope(x)
        length = np.sqrt(1 + slope**2)
        zero = np.zeros(len(x))
        tangent = np.c_[1 / length, zero, slope / length]
        h = 1e-3
        ahead = xyz + h * tangent
        behind = xyz - h * tangent
        rate = (self._c1(ahead[:, 0], ahead[:, 2])[0] - self._c1(behind[:, 0], behind[:, 2])[0]) / (
            2 * h
        )
        return (
            z - self.profile(x),
            self._c1(x, z)[0],
            rate,
            np.c_[-slope / length, zero, 1 / length],
            tangent,
        )

    def norm_target(self, xyz):
        """Exact ``|grad c1|`` along the tangent. An oracle, see ``CONSTRAINT_SETS``."""
        return self.truth(xyz)[2]

    def valid(self, xyz):
        """Points where the truth is defined: all those the trajectory table covers."""
        xyz = np.atleast_2d(xyz)
        x_ref = self._c1(xyz[:, 0], xyz[:, 2])[1]
        return (x_ref > self._table[0][0]) & (x_ref < self._table[0][-1])

    def isosurface(self, offset, n=400, xrange=(1.5, 10.5), y=1.0):
        """Points on the isosurface ``c0 = offset``, ordered by increasing c1."""
        x = np.linspace(*xrange, n)
        return np.c_[x, np.full(n, y), offset + self.profile(x)]


@dataclass(frozen=True)
class Method:
    """Interpolator and the points where its constraints are sampled."""

    name: str
    kind: str  # "P1" or "P2"
    rule: str  # "barycentre" (one point per tet), "centroid" or "gauss4"
    refine: int = 1  # multiplies the number of cells along each axis
    step: int = 1  # constrain every step-th tetrahedron, the builder default is 10


METHODS = {
    m.name: m
    for m in (
        Method("P1-step10", "P1", "barycentre", step=10),
        Method("P1", "P1", "barycentre"),
        Method("P1-fine", "P1", "barycentre", refine=2),
        Method("P2-centroid", "P2", "centroid"),
        Method("P2-gauss4", "P2", "gauss4"),
    )
}

# weights of (orthogonal to grad c0, orthogonal to fold axis, norm along tangent) and
# whether the norm target is the exact ``SimilarFold.norm_target`` (else 1). The exact
# target is an oracle: it needs the answer, so it shows what fixing the magnitude of
# grad c1 would buy, not a method that is available today.
CONSTRAINT_SETS = {
    "orthogonal": (1.0, 0.0, 0.0, False),
    "orthogonal+axis": (1.0, 1.0, 0.0, False),
    "oriented": (1.0, 1.0, 0.1, False),
    "oracle norm 0.1": (1.0, 1.0, 0.1, True),
    "oracle norm 1": (1.0, 1.0, 1.0, True),
}


def sample_data(fold, n, seed):
    """``n`` random value constraints ``[x, y, z, c1]`` from the truth."""
    rng = np.random.default_rng(seed)
    xyz = np.c_[rng.uniform(0.5, 11.5, n), rng.uniform(0.2, 1.8, n), rng.uniform(-1.2, 1.2, n)]
    return np.c_[xyz, fold.truth(xyz)[1]]


def _add_direction(interpolator, points, tets, vectors, weight, target, name):
    dn, _, inside = interpolator.support.evaluate_shape_derivatives(points, elements=tets)
    ok = inside & np.all(np.isfinite(dn), axis=(1, 2))
    rows = np.einsum("ij,ijk->ik", vectors[ok], dn[ok])
    target = np.broadcast_to(np.asarray(target, dtype=float), (len(points),))
    interpolator.add_constraints_to_least_squares(
        rows,
        target[ok],
        interpolator.support.elements[tets[ok]],
        w=weight[ok],
        name=name,
    )


def solve_coordinate1(fold, method, data, constraints="orthogonal", cells=(24, 3, 12)):
    """Build c1 from ``data`` and the known c0 with ``method``.

    ``constraints`` is a key of ``CONSTRAINT_SETS``. The constraint rows are
    weighted by ``w * sqrt(q)`` after row normalisation, as in
    ``P2FoldInterpolator``, so the energy per tetrahedron does not depend on
    the number of sample points.
    """
    method = METHODS[method] if isinstance(method, str) else method
    cells = np.asarray(cells) * method.refine
    step = (BOX[1] - BOX[0]) / cells
    if method.kind == "P1":
        support = TetMesh(origin=BOX[0], nsteps=cells + 1, step_vector=step)
        interpolator = P1Interpolator(support)
    else:
        support = P2UnstructuredTetMesh(origin=BOX[0], step_vector=step, nsteps_cells=cells)
        interpolator = P2Interpolator(support)
    interpolator.set_value_constraints(data)
    interpolator.setup_interpolator()
    if method.rule == "barycentre":
        points = support.barycentre
        tets = np.arange(len(points))
        q = np.ones(len(points))
    else:
        points, tets, q = tet_quadrature_points(support.nodes, support.elements, method.rule)
    keep = np.isin(tets, np.arange(0, len(support.elements), method.step))
    points, tets, q = points[keep], tets[keep], q[keep]
    sq = np.sqrt(q)
    _, _, _, grad0, tangent = fold.truth(points)
    axis = np.tile([0.0, 1.0, 0.0], (len(points), 1))
    w_orth, w_axis, w_norm, curved = CONSTRAINT_SETS[constraints]
    if w_orth > 0:
        _add_direction(interpolator, points, tets, grad0, w_orth * sq, 0.0, "orthogonal to c0")
    if w_axis > 0:
        _add_direction(interpolator, points, tets, axis, w_axis * sq, 0.0, "orthogonal to axis")
    if w_norm > 0:
        target = fold.norm_target(points) if curved else 1.0
        _add_direction(interpolator, points, tets, tangent, w_norm * sq, target, "norm along tangent")
    interpolator.solve_system("cg", solver_kwargs=dict(SOLVER))
    return interpolator


# (interpolator type, orthogonal_sample_rule, step) as passed to StructuralFrameBuilder
BUILDER_METHODS = {
    "FDI step10 (default)": ("FDI", None, 10),
    "PLI step10": ("PLI", None, 10),
    "P2 barycentre step10": ("P2", None, 10),
    "P2 centroid step1": ("P2", "centroid", 1),
    "P2 gauss4 step10": ("P2", "gauss4", 10),
    "P2 gauss4 step1": ("P2", "gauss4", 1),
}


class _Feature:
    """Adapts a model feature to the interpolator interface used by ``evaluate``."""

    def __init__(self, feature):
        self.feature = feature
        self.dof = int(getattr(feature.builder.interpolator, "dof", 0))

    def evaluate_value(self, xyz):
        return self.feature.evaluate_value(xyz)

    def evaluate_gradient(self, xyz):
        return self.feature.evaluate_gradient(xyz)


def build_coordinate1_with_builder(fold, method, data, nelements=800):
    """Build c1 with ``StructuralFrameBuilder`` and the known c0.

    ``method`` is a key of ``BUILDER_METHODS``. This is the path a model takes:
    the same orthogonality constraint, regularisation and solver as
    ``create_and_add_fold_frame``. Returns the second coordinate feature.
    """
    import pandas as pd
    from LoopStructural.geometry import BoundingBox
    from LoopStructural.modelling.features import LambdaGeologicalFeature
    from LoopStructural.modelling.features.builders import StructuralFrameBuilder
    from LoopStructural.utils.helper import all_heading

    interpolatortype, rule, step = BUILDER_METHODS[method]
    builder = StructuralFrameBuilder(
        interpolatortype=["FDI", interpolatortype, interpolatortype],
        bounding_box=BoundingBox(origin=BOX[0], maximum=BOX[1]),
        nelements=nelements,
        name="frame",
        buffer=0.0,
    )
    c0 = LambdaGeologicalFeature(
        lambda xyz: fold.truth(xyz)[0],
        name="c0",
        gradient_function=lambda xyz: fold.truth(xyz)[3],
    )
    builder.builders[0]._feature = c0
    builder.frame.features[0] = c0
    frame_data = pd.DataFrame(np.nan, index=range(len(data)), columns=all_heading())
    frame_data[["X", "Y", "Z", "val"]] = data
    frame_data["feature_name"], frame_data["coord"], frame_data["w"] = "frame", 1, 1.0
    builder.builders[1].add_data_from_data_frame(frame_data)
    builder.builders[2].add_data_from_data_frame(frame_data.iloc[:0])
    builder.build(w1=1.0, w3=0.0, orthogonal_sample_rule=rule, step=step)
    return builder.frame[1]


class Ramp:
    """The field ``c1 = x``, which ignores the folds. A method must beat it."""

    dof = 0

    def evaluate_value(self, xyz):
        return np.asarray(xyz)[:, 0].copy()

    def evaluate_gradient(self, xyz):
        return np.tile([1.0, 0.0, 0.0], (len(xyz), 1))


def evaluate(fold, interpolator, n=4000, seed=1):
    """Compare c1 with the truth. The scale and sign of c1 are free.

    Returns
    -------
    dict
        skill
            ``1 - rmse / rmse_ramp``, where rmse is the error of the best affine
            fit of the true c1 to the field and rmse_ramp the same for the ramp
            ``c1 = x``. Zero is no better than ignoring the folds, one is exact.
        stretch_bad
            fraction of points where the rate of change of c1 along the surface of
            c0, divided by its true value and by the median of that ratio, is
            outside 0.5 to 2. Counts reversals, collapse and blow up.
        reversal
            fraction of points where c1 decreases along the surface of c0.
        collapse
            fraction of points where ``|grad c1|`` is below 10% of its median.
        cos_med, cos_p95
            ``|cos|`` between ``grad c1`` and ``grad c0``. Zero is c1 constant along
            the normals of the folded surfaces.
        iso_reversal
            fraction of steps along traced isosurfaces of c0 (offsets -0.5, 0, 0.5)
            where c1 goes in the wrong direction.
        iso_stretch_cv
            coefficient of variation of the rate of change of c1 with distance
            along the same traces, divided by its true value. Zero is exactly
            proportional to the true c1 along the isosurface.
    """
    xyz = np.random.default_rng(seed).uniform(EVAL_LO, EVAL_HI, (n, 3))
    _, t, rate, grad0, tangent = fold.truth(xyz)
    value = interpolator.evaluate_value(xyz)
    gradient = interpolator.evaluate_gradient(xyz)
    ok = np.isfinite(value) & np.isfinite(gradient).all(axis=1)
    xyz, t, rate, grad0, tangent = (a[ok] for a in (xyz, t, rate, grad0, tangent))
    value, gradient = value[ok], gradient[ok]

    def rmse(x):
        return np.sqrt(np.mean((t - np.polyval(np.polyfit(x, t, 1), x)) ** 2))

    sign = 1.0 if np.polyfit(t, value, 1)[0] >= 0 else -1.0
    norm = np.linalg.norm(gradient, axis=1)
    along = sign * np.einsum("ij,ij->i", gradient, tangent)
    ratio = along / rate
    ratio = ratio / np.median(ratio)
    cosine = np.abs(np.einsum("ij,ij->i", gradient, grad0)) / np.maximum(norm, 1e-12)

    reversals, cvs = [], []
    for offset in (-0.5, 0.0, 0.5):
        line = fold.isosurface(offset)
        distance = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(line, axis=0), axis=1))]
        step = np.diff(distance)
        estimated = np.diff(sign * interpolator.evaluate_value(line)) / step
        true = np.diff(fold.truth(line)[1]) / step
        ok = np.isfinite(estimated)
        reversals.append(estimated[ok] <= 0)
        cvs.append(np.std(estimated[ok] / true[ok]) / abs(np.mean(estimated[ok] / true[ok])))
    return {
        "n": len(xyz),
        "skill": float(1 - rmse(value) / rmse(xyz[:, 0])),
        "stretch_bad": float(np.mean((ratio < 0.5) | (ratio > 2.0))),
        "reversal": float(np.mean(along <= 0)),
        "collapse": float(np.mean(norm < 0.1 * np.median(norm))),
        "cos_med": float(np.median(cosine)),
        "cos_p95": float(np.percentile(cosine, 95)),
        "iso_reversal": float(np.mean(np.concatenate(reversals))),
        "iso_stretch_cv": float(np.mean(cvs)),
        "dof": int(interpolator.dof),
    }


def compare(fold, methods, constraints, n_data, seeds, cells=(24, 3, 12)):
    """Metrics averaged over ``seeds`` random draws of the data."""
    rows = []
    for method in methods:
        runs = []
        for seed in seeds:
            interpolator = solve_coordinate1(
                fold, method, sample_data(fold, n_data, seed), constraints, cells
            )
            runs.append(evaluate(fold, interpolator))
        rows.append(
            {
                "method": method,
                "constraints": constraints,
                "n_data": n_data,
                **{k: float(np.mean([r[k] for r in runs])) for k in runs[0]},
                "worst_stretch_bad": float(max(r["stretch_bad"] for r in runs)),
            }
        )
    return rows


def compare_builder(fold, n_data, seeds, nelements=800):
    """Metrics for the builder path, averaged over ``seeds`` random draws."""
    rows = []
    for method in BUILDER_METHODS:
        runs = []
        for seed in seeds:
            feature = build_coordinate1_with_builder(
                fold, method, sample_data(fold, n_data, seed), nelements
            )
            runs.append(evaluate(fold, _Feature(feature)))
        rows.append(
            {
                "method": method,
                "n_data": n_data,
                **{k: float(np.mean([r[k] for r in runs])) for k in runs[0]},
                "worst_stretch_bad": float(max(r["stretch_bad"] for r in runs)),
            }
        )
    return rows


def print_header():
    print(
        f"{'':16s}{'n':>4s} {'method':22s}{'dof':>6s}{'skill':>7s}{'bad':>7s}{'worst':>7s}"
        f"{'rev':>7s}{'collapse':>9s}{'cos50':>7s}{'cos95':>7s}{'isoCV':>7s}"
    )


def print_row(r, label):
    print(
        f"{label:16s}{r['n_data']:4d} {r['method']:22s}{int(r['dof']):6d}{r['skill']:7.2f}"
        f"{r['stretch_bad']:7.3f}{r['worst_stretch_bad']:7.3f}{r['reversal']:7.3f}"
        f"{r['collapse']:9.3f}{r['cos_med']:7.2f}{r['cos_p95']:7.2f}{r['iso_stretch_cv']:7.2f}"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path, default=Path("validation/fold_frame_coordinate1"))
    parser.add_argument("--seeds", type=int, default=5)
    parser.add_argument("--builder", action="store_true", help="run through StructuralFrameBuilder")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    logging.disable(logging.CRITICAL)
    args.out.mkdir(parents=True, exist_ok=True)
    fold = SimilarFold()
    rows = []
    started = time.perf_counter()
    ramp = evaluate(fold, Ramp())
    print_row({"method": "ramp c1=x", "n_data": 0, "worst_stretch_bad": ramp["stretch_bad"], **ramp}, "")
    if args.builder:
        print_header()
        for n_data in (4, 12, 40):
            for r in compare_builder(fold, n_data, range(args.seeds)):
                print_row(r, "builder")
        print(f"total {time.perf_counter() - started:.0f} s")
        return
    for constraints in CONSTRAINT_SETS:
        for n_data in (4, 12, 40):
            rows += compare(fold, list(METHODS), constraints, n_data, range(args.seeds))
    print_header()
    for r in rows:
        print_row(r, r["constraints"])
    print(f"total {time.perf_counter() - started:.0f} s")
    import csv

    with open(args.out / "metrics.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    if not args.no_figures:
        make_figures(fold, args.out)


def make_figures(fold, out):
    """Section through c1 for each method, and where it fails to increase."""
    import matplotlib.pyplot as plt

    data = sample_data(fold, 8, 0)
    x = np.linspace(0.2, 11.8, 240)
    z = np.linspace(-2.8, 2.8, 120)
    xx, zz = np.meshgrid(x, z)
    section = np.c_[xx.ravel(), np.full(xx.size, 1.0), zz.ravel()]
    c0, t, _, _, tangent = fold.truth(section)
    valid = fold.valid(section)
    names = ["P1", "P1-fine", "P2-centroid", "P2-gauss4"]
    fig, axes = plt.subplots(2, 5, figsize=(20, 7), layout="constrained", sharex=True, sharey=True)
    panels = [("truth: c1 = arc length", np.where(valid, t, np.nan), None)]
    for name in names:
        interpolator = solve_coordinate1(fold, name, data, "orthogonal")
        value = interpolator.evaluate_value(section)
        gradient = interpolator.evaluate_gradient(section)
        along = np.einsum("ij,ij->i", gradient, tangent)
        along = along * (1.0 if np.nanmedian(along) >= 0 else -1.0)
        panels.append((name, np.where(valid, value, np.nan), np.where(valid, along, np.nan)))
    for column, (title, value, along) in enumerate(panels):
        ax = axes[0, column]
        image = ax.contourf(xx, zz, value.reshape(xx.shape), levels=24, cmap="viridis")
        ax.contour(xx, zz, np.where(valid, c0, np.nan).reshape(xx.shape), levels=12,
                   colors="white", linewidths=0.4)
        ax.set_title(title)
        fig.colorbar(image, ax=ax, shrink=0.7)
        ax = axes[1, column]
        if along is None:
            ax.axis("off")
            continue
        limit = np.nanpercentile(np.abs(along), 98)
        image = ax.pcolormesh(xx, zz, along.reshape(xx.shape), cmap="RdBu", vmin=-limit, vmax=limit,
                              shading="auto")
        ax.contour(xx, zz, along.reshape(xx.shape), levels=[0.0], colors="k", linewidths=0.8)
        ax.set_title("derivative of c1 along surface (red < 0)")
        fig.colorbar(image, ax=ax, shrink=0.7)
    fig.savefig(out / "coordinate1_sections.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
