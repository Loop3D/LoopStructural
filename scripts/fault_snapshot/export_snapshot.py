"""Export a built GeologicalModel as a language-neutral "fault snapshot" (proof of concept).

The snapshot holds what a native oracle needs to answer the three kiridashi oracle
queries (``values``, ``fault_distance``, ``tip``) without Python:

* the node arrays and grid geometry of the three frame coordinates of every fault
  and of every exported series,
* the displacement profile of every fault as parameters, or as regular samples for a
  profile type without an exact form (``linearise``),
* the order of the faults of every series, and the constants of the fault-block
  extension (``LoopStructural/utils/_fault_block.py``).

It also writes reference results computed by the Python implementation, to test the
native one against.

Output layout (``out_dir``)::

    manifest.json        schema, constants, faults, series (see SNAPSHOT_SPEC in PROMPT.md)
    data/*.bin           little-endian float64 arrays, C order, no header
    reference/index.json what each reference file holds, with shapes
    reference/*.bin      reference results

Only the plain case is supported; anything else raises ``UnsupportedModel`` and names why:
structured-grid supports, no regions (no abutting faults), faults whose frames are
displaced only by faults that come earlier in the series restore order, and profiles
that are ``FaultDisplacement`` objects. Ones / Zeros / CubicFunction / Composite are exact;
any other profile is sampled at a regular interval and stored as a linear interpolant.

Run from the LoopStructural checkout (it uses ``tests/fixtures/fault_models.py``)::

    PYTHONPATH=. python scripts/fault_snapshot/export_snapshot.py OUT_DIR
"""

from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

SCHEMA = "loop-fault-snapshot"
VERSION = 1


class UnsupportedModel(Exception):
    """The model has a feature that the snapshot cannot represent."""


# --------------------------------------------------------------------------- profiles


# Samples of a linearised profile (see ``linearise``), from the low end of its range
# to the high end. Enough that the interval stays well below the frame stencils.
PROFILE_SAMPLES = 2001


def linearise(profile, span: tuple[float, float], samples: int = PROFILE_SAMPLES) -> dict:
    """Sample any profile at a regular interval over ``span``: a piecewise linear profile.

    ``span`` is the range of the frame coordinate that the profile is evaluated on (the
    range of the node values of its grid, which bounds every interpolated value). Beyond
    the span the value is held at the end sample (``np.interp``).
    """
    lo, hi = float(span[0]), float(span[1])
    if not (np.isfinite(lo) and np.isfinite(hi) and hi > lo):
        raise UnsupportedModel(f"profile {type(profile).__name__}: bad sampling span {span}")
    x = np.linspace(lo, hi, samples)
    y = np.asarray(profile(x), dtype=float).reshape(-1)
    if y.shape != x.shape or not np.all(np.isfinite(y)):
        raise UnsupportedModel(
            f"profile {type(profile).__name__}: sampling did not give finite values"
        )
    # The continuation across the fault (``continue_profile`` -> ``_taylor``) takes
    # finite differences of the profile on six points with a step finer than the
    # samples: keep the exact values of the profile there, by side (-1, +1).
    from LoopStructural.utils import _fault_block as fb

    def scalar(v):
        return float(np.asarray(profile(np.array([v], dtype=float)), dtype=float)[0])

    stencil = [
        [scalar(side * fb.BLOCK_SIDE_EPSILON + side * k * fb.BLOCK_EXTENSION_STEP) for k in range(6)]
        for side in (-1, 1)
    ]
    return {
        "type": "linear",
        "start": lo,
        "step": float((hi - lo) / (samples - 1)),
        "values": [float(v) for v in y],
        "stencil": stencil,  # [side -1, side +1] x f(x0 + side k h), k = 0..5
        "source": type(profile).__name__,  # informational: what was linearised
    }


def profile_to_dict(profile, span: tuple[float, float] | None = None) -> dict:
    """Serialise a fault profile function as parameters.

    A profile that is not one of the exact types is linearised over ``span`` (see
    ``linearise``); without a span it is unsupported.
    """
    from LoopStructural.modelling.features.fault._fault_function import (
        Composite,
        CubicFunction,
        Ones,
        Zeros,
    )

    if isinstance(profile, Ones):
        return {"type": "ones"}
    if isinstance(profile, Zeros):
        return {"type": "zeros"}
    if isinstance(profile, Composite):
        # each side is only evaluated on its own side of the fault
        pos = None if span is None else (0.0, max(float(span[1]), 0.0))
        neg = None if span is None else (min(float(span[0]), 0.0), 0.0)
        return {
            "type": "composite",
            "positive": profile_to_dict(profile.positive, pos),
            "negative": profile_to_dict(profile.negative, neg),
        }
    if isinstance(profile, CubicFunction):
        profile.solve()
        return {
            "type": "cubic",
            "w": [float(x) for x in profile.w],  # w0 v^3 + w1 v^2 + w2 v + w3
            "lim": [float(profile.lim[0]), float(profile.lim[1])],
            # beyond lim the value is the value at the limit
            "limit_value": [
                float(profile._limit_value(profile.lim[0])),
                float(profile._limit_value(profile.lim[1])),
            ],
            "min_v": float(profile.min_v),  # floor, applied last
            "max_v": float(profile.max_v),  # ceiling, applied last
        }
    if span is not None and callable(profile):
        return linearise(profile, span)
    raise UnsupportedModel(
        f"profile {type(profile).__name__} is not Ones, Zeros, CubicFunction or Composite"
    )


def displacement_to_dict(fault) -> dict:
    from LoopStructural.modelling.features.fault._fault_function import FaultDisplacement

    ff = fault.faultfunction
    if not isinstance(ff, FaultDisplacement):
        raise UnsupportedModel(f"fault {fault.name}: faultfunction is not a FaultDisplacement")

    def span(c: int) -> tuple[float, float]:
        # the frame coordinate is interpolated from its nodes: their range bounds it
        nodes = np.asarray(fault[c].interpolator.c, dtype=float)
        return float(np.nanmin(nodes)), float(np.nanmax(nodes))

    return {
        "scale": float(ff.scale),
        "gx": profile_to_dict(ff.gx, span(0)),
        "gy": profile_to_dict(ff.gy, span(1)),
        "gz": profile_to_dict(ff.gz, span(2)),
    }


# ------------------------------------------------------------------------------ grids


class Writer:
    """Writes float64 arrays as raw little-endian files and remembers their shapes."""

    def __init__(self, directory: Path):
        self.directory = directory
        directory.mkdir(parents=True, exist_ok=True)
        self.index: dict[str, dict] = {}

    def write(self, name: str, array, **info) -> str:
        a = np.ascontiguousarray(array, dtype="<f8")
        path = self.directory / f"{name}.bin"
        path.write_bytes(a.tobytes())
        self.index[name] = {"file": path.name, "shape": list(a.shape), **info}
        return path.name


def check_regions(feature) -> None:
    """Only erosional-unconformity regions are allowed.

    The snapshot evaluates a field without its regions: the field is continued above
    an unconformity. The mesher limits an older series (and the faults that the
    unconformity truncates) to below the unconformity itself, from the history.
    """
    from LoopStructural.modelling.features import FeatureType

    for r in feature.regions:
        if getattr(r, "type", None) != FeatureType.UNCONFORMITY:
            raise UnsupportedModel(
                f"{feature.name}: has a region {getattr(r, 'name', r)!r} that is not an "
                "erosional unconformity (abutting faults are not supported)"
            )


def grid_to_dict(feature, writer: Writer, name: str) -> dict:
    """Geometry and node values of the interpolator of one scalar field feature."""
    from loop_common.supports import StructuredGrid

    check_regions(feature)
    interpolator = feature.interpolator
    # the builder refines the grid lazily: finish it before reading anything
    feature.builder.up_to_date()
    interpolator.update()
    support = interpolator.support
    if type(support) is not StructuredGrid:
        raise UnsupportedModel(f"{feature.name}: support is {type(support).__name__}")
    if not np.allclose(support.rotation_xy, np.eye(3)):
        raise UnsupportedModel(f"{feature.name}: the grid is rotated")
    box = interpolator.bounding_box
    b = box.project(np.zeros((1, 3)))[0]
    a = (box.project(np.eye(3)) - b).T  # local = a @ world + b
    bvec = box.reproject_vectors(np.eye(3)).T  # world_vector = bvec @ local_vector
    nodes = np.asarray(interpolator.c, dtype=float)
    ncells = [int(k) for k in support.nsteps_cells]
    if nodes.shape[0] != (ncells[0] + 1) * (ncells[1] + 1) * (ncells[2] + 1):
        raise UnsupportedModel(f"{feature.name}: node array does not match the grid")
    return {
        "origin": [float(x) for x in support.origin],
        "step": [float(x) for x in support.step_vector],
        "ncells": ncells,
        "local_from_world": {"A": a.ravel().tolist(), "b": b.tolist()},
        "world_vector_from_local": bvec.ravel().tolist(),
        "nodes": "data/" + writer.write(name, nodes),
    }


# ---------------------------------------------------------------------------- export


def export_snapshot(model, out_dir, series: list[str]) -> dict:
    """Write ``manifest.json`` and ``data/`` for the given series. Returns the manifest."""
    from LoopStructural.utils import _fault_block as fb

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    data = Writer(out / "data")

    # make sure everything is built: evaluating once finishes the lazy builders
    probe = np.full((1, 3), 0.5)
    for s in series:
        model.evaluate_feature_value(s, probe)

    fault_names: list[str] = []
    series_entries = []
    for s in series:
        feature = model[s]
        check_regions(feature)
        names = [f.name for f in feature.faults]  # youngest first: the restore order
        for n in names:
            if n not in fault_names:
                fault_names.append(n)
        series_entries.append(
            {
                "name": s,
                "faults": names,
                "grid": grid_to_dict(feature, data, f"series_{s}"),
            }
        )

    faults = []
    for n in fault_names:
        fault = model[n]
        frame = []
        for c in range(3):
            frame.append(grid_to_dict(fault[c], data, f"fault_{n}_g{'xyz'[c]}"))
        frame_faults = sorted({f.name for c in range(3) for f in fault[c].faults})
        tip = fault.evaluate_tip_value(np.full((1, 3), 0.5))
        faults.append(
            {
                "name": n,
                "displacement": float(fault.displacement),  # already scaled to model units
                "steps": int(fault.steps),
                "profile": displacement_to_dict(fault),
                "frame": frame,  # gx, gy, gz
                "frame_faults": frame_faults,  # younger faults that displace this fault
                "has_tip": bool(not np.isinf(tip[0])),
                "tip_ends": {
                    str(k): [int(e) for e in v] for k, v in fault._tip_ends().items()
                },
            }
        )

    # A fault of a series is evaluated at points that younger faults of the series have
    # already restored. Its frame faults must all be among those (restored once).
    for entry in series_entries:
        done: list[str] = []
        for n in entry["faults"]:
            frame_faults = next(f for f in faults if f["name"] == n)["frame_faults"]
            if not set(frame_faults) <= set(done):
                raise UnsupportedModel(
                    f"series {entry['name']}: fault {n} is displaced by {frame_faults}, "
                    f"which do not all come before it in {entry['faults']}"
                )
            done.append(n)

    manifest = {
        "schema": SCHEMA,
        "version": VERSION,
        "node_order": "x fastest: index = i + (ncells[0]+1) * (j + (ncells[1]+1) * k)",
        "float_format": "little-endian float64, C order, no header",
        "constants": {
            "BLOCK_SIDE_EPSILON": fb.BLOCK_SIDE_EPSILON,
            "BLOCK_EXTENSION": fb.BLOCK_EXTENSION,
            "BLOCK_EXTENSION_EXACT": fb.BLOCK_EXTENSION_EXACT,
            "BLOCK_EXTENSION_LIMIT": fb.BLOCK_EXTENSION_LIMIT,
            "BLOCK_EXTENSION_STEP": fb.BLOCK_EXTENSION_STEP,
        },
        "gradient": {
            # evaluate_feature_block_gradient: P1 gradient of the field over a regular
            # tetrahedron centred on the point (LoopStructural.utils.maths), whose four
            # vertices are restored through the faults of the series like any point.
            "method": "regular_tetrahedron",
            "scale": float(np.min(model.bounding_box.step_vector) / 10),
        },
        "domain": {
            "min": [float(x) for x in model.bounding_box.origin],
            "max": [float(x) for x in model.bounding_box.maximum],
        },
        "faults": faults,
        "series": series_entries,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1))
    return manifest


# ------------------------------------------------------------------------- references


def test_points(model, fault_names, n=1500, n_surface=200, seed=11) -> np.ndarray:
    """Random points around the model box, and points on or very near each fault."""
    rng = np.random.default_rng(seed)
    lo = np.asarray(model.bounding_box.origin)
    hi = np.asarray(model.bounding_box.maximum)
    span = hi - lo
    pts = [rng.uniform(lo - 0.15 * span, hi + 0.15 * span, (n, 3))]
    for name in fault_names:
        fault = model[name]
        p = rng.uniform(lo + 0.1 * span, hi - 0.1 * span, (n_surface, 3))
        for _ in range(3):  # Newton steps onto gx = 0
            gx = fault.evaluate_frame_value(p, 0)
            g = fault.evaluate_frame_gradient(p, 0)
            p = p - (gx / np.einsum("ij,ij->i", g, g))[:, None] * g
        pts.append(p)  # exactly on the fault surface
        pts.append(p + rng.normal(0.0, 1e-4, p.shape))  # a hair either side
    return np.concatenate(pts)


def write_references(model, out_dir, series: list[str]) -> None:
    """Reference results of the Python implementation (the oracle of kiridashi/loop.py)."""
    from LoopStructural.utils._fault_block import restore_through

    out = Path(out_dir)
    ref = Writer(out / "reference")
    fault_names = [n for e in json.loads((out / "manifest.json").read_text())["faults"] for n in [e["name"]]]
    pts = test_points(model, fault_names)
    ref.write("points", pts, description="test points, world coordinates (N,3)")

    for n in fault_names:
        fault = model[n]
        for c in range(3):
            ref.write(
                f"fault_{n}_frame_g{'xyz'[c]}",
                fault.evaluate_frame_value(pts, c),
                description=f"evaluate_frame_value(points, {c}), NaN outside the grid",
            )
            ref.write(
                f"fault_{n}_frame_grad_g{'xyz'[c]}",
                fault.evaluate_frame_gradient(pts, c),
                description=f"evaluate_frame_gradient(points, {c}), world frame (N,3)",
            )
        grad = fault.evaluate_frame_gradient(pts, 0)
        norm = np.linalg.norm(grad, axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            dist = np.where(norm > 0.0, fault.evaluate_frame_value(pts, 0) / norm, np.nan)
        ref.write(f"fault_{n}_distance", dist, description="gx / |grad gx| (LoopOracle.fault_distance)")
        ref.write(f"fault_{n}_tip", fault.evaluate_tip_value(pts), description="evaluate_tip_value; inf = no tip")
        ref.write(
            f"fault_{n}_apply",
            fault.apply_to_points(pts),
            description="apply_to_points(points) with no block, no restored_through",
        )
        ref.write(
            f"fault_{n}_apply_reverse",
            fault.apply_to_points(pts, reverse=True),
            description="apply_to_points(points, reverse=True)",
        )

    for s in series:
        feature = model[s]
        faults = list(feature.faults)
        names = [f.name for f in faults]
        ref.write(f"series_{s}_value", model.evaluate_feature_value(s, pts), description="evaluate_feature_value, no block")
        # the points after each fault of the series restore order (restore_through semantics)
        for j in range(len(faults)):
            ref.write(
                f"series_{s}_restored_after_{j}",
                restore_through(pts, faults[: j + 1]),
                description=f"points restored through faults {names[: j + 1]} (youngest first)",
            )
        blocks = []
        for sides in itertools.product((-1, 0, 1), repeat=len(names)):
            block = dict(zip(names, sides))
            key = "_".join(f"{n}{'m' if v < 0 else 'p' if v > 0 else 'z'}" for n, v in block.items())
            ref.write(f"series_{s}_block_{key}_value", model.evaluate_feature_block_value(s, pts, block),
                      description=f"evaluate_feature_block_value, block {block}")
            ref.write(f"series_{s}_block_{key}_gradient", model.evaluate_feature_block_gradient(s, pts, block),
                      description=f"evaluate_feature_block_gradient, block {block} (N,3)")
            blocks.append({"key": key, "block": block})
        ref.index[f"series_{s}_blocks"] = {"blocks": blocks}
    (out / "reference" / "index.json").write_text(json.dumps(ref.index, indent=1))


# ---------------------------------------------------------------------------- demos


def build_demos():
    """Small models from the test fixtures: (name, model, series)."""
    from LoopStructural import GeologicalModel
    from tests.fixtures.fault_models import YOUNG, add_fault, add_layers
    from tests.unit.modelling.test_fault_block_extension import sloped_cubic_profile

    demos = []

    m = GeologicalModel(np.zeros(3), np.ones(3))
    f = add_fault(m, "f1", (0.5, 0.5, 0.5), 60, 90, 0.15, steps=10)
    add_layers(m, [f])
    demos.append(("single_constant", m, ["strat"]))

    m = GeologicalModel(np.zeros(3), np.ones(3))
    f = add_fault(m, "f1", **YOUNG, steps=10, profile=sloped_cubic_profile())
    add_layers(m, [f])
    demos.append(("single_cubic", m, ["strat"]))

    m = GeologicalModel(np.zeros(3), np.ones(3))
    main = add_fault(m, "main", (0.4, 0.5, 0.5), 60, 90, 0.1, steps=10)
    splay = add_fault(m, "splay", (0.75, 0.5, 0.6), 60, 270, 0.08, steps=10)
    add_layers(m, [splay, main])  # no abutting relation: two independent faults
    demos.append(("two_faults", m, ["strat"]))
    return demos


def main(argv: list[str]) -> int:
    import logging

    logging.disable(logging.WARNING)  # LoopStructural warns for every refined tetrahedron
    if len(argv) != 2:
        print(__doc__)
        return 2
    root = Path(argv[1])
    for name, model, series in build_demos():
        out = root / name
        t0 = time.perf_counter()
        export_snapshot(model, out, series)
        write_references(model, out, series)
        print(f"exported {name} -> {out} ({time.perf_counter() - t0:.0f} s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
