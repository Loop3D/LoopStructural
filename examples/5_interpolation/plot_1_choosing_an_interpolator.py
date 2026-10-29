"""
5a. Choosing an interpolator
============================

LoopStructural can solve for the same implicit function on three kinds of support:

* ``FiniteDifferenceInterpolator`` (FDI) on a regular grid,
* ``P1Interpolator`` (PLI) on a mesh of linear tetrahedra,
* ``P2Interpolator`` on a mesh of quadratic tetrahedra.

They minimise the same kind of misfit (fit the data, keep the second derivative small) but
discretise it differently, so the same data give slightly different surfaces. To see how much
that matters we need to know the right answer, so this example builds a synthetic geological
model with a known scalar field: a gaussian dome. The data are sampled from it, then each
interpolator is asked to recover it, and we measure the error.

The same comparison, on several synthetic models, is what the interpolation regression tests
(``packages/loop_interpolation/tests/golden``) run.
"""

import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq

from LoopStructural.interpolators import (
    FiniteDifferenceInterpolator,
    P1Interpolator,
    P2Interpolator,
    P2UnstructuredTetMesh,
    StructuredGrid,
    TetMesh,
)

# %%
# A synthetic model with a known answer
# -------------------------------------
# The stratigraphic scalar field is the height above a gaussian dome, ``s = z - dome(x, y)``.
# Horizons are the surfaces ``s = 0.3, 0.5, 0.7``. Everything lives in the unit cube.

RELIEF, CENTRE, WIDTH = 0.2, 0.5, 0.2
HORIZONS = (0.3, 0.5, 0.7)


def truth(xyz):
    r2 = (xyz[:, 0] - CENTRE) ** 2 + (xyz[:, 1] - CENTRE) ** 2
    return xyz[:, 2] - RELIEF * np.exp(-r2 / (2 * WIDTH**2))


def truth_gradient(xyz):
    r2 = (xyz[:, 0] - CENTRE) ** 2 + (xyz[:, 1] - CENTRE) ** 2
    e = RELIEF * np.exp(-r2 / (2 * WIDTH**2)) / WIDTH**2
    return np.c_[e * (xyz[:, 0] - CENTRE), e * (xyz[:, 1] - CENTRE), np.ones(len(xyz))]


# %%
# Sample data from it the way a geologist would: contact points where nine boreholes cross
# each horizon, and five orientation measurements.


def borehole_contacts(xy):
    points = []
    for x, y in xy:
        for level in HORIZONS:

            def misfit(z, x=x, y=y, level=level):
                return truth(np.array([[x, y, z]]))[0] - level

            if misfit(0.0) * misfit(1.0) < 0:
                points.append([x, y, brentq(misfit, 0.0, 1.0), level])
    return np.array(points)


rng = np.random.default_rng(4)
borehole_xy = np.vstack(
    [[0.5, 0.5], [0.35, 0.5], [0.65, 0.5], [0.5, 0.35], [0.5, 0.65], rng.uniform(0.1, 0.9, (4, 2))]
)
values = borehole_contacts(borehole_xy)  # x, y, z, value
orientation_xyz = np.array(
    [[0.4, 0.4, 0.5], [0.6, 0.6, 0.5], [0.3, 0.7, 0.6], [0.8, 0.3, 0.4], [0.5, 0.5, 0.5]]
)
normals = np.c_[orientation_xyz, truth_gradient(orientation_xyz)]  # x, y, z, nx, ny, nz
print(f"{len(values)} contact points and {len(normals)} orientations")

# %%
# Three supports
# --------------
# FDI and PLI use a grid of 16 cells (17 nodes) per side. The quadratic mesh has far more
# unknowns per cell, so 8 cells per side gives a similar number of unknowns.

CELLS = 16
origin = np.zeros(3)


def make_interpolator(kind):
    if kind == "FDI":
        grid = StructuredGrid(
            origin=origin, nsteps_cells=np.full(3, CELLS), step_vector=np.full(3, 1 / CELLS)
        )
        return FiniteDifferenceInterpolator(grid)
    if kind == "PLI":
        mesh = TetMesh(
            origin=origin, nsteps_cells=np.full(3, CELLS), step_vector=np.full(3, 1 / CELLS)
        )
        return P1Interpolator(mesh)
    mesh = P2UnstructuredTetMesh(
        origin=origin, nsteps_cells=np.full(3, 8), step_vector=np.full(3, 1 / 8)
    )
    return P2Interpolator(mesh)


def interpolate(kind):
    interpolator = make_interpolator(kind)
    interpolator.set_value_constraints(values)
    interpolator.set_normal_constraints(normals)
    start = time.perf_counter()
    interpolator.setup_interpolator()
    interpolator.solve_system("cg")
    return interpolator, time.perf_counter() - start


results = {kind: interpolate(kind) for kind in ("FDI", "PLI", "P2")}
interpolators = {kind: interpolator for kind, (interpolator, _) in results.items()}

# %%
# How good is each one?
# ---------------------
# Evaluate on a regular set of points well inside the cube and compare with the truth:
#
# * ``rms error`` of the scalar field over the whole interior,
# * ``crest error``, the same near the top of the dome, where smoothing hurts most,
# * ``data misfit``, how closely the surface passes through the contact points,
# * ``median gradient angle``, how far the gradient direction is from the true one,
# * ``time``, to build and solve the system (indicative only, it depends on the machine).

axis = np.linspace(0.1, 0.9, 9)
points = np.array(np.meshgrid(axis, axis, axis, indexing="ij")).reshape(3, -1).T
crest = points[(np.abs(points[:, 0] - CENTRE) < 0.15) & (np.abs(points[:, 1] - CENTRE) < 0.15)]


def angle_degrees(a, b):
    cosine = np.sum(a * b, axis=1) / np.linalg.norm(a, axis=1) / np.linalg.norm(b, axis=1)
    return np.rad2deg(np.arccos(np.clip(cosine, -1, 1)))


rows = {}
for kind, interpolator in interpolators.items():
    error = interpolator.evaluate_value(points) - truth(points)
    crest_error = interpolator.evaluate_value(crest) - truth(crest)
    misfit = interpolator.evaluate_value(values[:, :3]) - values[:, 3]
    angle = angle_degrees(interpolator.evaluate_gradient(points), truth_gradient(points))
    rows[kind] = {
        "unknowns": interpolator.dof,
        "rms error": np.sqrt(np.mean(error**2)),
        "crest error": np.sqrt(np.mean(crest_error**2)),
        "data misfit": np.sqrt(np.mean(misfit**2)),
        "median gradient angle": np.median(angle),
        "time (s)": results[kind][1],
    }
print(pd.DataFrame(rows).T.astype({"unknowns": int}).round(4).to_string())

# %%
# Looking at a section
# --------------------
# The section through the middle of the dome (``y = 0.5``). The dashed lines are the true
# horizons and the solid lines are each interpolator's.

x = np.linspace(0, 1, 120)
z = np.linspace(0, 1, 120)
xx, zz = np.meshgrid(x, z)
section = np.c_[xx.ravel(), np.full(xx.size, 0.5), zz.ravel()]
true_field = truth(section).reshape(xx.shape)
near_section = np.abs(values[:, 1] - 0.5) < 0.2

fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
for ax, (kind, interpolator) in zip(axes, interpolators.items()):
    field = interpolator.evaluate_value(section).reshape(xx.shape)
    ax.contour(xx, zz, true_field, levels=HORIZONS, colors="k", linestyles="--", linewidths=1)
    ax.contour(xx, zz, field, levels=HORIZONS, colors="tab:red", linewidths=1.5)
    ax.plot(values[near_section, 0], values[near_section, 2], "ko", ms=4)
    ax.set_title(f"{kind}: rms error {rows[kind]['rms error']:.3f}")
    ax.set_xlabel("x")
    ax.set_aspect("equal")
axes[0].set_ylabel("z")
fig.suptitle("True horizons (dashed) and interpolated horizons (red), section y = 0.5")
fig.tight_layout()
plt.show()

# %%
# What to take from this
# ----------------------
# No interpolator wins on every measure here:
#
# * Over the whole interior P2 has the lowest error and the lowest gradient-angle error, then
#   PLI, then FDI.
# * Near the crest the order changes: FDI has the lowest error and PLI the highest.
# * P2 passes through the contact points (misfit about zero) and FDI has the largest misfit,
#   so FDI is the one that trades the data for smoothness.
# * The quadratic mesh is the slowest to solve, FDI the fastest. FDI is also the default
#   used by ``GeologicalModel`` (``interpolatortype="FDI"``).
#
# These results come from one synthetic model, one random draw of data and one resolution,
# so treat the ranking as a property of this example, not a rule. Errors of a few percent of
# the box height are normal for sparse data: the interpolator is inventing the surface between
# the data. More data, or a finer support, reduce them. The next example looks at the other
# lever, the regularisation weight.
