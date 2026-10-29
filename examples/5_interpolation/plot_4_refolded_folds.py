"""
5d. Refolded folds: when the fold axis changes
==============================================

In a refolded fold the axis of the first fold is itself folded by a later event, so the fold
axis is not constant: it plunges one way here and the other way a little further along. A fold
constraint that assumes a straight axis then pushes the surface the wrong way wherever the axis
has rotated.

This example builds the simplest refolded geometry, a dome-and-basin interference pattern, and
asks what the fold interpolators need to recover it from three contacts and six orientations:

* the early fold has its axis along ``y`` and shape ``sin(x)``,
* a later fold with axis along ``x`` and shape ``sin(y)`` tilts that axis, so its plunge
  changes along ``y``.

The scalar field is ``s = z - A2 sin(k2 x) - B3 sin(k3 y)``. The fold event describes this with
two functions: the limb rotation as a function of the coordinate across the fold, and the
*axis rotation* as a function of the coordinate along it. We compare an event that includes the
axis rotation with one that assumes a constant axis. For the sparse-data fold with a constant
axis see the previous example.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from loop_interpolation import (
    DiscreteFoldInterpolator,
    FDFoldInterpolator,
    FiniteDifferenceInterpolator,
    FoldEvent,
    P1Interpolator,
    P2FoldInterpolator,
    P2Interpolator,
    P2UnstructuredTetMesh,
    StructuredGrid,
    TetMesh,
)
from scipy.optimize import brentq

# %%
# The refolded surface
# --------------------

A2, WAVELENGTH2 = 0.10, 0.5  # early fold: amplitude and wavelength (axis along y)
B3, WAVELENGTH3 = 0.08, 0.8  # later fold: amplitude and wavelength (axis along x)
K2, K3 = 2 * np.pi / WAVELENGTH2, 2 * np.pi / WAVELENGTH3


def truth(xyz):
    return xyz[:, 2] - A2 * np.sin(K2 * xyz[:, 0]) - B3 * np.sin(K3 * xyz[:, 1])


def truth_gradient(xyz):
    return np.c_[
        -A2 * K2 * np.cos(K2 * xyz[:, 0]),
        -B3 * K3 * np.cos(K3 * xyz[:, 1]),
        np.ones(len(xyz)),
    ]


# %%
# Data: three contact points on the 0.5 horizon and six orientations at random places.


def contact(x, y, level=0.5):
    z = brentq(lambda z: truth(np.array([[x, y, z]]))[0] - level, 0.0, 1.0)
    return [x, y, z, level]


values = np.array([contact(0.25, 0.25), contact(0.75, 0.5), contact(0.4, 0.8)])
orientation_xyz = np.random.default_rng(5).uniform(0.15, 0.85, (6, 3))
normals = np.c_[orientation_xyz, truth_gradient(orientation_xyz)]

# %%
# The fold frame and two fold events
# ----------------------------------
# As in the previous example the fold frame is just ``gx = x`` (across the early fold) and
# ``gy = y`` (along it). The limb rotation is the dip of the early fold's limbs.
#
# The axis rotation is the plunge of the early fold's axis, which the later fold varies with
# ``gy``. The ``+ 180`` is a sign convention of the library: the frame's third axis must point
# against the orientation vectors. Note: the library may warn that ``fold_norm`` may need sign inversion; that warning is expected here because of this deliberate sign convention (the ``+ 180``/negated angles) and can be ignored.


class Coordinate:
    def __init__(self, axis):
        self.axis = axis

    def evaluate_value(self, xyz):
        return np.asarray(xyz)[:, self.axis].astype(float)

    def evaluate_gradient(self, xyz):
        gradient = np.zeros((len(xyz), 3))
        gradient[:, self.axis] = 1.0
        return gradient


class FoldFrame:
    features = (Coordinate(0), Coordinate(1))


def limb_rotation(gx):
    return -np.rad2deg(np.arctan(A2 * K2 * np.cos(K2 * gx)))


def axis_rotation(gy):
    return np.rad2deg(np.arctan(B3 * K3 * np.cos(K3 * gy))) + 180.0


def fold_with_rotating_axis():
    return FoldEvent(FoldFrame(), fold_limb_rotation=limb_rotation, fold_axis_rotation=axis_rotation)


def fold_with_constant_axis():
    return FoldEvent(
        FoldFrame(), fold_limb_rotation=limb_rotation, fold_axis=np.array([0.0, -1.0, 0.0])
    )


# %%
# Interpolators
# -------------
# The three plain interpolators, the three fold interpolators with the rotating axis, and the
# same two of them (``DFI`` and ``P2Fold``) with a constant axis.

CELLS = 16


def grid():
    return StructuredGrid(
        origin=np.zeros(3), nsteps_cells=np.full(3, CELLS), step_vector=np.full(3, 1 / CELLS)
    )


def tetrahedra():
    return TetMesh(
        origin=np.zeros(3), nsteps_cells=np.full(3, CELLS), step_vector=np.full(3, 1 / CELLS)
    )


def quadratic():
    return P2UnstructuredTetMesh(
        origin=np.zeros(3), nsteps_cells=np.full(3, 8), step_vector=np.full(3, 1 / 8)
    )


FACTORIES = {
    "FDI": lambda: FiniteDifferenceInterpolator(grid()),
    "PLI": lambda: P1Interpolator(tetrahedra()),
    "P2": lambda: P2Interpolator(quadratic()),
    "FDFold": lambda: FDFoldInterpolator(grid(), fold=fold_with_rotating_axis()),
    "DFI": lambda: DiscreteFoldInterpolator(tetrahedra(), fold=fold_with_rotating_axis()),
    "P2Fold": lambda: P2FoldInterpolator(quadratic(), fold=fold_with_rotating_axis()),
    "DFI (constant axis)": lambda: DiscreteFoldInterpolator(
        tetrahedra(), fold=fold_with_constant_axis()
    ),
    "P2Fold (constant axis)": lambda: P2FoldInterpolator(
        quadratic(), fold=fold_with_constant_axis()
    ),
}


def interpolate(kind):
    interpolator = FACTORIES[kind]()
    interpolator.set_value_constraints(values)
    interpolator.set_normal_constraints(normals)
    interpolator.setup_interpolator()
    interpolator.solve_system("cg")
    return interpolator


interpolators = {kind: interpolate(kind) for kind in FACTORIES}

# %%
# Results
# -------
# * ``rms error``: error in the scalar field over the interior,
# * ``median`` and ``90th percentile gradient angle``: how far the dip direction is from the
#   truth (degrees),
# * ``worst angle at the data``: the largest angle between the solution's gradient and one of
#   the six orientation measurements. A solution that respects the data should be near zero.

axis = np.linspace(0.1, 0.9, 9)
points = np.array(np.meshgrid(axis, axis, axis, indexing="ij")).reshape(3, -1).T


def angle_degrees(a, b):
    cosine = np.sum(a * b, axis=1) / np.linalg.norm(a, axis=1) / np.linalg.norm(b, axis=1)
    return np.rad2deg(np.arccos(np.clip(cosine, -1, 1)))


rows = {}
for kind, interpolator in interpolators.items():
    error = interpolator.evaluate_value(points) - truth(points)
    angle = angle_degrees(interpolator.evaluate_gradient(points), truth_gradient(points))
    at_data = angle_degrees(interpolator.evaluate_gradient(normals[:, :3]), normals[:, 3:])
    rows[kind] = {
        "rms error": np.sqrt(np.mean(error**2)),
        "median angle": np.median(angle),
        "p90 angle": np.percentile(angle, 90),
        "worst angle at the data": np.max(at_data),
    }
print(pd.DataFrame(rows).T.round(3).to_string())

# %%
# The height of the 0.5 horizon
# -----------------------------
# A map of the height ``z`` of the horizon ``s = 0.5`` for the truth, a plain interpolator, and
# the fold interpolator with a constant and with a rotating axis.

x = np.linspace(0.05, 0.95, 40)
xx, yy = np.meshgrid(x, x)
z_levels = np.linspace(0.0, 1.0, 80)


def horizon_height(field):
    """Height z at which ``field(x, y, z) = 0.5``, found by linear interpolation along z."""
    heights = np.full(xx.shape, np.nan)
    lattice = np.array(np.meshgrid(x, x, z_levels, indexing="ij")).reshape(3, -1).T
    values_ = np.asarray(field(lattice)).reshape(len(x), len(x), len(z_levels)) - 0.5
    for i in range(len(x)):
        for j in range(len(x)):
            column = values_[i, j]
            crossing = np.where(np.sign(column[:-1]) * np.sign(column[1:]) < 0)[0]
            if len(crossing):
                k = crossing[0]
                heights[j, i] = z_levels[k] + (z_levels[k + 1] - z_levels[k]) * column[k] / (
                    column[k] - column[k + 1]
                )
    return heights


panels = {
    "truth": truth,
    "FDI (no fold constraint)": interpolators["FDI"].evaluate_value,
    "P2Fold, constant axis": interpolators["P2Fold (constant axis)"].evaluate_value,
    "P2Fold, rotating axis": interpolators["P2Fold"].evaluate_value,
}
fig, axes = plt.subplots(1, 4, figsize=(15, 3.8), sharey=True)
for ax, (title, field) in zip(axes, panels.items()):
    image = ax.contourf(xx, yy, horizon_height(field), levels=np.linspace(0.35, 0.65, 13), extend="both")
    ax.plot(values[:, 0], values[:, 1], "ko", ms=4)
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_aspect("equal")
axes[0].set_ylabel("y")
fig.colorbar(image, ax=axes, label="height of the 0.5 horizon", shrink=0.8)
plt.show()

# %%
# What to take from this
# ----------------------
# * The plain interpolators miss the dome-and-basin pattern: the gradient direction is about 40
#   degrees out and the horizon is a smooth ramp.
# * With the axis rotation included, all three fold interpolators get the median direction to
#   within 2 to 10 degrees, and ``P2Fold`` reproduces the pattern in the map.
# * With a constant axis the median direction error is about 18 to 19 degrees for both ``DFI``
#   and ``P2Fold``, and the map shows what has gone wrong: the horizon is straight stripes
#   along ``y``. The later fold, which varies the axis, has been lost entirely.
# * A single error number can hide this. ``P2Fold`` with a constant axis has a scalar-field
#   error (0.065) as small as with the rotating axis (0.062). Look at the gradient direction
#   and at maps or sections, not only at the error.
# * The fold constraint can override the data. ``FDFold`` ends 40 degrees from one of the six
#   orientation measurements and ``DFI`` 14 degrees, while ``P2Fold`` stays within 6 degrees.
#   This example does not examine why. When a fold constraint is used it is worth checking the
#   solution against the measurements.
# * A fold constraint improves the direction without necessarily improving the scalar-field
#   error: ``DFI`` has a larger error (0.16) than the plain interpolators (about 0.10).
#
# The fold event can describe the limb rotation only as a function of the coordinate across the
# fold, so for this surface the constraint direction is itself a few degrees off in places.
# Treat these numbers as one example, not a benchmark.
