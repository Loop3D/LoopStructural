"""
5c. Fold constraints with sparse data
=====================================

A fold can be very hard to model from a few measurements: two contact points and one
orientation say little about a surface that folds back on itself. A plain interpolator,
which only asks for a smooth surface, will return something close to a plane.

The fold interpolators add the knowledge that the surface is folded. Given a *fold frame* (a
coordinate system with one axis across the fold and one along it) and a description of how the
limb dip changes across the fold, they constrain the *direction* of the gradient everywhere,
not only where there are data.

This example measures what that adds. The truth is a cylindrical fold, and each of six
interpolators is given the same three data: plain ``FDI``, ``PLI`` and ``P2``, and their
fold-constrained versions ``FDFold``, ``DFI`` (the fold version of PLI) and ``P2Fold``.

Real models get the fold frame from ``GeologicalModel.create_and_add_fold_frame``; see the
fold examples. Here we write it down analytically so that the answer is exactly known.
"""

import time

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
# The folded surface
# ------------------
# A cylindrical fold with its axis along ``y``: the horizons are ``s = z - h(x)`` with
# ``h(x) = A sin(2 pi x / wavelength)``, in the unit cube.

AMPLITUDE, WAVELENGTH = 0.12, 0.6
K = 2 * np.pi / WAVELENGTH


def height(x):
    return AMPLITUDE * np.sin(K * x)


def height_slope(x):
    return AMPLITUDE * K * np.cos(K * x)


def truth(xyz):
    return xyz[:, 2] - height(xyz[:, 0])


def truth_gradient(xyz):
    return np.c_[-height_slope(xyz[:, 0]), np.zeros(len(xyz)), np.ones(len(xyz))]


# %%
# Three data: two contact points on the 0.5 horizon, far apart, and one orientation in the middle.


def contact(x, y, level=0.5):
    z = brentq(lambda z: truth(np.array([[x, y, z]]))[0] - level, 0.0, 1.0)
    return [x, y, z, level]


values = np.array([contact(0.2, 0.3), contact(0.8, 0.7)])
orientation_xyz = np.array([[0.5, 0.5, 0.5]])
normals = np.c_[orientation_xyz, truth_gradient(orientation_xyz)]

# %%
# The fold frame and the fold event
# ---------------------------------
# The fold frame's first coordinate ``gx`` runs across the fold and the second, ``gy``, along
# the axis. For this fold they are just ``x`` and ``y``. A ``FoldEvent`` combines the frame with
# the fold axis and a function giving the limb rotation angle (the dip of the limb) in terms of
# ``gx``. For a cylindrical fold the axis is constant. The signs follow the library's
# convention for the direction of the frame's third axis. Note: the library may warn that ``fold_norm`` may need sign inversion; that warning is expected here because of this deliberate sign convention (the ``+ 180``/negated angles) and can be ignored.


class Coordinate:
    """Scalar field equal to one coordinate; its gradient is the unit vector."""

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
    return -np.rad2deg(np.arctan(height_slope(gx)))


def make_fold_event():
    return FoldEvent(
        FoldFrame(), fold_limb_rotation=limb_rotation, fold_axis=np.array([0.0, -1.0, 0.0])
    )


# %%
# Six interpolators
# -----------------
# FDI/FDFold use a grid of 16 cells per side, PLI/DFI a tetrahedral mesh of the same grid, and
# P2/P2Fold a mesh of quadratic tetrahedra with 8 cells per side.

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
    "FDFold": lambda: FDFoldInterpolator(grid(), fold=make_fold_event()),
    "DFI": lambda: DiscreteFoldInterpolator(tetrahedra(), fold=make_fold_event()),
    "P2Fold": lambda: P2FoldInterpolator(quadratic(), fold=make_fold_event()),
}


def interpolate(kind):
    interpolator = FACTORIES[kind]()
    interpolator.set_value_constraints(values)
    interpolator.set_normal_constraints(normals)
    start = time.perf_counter()
    interpolator.setup_interpolator()
    interpolator.solve_system("cg")
    return interpolator, time.perf_counter() - start


results = {kind: interpolate(kind) for kind in FACTORIES}

# %%
# Measuring the result
# --------------------
# On a regular set of points inside the cube compare each interpolator with the truth:
#
# * ``rms error``, the error in the scalar field (the position of the horizons),
# * ``median gradient angle``, how far the gradient direction (the dip of the surface) is from
#   the true one, in degrees.

axis = np.linspace(0.1, 0.9, 9)
points = np.array(np.meshgrid(axis, axis, axis, indexing="ij")).reshape(3, -1).T


def angle_degrees(a, b):
    cosine = np.sum(a * b, axis=1) / np.linalg.norm(a, axis=1) / np.linalg.norm(b, axis=1)
    return np.rad2deg(np.arccos(np.clip(cosine, -1, 1)))


rows = {}
for kind, (interpolator, elapsed) in results.items():
    error = interpolator.evaluate_value(points) - truth(points)
    angle = angle_degrees(interpolator.evaluate_gradient(points), truth_gradient(points))
    rows[kind] = {
        "unknowns": interpolator.dof,
        "rms error": np.sqrt(np.mean(error**2)),
        "median gradient angle": np.median(angle),
        "time (s)": elapsed,
    }
print(pd.DataFrame(rows).T.astype({"unknowns": int}).round(3).to_string())

# %%
# The section through the middle of the fold
# ------------------------------------------
# Dashed lines are the true horizons, red lines are each interpolator's. Contact points near
# the section are black dots.

x = np.linspace(0, 1, 100)
xx, zz = np.meshgrid(x, x)
section = np.c_[xx.ravel(), np.full(xx.size, 0.5), zz.ravel()]
true_field = truth(section).reshape(xx.shape)

fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True, sharey=True)
for ax, (kind, (interpolator, _)) in zip(axes.ravel(), results.items()):
    field = interpolator.evaluate_value(section).reshape(xx.shape)
    ax.contour(xx, zz, true_field, levels=[0.3, 0.5, 0.7], colors="k", linestyles="--", linewidths=1)
    ax.contour(xx, zz, field, levels=[0.3, 0.5, 0.7], colors="tab:red", linewidths=1.5)
    ax.plot(values[:, 0], values[:, 2], "ko", ms=5)
    ax.set_title(f"{kind}: gradient angle {rows[kind]['median gradient angle']:.0f} deg")
    ax.set_aspect("equal")
fig.tight_layout()
plt.show()

# %%
# What to take from this
# ----------------------
# * The three plain interpolators return straight, gently tilted horizons, about 45 degrees
#   out in the direction of the surface. With three data nothing in the problem says the surface
#   is folded.
# * All three fold interpolators bring the direction to within about 5 degrees. This is what
#   the fold constraint adds: the dip of the limbs everywhere, not only near the data.
# * Getting the direction right does not by itself put the horizons in the right place.
#   ``FDFold`` and ``DFI`` recover the folded shape but exaggerate its amplitude, and their
#   scalar-field error (0.11 and 0.12) is only modestly better than the plain interpolators
#   (about 0.16). ``P2Fold`` is the only one close to the truth in both respects. This example
#   does not investigate why; it is one set of data on one fold.
#
# The fold event here is exact because we wrote the geometry down. In a real model it is
# fitted to the structural data, so the quality of the fold constraint is only as good as that
# fit: see the fold examples, and the S-plot checks in the debugging guide.
