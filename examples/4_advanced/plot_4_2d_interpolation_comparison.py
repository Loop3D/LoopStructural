"""
============================================================
4d. Comparing scipy's RBF interpolator to LoopStructural 2D
============================================================
LoopStructural's discrete interpolators (piecewise linear "P1" and
piecewise quadratic "P2") are usually used on 3D tetrahedral meshes, but
the same interpolator classes also work on 2D triangulated meshes built
directly from a 2D bounding box.

This example compares that 2D interpolation against
:class:`scipy.interpolate.RBFInterpolator`, a widely used method for
interpolating scattered data with a global radial basis function. Both
approaches take a set of scattered (x, y, value) observations and
produce a continuous scalar field - the classic scattered-data
interpolation problem - but they make very different trade-offs.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RBFInterpolator

from LoopStructural.geometry import BoundingBox
from LoopStructural.interpolators import InterpolatorFactory
from LoopStructural.utils import rng

##############################################################################
# Test function and scattered samples
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Franke's function is a standard benchmark for scattered-data
# interpolation: smooth almost everywhere but with enough local structure
# (two bumps and a dip) that no low-order interpolator reproduces it
# exactly from a handful of samples.


def franke(x, y):
    term1 = 0.75 * np.exp(-((9 * x - 2) ** 2 + (9 * y - 2) ** 2) / 4)
    term2 = 0.75 * np.exp(-((9 * x + 1) ** 2) / 49 - (9 * y + 1) / 10)
    term3 = 0.5 * np.exp(-((9 * x - 7) ** 2 + (9 * y - 3) ** 2) / 4)
    term4 = -0.2 * np.exp(-((9 * x - 4) ** 2) - (9 * y - 7) ** 2)
    return term1 + term2 + term3 + term4


n_samples = 60
sample_xy = rng.random((n_samples, 2))
sample_val = franke(sample_xy[:, 0], sample_xy[:, 1])

# fine regular grid to evaluate and compare all three interpolants on
nx = ny = 100
gx, gy = np.meshgrid(np.linspace(0, 1, nx), np.linspace(0, 1, ny))
grid_xy = np.array([gx.flatten(), gy.flatten()]).T
true_val = franke(grid_xy[:, 0], grid_xy[:, 1]).reshape(ny, nx)

##############################################################################
# scipy RBFInterpolator
# ~~~~~~~~~~~~~~~~~~~~~~
# RBFInterpolator fits a global radial basis function so that the surface
# passes exactly through every sample point. It has no concept of a mesh -
# every evaluation is a weighted sum over *all* of the sample points.

rbf = RBFInterpolator(sample_xy, sample_val, kernel="thin_plate_spline")
rbf_val = rbf(grid_xy).reshape(ny, nx)

##############################################################################
# LoopStructural P1 and P2 interpolators
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LoopStructural instead triangulates the bounding box and solves a
# (sparse) least-squares system for the coefficients at each mesh node,
# combining the value constraints with a smoothing regularisation term.
# Increasing ``nelements`` gives the mesh more freedom to follow the data.

# pad the mesh slightly beyond [0, 1] so that evaluation points sitting
# exactly on the domain edge are safely inside an element rather than
# right on the mesh boundary. BoundingBox.with_buffer() would normally do
# this (via the interpolator factory's buffer= argument) but it doesn't
# yet support 2D bounding boxes, so the padding is done directly here.
bounding_box = BoundingBox(
    origin=np.array([-0.05, -0.05]), maximum=np.array([1.05, 1.05]), dimensions=2
)

value_constraints = np.hstack(
    [sample_xy, sample_val[:, None], np.ones((n_samples, 1))]
)

p1_interpolator = InterpolatorFactory.create_interpolator("P1", bounding_box, nelements=2000)
p1_interpolator.set_value_constraints(value_constraints)
p1_interpolator.setup_interpolator(regularisation=0.1)
p1_interpolator.solve_system(solver="lsmr")
p1_val = p1_interpolator.evaluate_value(grid_xy).reshape(ny, nx)

p2_interpolator = InterpolatorFactory.create_interpolator("P2", bounding_box, nelements=1000)
p2_interpolator.set_value_constraints(value_constraints)
p2_interpolator.setup_interpolator(regularisation=0.1)
p2_interpolator.solve_system(solver="lsmr")
p2_val = p2_interpolator.evaluate_value(grid_xy).reshape(ny, nx)

##############################################################################
# Visual comparison
# ~~~~~~~~~~~~~~~~~~

fig, axs = plt.subplots(2, 3, figsize=(18, 11))
levels = np.linspace(true_val.min(), true_val.max(), 15)

for ax, values, title in zip(
    axs[0],
    [true_val, rbf_val, p1_val],
    ["Franke's function (truth)", "scipy RBFInterpolator", "LoopStructural P1 (2D)"],
):
    cf = ax.contourf(gx, gy, values, levels=levels, cmap="viridis")
    ax.scatter(sample_xy[:, 0], sample_xy[:, 1], c="k", s=8)
    ax.set_title(title)
    fig.colorbar(cf, ax=ax, shrink=0.8)

error_levels = np.linspace(0, 0.3, 13)
axs[1, 0].axis("off")
for ax, values, title in zip(
    axs[1, 1:],
    [rbf_val, p1_val],
    ["RBF error", "P1 error"],
):
    err = np.abs(values - true_val)
    cf = ax.contourf(gx, gy, err, levels=error_levels, cmap="magma")
    ax.set_title(f"{title} (RMSE={np.sqrt(np.mean(err**2)):.3f})")
    fig.colorbar(cf, ax=ax, shrink=0.8)

# P2 gets its own row-2 slot too, swap it in over the blank axis
axs[1, 0].axis("on")
cf = axs[1, 0].contourf(gx, gy, p2_val, levels=levels, cmap="viridis")
axs[1, 0].scatter(sample_xy[:, 0], sample_xy[:, 1], c="k", s=8)
axs[1, 0].set_title("LoopStructural P2 (2D)")
fig.colorbar(cf, ax=axs[1, 0], shrink=0.8)

plt.tight_layout()
plt.show()

print("RMSE against Franke's function:")
print(f"  scipy RBF (thin_plate_spline): {np.sqrt(np.mean((rbf_val - true_val) ** 2)):.4f}")
print(f"  LoopStructural P1:             {np.sqrt(np.mean((p1_val - true_val) ** 2)):.4f}")
print(f"  LoopStructural P2:             {np.sqrt(np.mean((p2_val - true_val) ** 2)):.4f}")

##############################################################################
# Discussion
# ~~~~~~~~~~
# **scipy's RBFInterpolator**
#
# * Solves a dense ``n_samples x n_samples`` linear system - exact through
#   every point, but that cost grows quickly and the system can become
#   ill-conditioned as the number of samples grows or points cluster
#   together.
# * No mesh is involved, so there's no meaningful way to add a smoothing/
#   regularisation term, or to constrain gradients or normals - only
#   point values.
# * Trivial to set up for a one-off scattered-data fit.
#
# **LoopStructural's P1/P2 interpolators**
#
# * Solve a sparse least-squares system over mesh nodes, so cost scales
#   with the *mesh* resolution rather than the number of data points -
#   this is what makes it practical to combine thousands of geological
#   observations with a fine model resolution.
# * Value constraints are blended with a regularisation term
#   (``regularisation=`` above) rather than honoured exactly, which is
#   useful when data is noisy but means the fit isn't forced through
#   every sample point.
# * Can also take gradient and gradient-norm constraints natively - the
#   feature LoopStructural actually needs this interpolation machinery
#   for, since geological observations (bedding orientations, fault
#   planes) are as often directional as they are point values.
# * P2's quadratic shape functions let it follow curved structure with a
#   coarser mesh than P1 needs. Its regularisation combines the same
#   edge-jump smoothing P1 uses with a curvature-minimising term
#   (``minimise_grad_steepness``) that P1 doesn't need - both are applied
#   automatically by ``setup_interpolator()`` above.
#
# In short: RBF is a convenient, exact fit for smallish scattered
# datasets with no directional information; LoopStructural's discrete
# interpolators trade exactness at the sample points for scalability and
# the ability to fold in the directional constraints that dominate real
# geological datasets.
