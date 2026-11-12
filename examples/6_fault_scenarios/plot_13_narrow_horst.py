"""
6m. A narrow horst between two near-vertical faults
===================================================
A horst is a block that stays high between two normal faults that dip away
from each other. Here the faults are almost vertical (dip 89 degrees) and
only 0.04 apart at mid-depth (0.022 at the top, 0.058 at the base), so the
horst is a thin slab, about one cell of the interpolation support wide.

This example shows two things:

* The **interpolation** has no problem with the narrow block. A faulted
  series is interpolated in the restored position, where the stratigraphy
  is smooth and the horst is not a separate block. The field in the horst is
  correct even with a coarse support.
* The **surfaces** are a problem. Isosurfaces are extracted on a regular
  grid, and a grid cell that the faults cut mixes the values of three
  blocks. The pieces of the horizons in the horst are smeared or lost
  unless the output grid is very fine.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel
from LoopStructural.modelling.features.fault._fault_function import (
    Composite,
    FaultDisplacement,
    Ones,
    Zeros,
)
from LoopStructural.visualisation import Loop3DView

rng = np.random.default_rng(20261111)

######################################################################
# Helper functions
# ~~~~~~~~~~~~~~~~
# The same helpers as in the Y junction example.


def dip_vectors(dip, dip_direction):
    """Return the upward unit normal and the down-dip unit vector of a plane."""
    d, a = np.radians(dip), np.radians(dip_direction)
    horizontal = np.array([np.sin(a), np.cos(a), 0.0])
    normal = np.sin(d) * horizontal + np.array([0.0, 0.0, np.cos(d)])
    down_dip = np.cos(d) * horizontal - np.array([0.0, 0.0, np.sin(d)])
    return normal, down_dip


def planar_fault_data(name, center, dip, dip_direction, n=7, half_size=0.9):
    """Points on a planar fault, each with the fault normal (``coord = 0``)."""
    normal, down_dip = dip_vectors(dip, dip_direction)
    strike = np.cross(normal, down_dip)
    u, v = np.meshgrid(np.linspace(-half_size, half_size, n), np.linspace(-half_size, half_size, n))
    points = np.asarray(center) + u.reshape(-1, 1) * strike + v.reshape(-1, 1) * down_dip
    points = points[np.all((points > -0.1) & (points < 1.1), axis=1)]
    data = pd.DataFrame(points, columns=["X", "Y", "Z"])
    data[["nx", "ny", "nz"]] = normal
    data["val"] = 0.0
    data["coord"] = 0
    data["feature_name"] = name
    return data


def restore_through(faults_oldest_first):
    """Return a function that restores points through the faults, youngest first."""

    def restore(points):
        for fault in reversed(faults_oldest_first):
            points = fault.apply_to_points(points)
        return points

    return restore


def synthetic_observations(restore, field, horizons, feature_name, n=8000, half_width=0.015):
    """Random points near the horizons, with the value of the restored field."""
    points = rng.uniform(0.0, 1.0, size=(n, 3))
    values = field(restore(points))
    near = np.min(np.abs(values[:, None] - np.array(horizons)[None, :]), axis=1) < half_width
    data = pd.DataFrame(points[near], columns=["X", "Y", "Z"])
    data["val"] = values[near]
    data["feature_name"] = feature_name
    return data


constant_profile = FaultDisplacement(gx=Composite(Ones(), Zeros()), gy=Ones(), gz=Ones(), scale=1.0)

######################################################################
# Build the faults
# ~~~~~~~~~~~~~~~~
# The west fault dips toward -x and the east fault toward +x, so their
# hanging walls are on the outside and the horst stays high. Each fault
# has a displacement of 0.08. The faults do not displace each other.

model = GeologicalModel(np.zeros(3), np.ones(3))
faults = []
for name, x0, dip_direction in [("west", 0.48, 270), ("east", 0.52, 90)]:
    center = np.array([x0, 0.5, 0.5])
    normal, down_dip = dip_vectors(89, dip_direction)
    fault = model.create_and_add_fault(
        name,
        0.08,
        data=planar_fault_data(name, center, 89, dip_direction),
        fault_slip_vector=-down_dip,
        fault_normal_vector=normal,
        fault_center=center,
        major_axis=4.0,
        intermediate_axis=4.0,
        minor_axis=4.0,
        faultfunction=constant_profile,
        nelements=40_000,
        fault_buffer=0.5,
    )
    fault.steps = 1
    faults.append(fault)

######################################################################
# Synthetic stratigraphy with a coarse support
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

horizons = [0.25, 0.45, 0.65, 0.8]
strat_data = synthetic_observations(
    restore_through(faults), lambda p: p[:, 2], horizons, "strat", n=20000
)
strat = model.create_and_add_foliation("strat", data=strat_data, faults=faults, nelements=20_000)
print("interpolation cell size:", np.round(strat.interpolator.support.step_vector, 3))

######################################################################
# Section
# ~~~~~~~
# The section shows only the central part of the model (x from 0.4 to
# 0.6), with the horizontal scale exaggerated. The horizons in the horst
# are 0.08 above the blocks on each side, although the horst is about one
# interpolation cell wide.

x, z = np.meshgrid(np.linspace(0.4, 0.6, 201), np.linspace(0, 1, 201))
xyz = np.column_stack([x.ravel(), np.full(x.size, 0.5), z.ravel()])
field = strat.evaluate_value(xyz).reshape(x.shape)
fig, ax = plt.subplots(figsize=(6, 6))
ax.contourf(x, z, field, levels=20, cmap="viridis", alpha=0.6)
ax.contour(x, z, field, levels=horizons, colors="k", linewidths=1.5, linestyles="solid")
for name in ["west", "east"]:
    gx = model[name][0].evaluate_value(xyz).reshape(x.shape)
    ax.contour(x, z, gx, levels=[0], colors="r", linewidths=2)
ax.set_xlabel("x")
ax.set_ylabel("z")
ax.set_title("Section at y = 0.5 (x from 0.4 to 0.6)")
plt.show()

######################################################################
# Error in the horst
# ~~~~~~~~~~~~~~~~~~
# The horst did not move, so the exact restored height in the horst is the
# present-day height. Along the centre line of the horst the error is less
# than 0.001.

z = np.linspace(0.05, 0.95, 91)
centre_line = np.column_stack([np.full_like(z, 0.5), np.full_like(z, 0.5), z])
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(z, strat.evaluate_value(centre_line) - z)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("z (x = 0.5, y = 0.5)")
ax.set_ylabel("model - exact")
ax.set_title("Error on the centre line of the horst")
plt.show()

######################################################################
# Surfaces on a coarse and a fine output grid
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The isosurfaces are extracted on the grid of ``model.bounding_box``. On a
# coarse grid the cells that contain the horst also contain both faults, so
# the horst pieces of the horizons are smeared into the neighbouring blocks.
# On a finer grid they are better, but the cost increases with the cube of
# the resolution, over the whole model.

for nelements in [1e4, 4e5]:
    model.bounding_box.nelements = nelements
    print("output grid cell:", np.round(model.bounding_box.step_vector, 3))
    view = Loop3DView(model)
    view.plot_surface(strat, value=horizons, paint_with=strat, cmap="viridis")
    view.camera_position = "xz"
    view.camera.azimuth = 15
    view.camera.elevation = 15
    view.display()
