"""
6f. Crossing faults: a younger fault offsets an older fault
==========================================================
When a younger fault cuts an older fault, the older fault surface is itself
offset. The model then has four fault blocks, and the older fault is two
separate patches.

In LoopStructural the older fault is built with ``faults=[younger]``. Its
frame is then interpolated in the restored position (before the younger
fault moved), exactly like a faulted stratigraphic series. The series lists
both faults, oldest first.

.. warning::
   This scenario needs two workarounds in the current version of
   LoopStructural (see :ref:`fault-scenario-future-work`):

   * **Double restoration.** When a series is restored, LoopStructural
     restores a point through the younger fault, and then evaluates the
     older fault frame, which restores the point through the younger fault
     again. Points in the hanging wall of the younger fault are then tested
     against the wrong position of the older fault (about 6 % of the points
     of this model).
   * **Speed.** The slip direction of a faulted fault is calculated with
     small tetrahedra that are moved through the younger fault, which costs
     some milliseconds for each point.

   The helper ``fault_network_workaround`` below fixes both for one fault.
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
# The same helpers as in the other fault scenario examples, and
# ``move_forward``, which moves restored points to their present-day
# position. It inverts the restoration by fixed-point iteration: it finds
# the point ``p`` whose restored position is the given point.


def dip_vectors(dip, dip_direction):
    """Return the upward unit normal and the down-dip unit vector of a plane."""
    d, a = np.radians(dip), np.radians(dip_direction)
    horizontal = np.array([np.sin(a), np.cos(a), 0.0])
    normal = np.sin(d) * horizontal + np.array([0.0, 0.0, np.cos(d)])
    down_dip = np.cos(d) * horizontal - np.array([0.0, 0.0, np.sin(d)])
    return normal, down_dip


def planar_fault_points(center, dip, dip_direction, n=9, half_size=0.9):
    """Points on a planar fault and the fault normal."""
    normal, down_dip = dip_vectors(dip, dip_direction)
    strike = np.cross(normal, down_dip)
    u, v = np.meshgrid(np.linspace(-half_size, half_size, n), np.linspace(-half_size, half_size, n))
    points = np.asarray(center) + u.reshape(-1, 1) * strike + v.reshape(-1, 1) * down_dip
    return points, normal


def fault_data(name, points, normal):
    """Fault observations (``coord = 0``) inside the model, with a normal each."""
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


def move_forward(points, faults_oldest_first, iterations=6):
    """Present-day position of restored points (inverse of the restoration)."""
    restore = restore_through(faults_oldest_first)
    present = np.array(points, dtype=float)
    for _ in range(iterations):
        present += points - restore(present)
    return present


def synthetic_observations(restore, field, horizons, feature_name, n=8000, half_width=0.015):
    """Random points near the horizons, with the value of the restored field."""
    points = rng.uniform(0.0, 1.0, size=(n, 3))
    values = field(restore(points))
    near = np.min(np.abs(values[:, None] - np.array(horizons)[None, :]), axis=1) < half_width
    data = pd.DataFrame(points[near], columns=["X", "Y", "Z"])
    data["val"] = values[near]
    data["feature_name"] = feature_name
    return data


def plot_section(model, series, faults, y=0.5, ax=None, title=None):
    """Vertical section at constant y. ``series`` maps a series name to its horizon values."""
    x, z = np.meshgrid(np.linspace(0, 1, 201), np.linspace(0, 1, 201))
    xyz = np.column_stack([x.ravel(), np.full(x.size, y), z.ravel()])
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    for name, values in series.items():
        field = model[name].evaluate_value(xyz).reshape(x.shape)
        ax.contourf(x, z, field, levels=20, cmap="viridis", alpha=0.6)
        ax.contour(
            x, z, field, levels=sorted(values), colors="k", linewidths=1.5, linestyles="solid"
        )
    for name in faults:
        gx = model[name][0].evaluate_value(xyz).reshape(x.shape)
        ax.contour(x, z, gx, levels=[0], colors="r", linewidths=2)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.set_title(title or f"Section at y = {y}")
    return ax


######################################################################
# Workaround for faults that other faults displace
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ``restore_once``: when a series restores a point through this fault, the
# point is already restored through the younger faults, so the frame of
# this fault is evaluated with its own faults disabled.
#
# ``restored_slip_direction``: the slip direction is the gradient of ``gy``
# at the restored point. This is exact where the younger faults move rigidly
# (as here) and an approximation in their drag zone.


def fault_network_workaround(fault):
    frame = [fault[i] for i in range(3)]

    def restore_once(points, reverse=False, precomputed_gx=None):
        for feature in frame:
            feature.faults_enabled = False
        try:
            return type(fault).apply_to_points(fault, points, reverse=reverse)
        finally:
            for feature in frame:
                feature.faults_enabled = True

    def restored_slip_direction(points, ignore_regions=False, element_scale_parameter=None):
        gy = frame[1]
        gy.builder.up_to_date()
        return gy.interpolator.evaluate_gradient(gy._apply_faults(np.asarray(points, dtype=float)))

    fault.apply_to_points = restore_once
    frame[1].evaluate_gradient = restored_slip_direction


constant_profile = FaultDisplacement(gx=Composite(Ones(), Zeros()), gy=Ones(), gz=Ones(), scale=1.0)

######################################################################
# Build the younger fault
# ~~~~~~~~~~~~~~~~~~~~~~~
# The faults are built youngest first, because the older fault needs the
# younger fault to restore its observations. The younger fault dips 55
# degrees toward -x.

model = GeologicalModel(np.zeros(3), np.ones(3))

young_center = np.array([0.6, 0.5, 0.5])
young_points, young_normal = planar_fault_points(young_center, 55, 270)
_, young_down_dip = dip_vectors(55, 270)
young = model.create_and_add_fault(
    "young",
    0.08,
    data=fault_data("young", young_points, young_normal),
    fault_slip_vector=-young_down_dip,
    fault_normal_vector=young_normal,
    fault_center=young_center,
    major_axis=4.0,
    intermediate_axis=4.0,
    minor_axis=4.0,
    faultfunction=constant_profile,
    nelements=20_000,
    fault_buffer=0.5,
)
young.steps = 1

######################################################################
# Observations of the older fault
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Before the younger fault moved, the older fault was a plane dipping 60
# degrees toward +x. Its present-day observations are the plane moved
# forward through the younger fault. Points within 0.04 of the younger
# fault are removed, because their side of the younger fault is uncertain.

old_center = np.array([0.4, 0.5, 0.5])
old_points, old_normal = planar_fault_points(old_center, 60, 90)
_, old_down_dip = dip_vectors(60, 90)
old_points = move_forward(old_points, [young])
far_from_young = np.abs((old_points - young_center) @ young_normal) > 0.04
old_points = old_points[far_from_young]
old_center = move_forward(old_center[None, :], [young])[0]

######################################################################
# Build the older fault
# ~~~~~~~~~~~~~~~~~~~~~
# ``faults=[young]`` makes the younger fault displace the older fault.

old = model.create_and_add_fault(
    "old",
    0.10,
    data=fault_data("old", old_points, old_normal),
    fault_slip_vector=-old_down_dip,
    fault_normal_vector=old_normal,
    fault_center=old_center,
    major_axis=4.0,
    intermediate_axis=4.0,
    minor_axis=4.0,
    faultfunction=constant_profile,
    faults=[young],
    nelements=20_000,
    fault_buffer=0.5,
)
old.steps = 1
fault_network_workaround(old)

######################################################################
# Synthetic stratigraphy
# ~~~~~~~~~~~~~~~~~~~~~~
# The layers were flat before both faults. The series lists the faults
# oldest first.

horizons = [0.25, 0.45, 0.65, 0.8]
restore = restore_through([old, young])
strat_data = synthetic_observations(restore, lambda p: p[:, 2], horizons, "strat")
strat = model.create_and_add_foliation(
    "strat", data=strat_data, faults=[old, young], nelements=25_000
)

######################################################################
# 3D view
# ~~~~~~~

view = Loop3DView(model)
view.plot_surface(young, value=[0], colour="darkred", opacity=0.6)
view.plot_surface(old, value=[0], colour="grey", opacity=0.6)
view.plot_surface(strat, value=horizons, paint_with=strat, cmap="viridis")
view.camera_position = "xz"
view.camera.azimuth = 30
view.camera.elevation = 20
view.display()

######################################################################
# Section
# ~~~~~~~
# The trace of the older fault is offset where it crosses the younger
# fault. Each horizon has four pieces: one in each fault block.

plot_section(model, {"strat": horizons}, ["old", "young"])
plt.show()

######################################################################
# Check against the exact solution
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Both faults are planar and move rigidly, so the exact restored height of
# any point is easy to calculate. The error of the model is small
# everywhere except within a few cells of the faults.


def exact_restored_height(points):
    p = points.copy()
    in_young_hanging_wall = (p - young_center) @ young_normal > 0
    p[in_young_hanging_wall] -= 0.08 * young_down_dip
    in_old_hanging_wall = (p - np.array([0.4, 0.5, 0.5])) @ old_normal > 0
    p[in_old_hanging_wall] -= 0.10 * old_down_dip
    return p[:, 2]


x, z = np.meshgrid(np.linspace(0, 1, 201), np.linspace(0, 1, 201))
xyz = np.column_stack([x.ravel(), np.full(x.size, 0.5), z.ravel()])
error = np.abs(strat.evaluate_value(xyz) - exact_restored_height(xyz))
print(
    f"median error {np.nanmedian(error):.1e}; points with error > 0.02: {np.mean(error > 0.02):.2%}"
)

fig, ax = plt.subplots(figsize=(6, 6))
c = ax.pcolormesh(x, z, error.reshape(x.shape), vmin=0, vmax=0.05, cmap="Reds")
fig.colorbar(c, ax=ax, label="|model - exact|")
ax.set_aspect("equal")
ax.set_title("Error of the stratigraphic field")
plt.show()
