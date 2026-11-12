"""
6g. A graben: two facing normal faults
======================================
A graben is a block that dropped down between two normal faults that dip
toward each other. In this model the two faults dip 60 degrees and meet a
little below the base of the model, so the floor of the graben is only
about 0.22 wide at z = 0.

Both faults have a constant slip vector, so each fault moves its whole
hanging-wall half-space. The west fault is in the hanging wall of the east
fault, so the east fault also moves the west fault. The faults must
therefore have an age relation, even if they are of the same age in
nature: here the east fault is the younger fault, and the west fault is
built with ``faults=[east]``.

.. warning::
   Like the crossing faults example, this model uses
   ``fault_network_workaround`` for the fault that the other fault
   displaces (see :ref:`fault-scenario-future-work`). With fault profiles
   that decay away from the fault (``BaseFault``), the faults do not need an
   age relation if the hanging wall of one fault does not reach the other
   fault. LoopStructural has no way to declare two faults as the same age
   (conjugate faults that both stop at their intersection).
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
# The same helpers as in the crossing faults example.


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


def fault_network_workaround(fault):
    """Restore once through a displaced fault, with a fast slip direction."""
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
# Build the east fault (the younger fault)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

model = GeologicalModel(np.zeros(3), np.ones(3))

east_center = np.array([0.85, 0.5, 0.5])
east_points, east_normal = planar_fault_points(east_center, 60, 270)
_, east_down_dip = dip_vectors(60, 270)
east = model.create_and_add_fault(
    "east",
    0.10,
    data=fault_data("east", east_points, east_normal),
    fault_slip_vector=-east_down_dip,
    fault_normal_vector=east_normal,
    fault_center=east_center,
    major_axis=4.0,
    intermediate_axis=4.0,
    minor_axis=4.0,
    faultfunction=constant_profile,
    nelements=20_000,
    fault_buffer=0.5,
)
east.steps = 1

######################################################################
# Build the west fault (displaced by the east fault)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The west fault is a plane in its restored position. It is moved forward
# through the east fault to give the present-day observations.

west_center = np.array([0.15, 0.5, 0.5])
west_points, west_normal = planar_fault_points(west_center, 60, 90)
_, west_down_dip = dip_vectors(60, 90)
west_points = move_forward(west_points, [east])
west_points = west_points[np.abs((west_points - east_center) @ east_normal) > 0.04]
west = model.create_and_add_fault(
    "west",
    0.10,
    data=fault_data("west", west_points, west_normal),
    fault_slip_vector=-west_down_dip,
    fault_normal_vector=west_normal,
    fault_center=move_forward(west_center[None, :], [east])[0],
    major_axis=4.0,
    intermediate_axis=4.0,
    minor_axis=4.0,
    faultfunction=constant_profile,
    faults=[east],
    nelements=20_000,
    fault_buffer=0.5,
)
west.steps = 1
fault_network_workaround(west)

######################################################################
# Synthetic stratigraphy
# ~~~~~~~~~~~~~~~~~~~~~~

horizons = [0.25, 0.45, 0.65, 0.8]
restore = restore_through([west, east])
strat_data = synthetic_observations(restore, lambda p: p[:, 2], horizons, "strat")
strat = model.create_and_add_foliation(
    "strat", data=strat_data, faults=[west, east], nelements=25_000
)

######################################################################
# 3D view
# ~~~~~~~

view = Loop3DView(model)
view.plot_surface(west, value=[0], colour="grey", opacity=0.6)
view.plot_surface(east, value=[0], colour="darkred", opacity=0.6)
view.plot_surface(strat, value=horizons, paint_with=strat, cmap="viridis")
view.camera_position = "xz"
view.camera.azimuth = 30
view.camera.elevation = 20
view.display()

######################################################################
# Section
# ~~~~~~~
# Most of the section is in the graben: the faults leave the model through
# its sides above z = 0.6 to 0.75. The throw of each fault is
# :math:`0.10 \sin 60^\circ = 0.087`.
#
# A constant profile moves the **whole** hanging-wall half-space. The east
# fault moved everything to the west of it (the graben and the west flank),
# and the west fault moved everything to the east of it (the graben and the
# east flank). So the graben dropped 0.17, and each flank dropped 0.087.
# The relative offsets are those of a graben, but the flanks are not fixed.
# A ``BaseFault`` profile with a decay zone that is narrower than the
# graben keeps the flanks fixed.

plot_section(model, {"strat": horizons}, ["west", "east"])
plt.show()

######################################################################
# Check against the exact solution
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Both faults are planar in their restored position and move rigidly, so
# the exact restored height is easy to calculate: restore through the east
# fault, then through the west fault in its restored position.


def exact_restored_height(points):
    p = points.copy()
    in_east_hanging_wall = (p - east_center) @ east_normal > 0
    p[in_east_hanging_wall] -= 0.10 * east_down_dip
    in_west_hanging_wall = (p - west_center) @ west_normal > 0
    p[in_west_hanging_wall] -= 0.10 * west_down_dip
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

######################################################################
# Graben width
# ~~~~~~~~~~~~
# The distance between the two faults decreases with depth. A block that
# narrows to a few cells of the interpolation support is difficult to model
# and to mesh (see also the narrow horst example).

x = np.linspace(0, 1, 2001)
for z0 in [1.0, 0.5, 0.0]:
    line = np.column_stack([x, np.full_like(x, 0.5), np.full_like(x, z0)])
    inside = (west[0].evaluate_value(line) > 0) & (east[0].evaluate_value(line) > 0)
    print(f"graben width at z = {z0}: {inside.sum() * (x[1] - x[0]):.3f}")
