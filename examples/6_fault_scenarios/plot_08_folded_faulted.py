"""
6h. Folded strata cut by a normal fault and a reverse fault
===========================================================
A more realistic structural framework: upright folds, cut after folding by
a normal fault and a reverse fault. The folds trend 30 degrees oblique to
the faults, so the cut-off lines of the horizons on the faults are curved,
and the angle between a horizon and a fault changes along the fault.

Both faults use the ``BaseFault`` profile: the throw varies along strike
(``major_axis``), and both walls are dragged over 0.3 from each fault
(``minor_axis = 0.6``). The decay zone of the reverse fault reaches the
normal fault, so the reverse fault is younger and displaces the normal
fault.

The folds are given directly by the synthetic data (a sinusoidal folded
surface), and the stratigraphy is interpolated as a faulted foliation. To
model the folds with a fold frame, see the fold examples (section 2).

.. warning::
   This model uses ``fault_network_workaround`` for the normal fault (see
   :ref:`fault-scenario-future-work`). In the drag zone of the reverse
   fault, the fast slip direction of the workaround is an approximation.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel
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


######################################################################
# Folded stratigraphy before faulting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The stratigraphic value is the height above a folded reference surface:
# wavelength 0.5, amplitude 0.07, fold axis trending 30 degrees from y.


def folded_field(points):
    s = points[:, 0] * np.cos(np.radians(30)) + points[:, 1] * np.sin(np.radians(30))
    return points[:, 2] - 0.5 - 0.07 * np.sin(2 * np.pi * s / 0.5)


horizons = [-0.3, -0.15, 0.0, 0.15, 0.3]

######################################################################
# Build the reverse fault (the younger fault)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The reverse fault dips 45 degrees toward +x. Its hanging wall moves
# up-dip, so the restoration direction (``fault_slip_vector``) is down-dip.

model = GeologicalModel(np.zeros(3), np.ones(3))

reverse_center = np.array([0.8, 0.5, 0.5])
reverse_points, reverse_normal = planar_fault_points(reverse_center, 45, 90)
_, reverse_down_dip = dip_vectors(45, 90)
reverse = model.create_and_add_fault(
    "reverse",
    0.10,
    data=fault_data("reverse", reverse_points, reverse_normal),
    fault_slip_vector=reverse_down_dip,
    fault_normal_vector=reverse_normal,
    fault_center=reverse_center,
    major_axis=0.9,
    intermediate_axis=4.0,
    minor_axis=0.6,
    faultfunction="BaseFault",
    nelements=20_000,
    fault_buffer=0.5,
)
reverse.steps = 10

######################################################################
# Build the normal fault (displaced by the reverse fault)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The normal fault dips 60 degrees toward +x. Its observations are moved
# forward through the reverse fault. In the drag zone of the reverse fault
# the normal fault surface is bent.

normal_center = np.array([0.35, 0.5, 0.5])
normal_points, normal_normal = planar_fault_points(normal_center, 60, 90)
_, normal_down_dip = dip_vectors(60, 90)
normal_points = move_forward(normal_points, [reverse])
normal_points = normal_points[np.abs((normal_points - reverse_center) @ reverse_normal) > 0.04]
normal = model.create_and_add_fault(
    "normal",
    0.12,
    data=fault_data("normal", normal_points, normal_normal),
    fault_slip_vector=-normal_down_dip,
    fault_normal_vector=normal_normal,
    fault_center=move_forward(normal_center[None, :], [reverse])[0],
    major_axis=0.9,
    intermediate_axis=4.0,
    minor_axis=0.6,
    faultfunction="BaseFault",
    faults=[reverse],
    nelements=20_000,
    fault_buffer=0.5,
)
normal.steps = 10
fault_network_workaround(normal)

######################################################################
# Synthetic stratigraphy
# ~~~~~~~~~~~~~~~~~~~~~~

restore = restore_through([normal, reverse])
strat_data = synthetic_observations(restore, folded_field, horizons, "strat", n=10000)
print(f"{len(strat_data)} stratigraphic observations")
strat = model.create_and_add_foliation(
    "strat", data=strat_data, faults=[normal, reverse], nelements=30_000
)

######################################################################
# 3D view
# ~~~~~~~

view = Loop3DView(model)
view.plot_surface(normal, value=[0], colour="grey", opacity=0.6)
view.plot_surface(reverse, value=[0], colour="darkred", opacity=0.6)
view.plot_surface(strat, value=horizons, paint_with=strat, cmap="viridis")
view.camera_position = "xz"
view.camera.azimuth = 30
view.camera.elevation = 20
view.display()

######################################################################
# Sections
# ~~~~~~~~
# The normal fault drops its hanging wall, the reverse fault lifts its
# hanging wall, and the horizons on both sides of each fault are dragged.
# The folds are oblique to the sections, so they move along the section
# from one section to the next. The throw is largest at y = 0.5 (the fault
# centres) and smaller toward y = 0 and y = 1.

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for ax, y in zip(axes, [0.2, 0.5, 0.8]):
    plot_section(model, {"strat": horizons}, ["normal", "reverse"], y=y, ax=ax)
plt.show()

######################################################################
# Map of a horizon
# ~~~~~~~~~~~~~~~~
# The height of the horizon ``strat = 0`` shows the curved cut-off lines.
# Where the reverse fault repeats the horizon, the map shows the highest
# piece.

x, y, z = np.meshgrid(
    np.linspace(0, 1, 61), np.linspace(0, 1, 61), np.linspace(0, 1, 61), indexing="ij"
)
grid = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
above = (strat.evaluate_value(grid) >= 0.0).reshape(x.shape)
crosses = above[:, :, 1:] != above[:, :, :-1]
top = crosses.shape[2] - 1 - np.argmax(crosses[:, :, ::-1], axis=2)
height = np.where(crosses.any(axis=2), z[0, 0, top], np.nan)
fig, ax = plt.subplots(figsize=(6, 5))
c = ax.pcolormesh(x[:, :, 0], y[:, :, 0], height, cmap="terrain", shading="auto")
fig.colorbar(c, ax=ax, label="height of strat = 0")
ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Map of the horizon strat = 0")
plt.show()
