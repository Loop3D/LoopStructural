"""
6d. A listric fault and hanging-wall rollover
=============================================
A listric fault is curved: its dip decreases with depth. When the hanging
wall slides down a curved fault, it rotates, and the hanging-wall layers
bend into a **rollover** anticline next to the fault.

LoopStructural moves points along the slip direction of the fault frame
(the gradient of ``gy``), which follows the curved fault surface. A curved
fault therefore needs no special treatment: give observations of the curved
surface and the displacement follows it.

In this example:

* the fault dip decreases from 72 degrees at the top to 22 degrees at the
  base;
* the fault trace is bowed in map view;
* the ``BaseFault3D`` profile makes the throw vary along strike and dip.
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
# The same helpers as in the other fault scenario examples.


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
# Listric fault geometry
# ~~~~~~~~~~~~~~~~~~~~~~
# The fault dips toward +x. Its dip changes linearly from 72 degrees at
# z = 1.05 to 22 degrees at z = -0.05. We integrate the dip to get the
# horizontal position as a function of depth, and add a bow of 0.08 to the
# trace in map view.

x_top, z_top, depth = 0.25, 1.05, 1.1
depths = np.linspace(0.0, depth, 401)
dips = np.radians(72 + (22 - 72) * depths / depth)
heave = np.concatenate([[0.0], np.cumsum(np.diff(depths) / np.tan(dips[1:]))])


def listric_point(y, s):
    """Position on the fault at map coordinate y and depth s below the top."""
    x = x_top + np.interp(s, depths, heave) + 0.08 * (1.0 - (2.0 * y - 1.0) ** 2)
    return np.column_stack([x, y, z_top - s])


def listric_normal(y, s, h=1e-4):
    """Unit normal of the fault, pointing up into the hanging wall."""
    along_strike = listric_point(y + h, s) - listric_point(y - h, s)
    along_dip = listric_point(y, s + h) - listric_point(y, s - h)
    n = np.cross(along_strike, along_dip)
    n /= np.linalg.norm(n, axis=1)[:, None]
    n[n[:, 2] < 0] *= -1.0
    return n


y, s = np.meshgrid(np.linspace(-0.1, 1.1, 11), np.linspace(0.0, depth, 15))
y, s = y.ravel(), s.ravel()
fault_points = listric_point(y, s)
fault_normals = listric_normal(y, s)
fault_data = pd.DataFrame(fault_points, columns=["X", "Y", "Z"])
fault_data[["nx", "ny", "nz"]] = fault_normals
fault_data["val"] = 0.0
fault_data["coord"] = 0
fault_data["feature_name"] = "fault"

######################################################################
# The fault centre, normal and slip direction are taken at mid-depth on the
# centre line. The slip direction is down the dip of the fault there; away
# from the centre LoopStructural follows the interpolated fault frame.

center = listric_point(np.array([0.5]), np.array([depth / 2]))[0]
center_normal = listric_normal(np.array([0.5]), np.array([depth / 2]))[0]
below = listric_point(np.array([0.5]), np.array([depth / 2 + 1e-3]))[0]
center_down_dip = (below - center) / np.linalg.norm(below - center)
print("fault centre:", np.round(center, 3))

######################################################################
# Build the model
# ~~~~~~~~~~~~~~~

model = GeologicalModel(np.zeros(3), np.ones(3))
fault = model.create_and_add_fault(
    "fault",
    0.15,
    data=fault_data,
    fault_slip_vector=-center_down_dip,
    fault_normal_vector=center_normal,
    fault_center=center,
    major_axis=0.9,
    intermediate_axis=1.0,
    minor_axis=0.8,
    faultfunction="BaseFault3D",
    nelements=30_000,
    fault_buffer=0.5,
)
fault.steps = 10

horizons = [0.3, 0.5, 0.7, 0.85]
strat_data = synthetic_observations(fault.apply_to_points, lambda p: p[:, 2], horizons, "strat")
strat = model.create_and_add_foliation("strat", data=strat_data, faults=[fault], nelements=30_000)

######################################################################
# 3D view
# ~~~~~~~

view = Loop3DView(model)
view.plot_surface(fault, value=[0], colour="grey", opacity=0.5)
view.plot_surface(strat, value=horizons, paint_with=strat, cmap="viridis")
view.camera_position = "xz"
view.camera.azimuth = 30
view.camera.elevation = 20
view.display()

######################################################################
# Section
# ~~~~~~~
# The fault flattens with depth, and the hanging-wall horizons bend down
# toward the fault (rollover). The movement follows the fault, so the
# deepest hanging-wall horizon moves almost horizontally.

ax = plot_section(model, {"strat": horizons}, ["fault"])
ax.plot(*listric_point(np.full(50, 0.5), np.linspace(0, depth, 50))[:, [0, 2]].T, "w--", lw=1)
plt.show()

######################################################################
# Slip direction
# ~~~~~~~~~~~~~~
# The arrows show the restoration direction of the hanging wall (the
# gradient of ``gy``) on the section. They are parallel to the fault and
# turn with it.

x, z = np.meshgrid(np.linspace(0.05, 0.95, 13), np.linspace(0.05, 0.95, 13))
xyz = np.column_stack([x.ravel(), np.full(x.size, 0.5), z.ravel()])
slip = fault[1].evaluate_gradient(xyz)
slip /= np.linalg.norm(slip, axis=1)[:, None]
ax = plot_section(model, {"strat": horizons}, ["fault"], title="Slip direction (gy gradient)")
ax.quiver(xyz[:, 0], xyz[:, 2], slip[:, 0], slip[:, 2], color="w", scale=25)
plt.show()
