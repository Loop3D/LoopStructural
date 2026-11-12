"""
6j. Displacement gradient: sharp, medium and broad drag
=======================================================
Two faults can have the same offset at the fault but very different
deformation near the fault. If the displacement decays to zero over a short
distance, the horizons next to the fault are strongly bent (drag). If it
decays over a long distance, the walls are almost rigid near the fault and
the deformation is spread over a wide zone.

With the ``BaseFault`` profile, the decay distance is half of the
``minor_axis``. This example builds the same fault three times with
``minor_axis`` = 0.12, 0.4 and 1.2 (decay over 0.06, 0.2 and 0.6) and
compares the results.

.. note::
   With ``BaseFault`` each wall moves half of the displacement, so the
   offset at the fault is the same for the three models. The decay is the
   same on both sides of the fault. To give each wall a different decay
   distance, use a custom profile (see example 3c, *Defining the fault
   displacement function*).
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
# One function builds the model for a given minor axis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The fault dips 60 degrees toward +x, with a displacement of 0.15. The
# major axis (1.6) makes the throw decrease a little toward y = 0 and
# y = 1. A narrow decay zone needs a fine interpolation support, so the
# number of elements of the stratigraphy is larger than in the other
# examples.

center = np.array([0.5, 0.5, 0.5])
normal, down_dip = dip_vectors(60, 90)
horizons = [0.25, 0.45, 0.65, 0.8]


def build_model(minor_axis):
    model = GeologicalModel(np.zeros(3), np.ones(3))
    fault = model.create_and_add_fault(
        "fault",
        0.15,
        data=planar_fault_data("fault", center, 60, 90),
        fault_slip_vector=-down_dip,
        fault_normal_vector=normal,
        fault_center=center,
        major_axis=1.6,
        intermediate_axis=4.0,
        minor_axis=minor_axis,
        faultfunction="BaseFault",
        nelements=20_000,
        fault_buffer=0.5,
    )
    fault.steps = 10
    data = synthetic_observations(fault.apply_to_points, lambda p: p[:, 2], horizons, "strat")
    model.create_and_add_foliation("strat", data=data, faults=[fault], nelements=40_000)
    return model


minor_axes = {"sharp": 0.12, "medium": 0.4, "broad": 1.2}
models = {label: build_model(minor) for label, minor in minor_axes.items()}

######################################################################
# Sections
# ~~~~~~~~
# The offset at the fault is the same in the three models. With the sharp
# profile the horizons are almost vertical next to the fault: the decay
# width (0.06) is less than the movement of each wall (0.075), so the drag
# zone is strongly sheared. With the broad profile the horizons are tilted
# over most of the section.

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for ax, (label, model) in zip(axes, models.items()):
    plot_section(
        model,
        {"strat": horizons},
        ["fault"],
        ax=ax,
        title=f"{label}: minor_axis = {minor_axes[label]}",
    )
plt.show()

######################################################################
# Movement along a horizontal line
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The vertical component of the movement on the line z = 0.5 through the
# fault centre. Each wall moves half of the displacement at the fault. The
# jump at the fault is the same; the gradient away from it is different.

x = np.linspace(0.0, 1.0, 801)
line = np.column_stack([x, np.full_like(x, 0.5), np.full_like(x, 0.5)])
fig, ax = plt.subplots(figsize=(8, 4))
for label, model in models.items():
    movement = line - model["fault"].apply_to_points(line)
    ax.plot(x, movement[:, 2], label=f"{label} ({minor_axes[label]})")
ax.axvline(center[0], color="r", lw=1)
ax.set_xlabel("x (y = 0.5, z = 0.5)")
ax.set_ylabel("vertical movement")
ax.legend(title="minor_axis")
plt.show()

######################################################################
# 3D view of the sharp model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# A narrow drag zone is a hard case for isosurface extraction: the cells
# that the fault cuts contain a large change of value, and the surfaces
# are smeared across the fault.

model = models["sharp"]
view = Loop3DView(model)
view.plot_surface(model["fault"], value=[0], colour="grey", opacity=0.5)
view.plot_surface(model["strat"], value=horizons, paint_with=model["strat"], cmap="viridis")
view.camera_position = "xz"
view.camera.azimuth = 30
view.camera.elevation = 20
view.display()
