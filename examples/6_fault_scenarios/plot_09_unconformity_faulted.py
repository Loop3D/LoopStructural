"""
6i. Faults before and after an angular unconformity
===================================================
The relative age of a fault and an unconformity controls which rocks the
fault cuts:

1. Basement layers were tilted 20 degrees and cut by an **early fault**.
2. The basement was eroded, and flat cover layers were deposited on the
   erosion surface (an **angular unconformity**).
3. A **late fault** cut the cover, the basement and the early fault.

In LoopStructural, the faults that move a series are the faults in its
``faults`` list. The cover is moved only by the late fault; the basement is
moved by both faults. The unconformity is the zero value of the cover
series (``model.add_unconformity``): the basement exists only below it.

.. warning::
   * The early fault is displaced by the late fault, so this model uses
     ``fault_network_workaround`` (see :ref:`fault-scenario-future-work`).
   * LoopStructural does not cut a fault surface at an unconformity: the
     surface of the early fault continues up into the cover. It does not
     move the cover, because the cover series does not list it.
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
# ``synthetic_observations`` has one more option here, ``keep``, to keep
# only the points where a rock unit exists.


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


def synthetic_observations(
    restore, field, horizons, feature_name, n=8000, half_width=0.015, keep=None
):
    """Random points near the horizons, with the value of the restored field."""
    points = rng.uniform(0.0, 1.0, size=(n, 3))
    if keep is not None:
        points = points[keep(points)]
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
    for (name, values), cmap in zip(series.items(), ["viridis", "copper", "cividis"]):
        field = model[name].evaluate_value(xyz).reshape(x.shape)
        ax.contourf(x, z, field, levels=20, cmap=cmap, alpha=0.6)
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
# Stratigraphy before faulting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The cover value is the height above z = 0.6, and the unconformity is
# cover = 0. The basement value is the distance above a plane that dips 20
# degrees toward -x.


def cover_field(points):
    return points[:, 2] - 0.6


def basement_field(points):
    tilt = np.radians(20)
    return -np.sin(tilt) * points[:, 0] + np.cos(tilt) * points[:, 2] - 0.1


cover_horizons = [0.0, 0.15, 0.3]
basement_horizons = [0.05, 0.2, 0.35, 0.5]

######################################################################
# Build the late fault
# ~~~~~~~~~~~~~~~~~~~~

model = GeologicalModel(np.zeros(3), np.ones(3))

late_center = np.array([0.65, 0.5, 0.5])
late_points, late_normal = planar_fault_points(late_center, 65, 270)
_, late_down_dip = dip_vectors(65, 270)
late = model.create_and_add_fault(
    "late_fault",
    0.12,
    data=fault_data("late_fault", late_points, late_normal),
    fault_slip_vector=-late_down_dip,
    fault_normal_vector=late_normal,
    fault_center=late_center,
    major_axis=4.0,
    intermediate_axis=4.0,
    minor_axis=4.0,
    faultfunction=constant_profile,
    nelements=20_000,
    fault_buffer=0.5,
)
late.steps = 1

######################################################################
# Build the early fault
# ~~~~~~~~~~~~~~~~~~~~~
# The early fault was a plane before the late fault moved it. We use the
# observations below the unconformity only, as in a real data set (the
# early fault is not seen in the cover).

early_center = np.array([0.3, 0.5, 0.4])
early_points, early_normal = planar_fault_points(early_center, 55, 90, n=11)
_, early_down_dip = dip_vectors(55, 90)
early_points = early_points[cover_field(early_points) < 0.0]
early_points = move_forward(early_points, [late])
early_points = early_points[np.abs((early_points - late_center) @ late_normal) > 0.04]
early = model.create_and_add_fault(
    "early_fault",
    0.10,
    data=fault_data("early_fault", early_points, early_normal),
    fault_slip_vector=-early_down_dip,
    fault_normal_vector=early_normal,
    fault_center=move_forward(early_center[None, :], [late])[0],
    major_axis=4.0,
    intermediate_axis=4.0,
    minor_axis=4.0,
    faultfunction=constant_profile,
    faults=[late],
    nelements=20_000,
    fault_buffer=0.5,
)
early.steps = 1
fault_network_workaround(early)

######################################################################
# Cover series and the unconformity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Series are built youngest first. ``add_unconformity`` makes every series
# that is built after the cover exist only where ``cover < 0``.

cover_data = synthetic_observations(restore_through([late]), cover_field, cover_horizons, "cover")
cover = model.create_and_add_foliation("cover", data=cover_data, faults=[late], nelements=20_000)
model.add_unconformity(cover, 0.0)

######################################################################
# Basement series
# ~~~~~~~~~~~~~~~
# Basement observations exist only below the unconformity.


def below_unconformity(points):
    return cover_field(late.apply_to_points(points)) < 0.0


basement_data = synthetic_observations(
    restore_through([early, late]),
    basement_field,
    basement_horizons,
    "basement",
    n=12000,
    keep=below_unconformity,
)
basement = model.create_and_add_foliation(
    "basement", data=basement_data, faults=[early, late], nelements=25_000
)

######################################################################
# 3D view
# ~~~~~~~

view = Loop3DView(model)
view.plot_surface(late, value=[0], colour="darkred", opacity=0.6)
view.plot_surface(early, value=[0], colour="grey", opacity=0.6)
view.plot_surface(cover, value=cover_horizons, colour="tan")
view.plot_surface(basement, value=basement_horizons, paint_with=basement, cmap="viridis")
view.camera_position = "xz"
view.camera.azimuth = 30
view.camera.elevation = 20
view.display()

######################################################################
# Section
# ~~~~~~~
# * The cover (above the unconformity) is offset only by the late fault.
# * The basement is offset by both faults. The early fault is offset by the
#   late fault.
# * Erosion happened after the early fault. On the upthrown side of the
#   early fault (its footwall) more basement was eroded, so the upper
#   basement horizons are missing there.
# * The trace of the early fault continues above the unconformity. In
#   nature it ends at the unconformity.

plot_section(
    model, {"cover": cover_horizons, "basement": basement_horizons}, ["early_fault", "late_fault"]
)
plt.show()
