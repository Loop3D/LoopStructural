"""
6e. An antithetic fault that abuts a main fault (Y junction)
============================================================
Faults often end on other faults. Here an antithetic fault in the hanging
wall of a main normal fault meets the main fault at depth. Below the
junction the antithetic fault does not exist, so the two faults make a
**Y junction**, and the model has three fault blocks.

LoopStructural models this with :code:`FaultSegment.add_abutting_fault`.
The younger (abutting) fault is limited to one side of the older fault: its
frame coordinates are only defined where the older fault's ``gx`` has the
given sign, so it neither exists nor moves points on the other side.
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
# ``restore_through``, which restores points through several faults,
# youngest first (the order that LoopStructural uses).


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


constant_profile = FaultDisplacement(gx=Composite(Ones(), Zeros()), gy=Ones(), gz=Ones(), scale=1.0)

######################################################################
# Build the main fault
# ~~~~~~~~~~~~~~~~~~~~
# Both faults have a constant slip vector (rigid hanging walls), so all of
# the offset is at the faults.

model = GeologicalModel(np.zeros(3), np.ones(3))

main_center = np.array([0.4, 0.5, 0.5])
main_normal, main_down_dip = dip_vectors(60, 90)
main = model.create_and_add_fault(
    "main",
    0.12,
    data=planar_fault_data("main", main_center, 60, 90),
    fault_slip_vector=-main_down_dip,
    fault_normal_vector=main_normal,
    fault_center=main_center,
    major_axis=4.0,
    intermediate_axis=4.0,
    minor_axis=4.0,
    faultfunction=constant_profile,
    nelements=20_000,
    fault_buffer=0.5,
)
main.steps = 1

######################################################################
# Build the antithetic fault and make it abut the main fault
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The antithetic fault dips 70 degrees toward -x, toward the main fault.
# ``positive=True`` keeps it on the side of the main fault where
# ``gx > 0`` (the hanging wall). If ``positive`` is not given,
# LoopStructural uses the side where most of the fault observations are.

branch_center = np.array([0.75, 0.5, 0.6])
branch_normal, branch_down_dip = dip_vectors(70, 270)
branch = model.create_and_add_fault(
    "branch",
    0.06,
    data=planar_fault_data("branch", branch_center, 70, 270),
    fault_slip_vector=-branch_down_dip,
    fault_normal_vector=branch_normal,
    fault_center=branch_center,
    major_axis=4.0,
    intermediate_axis=4.0,
    minor_axis=4.0,
    faultfunction=constant_profile,
    nelements=20_000,
    fault_buffer=0.5,
)
branch.steps = 1
branch.add_abutting_fault(main, positive=True)

######################################################################
# Synthetic stratigraphy
# ~~~~~~~~~~~~~~~~~~~~~~
# The series lists its faults oldest first. LoopStructural restores the
# points through the youngest fault first.

horizons = [0.25, 0.45, 0.65, 0.8]
restore = restore_through([main, branch])
strat_data = synthetic_observations(restore, lambda p: p[:, 2], horizons, "strat")
strat = model.create_and_add_foliation(
    "strat", data=strat_data, faults=[main, branch], nelements=25_000
)

######################################################################
# 3D view
# ~~~~~~~
# The antithetic fault surface ends on the main fault.

view = Loop3DView(model)
view.plot_surface(main, value=[0], colour="grey", opacity=0.5)
view.plot_surface(branch, value=[0], colour="darkred", opacity=0.6)
view.plot_surface(strat, value=horizons, paint_with=strat, cmap="viridis")
view.camera_position = "xz"
view.camera.azimuth = 30
view.camera.elevation = 20
view.display()

######################################################################
# Section
# ~~~~~~~
# The trace of the antithetic fault is only drawn where it exists. The
# junction is near z = 0.17. The east block is in the hanging wall of the
# main fault only: it moved down 0.10. The wedge between the two faults is
# in the hanging wall of both faults: it moved down 0.10 with the main
# fault and then 0.056 more on the antithetic fault, so it is a small
# graben.

plot_section(model, {"strat": horizons}, ["main", "branch"])
plt.show()

######################################################################
# Fault blocks
# ~~~~~~~~~~~~
# The side of each fault (the sign of ``gx``) gives the fault block. The
# antithetic fault has no value (NaN) in the footwall of the main fault.

x, z = np.meshgrid(np.linspace(0, 1, 201), np.linspace(0, 1, 201))
xyz = np.column_stack([x.ravel(), np.full(x.size, 0.5), z.ravel()])
main_side = main[0].evaluate_value(xyz) > 0
branch_side = np.nan_to_num(branch[0].evaluate_value(xyz), nan=-1.0) > 0
block = np.where(~main_side, 0, np.where(branch_side, 2, 1)).reshape(x.shape)
fig, ax = plt.subplots(figsize=(6, 6))
ax.contourf(x, z, block, levels=[-0.5, 0.5, 1.5, 2.5], colors=["#c9d6df", "#f3d9a4", "#e8a598"])
ax.text(0.12, 0.5, "footwall\nof main", ha="center")
ax.text(0.5, 0.75, "hanging wall\nof both faults", ha="center")
ax.text(0.86, 0.4, "hanging wall\nof main,\nfootwall\nof branch", ha="center")
ax.set_aspect("equal")
ax.set_title("Three fault blocks")
plt.show()
