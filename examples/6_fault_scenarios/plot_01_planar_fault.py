"""
6a. A planar fault with a constant slip vector
==============================================
The simplest faulted model: one planar normal fault that moves its hanging
wall rigidly, with no decay of displacement away from the fault. Every
horizon is offset by the same amount at every position on the fault, so the
fault is a pure step.

This example also introduces the workflow that all of the fault scenario
examples use:

1. Define the fault geometry (a few points on the fault surface with their
   normals) and the kinematics (displacement, slip vector, axes and
   displacement profile).
2. Build the fault first.
3. Make *synthetic* stratigraphic observations that agree with the fault
   kinematics: take random present-day points, restore them through the
   fault with :code:`fault.apply_to_points`, and give each point the value
   of a known (here flat) stratigraphy in the restored position.
4. Interpolate the faulted stratigraphy from these observations, and look at
   the result in 3D and on a vertical section.

Because the stratigraphy is known in the restored position, the only source
of offset in the model is the fault. This makes the effect of each fault
parameter easy to see.

All of the scenario examples use the unit cube :math:`[0, 1]^3` with z up.
The scenarios are taken from the fault meshing benchmark suite
(``fault_mesh_benchmarks``).
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
# Fault geometry
# ~~~~~~~~~~~~~~
# The fault dips 60 degrees toward +x and strikes along y. We calculate the
# unit normal (up, into the hanging wall) and the down-dip direction from the
# dip and the dip direction (measured clockwise from north, +y).


def dip_vectors(dip, dip_direction):
    """Return the upward unit normal and the down-dip unit vector of a plane."""
    d, a = np.radians(dip), np.radians(dip_direction)
    horizontal = np.array([np.sin(a), np.cos(a), 0.0])
    normal = np.sin(d) * horizontal + np.array([0.0, 0.0, np.cos(d)])
    down_dip = np.cos(d) * horizontal - np.array([0.0, 0.0, np.sin(d)])
    return normal, down_dip


fault_center = np.array([0.5, 0.5, 0.5])
fault_normal, fault_down_dip = dip_vectors(60, 90)
fault_strike = np.cross(fault_normal, fault_down_dip)

######################################################################
# Fault observations are a 7 x 7 grid of points on the fault plane, each with
# the fault normal. ``coord = 0`` marks them as observations of the fault
# surface (the first coordinate of the fault frame).

u, v = np.meshgrid(np.linspace(-0.9, 0.9, 7), np.linspace(-0.9, 0.9, 7))
fault_points = fault_center + u.reshape(-1, 1) * fault_strike + v.reshape(-1, 1) * fault_down_dip
fault_points = fault_points[np.all((fault_points > -0.1) & (fault_points < 1.1), axis=1)]
fault_data = pd.DataFrame(fault_points, columns=["X", "Y", "Z"])
fault_data[["nx", "ny", "nz"]] = fault_normal
fault_data["val"] = 0.0
fault_data["coord"] = 0
fault_data["feature_name"] = "fault"

######################################################################
# A constant displacement profile
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The displacement of a point is the product of three 1D profiles, one for
# each fault frame coordinate (``gx`` normal to the fault, ``gy`` along the
# slip direction, ``gz`` along strike). A rigid hanging wall is a step in
# ``gx`` (1 on the hanging wall, 0 on the footwall) and a constant 1 along
# ``gy`` and ``gz``.
#
# .. note::
#    ``FaultDisplacement`` multiplies the profiles by ``scale``, which is
#    0.5 by default (each wall of the default profile moves half of the
#    displacement). A one-sided profile must set ``scale=1.0``, or the
#    hanging wall moves only half of the displacement. LoopStructural has no
#    named constant (or "step") profile, and the profile classes are in a
#    private module. See :ref:`fault-scenario-future-work`.

constant_profile = FaultDisplacement(gx=Composite(Ones(), Zeros()), gy=Ones(), gz=Ones(), scale=1.0)
constant_profile.plot()

######################################################################
# Build the fault
# ~~~~~~~~~~~~~~~
# ``fault_slip_vector`` is the direction that *restores* the hanging wall
# (LoopStructural convention). For a normal fault the hanging wall moves
# down-dip, so the restoration direction is up-dip. ``displacement`` is the
# full offset along the slip vector.

model = GeologicalModel(np.zeros(3), np.ones(3))
fault = model.create_and_add_fault(
    "fault",
    0.15,
    data=fault_data,
    fault_slip_vector=-fault_down_dip,
    fault_normal_vector=fault_normal,
    fault_center=fault_center,
    major_axis=4.0,
    intermediate_axis=4.0,
    minor_axis=4.0,
    faultfunction=constant_profile,
    nelements=20_000,
    fault_buffer=0.5,
)
# A rigid translation needs only one integration step.
fault.steps = 1

######################################################################
# Synthetic stratigraphic observations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Before faulting the layers were flat, so the stratigraphic value is the
# height of a point in its restored position. We keep the points that are
# within 0.015 of four horizons, so the observations look like contacts
# picked on each side of the fault.

horizons = [0.25, 0.45, 0.65, 0.8]


def synthetic_observations(restore, field, feature_name, n=8000, half_width=0.015):
    points = rng.uniform(0.0, 1.0, size=(n, 3))
    values = field(restore(points))
    near = np.min(np.abs(values[:, None] - np.array(horizons)[None, :]), axis=1) < half_width
    data = pd.DataFrame(points[near], columns=["X", "Y", "Z"])
    data["val"] = values[near]
    data["feature_name"] = feature_name
    return data


strat_data = synthetic_observations(fault.apply_to_points, lambda p: p[:, 2], "strat")
print(f"{len(strat_data)} stratigraphic observations")

######################################################################
# Interpolate the faulted stratigraphy
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ``faults=[fault]`` tells LoopStructural to interpolate the stratigraphy in
# the restored position, so the interpolated field is smooth and the fault
# makes the offset.

strat = model.create_and_add_foliation("strat", data=strat_data, faults=[fault], nelements=20_000)

######################################################################
# 3D view
# ~~~~~~~

view = Loop3DView(model)
view.plot_surface(fault, value=[0], colour="grey", opacity=0.6)
view.plot_surface(strat, value=horizons, paint_with=strat, cmap="viridis")
view.plot_data(strat, vector=False)
view.camera_position = "xz"
view.camera.azimuth = 30
view.camera.elevation = 20
view.display()

######################################################################
# Section at y = 0.5
# ~~~~~~~~~~~~~~~~~~
# The section shows the stratigraphic field, the four horizons and the fault
# trace (the zero value of the first fault frame coordinate). Every horizon
# has the same throw: :math:`0.15 \sin 60^\circ \approx 0.13`.


def plot_section(model, series, faults, values, y=0.5, ax=None, title=None):
    """Plot a vertical section of a model at constant y."""
    x, z = np.meshgrid(np.linspace(0, 1, 201), np.linspace(0, 1, 201))
    xyz = np.column_stack([x.ravel(), np.full(x.size, y), z.ravel()])
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    field = model.evaluate_feature_value(series, xyz).reshape(x.shape)
    ax.contourf(x, z, field, levels=20, cmap="viridis", alpha=0.6)
    ax.contour(x, z, field, levels=sorted(values), colors="k", linewidths=1.5, linestyles="solid")
    for name in faults:
        gx = model[name][0].evaluate_value(xyz).reshape(x.shape)
        ax.contour(x, z, gx, levels=[0], colors="r", linewidths=2)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.set_title(title or f"Section at y = {y}")
    return ax


plot_section(model, "strat", ["fault"], horizons)
plt.show()

######################################################################
# Check the offset
# ~~~~~~~~~~~~~~~~
# ``fault.apply_to_points`` restores points through the fault. The
# difference between a present-day point and its restored position is the
# fault movement at that point: 0 in the footwall, and 0.15 down-dip
# everywhere in the hanging wall.

points = np.array([[0.2, 0.5, 0.5], [0.8, 0.5, 0.2], [0.8, 0.5, 0.9]])
movement = points - fault.apply_to_points(points)
print("movement (x, y, z):")
print(np.round(movement, 4))
print("expected hanging-wall movement:", np.round(0.15 * fault_down_dip, 4))
