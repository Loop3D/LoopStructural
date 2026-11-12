"""
6c. A one-sided fault with a tip line
=====================================
In the previous example both walls of the fault move, in opposite
directions. Here only the **hanging wall** moves and the footwall is fixed.
This is a common assumption when a footwall block is rigid basement, and it
gives a different shape: the full offset is a step on one side of the fault.

The profile is built from the profile functions of LoopStructural:

* ``gx``: a :code:`Composite` of the hanging-wall cubic of ``BaseFault``
  (1 at the fault, 0 at ``minor_axis / 2``) and :code:`Zeros` on the
  footwall;
* ``gy`` and ``gz``: the same smooth peak as ``BaseFault3D``, so the tip line
  is in the same position as in the previous example.

.. note::
   There is no named one-sided profile in LoopStructural. The profile
   classes and ``smooth_peak`` are imported from a private module, and the
   ``scale`` of :code:`FaultDisplacement` must be set to 1 (the default 0.5
   is for profiles where each wall moves half of the displacement).
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel
from LoopStructural.modelling.features.fault._fault_function import (
    BaseFault,
    Composite,
    FaultDisplacement,
    Zeros,
    smooth_peak,
)
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
# A hanging-wall-only profile
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

hanging_wall_only = FaultDisplacement(
    gx=Composite(BaseFault.hw, Zeros()), gy=smooth_peak, gz=smooth_peak, scale=1.0
)
hanging_wall_only.plot()

######################################################################
# Build the fault
# ~~~~~~~~~~~~~~~
# The plane and the tip line are the same as in the previous example. The
# minor axis is 1.0, so the hanging-wall movement decays over 0.5 from the
# fault.

center = np.array([0.5, 0.5, 0.5])
normal, down_dip = dip_vectors(70, 90)
strike = np.cross(normal, down_dip)

model = GeologicalModel(np.zeros(3), np.ones(3))
fault = model.create_and_add_fault(
    "fault",
    0.2,
    data=planar_fault_data("fault", center, 70, 90),
    fault_slip_vector=-down_dip,
    fault_normal_vector=normal,
    fault_center=center,
    major_axis=0.35,
    intermediate_axis=0.3,
    minor_axis=1.0,
    faultfunction=hanging_wall_only,
    nelements=20_000,
    fault_buffer=0.5,
)
fault.steps = 10

horizons = [0.25, 0.45, 0.65, 0.8]
strat_data = synthetic_observations(fault.apply_to_points, lambda p: p[:, 2], horizons, "strat")
strat = model.create_and_add_foliation("strat", data=strat_data, faults=[fault], nelements=25_000)

######################################################################
# 3D view
# ~~~~~~~

view = Loop3DView(model)
view.plot_surface(fault, value=[0], colour="grey", opacity=0.4)
view.plot_surface(strat, value=horizons, paint_with=strat, cmap="viridis")
view.camera_position = "xz"
view.camera.azimuth = 30
view.camera.elevation = 20
view.display()

######################################################################
# Section through the fault centre
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The footwall horizons are flat right up to the fault. All of the offset
# is in the hanging wall, which sags toward the fault.

plot_section(model, {"strat": horizons}, ["fault"], y=0.5)
plt.show()

######################################################################
# Movement on each side of the fault
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Along a horizontal line through the fault centre, the footwall movement
# is zero and the hanging-wall movement is a step at the fault that decays
# to zero at 0.5 from the fault. The step is about 0.17, not 0.2: a point
# moves in 10 steps, and the ``gy`` profile decreases as the point moves
# toward the upper tip line.

x = np.linspace(0.0, 1.0, 401)
line = np.column_stack([x, np.full_like(x, 0.5), np.full_like(x, 0.5)])
movement = np.linalg.norm(line - fault.apply_to_points(line), axis=1)
fig, ax = plt.subplots(figsize=(7, 3))
ax.plot(x, movement)
ax.axvline(center[0], color="r", lw=1)
ax.set_xlabel("x (y = 0.5, z = 0.5)")
ax.set_ylabel("movement")
ax.set_title("Hanging wall moves, footwall fixed")
plt.show()
