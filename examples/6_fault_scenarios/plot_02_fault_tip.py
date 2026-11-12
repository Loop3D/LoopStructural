"""
6b. A finite fault with a tip line inside the model
===================================================
Real faults end. The displacement is largest near the fault centre and goes
to zero on the **tip line**, so a horizon that the fault cuts has an offset
that dies out along strike and along dip.

This example uses the built-in ``BaseFault3D`` profile. The displacement is
the product of three 1D profiles of the fault frame coordinates:

* ``gx`` (normal to the fault): both walls move, in opposite directions, and
  the movement decays to zero at ``minor_axis / 2`` from the fault;
* ``gy`` (along the slip direction) and ``gz`` (along strike): a smooth peak
  that is 1 at the centre and 0 at ``intermediate_axis`` and ``major_axis``
  from the centre.

The fault is fully inside the model, so each horizon is one surface: a
horizon that the fault cuts has a slit, and a horizon above or below the tip
line is only folded (drag) or not changed.
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
# The same helpers are used in every fault scenario example: plane vectors
# from dip and dip direction, fault observations on a plane, synthetic
# stratigraphic observations, and a vertical section.


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
# Build the fault
# ~~~~~~~~~~~~~~~
# The fault dips 70 degrees toward +x. The axes set the size of the fault:
# the tip is 0.35 from the centre along strike (``major_axis``) and 0.3
# along the slip direction (``intermediate_axis``). The movement decays to
# zero at 0.25 from the fault (``minor_axis / 2``).
#
# .. note::
#    The fault frame coordinates are scaled differently: ``gz`` is the
#    distance divided by ``major_axis`` and ``gy`` the distance divided by
#    ``intermediate_axis`` (these two are radii), but ``gx`` is the distance
#    divided by ``minor_axis / 2`` (a full width).

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
    minor_axis=0.5,
    faultfunction="BaseFault3D",
    nelements=20_000,
    fault_buffer=0.5,
)
# The movement varies in space, so the fault is applied in 10 small steps
# that follow the curved fault frame.
fault.steps = 10

######################################################################
# The three profiles of ``BaseFault3D``:

fault.faultfunction.plot()

######################################################################
# Synthetic stratigraphy and interpolation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The layers were flat before faulting. The horizon at 0.8 is above the tip
# line (the upper tip is at about z = 0.78), so the fault does not cut it.

horizons = [0.25, 0.45, 0.65, 0.8]
strat_data = synthetic_observations(fault.apply_to_points, lambda p: p[:, 2], horizons, "strat")
strat = model.create_and_add_foliation("strat", data=strat_data, faults=[fault], nelements=25_000)

######################################################################
# 3D view
# ~~~~~~~
# The fault surface of the model (``gx = 0``) extends to the model
# boundary; the fault only moves points inside the tip line.

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
# On the section at y = 0.5 the offset is largest at mid-depth and goes to
# zero at the upper and lower tips. The hanging wall moves down and the
# footwall moves up, and both walls are dragged toward the fault. On the
# section at y = 0.9 (beyond the lateral tip at y = 0.85) the horizons are
# not offset.

fig, axes = plt.subplots(1, 2, figsize=(12, 6))
plot_section(model, {"strat": horizons}, ["fault"], y=0.5, ax=axes[0])
plot_section(model, {"strat": horizons}, ["fault"], y=0.9, ax=axes[1])
plt.show()

######################################################################
# Offset on the fault surface
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The offset at a point on the fault is the movement of the hanging wall
# relative to the footwall. We evaluate the movement just on each side of
# the fault with ``fault.apply_to_points``.
#
# The tip line is a **rectangle** in the fault plane (0.7 along strike by
# 0.6 along dip), because the profile is a product of a ``gy`` profile and a
# ``gz`` profile. Natural tip lines are closer to ellipses.

s, d = np.meshgrid(np.linspace(-0.45, 0.45, 61), np.linspace(-0.45, 0.45, 61))
on_fault = center + s.reshape(-1, 1) * strike + d.reshape(-1, 1) * down_dip
eps = 0.005 * normal
hanging_wall = on_fault + eps - fault.apply_to_points(on_fault + eps)
footwall = on_fault - eps - fault.apply_to_points(on_fault - eps)
offset = np.linalg.norm(hanging_wall - footwall, axis=1).reshape(s.shape)

fig, ax = plt.subplots(figsize=(6, 5))
c = ax.contourf(s, d, offset, levels=np.linspace(0, 0.2, 21), cmap="magma")
ax.add_patch(plt.Rectangle((-0.35, -0.3), 0.7, 0.6, fill=False, ec="w", ls="--", lw=1.5))
fig.colorbar(c, ax=ax, label="offset")
ax.set_xlabel("distance along strike")
ax.set_ylabel("distance down dip")
ax.set_aspect("equal")
ax.invert_yaxis()
ax.set_title("Offset on the fault plane (dashed: tip line)")
plt.show()
