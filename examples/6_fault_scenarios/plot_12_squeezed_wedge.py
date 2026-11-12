"""
6l. A thin wedge between two faults that meet at a small angle
==============================================================
Two synthetic normal faults (both dipping toward +x) meet at 20 degrees.
The steeper fault abuts the shallower fault, so the block between them is a
**wedge** that thins from 0.27 at the top of the model to zero at the
junction line near z = 0.38.

This is a hard case for any modelling method: near the junction the wedge
is thinner than the interpolation support, the angle between the faults is
small, and the horizons in the wedge are cut by two faults within a short
distance.
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
# The same helpers as in the Y junction example.


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
# Build the two faults
# ~~~~~~~~~~~~~~~~~~~~
# Fault ``a`` dips 55 degrees and has a displacement of 0.10. Fault ``b``
# dips 75 degrees, has a displacement of 0.05, and abuts fault ``a`` on its
# hanging-wall side.

model = GeologicalModel(np.zeros(3), np.ones(3))

a_center = np.array([0.45, 0.5, 0.5])
a_normal, a_down_dip = dip_vectors(55, 90)
fault_a = model.create_and_add_fault(
    "fault_a",
    0.10,
    data=planar_fault_data("fault_a", a_center, 55, 90),
    fault_slip_vector=-a_down_dip,
    fault_normal_vector=a_normal,
    fault_center=a_center,
    major_axis=4.0,
    intermediate_axis=4.0,
    minor_axis=4.0,
    faultfunction=constant_profile,
    nelements=20_000,
    fault_buffer=0.5,
)
fault_a.steps = 1

b_center = np.array([0.5, 0.5, 0.5])
b_normal, b_down_dip = dip_vectors(75, 90)
fault_b = model.create_and_add_fault(
    "fault_b",
    0.05,
    data=planar_fault_data("fault_b", b_center, 75, 90),
    fault_slip_vector=-b_down_dip,
    fault_normal_vector=b_normal,
    fault_center=b_center,
    major_axis=4.0,
    intermediate_axis=4.0,
    minor_axis=4.0,
    faultfunction=constant_profile,
    nelements=20_000,
    fault_buffer=0.5,
)
fault_b.steps = 1
fault_b.add_abutting_fault(fault_a, positive=True)

######################################################################
# Synthetic stratigraphy
# ~~~~~~~~~~~~~~~~~~~~~~

horizons = [0.2, 0.45, 0.6, 0.75, 0.9]
restore = restore_through([fault_a, fault_b])
strat_data = synthetic_observations(restore, lambda p: p[:, 2], horizons, "strat")
strat = model.create_and_add_foliation(
    "strat", data=strat_data, faults=[fault_a, fault_b], nelements=40_000
)

######################################################################
# 3D view
# ~~~~~~~

view = Loop3DView(model)
view.plot_surface(fault_a, value=[0], colour="grey", opacity=0.6)
view.plot_surface(fault_b, value=[0], colour="darkred", opacity=0.6)
view.plot_surface(strat, value=horizons, paint_with=strat, cmap="viridis")
view.camera_position = "xz"
view.camera.azimuth = 30
view.camera.elevation = 20
view.display()

######################################################################
# Section
# ~~~~~~~
# The wedge is between the two red lines above the junction. The horizons
# in the wedge are short pieces, and near the junction they are only a few
# interpolation cells long.

ax = plot_section(model, {"strat": horizons}, ["fault_a", "fault_b"])
plt.show()

######################################################################
# Wedge thickness
# ~~~~~~~~~~~~~~~
# The width of the wedge along x, at several heights. Compare it with the
# size of an interpolation cell of the stratigraphy.

x = np.linspace(0, 1, 4001)
for z0 in [1.0, 0.8, 0.6, 0.45, 0.4]:
    line = np.column_stack([x, np.full_like(x, 0.5), np.full_like(x, z0)])
    in_wedge = (fault_a[0].evaluate_value(line) > 0) & (
        np.nan_to_num(fault_b[0].evaluate_value(line), nan=1.0) < 0
    )
    print(f"wedge width at z = {z0}: {in_wedge.sum() * (x[1] - x[0]):.3f}")
print("interpolation cell size of strat:", np.round(strat.interpolator.support.step_vector, 3))
