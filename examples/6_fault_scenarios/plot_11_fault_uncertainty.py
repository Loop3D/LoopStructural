"""
6k. Fault uncertainty: realisations of an uncertain fault
=========================================================
The geometry and kinematics of a fault are rarely known exactly. One way to
show the effect of this uncertainty is to build many models, each with
fault parameters drawn from a probability distribution, and to compare
them.

This example draws four parameters of a planar normal fault:

* dip: normal distribution, mean 60, standard deviation 6 (clipped to 45 to
  75 degrees);
* displacement: normal, mean 0.12, standard deviation 0.03 (clipped to 0.05
  to 0.2);
* minor axis (width of the drag zone): log-normal around 0.4 (clipped to
  0.1 to 1.5);
* position of the fault centre along x: uniform between 0.4 and 0.6.

Each realisation uses its own random seed, so the family of models is
reproducible.

.. note::
   LoopStructural has no built-in sampler for fault parameters, so the loop
   is written here by hand. Each realisation interpolates the fault and the
   stratigraphy again, which is the main cost of a stochastic study.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel

######################################################################
# Helper functions
# ~~~~~~~~~~~~~~~~
# The same helpers as in the other fault scenario examples, with the random
# generator passed as an argument.


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


def synthetic_observations(rng, restore, field, horizons, feature_name, n=8000, half_width=0.015):
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
# Draw the fault parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~


def draw_parameters(seed):
    rng = np.random.default_rng(seed)
    return {
        "seed": seed,
        "dip": float(np.clip(rng.normal(60, 6), 45, 75)),
        "displacement": float(np.clip(rng.normal(0.12, 0.03), 0.05, 0.2)),
        "minor_axis": float(np.clip(rng.lognormal(np.log(0.4), 0.5), 0.1, 1.5)),
        "center_x": float(rng.uniform(0.4, 0.6)),
    }


realisations = [draw_parameters(1000 + k) for k in range(4)]
print(pd.DataFrame(realisations).round(3).to_string(index=False))

######################################################################
# Build one model for each realisation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The synthetic observations are made again for each realisation, so each
# model is a possible model of the same (flat, faulted) stratigraphy.

horizons = [0.25, 0.45, 0.65, 0.8]


def build_model(parameters):
    rng = np.random.default_rng(parameters["seed"])
    center = np.array([parameters["center_x"], 0.5, 0.5])
    normal, down_dip = dip_vectors(parameters["dip"], 90)
    model = GeologicalModel(np.zeros(3), np.ones(3))
    fault = model.create_and_add_fault(
        "fault",
        parameters["displacement"],
        data=planar_fault_data("fault", center, parameters["dip"], 90),
        fault_slip_vector=-down_dip,
        fault_normal_vector=normal,
        fault_center=center,
        major_axis=1.6,
        intermediate_axis=4.0,
        minor_axis=parameters["minor_axis"],
        faultfunction="BaseFault",
        nelements=15_000,
        fault_buffer=0.5,
    )
    fault.steps = 10
    data = synthetic_observations(rng, fault.apply_to_points, lambda p: p[:, 2], horizons, "strat")
    model.create_and_add_foliation("strat", data=data, faults=[fault], nelements=25_000)
    return model


models = [build_model(p) for p in realisations]

######################################################################
# Sections of each realisation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig, axes = plt.subplots(1, len(models), figsize=(5 * len(models), 5))
for ax, model, p in zip(axes, models, realisations):
    title = f"dip {p['dip']:.0f}, d {p['displacement']:.2f}, minor {p['minor_axis']:.2f}"
    plot_section(model, {"strat": horizons}, ["fault"], ax=ax, title=title)
plt.show()

######################################################################
# Where is the model uncertain?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A simple measure: the standard deviation of the stratigraphic value over
# the realisations, at each point of the section. It is largest between
# the fault traces of the realisations (a point is in the hanging wall in
# some realisations and in the footwall in others). It decreases away from
# the faults, and is not zero where a broad drag zone reaches.

x, z = np.meshgrid(np.linspace(0, 1, 201), np.linspace(0, 1, 201))
xyz = np.column_stack([x.ravel(), np.full(x.size, 0.5), z.ravel()])
values = np.array([model["strat"].evaluate_value(xyz) for model in models])
spread = values.std(axis=0).reshape(x.shape)

fig, ax = plt.subplots(figsize=(6, 6))
c = ax.pcolormesh(x, z, spread, cmap="Reds")
for model in models:
    gx = model["fault"][0].evaluate_value(xyz).reshape(x.shape)
    ax.contour(x, z, gx, levels=[0], colors="k", linewidths=1)
fig.colorbar(c, ax=ax, label="standard deviation of strat")
ax.set_aspect("equal")
ax.set_title("Spread of the realisations (black: fault traces)")
plt.show()
