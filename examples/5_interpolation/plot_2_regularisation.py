"""
5b. Regularisation: smoothing versus fitting
============================================

An implicit surface has to do two things at once: pass through the data, and behave sensibly
where there are no data. The interpolators do both by minimising a weighted sum of terms,

* data terms, which penalise the distance from the contact points and the angle to the
  orientation measurements, and
* a regularisation term, which penalises the second derivative of the scalar field, so that
  the surface is as smooth as the data allow.

The ``regularisation`` weight sets how much the smoothness term counts. Too small and the
surface chases every error in the data; too large and it flattens real structure. This example
measures that trade-off on a synthetic dome whose true shape is known, with clean and with noisy
data. For noise-free data the best weight is the smallest one that still solves, while for noisy
data the best weight is larger.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from LoopStructural.interpolators import FiniteDifferenceInterpolator, StructuredGrid

# %%
# The synthetic model and its data
# --------------------------------
# A gaussian dome in the unit cube, ``s = z - dome(x, y)``, with horizons at ``s = 0.3, 0.5, 0.7``
# (the same model as the previous example).

RELIEF, CENTRE, WIDTH = 0.2, 0.5, 0.2
HORIZONS = (0.3, 0.5, 0.7)


def truth(xyz):
    r2 = (xyz[:, 0] - CENTRE) ** 2 + (xyz[:, 1] - CENTRE) ** 2
    return xyz[:, 2] - RELIEF * np.exp(-r2 / (2 * WIDTH**2))


def truth_gradient(xyz):
    r2 = (xyz[:, 0] - CENTRE) ** 2 + (xyz[:, 1] - CENTRE) ** 2
    e = RELIEF * np.exp(-r2 / (2 * WIDTH**2)) / WIDTH**2
    return np.c_[e * (xyz[:, 0] - CENTRE), e * (xyz[:, 1] - CENTRE), np.ones(len(xyz))]


def borehole_contacts(xy):
    points = []
    for x, y in xy:
        for level in HORIZONS:

            def misfit(z, x=x, y=y, level=level):
                return truth(np.array([[x, y, z]]))[0] - level

            if misfit(0.0) * misfit(1.0) < 0:
                points.append([x, y, brentq(misfit, 0.0, 1.0), level])
    return np.array(points)


rng = np.random.default_rng(4)
borehole_xy = np.vstack(
    [[0.5, 0.5], [0.35, 0.5], [0.65, 0.5], [0.5, 0.35], [0.5, 0.65], rng.uniform(0.1, 0.9, (4, 2))]
)
exact_values = borehole_contacts(borehole_xy)
orientation_xyz = np.array(
    [[0.4, 0.4, 0.5], [0.6, 0.6, 0.5], [0.3, 0.7, 0.6], [0.8, 0.3, 0.4], [0.5, 0.5, 0.5]]
)
exact_normals = np.c_[orientation_xyz, truth_gradient(orientation_xyz)]


def add_noise(noise, seed):
    """Perturb the contacts' height by ``noise`` and the orientation vectors by ``2.5 * noise``."""
    generator = np.random.default_rng(seed)
    values = exact_values.copy()
    values[:, 2] += generator.normal(0, noise, len(values))
    normals = exact_normals.copy()
    normals[:, 3:] += generator.normal(0, 2.5 * noise, (len(normals), 3))
    return values, normals


# %%
# Solving with a given regularisation weight
# ------------------------------------------
# ``regularisation`` is passed to ``setup_interpolator``. Its default is 1, which
# gives each second-derivative term a weight of 1 against data weights of 1 for contacts and 10
# for orientation vectors. We use a grid of 16 cells per side.


def solve(values, normals, regularisation):
    grid = StructuredGrid(
        origin=np.zeros(3), nsteps_cells=np.full(3, 16), step_vector=np.full(3, 1 / 16)
    )
    interpolator = FiniteDifferenceInterpolator(grid)
    interpolator.set_value_constraints(values)
    interpolator.set_normal_constraints(normals)
    interpolator.setup_interpolator(regularisation=regularisation)
    solved = interpolator.solve_system("cg")
    return interpolator, bool(solved)


axis = np.linspace(0.1, 0.9, 9)
points = np.array(np.meshgrid(axis, axis, axis, indexing="ij")).reshape(3, -1).T


def errors(interpolator, values):
    """RMS error against the truth over the interior, and RMS distance from the contact points.

    Noise can move a contact outside the model, where there is no value: those are skipped.
    """
    field_error = interpolator.evaluate_value(points) - truth(points)
    data_error = interpolator.evaluate_value(values[:, :3]) - values[:, 3]
    return np.sqrt(np.mean(field_error**2)), np.sqrt(np.nanmean(data_error**2))


# %%
# Sweeping the weight
# -------------------
# For each weight, solve with exact data and with noisy data (height noise 0.1, a tenth of the
# box, averaged over five random draws so that the curve does not depend on one lucky sample),
# and compare the result with the truth.

WEIGHTS = (0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)
NOISE = 0.1
DRAWS = 5

curves = {"exact data": [], "noisy data": []}
misfits = {"exact data": [], "noisy data": []}
not_solved = {"exact data": [], "noisy data": []}
for weight in WEIGHTS:
    interpolator, solved = solve(exact_values, exact_normals, weight)
    error, misfit = errors(interpolator, exact_values)
    curves["exact data"].append(error)
    misfits["exact data"].append(misfit)
    not_solved["exact data"].append(not solved)

    draw_errors, draw_misfits, draw_unsolved = [], [], 0
    for draw in range(DRAWS):
        values, normals = add_noise(NOISE, seed=100 + draw)
        interpolator, solved = solve(values, normals, weight)
        error, misfit = errors(interpolator, values)
        draw_errors.append(error)
        draw_misfits.append(misfit)
        draw_unsolved += not solved
    curves["noisy data"].append(np.mean(draw_errors))
    misfits["noisy data"].append(np.mean(draw_misfits))
    not_solved["noisy data"].append(draw_unsolved > 0)

print(f"{'weight':>8} | {'exact: error':>12} {'misfit':>8} | {'noisy: error':>12} {'misfit':>8}")
for i, weight in enumerate(WEIGHTS):
    flag = "  <- solver did not converge" if not_solved["exact data"][i] or not_solved["noisy data"][i] else ""
    print(
        f"{weight:>8} | {curves['exact data'][i]:12.4f} {misfits['exact data'][i]:8.4f} |"
        f" {curves['noisy data'][i]:12.4f} {misfits['noisy data'][i]:8.4f}{flag}"
    )

# %%
# The error against the truth, and the sections
# ---------------------------------------------
# Left: error against the true field for each weight. The default weight (1) is marked. Right:
# the middle section through the dome for one noisy data set at a weak, a moderate and a strong
# weight, against the true horizons (dashed).

fig = plt.figure(figsize=(13, 4.2))
ax = fig.add_subplot(1, 4, 1)
for label, colour in (("exact data", "tab:blue"), ("noisy data", "tab:red")):
    ax.plot(WEIGHTS, curves[label], "o-", color=colour, label=label)
ax.axvline(1, color="grey", linestyle=":")
ax.set_xscale("log")
ax.set_xlabel("regularisation weight")
ax.set_ylabel("rms error against the truth")
ax.legend()

noisy_values, noisy_normals = add_noise(NOISE, seed=100)
x = np.linspace(0, 1, 100)
xx, zz = np.meshgrid(x, x)
section = np.c_[xx.ravel(), np.full(xx.size, 0.5), zz.ravel()]
for i, weight in enumerate((0.03, 1, 30), start=2):
    interpolator, _ = solve(noisy_values, noisy_normals, weight)
    sub = fig.add_subplot(1, 4, i)
    sub.contour(xx, zz, truth(section).reshape(xx.shape), levels=HORIZONS, colors="k", linestyles="--")
    sub.contour(
        xx, zz, interpolator.evaluate_value(section).reshape(xx.shape), levels=HORIZONS, colors="tab:red"
    )
    near = np.abs(noisy_values[:, 1] - 0.5) < 0.2
    sub.plot(noisy_values[near, 0], noisy_values[near, 2], "ko", ms=3)
    sub.set_title(f"weight {weight}")
    sub.set_aspect("equal")
fig.tight_layout()
plt.show()

# %%
# Reading the results
# -------------------
# * With exact data the error rises with the weight, because any smoothing moves the surface
#   away from the truth, which the data describe. The weakest weights fit the data almost
#   exactly (misfit near zero).
# * With noisy data the picture changes. A very small weight fits the noise: the error is
#   several times larger than at the default, and the horizons in the left section wiggle. A
#   very large weight flattens the dome: the error levels off near 0.12, which is roughly the
#   error of a flat horizon. The lowest error is in between, near the default weight of 1.
# * The data misfit always grows with the weight. It tells you how far the surface has been
#   pulled from the data, but a small misfit is not the same as a small error: the noisy data
#   are matched best by the weights that are furthest from the truth.
#
# The right weight depends on how much you trust the data, which is a geological judgement.
# A practical approach, as in the debugging guide, is to lower the weight until the surface
# starts to follow noise, then raise it until the structure starts to flatten, and pick a value
# between the two. The data weights (``cpw`` for contacts, ``npw`` for orientation vectors,
# ``gpw`` and ``tpw`` for the other constraint types) are the same trade-off from the other side.
