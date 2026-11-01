"""
5f. The native solver
=====================

The native solver is a Rust implementation of the same conjugate gradient solve that the Python
interpolator does. It solves the same discrete system, uses several threads, and, when it is chosen
*before* the interpolator is set up, the Python side never builds the sparse matrix, which is where
most of the memory goes.

This example checks the claims on a problem big enough to matter (about 55 000 unknowns) and looks
at the one thing that catches people out: the native solver stops at a default iteration limit, and
what it returns when it stops early is not the answer. The native solver needs the compiled
package ``loop_interpolation_core``; install it with ``pip install loop-interpolation[native]``.
It works on uniform 3D grids, which is what ``FDI`` uses (see the user guide for the full list of
limits).
"""

import importlib.util
import time
import tracemalloc

import numpy as np
import pandas as pd
from loop_common.geometry import BoundingBox
from scipy.optimize import brentq

from LoopStructural import InterpolatorBuilder

if importlib.util.find_spec("loop_interpolation_core") is None:
    raise ImportError(
        "This example needs the compiled native core: pip install loop-interpolation[native]"
    )

# %%
# The problem
# -----------
# The gaussian dome used in the earlier interpolation examples, in the unit cube, sampled with nine
# boreholes and five orientation measurements.

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
values = borehole_contacts(borehole_xy)
orientation_xyz = np.array(
    [[0.4, 0.4, 0.5], [0.6, 0.6, 0.5], [0.3, 0.7, 0.6], [0.8, 0.3, 0.4], [0.5, 0.5, 0.5]]
)
normals = np.c_[orientation_xyz, truth_gradient(orientation_xyz)]

axis = np.linspace(0.1, 0.9, 9)
points = np.array(np.meshgrid(axis, axis, axis, indexing="ij")).reshape(3, -1).T
bounding_box = BoundingBox(origin=[0.0, 0.0, 0.0], maximum=[1.0, 1.0, 1.0])

# %%
# A helper that builds, solves and measures
# -----------------------------------------
# ``InterpolatorBuilder`` chooses the solver with ``use_solver``. Whether that happens before or
# after ``setup_interpolator`` matters for the native solver, so the helper takes it as a
# parameter. ``tracemalloc`` measures the *Python* memory used while setting up and solving; the
# native core allocates its own memory outside of Python, which this does not see.


def run(solver, *, choose_solver="before", **solver_kwargs):
    builder = InterpolatorBuilder(
        "FDI", bounding_box=bounding_box, nelements=50000, buffer=0.1
    )
    builder.add_value_constraints(values).add_normal_constraints(normals)
    if choose_solver == "before":
        builder.use_solver(solver, **solver_kwargs)

    tracemalloc.start()
    start = time.perf_counter()
    builder.setup_interpolator()
    if choose_solver == "after":
        builder.use_solver(solver, **solver_kwargs)
    builder.solve()
    elapsed = time.perf_counter() - start
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()

    interpolator = builder.build()
    timing = getattr(interpolator, "latest_solve_timing", None) or {}
    return {
        "builder": builder,
        "interpolator": interpolator,
        "field": interpolator.evaluate_value(points),
        "seconds": elapsed,
        "python peak MB": peak / 1e6,
        "converged": bool(interpolator.up_to_date),
        "iterations": timing.get("iterations"),
        "unknowns": interpolator.dof,
    }


# %%
# The Python solver as the reference
# ----------------------------------

python = run("cg")
print(f"{python['unknowns']} unknowns; Python CG: {python['seconds']:.1f} s")

# %%
# The native solver with its defaults
# -----------------------------------
# Same problem, native solver. It is faster, but look at ``converged`` and at the error.

native_default = run("native")

error = lambda result: np.sqrt(np.mean((result["field"] - truth(points)) ** 2))
print(f"converged: {native_default['converged']} after {native_default['iterations']} iterations")
print(f"rms error against the truth: Python {error(python):.4f}, native {error(native_default):.4f}")
print(f"largest difference from the Python field: {np.max(np.abs(native_default['field'] - python['field'])):.2f}")

# %%
# The native solver stopped at its default limit of 5000 iterations without converging, logged a
# warning, and left the interpolator marked as *not* up to date. The field it returned is an
# unfinished solution and is far from the Python one, even though the residual in the warning
# looks close to the tolerance (about 3e-8 against 1e-8): on this problem a small residual does not
# mean a small error in the field. Nothing raised, so the way to catch this is to read the warning
# or check ``interpolator.up_to_date``, not to judge by the residual.
#
# The limit is a keyword, ``maxiter``, passed through ``use_solver``. This problem needs a little
# more than 5000 iterations, so raising the limit is enough:

native = run("native", maxiter=50000)
print(f"converged: {native['converged']} after {native['iterations']} iterations")
print(f"largest difference from the Python field: {np.max(np.abs(native['field'] - python['field'])):.1e}")
print(f"time: Python {python['seconds']:.1f} s, native {native['seconds']:.1f} s")

# %%
# Now the native field agrees with the Python field to about 1e-8, which is the solver tolerance,
# and it is faster. How much faster depends on the machine and the problem size. In a separate run
# on a problem twice this size the gap was larger (about 3 times, against about 1.7 times here).

# %%
# Memory: choose the solver before the setup
# ------------------------------------------
# The Python memory used while setting up and solving, for the Python solver, for the native
# solver chosen before ``setup_interpolator``, and for the native solver chosen after it.

after = run("native", choose_solver="after", maxiter=50000)
memory = pd.DataFrame(
    {
        "Python CG": python,
        "native, chosen before setup": native,
        "native, chosen after setup": after,
    }
).T[["python peak MB", "seconds"]]
print(memory.astype(float).round(1).to_string())
print(
    "same field whichever way: largest difference between the two native runs "
    f"{np.max(np.abs(after['field'] - native['field'])):.1e}"
)

# %%
# Chosen before setup, the Python side never builds the sparse matrix and its peak memory is a
# small fraction of the Python solver's. Chosen after setup, the matrix has already been built, so
# the peak is many times higher than that (the solver logs a warning to say so), although still
# below the Python solver's here. The answer is the same either way.

# %%
# Starting from the previous solution
# -----------------------------------
# After a small change, for example a new regularisation weight, the previous field is a good
# starting point. ``warm_start=True`` starts the solver from it instead of from zero.


def resolve(warm_start):
    builder = InterpolatorBuilder("FDI", bounding_box=bounding_box, nelements=50000, buffer=0.1)
    builder.use_solver("native", maxiter=50000)
    builder.add_value_constraints(values).add_normal_constraints(normals)
    builder.setup_interpolator().solve()
    first = builder.interpolator.latest_solve_timing["iterations"]

    builder.setup_interpolator(regularisation=1.1)  # a slightly stronger smoothing weight
    start = time.perf_counter()
    builder.solve(warm_start=warm_start)
    seconds = time.perf_counter() - start
    return first, builder.interpolator.latest_solve_timing["iterations"], seconds


first, cold_iterations, cold_seconds = resolve(False)
_, warm_iterations, warm_seconds = resolve(True)
print(f"first solve: {first} iterations")
print(f"new weight, starting from zero:     {cold_iterations} iterations, {cold_seconds:.2f} s")
print(f"new weight, warm start:             {warm_iterations} iterations, {warm_seconds:.2f} s")

# %%
# What to take from this
# ----------------------
# * The native solver gives the same answer as the Python solver, to the solver tolerance, and
#   is faster. The saving in Python memory only appears if the solver is chosen **before**
#   ``setup_interpolator``.
# * **Check that it converged.** The default limit of 5000 iterations is not enough for a problem
#   of this size, and what it returns when it stops is not the answer. Look for the warning, check
#   ``interpolator.up_to_date``, and raise ``maxiter`` if needed.
# * ``warm_start=True`` saves iterations when you solve again after a small change.
#
# The timings and memory figures are from one machine and one problem; they show the shape of the
# behaviour, not a benchmark. For the full list of options and unsupported settings see the native
# solver page of the user guide.
