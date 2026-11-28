"""Warm start and ridge factor with solver="native".

solve_system gives ``x0``/``warm_start`` and ``ridge_factor`` to the core.
A warm start must give the same field as a cold start. A large ridge
factor must pull the field to zero.
All tests skip while the native core (loop_interpolation_core) is not installed.
"""

import numpy as np
import pytest

# Skip the whole file when the native core is not installed.
pytest.importorskip("loop_interpolation_core")

from loop_interpolation import FiniteDifferenceInterpolator as FDI
from loop_interpolation import StructuredGrid

FIELD_TOL = 1e-6  # fraction of the field range


def _make_interpolator(bounds=False):
    """Make a 20^3 grid with value data, set up for the native solver."""
    origin = np.array([-0.1, -0.1, -0.1])
    maximum = np.array([1.1, 1.1, 1.1])
    nsteps = np.array([20, 20, 20])
    grid = StructuredGrid(origin=origin, nsteps=nsteps, step_vector=(maximum - origin) / nsteps)
    interpolator = FDI(grid)
    interpolator.solver = "native"
    rng = np.random.default_rng(0)
    xyz = rng.uniform(0.1, 0.9, size=(60, 3))
    values = xyz[:, 2] + 0.1 * np.sin(2.0 * xyz[:, 0]) * xyz[:, 1]
    interpolator.set_value_constraints(np.column_stack([xyz, values, np.ones(60)]))
    if bounds:
        bound_xyz = np.array([[0.5, 0.5, 0.7], [0.4, 0.6, 0.8]])
        interpolator.set_value_inequality_constraints(
            np.column_stack([bound_xyz, np.full(2, -5.0), np.full(2, 0.5), np.ones(2)])
        )
    interpolator.setup()
    return interpolator


def _range(c):
    return float(np.max(c) - np.min(c))


def test_native_warm_start_needs_fewer_iterations():
    interpolator = _make_interpolator()
    assert interpolator.solve_system("native")
    cold = interpolator.c.copy()
    cold_iterations = interpolator.latest_native_report["iterations"]

    assert interpolator.solve_system("native", solver_kwargs={"warm_start": True})
    warm_iterations = interpolator.latest_native_report["iterations"]
    assert warm_iterations * 10 < cold_iterations, (warm_iterations, cold_iterations)
    assert np.max(np.abs(interpolator.c - cold)) <= FIELD_TOL * _range(cold)


def test_native_x0_gives_the_same_field():
    interpolator = _make_interpolator()
    assert interpolator.solve_system("native")
    cold = interpolator.c.copy()

    rough = np.sin(np.arange(interpolator.support.n_nodes, dtype=float))
    assert interpolator.solve_system("native", solver_kwargs={"x0": rough})
    assert np.max(np.abs(interpolator.c - cold)) <= FIELD_TOL * _range(cold)


def test_native_warm_start_with_bounds():
    interpolator = _make_interpolator(bounds=True)
    kwargs = {"maxiter": 50000, "admm_abs_tol": 1e-7, "admm_rel_tol": 1e-7}
    assert interpolator.solve_system("native", solver_kwargs=dict(kwargs))
    cold = interpolator.c.copy()

    assert interpolator.solve_system("native", solver_kwargs={**kwargs, "warm_start": True})
    assert interpolator.latest_native_report["converged"]
    # ADMM stops at a tolerance, so the two fields agree less closely.
    assert np.max(np.abs(interpolator.c - cold)) <= 1e-3 * _range(cold)


def test_native_ridge_factor_pulls_the_field_to_zero():
    interpolator = _make_interpolator()
    assert interpolator.solve_system("native", solver_kwargs={"ridge_factor": 0.0})
    no_ridge = np.linalg.norm(interpolator.c)
    assert interpolator.solve_system("native", solver_kwargs={"ridge_factor": 1.0})
    assert np.linalg.norm(interpolator.c) < 0.9 * no_ridge


def test_native_ridge_off_is_the_same_as_ridge_factor_zero():
    interpolator = _make_interpolator()
    assert interpolator.solve_system("native", solver_kwargs={"ridge_factor": 0.0})
    zero = interpolator.c.copy()
    assert interpolator.solve_system(
        "native", solver_kwargs={"add_ridge_regularisation": False, "ridge_factor": 1.0}
    )
    assert np.allclose(interpolator.c, zero)


def test_native_matches_cg_with_the_same_ridge_factor():
    """The core and the matrix cg solver use the same ridge term."""
    native = _make_interpolator()
    assert native.solve_system("native", solver_kwargs={"ridge_factor": 0.1})

    origin = np.array([-0.1, -0.1, -0.1])
    maximum = np.array([1.1, 1.1, 1.1])
    nsteps = np.array([20, 20, 20])
    grid = StructuredGrid(origin=origin, nsteps=nsteps, step_vector=(maximum - origin) / nsteps)
    matrix = FDI(grid)
    matrix.apply_scaling_matrix = False
    rng = np.random.default_rng(0)
    xyz = rng.uniform(0.1, 0.9, size=(60, 3))
    values = xyz[:, 2] + 0.1 * np.sin(2.0 * xyz[:, 0]) * xyz[:, 1]
    matrix.set_value_constraints(np.column_stack([xyz, values, np.ones(60)]))
    matrix.setup()
    assert matrix.solve_system("cg", solver_kwargs={"ridge_factor": 0.1, "rtol": 1e-10})

    assert np.max(np.abs(native.c - matrix.c)) <= 1e-5 * _range(matrix.c)
