"""Gradient orthogonal constraints, the iteration limit and a solve that does
not converge, with solver="native".

All tests skip while the native core (loop_interpolation_core) is not installed.
"""

import numpy as np
import pytest

# Skip the whole file when the native core is not installed.
pytest.importorskip("loop_interpolation_core")

from loop_interpolation import FiniteDifferenceInterpolator as FDI
from loop_interpolation import StructuredGrid, _native

FIELD_TOL = 1e-4  # fraction of the field range


def _make_interpolator(solver, b=0.0):
    """Make a 15^3 grid with value data and a gradient orthogonal constraint."""
    origin = np.array([-0.1, -0.1, -0.1])
    maximum = np.array([1.1, 1.1, 1.1])
    nsteps = np.array([15, 15, 15])
    grid = StructuredGrid(origin=origin, nsteps=nsteps, step_vector=(maximum - origin) / nsteps)
    interpolator = FDI(grid)
    interpolator.solver = solver
    rng = np.random.default_rng(0)
    xyz = rng.uniform(0.1, 0.9, size=(60, 3))
    values = xyz[:, 2] + 0.3 * xyz[:, 0]
    interpolator.set_value_constraints(np.column_stack([xyz, values, np.ones(60)]))
    interpolator.setup()
    # as a structural frame: the gradient is orthogonal to x at the barycentres
    points = grid.barycentre[::3]
    vectors = np.tile([1.0, 0.0, 0.0], (points.shape[0], 1))
    interpolator.add_gradient_orthogonal_constraints(
        points, vectors, w=1.0, b=b, name="frame_orthogonal"
    )
    return interpolator


def _range(c):
    return float(np.max(c) - np.min(c))


def test_a_native_setup_keeps_orthogonal_constraints_as_raw_data():
    interpolator = _make_interpolator("native")
    assert "frame_orthogonal" in interpolator._native_orthogonal
    assert "frame_orthogonal" not in interpolator.constraints
    # the core can solve it: no NotImplementedError
    _native.check_supported(interpolator, interpolator._setup_kwargs, {})


def test_orthogonal_constraints_give_the_same_field_as_cg():
    native = _make_interpolator("native")
    assert native.solve_system("native")
    matrix = _make_interpolator("cg")
    assert "frame_orthogonal" in matrix.constraints
    assert matrix.solve_system("cg", solver_kwargs={"rtol": 1e-12, "atol": 0.0})
    assert np.max(np.abs(native.c - matrix.c)) <= FIELD_TOL * _range(matrix.c)
    # the constraint has an effect: the gradient in x is small
    gradient = native.support.evaluate_gradient(native.support.barycentre[::3], native.c)
    assert np.median(np.abs(gradient[:, 0])) < 0.1 * np.median(np.abs(gradient[:, 2]))


def test_an_orthogonal_constraint_with_b_is_not_supported():
    interpolator = _make_interpolator("native", b=0.5)
    assert "frame_orthogonal" in interpolator.constraints
    with pytest.raises(NotImplementedError, match="add_constraints_to_least_squares"):
        _native.check_supported(interpolator, interpolator._setup_kwargs, {})


def test_the_default_iteration_limit_follows_the_grid():
    assert _native._core_settings(None, {}, 1000)["max_iter"] == 10 * 1000
    assert _native._core_settings(None, {"maxiter": 7}, 1000)["max_iter"] == 7
    assert "max_iter" not in _native._core_settings(None, {}, None)


def test_a_solve_that_does_not_converge_replaces_a_field_of_another_grid():
    interpolator = _make_interpolator("native")
    # the field of an old, smaller grid
    interpolator.c = np.zeros(10)
    assert not interpolator.solve_system("native", solver_kwargs={"maxiter": 1})
    assert interpolator.c.shape == (interpolator.support.n_nodes,)
    assert "field" in interpolator.latest_native_report


def test_a_solve_that_does_not_converge_keeps_a_field_of_this_grid():
    interpolator = _make_interpolator("native")
    assert interpolator.solve_system("native")
    good = interpolator.c.copy()
    assert not interpolator.solve_system("native", solver_kwargs={"maxiter": 1})
    assert np.array_equal(interpolator.c, good)
