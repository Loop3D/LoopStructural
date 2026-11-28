"""The cg stop on an estimate of the position error (position_tolerance)."""

import logging

import numpy as np
import pytest
from loop_interpolation import FiniteDifferenceInterpolator, StructuredGrid, _solver_strategy
from loop_interpolation._solver_strategy import (
    CG_POSITION_CHECK,
    CG_POSITION_DELAY,
    cg_to_position_tolerance,
)
from scipy import sparse
from scipy.sparse.linalg import cg


def _spd(n=400, seed=0, ridge=1e-6):
    """A 1D bi-Laplacian plus a small ridge: cg stagnates, the solution is known.

    The change of the solution over 50 iterations is here up to 500 times
    smaller than the error, so a stop on the change alone fails on it.
    """
    lap = sparse.diags([-1.0, 2.0, -1.0], [-1, 0, 1], shape=(n, n), format="csr")
    N = (lap.T @ lap + ridge * sparse.eye(n)).tocsr()
    x = np.sin(np.linspace(0.0, 3.0, n)) + np.random.default_rng(seed).normal(0.0, 0.01, n)
    return N, N @ x, x


def _never(x, dx):
    return float("inf")


def test_iterates_are_those_of_scipy():
    N, rhs, _ = _spd()
    for maxiter in (1, 7, 60, 135):
        ours, info, report = cg_to_position_tolerance(N, rhs, 1.0, _never, maxiter=maxiter)
        theirs, _ = cg(N, rhs, rtol=0.0, atol=0.0, maxiter=maxiter)
        assert info == maxiter and report["stop"] == "maxiter"
        assert report["iterations"] == maxiter
        np.testing.assert_array_equal(ours, theirs)


def test_zero_right_hand_side():
    N, rhs, _ = _spd(50)
    x, info, report = cg_to_position_tolerance(N, np.zeros_like(rhs), 1.0, _never)
    assert info == 0 and report["stop"] == "zero_rhs"
    assert np.all(x == 0.0)


def test_stops_on_the_position_estimate():
    N, rhs, exact = _spd()
    gradient = 2.0
    calls = []

    def position_error(x, dx):
        calls.append(len(calls))
        return float(np.max(np.abs(dx)) / gradient)

    tolerance = 1e-4
    x, info, report = cg_to_position_tolerance(N, rhs, tolerance, position_error)
    assert info == 0 and report["stop"] == "position"
    assert report["position_error_estimate"] <= tolerance
    # the estimate is the change scaled up, never less than the change
    assert report["position_error_change"] <= report["position_error_estimate"]
    assert report["smallest_ritz_value"] == pytest.approx(1e-6, rel=0.05)
    assert report["iterations"] % CG_POSITION_CHECK == 0
    assert report["iterations"] >= CG_POSITION_DELAY
    # the true position error is within the tolerance
    assert np.max(np.abs(x - exact)) / gradient <= tolerance
    # a tighter tolerance takes more iterations and is more accurate
    x2, _, report2 = cg_to_position_tolerance(N, rhs, 1e-6, position_error)
    assert report2["iterations"] > report["iterations"]
    assert np.max(np.abs(x2 - exact)) < np.max(np.abs(x - exact))
    assert np.max(np.abs(x2 - exact)) / gradient <= 1e-6


@pytest.mark.parametrize("ridge", [1e-6, 1e-9])
@pytest.mark.parametrize("tolerance", [1e-2, 1e-4])
def test_the_stop_holds_on_a_stagnating_system(ridge, tolerance):
    """Where the change of the solution alone is 10 to 500 times too small."""
    N, rhs, exact = _spd(ridge=ridge)
    fn = lambda x, dx: float(np.max(np.abs(dx)))  # noqa: E731
    x, info, report = cg_to_position_tolerance(N, rhs, tolerance, fn, maxiter=40 * len(rhs))
    assert info == 0 and report["stop"] == "position"
    assert np.max(np.abs(x - exact)) <= tolerance


def test_the_estimate_is_not_tested_before_the_residual_guard():
    N, rhs, _ = _spd()
    seen = []

    def position_error(x, dx):
        seen.append(np.linalg.norm(rhs - N @ x) / np.linalg.norm(rhs))
        return 1e-30  # would stop at once

    _, info, report = cg_to_position_tolerance(N, rhs, 1e30, position_error)
    assert info == 0 and report["stop"] == "position"
    assert len(seen) == 1
    assert seen[0] < _solver_strategy.CG_POSITION_GUARD_RTOL * (1 + 1e-6)


def test_stops_at_the_residual_floor_without_an_estimate():
    N, rhs, exact = _spd(60)
    x, info, report = cg_to_position_tolerance(N, rhs, 1e-30, _never)
    assert info == 0 and report["stop"] == "residual"
    assert report["relative_residual"] < 1e-12
    np.testing.assert_allclose(x, exact, atol=1e-6)


def test_repeatable():
    N, rhs, _ = _spd()
    fn = lambda x, dx: float(np.max(np.abs(dx)))  # noqa: E731
    a, _, _ = cg_to_position_tolerance(N, rhs, 1e-5, fn)
    b, _, _ = cg_to_position_tolerance(N, rhs, 1e-5, fn)
    np.testing.assert_array_equal(a, b)


def _interpolator(scale=1.0, matrix_free=False):
    """A dipping plane with noise in a box of side 20 * scale."""
    grid = StructuredGrid(
        origin=np.zeros(3), nsteps=np.array([21, 21, 21]), step_vector=np.full(3, scale)
    )
    interp = FiniteDifferenceInterpolator(grid)
    rng = np.random.default_rng(3)
    xyz = rng.uniform(2.0, 18.0, size=(60, 3)) * scale
    values = xyz[:, 2] + 0.3 * xyz[:, 0] + scale * np.sin(xyz[:, 1] / scale)
    interp.set_value_constraints(np.column_stack([xyz, values, np.ones(len(xyz))]))
    if matrix_free:
        interp.apply_scaling_matrix = False
    interp.setup_interpolator(regularisation_matrix_free=matrix_free)
    return interp


def _position_error(interp, reference):
    nodes = interp.support.nodes
    gradient = np.linalg.norm(interp.support.evaluate_gradient(nodes, reference), axis=1)
    return np.max(np.abs(interp.c - reference) / gradient)


@pytest.mark.parametrize("scale", [1.0, 500.0])
@pytest.mark.parametrize("matrix_free", [False, True])
def test_interpolator_position_tolerance(scale, matrix_free):
    """The position error is within the tolerance, in the units of the support."""
    reference = _interpolator(scale, matrix_free)
    assert reference.solve_system(solver="cg", solver_kwargs={"rtol": 1e-13, "atol": 0.0})
    extent = 20.0 * scale
    errors = []
    for fraction in (1e-3, 1e-5):
        tolerance = fraction * extent
        interp = _interpolator(scale, matrix_free)
        assert interp.solve_system(solver="cg", position_tolerance=tolerance)
        timing = interp.get_last_solve_timing()
        assert timing["cg_stop"] == "position"
        assert timing["position_tolerance"] == tolerance
        assert timing["position_error_estimate"] <= tolerance
        error = _position_error(interp, reference.c)
        assert error <= tolerance
        errors.append(error / extent)
    assert errors[1] < errors[0]


def test_position_tolerance_does_not_depend_on_the_length_unit():
    """The same model in metres and in kilometres stops at the same iteration."""
    iterations = []
    for scale in (1.0, 1000.0):
        interp = _interpolator(scale)
        interp.solve_system(solver="cg", position_tolerance=1e-4 * 20.0 * scale)
        iterations.append(interp.get_last_solve_timing()["cg_iterations"])
    assert abs(iterations[0] - iterations[1]) <= CG_POSITION_CHECK


def test_without_a_position_tolerance_nothing_changes():
    a = _interpolator()
    a.solve_system(solver="cg", tol=1e-9)
    b = _interpolator()
    b.solve_system(solver="cg", tol=1e-9, position_tolerance=None)
    np.testing.assert_array_equal(a.c, b.c)
    assert "cg_stop" not in a.get_last_solve_timing()
    assert "cg_stop" not in b.get_last_solve_timing()


def test_tol_and_rtol_do_not_stop_a_position_solve(caplog):
    reference = _interpolator()
    reference.solve_system(solver="cg", position_tolerance=1e-3)
    interp = _interpolator()
    with caplog.at_level(logging.WARNING):
        # this tol and rtol alone would stop cg after a few iterations
        interp.solve_system(
            solver="cg", tol=1e3, solver_kwargs={"rtol": 0.5}, position_tolerance=1e-3
        )
    np.testing.assert_array_equal(interp.c, reference.c)
    assert "ignored" in caplog.text


def test_position_tolerance_must_be_positive():
    interp = _interpolator()
    with pytest.raises(ValueError, match="position_tolerance"):
        interp.solve_system(solver="cg", position_tolerance=0.0)


def test_other_solvers_ignore_the_position_tolerance(caplog):
    a = _interpolator()
    a.solve_system(solver="lsmr", tol=1e-10)
    b = _interpolator()
    with caplog.at_level(logging.WARNING):
        b.solve_system(solver="lsmr", tol=1e-10, position_tolerance=1e-3)
    np.testing.assert_array_equal(a.c, b.c)
    assert "position_tolerance is used only by solver='cg'" in caplog.text
