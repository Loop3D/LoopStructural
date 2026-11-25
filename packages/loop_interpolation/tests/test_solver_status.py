"""The solve status must be true, a failed solve must not change solver, and NaN in gives NaN out."""

import numpy as np
import pytest
from loop_interpolation import FiniteDifferenceInterpolator, StructuredGrid


def _interpolator():
    grid = StructuredGrid(origin=np.zeros(3), nsteps=np.array([10, 10, 10]), step_vector=np.ones(3))
    interpolator = FiniteDifferenceInterpolator(grid)
    interpolator.set_value_constraints(
        np.array([[2.0, 2.0, 2.0, 0.0, 1.0], [7.0, 7.0, 7.0, 1.0, 1.0]])
    )
    interpolator.setup_interpolator()
    return interpolator


def test_lsmr_at_iteration_limit_is_not_a_success():
    interpolator = _interpolator()
    assert interpolator.solve_system("lsmr", solver_kwargs={"maxiter": 2}) is False
    assert not interpolator.up_to_date


def test_admm_without_inequalities_falls_back_to_lsmr():
    lsmr = _interpolator()
    assert lsmr.solve_system("lsmr") is True
    for kwargs in ({}, {"admm_abs_tol": 1e-5}):
        admm = _interpolator()
        assert admm.solve_system("admm", solver_kwargs=kwargs) is True
        assert np.allclose(admm.c, lsmr.c)


def test_update_reuses_the_solver_and_kwargs_of_the_last_solve(monkeypatch):
    interpolator = _interpolator()
    interpolator.solve_system("lsmr", tol=1e-9, solver_kwargs={"maxiter": 5000})
    calls = []
    original = interpolator.solve_system

    def spy(solver=None, **kwargs):
        calls.append((solver, kwargs))
        return original(solver, **kwargs)

    monkeypatch.setattr(interpolator, "solve_system", spy)
    interpolator.up_to_date = False
    interpolator.evaluate_value(np.array([[5.0, 5.0, 5.0]]))
    assert calls == [("lsmr", {"tol": 1e-9, "solver_kwargs": {"maxiter": 5000}})]


def test_failed_solve_is_not_repeated_at_each_evaluation(monkeypatch):
    interpolator = _interpolator()
    assert interpolator.solve_system("lsmr", solver_kwargs={"maxiter": 2}) is False
    calls = []
    original = interpolator.solve_system

    def spy(*args, **kwargs):
        calls.append(args)
        return original(*args, **kwargs)

    monkeypatch.setattr(interpolator, "solve_system", spy)
    points = np.array([[5.0, 5.0, 5.0]])
    for _ in range(3):
        interpolator.evaluate_value(points)
    assert len(calls) == 1
    # a change to the constraints lets update() solve again
    interpolator.set_value_constraints(np.array([[2.0, 2.0, 2.0, 0.0, 1.0]]))
    interpolator.evaluate_value(points)
    assert len(calls) == 2


def test_nan_point_evaluates_to_nan():
    interpolator = _interpolator()
    interpolator.solve_system("lsmr")
    values = interpolator.evaluate_value(np.array([[np.nan, 5.0, 5.0], [5.0, 5.0, 5.0]]))
    assert np.isnan(values[0])
    assert np.isfinite(values[1])
    gradient = interpolator.evaluate_gradient(np.array([[np.nan, 5.0, 5.0], [5.0, 5.0, 5.0]]))
    assert np.all(np.isnan(gradient[0]))
    assert np.all(np.isfinite(gradient[1]))


def test_fdi_does_not_build_kdtree_when_weight_scale_is_off(monkeypatch):
    from loop_interpolation import _finite_difference_interpolator as fdi

    def fail(*args, **kwargs):
        raise AssertionError("KD-tree built with use_regularisation_weight_scale=False")

    monkeypatch.setattr(fdi, "KDTree", fail)
    interpolator = FiniteDifferenceInterpolator(
        StructuredGrid(origin=np.zeros(3), nsteps=np.array([5, 5, 5]), step_vector=np.ones(3))
    )
    interpolator.set_normal_constraints(np.array([[2.0, 2.0, 2.0, 0.0, 0.0, 1.0, 1.0]]))
    interpolator.setup_interpolator()
    assert interpolator.n_n == 1


if __name__ == "__main__":
    pytest.main([__file__])
