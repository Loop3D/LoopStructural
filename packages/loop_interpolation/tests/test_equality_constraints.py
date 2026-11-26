"""Equality constraints must reach the solve (matrix path and native path)."""

import numpy as np
from loop_interpolation import FiniteDifferenceInterpolator, StructuredGrid, _native


def _make_interpolator(solver="cg"):
    grid = StructuredGrid(origin=np.zeros(3), nsteps=np.array([6, 6, 6]), step_vector=np.ones(3))
    interpolator = FiniteDifferenceInterpolator(grid)
    interpolator.solver = solver
    interpolator.set_value_constraints(np.array([[2.5, 2.5, 2.5, 1.0, 1.0]]))
    interpolator.setup()
    return interpolator


def test_equality_constraint_fixes_the_node_value():
    interpolator = _make_interpolator()
    interpolator.add_equality_constraints([0, 5], [100.0, 50.0], name="fixed")
    interpolator.solve_system("cg")
    np.testing.assert_allclose(interpolator.c[[0, 5]], [100.0, 50.0], rtol=1e-2)
    # a new setup removes the constraint
    interpolator.setup()
    assert "fixed" not in interpolator.constraints


def test_native_setup_sends_equality_as_node_values():
    interpolator = _make_interpolator("native")
    interpolator.add_equality_constraints([0, 5], [100.0, 50.0], name="fixed")
    assert "fixed" not in interpolator.constraints
    _native.check_supported(interpolator, interpolator._setup_kwargs, {})
    calls = []

    class Problem:
        def add_value(self, points, values, weights):
            calls.append((points, values, weights))

        def __getattr__(self, name):
            return lambda *args, **kwargs: None

    _native._add_data(Problem(), interpolator)
    points, values, weights = calls[-1]
    np.testing.assert_allclose(points.reshape(-1, 3), interpolator.support.nodes[[0, 5]])
    np.testing.assert_allclose(values, [100.0, 50.0])
    assert np.all(weights > 1.0)
