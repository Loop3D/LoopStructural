"""Inequality problems: the native solver must agree with the Python ADMM solver.

Story 6 scenario 1, SC-004. Two cases on a 20^3 grid: value bounds and pair bounds.
For each case, two things are checked:
  1. Each solver keeps the constraint violation below 1e-4 of the field range.
  2. The two fields differ by less than 1e-3 of the field range.
All tests skip while the native core (loop_interpolation_core) is not installed.
"""

import numpy as np
import pytest

# Skip the whole file when the native core is not installed.
pytest.importorskip("loop_interpolation_core")

from loop_interpolation import FiniteDifferenceInterpolator as FDI  # noqa: E402
from loop_interpolation import StructuredGrid  # noqa: E402

VIOLATION_TOL = 1e-4  # fraction of the field range
FIELD_TOL = 1e-3  # fraction of the field range
PAIR_UPPER_BOUND = 1e-6  # value for setup kwarg inequality_pair_upper_bound


# Settings for each solver. Both use tight ADMM tolerances,
# so that the bounds are met to the limit of this test.
SOLVER_KWARGS = {
    "native": {"maxiter": 50000, "admm_abs_tol": 1e-7, "admm_rel_tol": 1e-7},
    # Weight 1.0 makes the bound penalty strong enough for the Python ADMM.
    "admm": {
        "nmajor": 500,
        "admm_weight": 1.0,
        "admm_abs_tol": 1e-7,
        "admm_rel_tol": 1e-7,
        "atol": 1e-10,
        "btol": 1e-10,
        "maxiter": 5000,
    },
}


def _make_interpolator(solver):
    """Make a 20^3 structured grid interpolator that uses the given solver."""
    origin = np.array([-0.1, -0.1, -0.1])
    maximum = np.array([1.1, 1.1, 1.1])
    nsteps = np.array([20, 20, 20])
    grid = StructuredGrid(origin=origin, nsteps=nsteps, step_vector=(maximum - origin) / nsteps)
    interpolator = FDI(grid)
    # Set the solver before setup, so the native path can skip the sparse rows.
    interpolator.solver = solver
    # The Python path scales the columns of the main matrix but not the bound matrix,
    # so the bounds would apply to scaled values. Switch the scaling off for a fair oracle.
    interpolator.apply_scaling_matrix = False
    return interpolator


def _add_ramp_data(interpolator):
    """Add 60 value rows from a smooth field that rises with z (about 0 to 1).

    The many rows make the unbounded field well determined,
    so both solvers can reach the bounds.
    """
    rng = np.random.default_rng(0)
    xyz = rng.uniform(0.1, 0.9, size=(60, 3))
    values = xyz[:, 2] + 0.1 * np.sin(2.0 * xyz[:, 0]) * xyz[:, 1]
    interpolator.set_value_constraints(np.column_stack([xyz, values, np.ones(60)]))


def _evaluation_points():
    """Points on a regular grid in the unit cube, used to measure the field."""
    axis = np.linspace(0.05, 0.95, 12)
    return np.array(np.meshgrid(axis, axis, axis, indexing="ij")).reshape(3, -1).T


def _check_converged(interpolator, solver):
    """The native core must report convergence."""
    if solver == "native":
        report = interpolator.latest_native_report
        assert report["converged"], f"native did not converge: {report}"


def _solve_value_bounds(solver):
    """Ramp field with value bounds that cut the field in the middle."""
    interpolator = _make_interpolator(solver)
    _add_ramp_data(interpolator)
    # Columns: X, Y, Z, lower, upper, weight. The upper bound is below the ramp value.
    bound_xyz = np.array(
        [[0.3, 0.3, 0.7], [0.5, 0.5, 0.7], [0.7, 0.7, 0.7], [0.4, 0.6, 0.8], [0.6, 0.4, 0.8]]
    )
    bounds = np.column_stack(
        [bound_xyz, np.full(5, -5.0), np.full(5, 0.5), np.ones(5)]
    )
    interpolator.set_value_inequality_constraints(bounds)
    interpolator.setup()
    interpolator.solve_system(solver=solver, solver_kwargs=dict(SOLVER_KWARGS[solver]))
    _check_converged(interpolator, solver)
    return interpolator, bound_xyz


def _solve_pair_bounds(solver):
    """Ramp field with pair bounds: every point of rock 0 must be below rock 1."""
    interpolator = _make_interpolator(solver)
    _add_ramp_data(interpolator)
    # Columns: X, Y, Z, rock_id, weight. Rock 0 is high in z and rock 1 is low in z,
    # so the ramp breaks the bound and the solver must act.
    rock0 = np.array([[0.3, 0.3, 0.8], [0.5, 0.5, 0.8], [0.7, 0.7, 0.8]])
    rock1 = np.array([[0.3, 0.7, 0.3], [0.5, 0.5, 0.3], [0.7, 0.3, 0.3]])
    pairs = np.vstack(
        [
            np.column_stack([rock0, np.zeros(3), np.ones(3)]),
            np.column_stack([rock1, np.ones(3), np.ones(3)]),
        ]
    )
    interpolator.set_inequality_pairs_constraints(pairs)
    interpolator.setup(
        inequality_pairs=[(0, 1)],
        inequality_pair_upper_bound=PAIR_UPPER_BOUND,
        inequality_pair_lower_bound=-np.inf,
    )
    interpolator.solve_system(solver=solver, solver_kwargs=dict(SOLVER_KWARGS[solver]))
    _check_converged(interpolator, solver)
    return interpolator, rock0, rock1


def _field_range(field):
    return float(np.nanmax(field) - np.nanmin(field))


def test_native_matches_admm_with_value_bounds():
    native, bound_xyz = _solve_value_bounds("native")
    admm, _ = _solve_value_bounds("admm")
    points = _evaluation_points()
    field_native = native.evaluate_value(points)
    field_admm = admm.evaluate_value(points)
    value_range = _field_range(field_admm)
    assert value_range > 0

    for name, solved in (("native", native), ("admm", admm)):
        at_bounds = solved.evaluate_value(bound_xyz)
        # The upper bound is 0.5 and the lower bound is -5 (far below the field).
        violation = max(float(np.max(at_bounds - 0.5)), 0.0)
        assert violation <= VIOLATION_TOL * value_range, f"{name}: violation {violation}"

    difference = float(np.nanmax(np.abs(field_native - field_admm)))
    assert difference <= FIELD_TOL * value_range, f"max field difference {difference}"


def test_native_matches_admm_with_pair_bounds():
    native, rock0, rock1 = _solve_pair_bounds("native")
    admm, _, _ = _solve_pair_bounds("admm")
    points = _evaluation_points()
    field_native = native.evaluate_value(points)
    field_admm = admm.evaluate_value(points)
    value_range = _field_range(field_admm)
    assert value_range > 0

    for name, solved in (("native", native), ("admm", admm)):
        upper = solved.evaluate_value(rock0)
        lower = solved.evaluate_value(rock1)
        # Each row needs f(rock 0 point) - f(rock 1 point) <= upper bound,
        # for all pairs of points.
        violation = max(float(np.max(upper) - np.min(lower)) - PAIR_UPPER_BOUND, 0.0)
        assert violation <= VIOLATION_TOL * value_range, f"{name}: violation {violation}"

    difference = float(np.nanmax(np.abs(field_native - field_admm)))
    assert difference <= FIELD_TOL * value_range, f"max field difference {difference}"


# --- Cap on the pair-bound rows (finding 5, decision D-R1) ---------------------


class _CountingProblem:
    """Stub of the core problem. It keeps the pair rows and the largest chunk."""

    def __init__(self, keep_rows=False):
        self.keep_rows = keep_rows
        self.n_rows = 0
        self.max_chunk = 0
        self.calls = 0
        self.a = []
        self.b = []

    def add_value_bound(self, *args):
        pass

    def add_pair_bound(self, points_a, points_b, lower, upper, weights):
        self.calls += 1
        self.n_rows += lower.size
        self.max_chunk = max(self.max_chunk, lower.size)
        if self.keep_rows:
            self.a.append(points_a.reshape(-1, 3).copy())
            self.b.append(points_b.reshape(-1, 3).copy())


class _StubInterpolator:
    def __init__(self, pairs):
        self._pairs = pairs

    def get_inequality_value_constraints(self):
        return None

    def get_inequality_pairs_constraints(self):
        return self._pairs


class _RecordingLogger:
    """Replace the module logger. The loop_common logger does not propagate to caplog."""

    def __init__(self):
        self.warnings = []

    def warning(self, message, *args):
        self.warnings.append(message % args)

    def debug(self, *args, **kwargs):
        pass

    info = debug


def _two_rock_pairs(n_each, seed=1):
    rng = np.random.default_rng(seed)
    a = np.column_stack([rng.uniform(0, 1, (n_each, 3)), np.zeros(n_each), np.ones(n_each)])
    b = np.column_stack([rng.uniform(0, 1, (n_each, 3)), np.ones(n_each), np.ones(n_each)])
    return np.vstack([a, b])


def test_pair_bounds_cap_keeps_memory_low_and_warns(monkeypatch):
    import tracemalloc

    from loop_interpolation import _native

    recorder = _RecordingLogger()
    monkeypatch.setattr(_native, "logger", recorder)
    interpolator = _StubInterpolator(_two_rock_pairs(5000))  # 25 million combinations
    problem = _CountingProblem()
    tracemalloc.start()
    try:
        _native._add_inequalities(problem, interpolator, {})
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert peak < 200 * 1024**2, f"peak {peak / 1024**2:.0f} MB"
    assert problem.n_rows == 100_000
    assert problem.max_chunk <= 100_000
    assert len(recorder.warnings) == 1
    text = recorder.warnings[0]
    assert "25000000" in text and "100000" in text


def test_pair_bounds_below_cap_send_all_combinations(monkeypatch):
    from loop_interpolation import _native

    recorder = _RecordingLogger()
    monkeypatch.setattr(_native, "logger", recorder)
    data = _two_rock_pairs(7)
    problem = _CountingProblem(keep_rows=True)
    _native._add_inequalities(problem, _StubInterpolator(data), {})
    assert not recorder.warnings
    a = np.vstack(problem.a)
    b = np.vstack(problem.b)
    # The matrix path pairs every point of the first id with every point of the second.
    expected_a = np.repeat(data[:7, :3], 7, axis=0)
    expected_b = np.tile(data[7:, :3], (7, 1))
    np.testing.assert_array_equal(a, expected_a)
    np.testing.assert_array_equal(b, expected_b)


def test_pair_bounds_sample_is_repeatable_and_inside_the_product():
    from loop_interpolation import _native

    def rows(seed):
        problem = _CountingProblem(keep_rows=True)
        setup = {"inequality_pair_max_rows": 50}
        if seed is not None:
            setup["inequality_pair_seed"] = seed
        _native._add_inequalities(problem, _StubInterpolator(_two_rock_pairs(30)), setup)
        return np.vstack(problem.a), np.vstack(problem.b)

    first, second = rows(None), rows(None)
    np.testing.assert_array_equal(first[0], second[0])
    np.testing.assert_array_equal(first[1], second[1])
    assert first[0].shape[0] == 50
    # No combination appears twice.
    assert np.unique(np.hstack(first), axis=0).shape[0] == 50
    other = rows(5)
    assert not np.array_equal(first[0], other[0]) or not np.array_equal(first[1], other[1])


def test_pair_bounds_bad_cap_is_rejected():
    from loop_interpolation import _native

    with pytest.raises(ValueError, match="inequality_pair_max_rows"):
        _native._add_inequalities(
            _CountingProblem(),
            _StubInterpolator(_two_rock_pairs(3)),
            {"inequality_pair_max_rows": 0},
        )


def test_pair_bounds_nan_rows_are_dropped():
    from loop_interpolation import _native

    data = _two_rock_pairs(4)
    data[0, 0] = np.nan
    problem = _CountingProblem()
    _native._add_inequalities(problem, _StubInterpolator(data), {})
    assert problem.n_rows == 3 * 4


def test_native_with_cap_above_size_matches_admm_and_repeats():
    # Cap above the size of the product: same field as the ADMM oracle.
    def solve(solver, **extra):
        interpolator = _make_interpolator(solver)
        _add_ramp_data(interpolator)
        rock0 = np.array([[0.3, 0.3, 0.8], [0.5, 0.5, 0.8], [0.7, 0.7, 0.8]])
        rock1 = np.array([[0.3, 0.7, 0.3], [0.5, 0.5, 0.3], [0.7, 0.3, 0.3]])
        pairs = np.vstack(
            [
                np.column_stack([rock0, np.zeros(3), np.ones(3)]),
                np.column_stack([rock1, np.ones(3), np.ones(3)]),
            ]
        )
        interpolator.set_inequality_pairs_constraints(pairs)
        interpolator.setup(
            inequality_pairs=[(0, 1)],
            inequality_pair_upper_bound=PAIR_UPPER_BOUND,
            inequality_pair_lower_bound=-np.inf,
            **extra,
        )
        interpolator.solve_system(solver=solver, solver_kwargs=dict(SOLVER_KWARGS[solver]))
        return interpolator.evaluate_value(_evaluation_points())

    first = solve("native", inequality_pair_max_rows=1000)
    second = solve("native", inequality_pair_max_rows=1000)
    admm = solve("admm")
    scale = _field_range(admm)
    assert float(np.max(np.abs(first - admm))) <= FIELD_TOL * scale
    # The same seed gives the same rows. The parallel sums of the core do
    # not add in a fixed order, so two runs can differ in the last bits
    # (SC-003). On macOS the difference was 1.4e-10.
    np.testing.assert_allclose(first, second, rtol=0, atol=1e-9 * scale)

    # A cap below the size of the product (4 of 9 rows) is repeatable with the same seed.
    capped_1 = solve("native", inequality_pair_max_rows=4)
    capped_2 = solve("native", inequality_pair_max_rows=4)
    np.testing.assert_allclose(capped_1, capped_2, rtol=0, atol=1e-9 * scale)
