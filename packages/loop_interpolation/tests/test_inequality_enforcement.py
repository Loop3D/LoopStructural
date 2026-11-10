"""Inequality constraints must be enforced, or at least not silently ignored.

``InterpolatorBuilder.add_inequality_constraints`` takes rows ``x, y, z, lower, upper`` and the
documented meaning is ``lower <= f(x, y, z) <= upper``. Inequality constraints are kept apart
from the equality system and are only used by the ADMM solver.

The scenario is a sparse model on a gaussian dome whose true field is known: one borehole (three
contact points on the 0.3, 0.5 and 0.7 horizons) and two orientations. Thirty random points whose
true value is clearly above or below the 0.5 horizon become one-sided bounds on it:
``0.5 <= f`` above, ``f <= 0.5`` below. Without the bounds this data leaves 4 of the 24 bounded
points on the wrong side of 0.5.

Diagnosis:

* ``solve_system`` column-scales the equality system (``apply_scaling_matrix``, on by default) and
  solves for ``y`` where ``c = S @ y``. The inequality matrix ``Q`` was passed to ADMM unscaled,
  so the bounds were applied to ``y`` instead of to the field ``Q @ c``. With scaling on, even a
  large penalty left bounds violated (``admm_weight=10``: 10 of 24). Fixed by passing ``Q @ S``.
* With that fixed ADMM behaves like ordinary ADMM. Its default penalty ``admm_weight=0.01`` is too
  weak for this problem: after the default 200 outer iterations the primal residual is 0.056
  against a tolerance of 0.0025 (still falling after 5000 iterations: 0.0145), so 4 of 24 bounds
  are still violated. ``admm_weight=0.1`` stops at iteration 150, ``admm_weight=1`` at iteration
  90 with every bound satisfied. The solver reaches ``nmajor`` without meeting its tolerance and
  returns without a warning.
* ``nmajor=10, maxiter=50`` (the arguments of the existing builder test) enforces the bounds
  only by accident: it is an unconverged, truncated solve, and ``maxiter=100`` with the same
  ``nmajor`` violates them again.
* With the ``cg`` solver the inequality data have no effect, and nothing is logged.

The default-penalty and ``cg`` cases are strict xfails. When they start to pass they XPASS, which
fails the run: remove the markers then.
"""
import logging

import numpy as np
import pytest
from loop_common.geometry import BoundingBox
from loop_interpolation import InterpolatorBuilder
from scipy.optimize import brentq

TOLERANCE = 1e-3
LEVEL = 0.5


def truth(xyz):
    r2 = (xyz[:, 0] - 0.5) ** 2 + (xyz[:, 1] - 0.5) ** 2
    return xyz[:, 2] - 0.2 * np.exp(-r2 / (2 * 0.2**2))


def truth_gradient(xyz):
    r2 = (xyz[:, 0] - 0.5) ** 2 + (xyz[:, 1] - 0.5) ** 2
    e = 0.2 * np.exp(-r2 / (2 * 0.2**2)) / 0.2**2
    return np.c_[e * (xyz[:, 0] - 0.5), e * (xyz[:, 1] - 0.5), np.ones(len(xyz))]


def sparse_data():
    contacts = []
    for level in (0.3, 0.5, 0.7):
        z = brentq(lambda z, level=level: truth(np.array([[0.5, 0.5, z]]))[0] - level, 0.0, 1.0)
        contacts.append([0.5, 0.5, z, level])
    orientation_xyz = np.array([[0.4, 0.4, 0.5], [0.6, 0.6, 0.5]])
    normals = np.c_[orientation_xyz, truth_gradient(orientation_xyz)]

    rng = np.random.default_rng(9)
    candidates = rng.uniform(0.1, 0.9, (30, 3))
    rows = []
    for point, value in zip(candidates, truth(candidates)):
        if value > LEVEL + 0.05:
            rows.append([*point, LEVEL, 10.0])  # f >= 0.5
        elif value < LEVEL - 0.05:
            rows.append([*point, -10.0, LEVEL])  # f <= 0.5
    return np.array(contacts), normals, np.array(rows)


def builder(with_inequality=True):
    contacts, normals, inequality = sparse_data()
    bounding_box = BoundingBox(np.array([0.0, 0.0, 0.0]), np.array([1.0, 1.0, 1.0]))
    built = InterpolatorBuilder(
        "FDI", bounding_box=bounding_box, nelements=4000, buffer=0.1
    )
    built.add_value_constraints(contacts).add_normal_constraints(normals)
    if with_inequality:
        built.add_inequality_constraints(inequality)
    return built


def violations(interpolator):
    """How far each bounded point is outside its bounds (0 where satisfied)."""
    _, _, inequality = sparse_data()
    value = interpolator.evaluate_value(inequality[:, :3])
    below = np.maximum(inequality[:, 3] - value, 0.0)
    above = np.maximum(value - inequality[:, 4], 0.0)
    return below + above


def test_the_unconstrained_solution_violates_some_bounds():
    """The scenario is only meaningful if the data alone get some bounded points wrong."""
    interpolator = builder(with_inequality=False).setup_interpolator().solve().build()
    assert np.count_nonzero(violations(interpolator) > TOLERANCE) >= 3


def test_admm_with_the_arguments_of_the_existing_builder_test_enforces_the_bounds():
    interpolator = (
        builder()
        .setup_interpolator()
        .use_solver("admm", nmajor=10, admm_weight=0.01, maxiter=50)
        .solve()
        .build()
    )
    assert violations(interpolator).max() <= TOLERANCE


@pytest.mark.xfail(
    strict=True,
    reason="the default penalty admm_weight=0.01 is too weak here: after 200 iterations the primal "
    "residual (0.056) is far above the tolerance (0.0025) and 4 of 24 bounds are violated (the "
    "solver now warns about this, see the non-convergence tests below)",
)
def test_admm_with_the_default_penalty_enforces_the_bounds():
    interpolator = builder().setup_interpolator().use_solver("admm").solve().build()
    assert violations(interpolator).max() <= TOLERANCE


@pytest.mark.parametrize("apply_scaling_matrix", [True, False])
@pytest.mark.parametrize("admm_weight", [1.0, 10.0])
def test_admm_enforces_the_bounds_with_and_without_column_scaling(apply_scaling_matrix, admm_weight):
    """The inequality rows must be scaled with the equality system (they were not)."""
    built = builder()
    built.interpolator.apply_scaling_matrix = apply_scaling_matrix
    interpolator = built.setup_interpolator().use_solver("admm", admm_weight=admm_weight).solve().build()
    assert violations(interpolator).max() <= TOLERANCE


@pytest.mark.xfail(
    strict=True,
    reason="expectation, not a documented contract: with solver 'cg' the inequality data are "
    "dropped without a warning",
)
def test_cg_with_inequality_constraints_warns_that_they_are_not_used(caplog):
    with caplog.at_level(logging.WARNING):
        builder().setup_interpolator().use_solver("cg").solve().build()
    assert any("inequality" in record.getMessage().lower() for record in caplog.records)


# ------------------------------------------------------------------ non-convergence is reported
def not_converged_warnings(caplog):
    return [r for r in caplog.records if "without converging" in r.getMessage()]


def test_the_default_penalty_that_does_not_converge_says_so(caplog):
    with caplog.at_level(logging.WARNING):
        builder().setup_interpolator().use_solver("admm").solve().build()
    (warning,) = not_converged_warnings(caplog)
    message = warning.getMessage()
    assert "ADMM used all 200 iterations" in message
    assert "primal residual" in message and "tolerance" in message
    assert "admm_weight" in message
    # the primal residual missed its tolerance here, so the bounds really may be violated
    assert "may not be fully satisfied" in message


def test_the_constant_norm_polish_warns_once_not_once_per_iteration(caplog):
    """The polish re-solves with ADMM each iteration; the warning is for the whole solve."""
    with caplog.at_level(logging.WARNING):
        (
            builder()
            .setup_interpolator()
            .use_solver("admm", constant_norm_iterations=3, constant_norm_weight=0.1)
            .solve()
            .build()
        )
    assert len(not_converged_warnings(caplog)) == 1


def test_a_penalty_that_converges_does_not_warn(caplog):
    with caplog.at_level(logging.WARNING):
        interpolator = (
            builder().setup_interpolator().use_solver("admm", admm_weight=1.0).solve().build()
        )
    assert not not_converged_warnings(caplog)
    assert violations(interpolator).max() <= TOLERANCE


def small_problem():
    """Data say the unknowns sum to 3; the first two are bounded to [0, 0.2]."""
    from scipy import sparse

    n = 6
    a = sparse.vstack(
        [sparse.identity(n, format="csr") * 0.2, sparse.csr_matrix(np.ones((1, n)))]
    ).tocsr()
    b = np.r_[np.zeros(n), 3.0]
    q = sparse.csr_matrix(np.eye(n)[:2])
    bounds = np.array([[0.0, 0.2, 1.0], [0.0, 0.2, 1.0]])
    return a, b, q, bounds, np.zeros(n)


def solve_small(**kwargs):
    from loop_interpolation.loopsolver import admm_solve

    a, b, q, bounds, x0 = small_problem()
    return admm_solve(a, b, q, bounds, x0, linsys_solver="lsmr", **kwargs)


def test_solver_converged_in_time_does_not_warn(caplog):
    with caplog.at_level(logging.WARNING):
        x = solve_small(admm_weight=1.0, nmajor=200)
    assert not not_converged_warnings(caplog)
    np.testing.assert_allclose(x[:2], 0.2, atol=1e-3)


def test_solver_out_of_iterations_warns_and_reports_the_residuals(caplog):
    with caplog.at_level(logging.WARNING):
        x = solve_small(admm_weight=1e-4, nmajor=20)
    (warning,) = not_converged_warnings(caplog)
    assert "all 20 iterations" in warning.getMessage()
    assert "the final value used was 0.0001" in warning.getMessage()
    # and the bounds really are violated, so the warning is true
    assert x[:2].max() > 0.2 + 0.1


def test_dual_only_miss_does_not_claim_the_bounds_are_violated(caplog):
    """admm_weight=1 with 20 iterations meets the primal tolerance but not the dual one."""
    with caplog.at_level(logging.WARNING):
        interpolator = (
            builder()
            .setup_interpolator()
            .use_solver("admm", admm_weight=1.0, nmajor=20)
            .solve()
            .build()
        )
    (warning,) = not_converged_warnings(caplog)
    message = warning.getMessage()
    assert "satisfied to tolerance" in message
    assert "may not be fully satisfied" not in message
    assert violations(interpolator).max() <= TOLERANCE  # and the message is true


def test_warning_can_be_silenced_per_solve(caplog):
    with caplog.at_level(logging.WARNING):
        solve_small(admm_weight=1e-4, nmajor=20, warn_on_nonconvergence=False)
    assert not not_converged_warnings(caplog)


def test_solver_stopped_by_model_update_tolerance_does_not_warn(caplog):
    """Stopping because the model stopped changing was asked for; it is not a failure."""
    with caplog.at_level(logging.WARNING):
        solve_small(admm_weight=1e-4, nmajor=20, model_update_tol=10.0)
    assert not not_converged_warnings(caplog)


def test_a_solver_with_fewer_iterations_than_the_minimum_warns(caplog):
    """It cannot be declared converged before min_iterations, so it used all of its iterations."""
    with caplog.at_level(logging.WARNING):
        solve_small(admm_weight=1.0, nmajor=3)
    assert len(not_converged_warnings(caplog)) == 1
