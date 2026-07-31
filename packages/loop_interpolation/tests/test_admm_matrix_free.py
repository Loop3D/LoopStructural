"""Tests for the matrix-free CG normal-equations path in the ADMM solver.

These tests exercise ``admm_solve`` directly (see
``loop_interpolation.loopsolver.admm_solver``) with a small inequality-constrained
least-squares problem, comparing the default (materialized Gram matrix) path
against the opt-in ``matrix_free=True`` path.
"""

import numpy as np
from loop_interpolation.loopsolver import admm_solve
from loop_interpolation.loopsolver import admm_solver as admm_solver_module
from loop_interpolation.loopsolver.admm_solver import (
    _normal_equations_diagonal,
    _normal_equations_operator,
)
from scipy import sparse
from scipy.sparse.linalg import LinearOperator


def _build_inequality_problem(seed=42):
    """Small, well-conditioned inequality-constrained least-squares problem.

    A stacks random data-fit rows with a ridge (identity) block so the
    normal equations are well conditioned for CG. Q is a random sparse
    matrix (not identity) so matrix.T @ matrix has real fill-in, matching
    the scenario the matrix-free path is meant to help with.
    """
    rng = np.random.default_rng(seed)
    n = 15

    A_data = sparse.random(20, n, density=0.3, format="csr", random_state=rng)
    A_ridge = sparse.identity(n, format="csr") * 0.2
    A = sparse.vstack([A_data, A_ridge]).tocsr()

    x_true = rng.uniform(-1.0, 1.0, size=n)
    b = np.asarray(A @ x_true, dtype=float)

    Q = sparse.random(10, n, density=0.4, format="csr", random_state=rng)
    qx_true = np.asarray(Q @ x_true).ravel()
    bounds = np.column_stack(
        [qx_true - 0.05, qx_true + 0.05, np.ones_like(qx_true)]
    )

    x0 = np.zeros(n)
    return A, b, Q, bounds, x0


COMMON_KWARGS = dict(
    admm_weight=0.05,
    nmajor=40,
    linsys_solver="cg",
    cg_preconditioner=True,
    active_set=False,
    reuse_inner_solve=True,
)


def test_matrix_free_matches_dense_cg_solution():
    A, b, Q, bounds, x0 = _build_inequality_problem()

    sol_dense = admm_solve(A, b, Q, bounds, x0.copy(), matrix_free=False, **COMMON_KWARGS)
    sol_matrix_free = admm_solve(A, b, Q, bounds, x0.copy(), matrix_free=True, **COMMON_KWARGS)

    max_abs_diff = float(np.max(np.abs(sol_dense - sol_matrix_free)))
    max_rel_diff = float(
        max_abs_diff / max(np.max(np.abs(sol_dense)), 1e-12)
    )

    assert np.allclose(sol_dense, sol_matrix_free, atol=1e-5, rtol=1e-5), (
        f"matrix_free solution diverged from dense solution: "
        f"max_abs_diff={max_abs_diff}, max_rel_diff={max_rel_diff}"
    )


def test_matrix_free_default_is_false_and_preserves_behavior():
    # matrix_free defaults to False; omitting it should behave identically to
    # explicitly passing matrix_free=False.
    A, b, Q, bounds, x0 = _build_inequality_problem()

    sol_default = admm_solve(A, b, Q, bounds, x0.copy(), **COMMON_KWARGS)
    sol_explicit_false = admm_solve(A, b, Q, bounds, x0.copy(), matrix_free=False, **COMMON_KWARGS)

    assert np.array_equal(sol_default, sol_explicit_false)


def test_matrix_free_uses_linear_operator_not_materialized_gram(monkeypatch):
    """The CG solver should receive a LinearOperator (never a materialized
    sparse Gram matrix) as its system operator when matrix_free=True, and the
    opposite (a materialized sparse matrix) when matrix_free=False.
    """
    A, b, Q, bounds, x0 = _build_inequality_problem()

    captured = {}
    real_cg = admm_solver_module.cg

    def spy_cg(lhs, rhs, **kwargs):
        captured["lhs"] = lhs
        captured["lhs_type"] = type(lhs)
        return real_cg(lhs, rhs, **kwargs)

    monkeypatch.setattr(admm_solver_module, "cg", spy_cg)

    small_kwargs = dict(COMMON_KWARGS)
    small_kwargs["nmajor"] = 3

    admm_solve(A, b, Q, bounds, x0.copy(), matrix_free=True, **small_kwargs)
    assert isinstance(captured["lhs"], LinearOperator), (
        f"expected LinearOperator with matrix_free=True, got {captured['lhs_type']}"
    )

    captured.clear()
    admm_solve(A, b, Q, bounds, x0.copy(), matrix_free=False, **small_kwargs)
    assert not isinstance(captured["lhs"], LinearOperator), (
        "matrix_free=False unexpectedly produced a LinearOperator "
        "(materialized Gram matrix expected)"
    )
    assert sparse.issparse(captured["lhs"]), (
        f"expected a materialized sparse Gram matrix, got {captured['lhs_type']}"
    )


def test_normal_equations_operator_matches_explicit_product():
    rng = np.random.default_rng(7)
    M = sparse.random(12, 6, density=0.5, format="csr", random_state=rng)

    op = _normal_equations_operator(M)
    explicit = (M.T @ M).toarray()

    x = rng.standard_normal(M.shape[1])
    assert np.allclose(op @ x, explicit @ x, atol=1e-10)

    # symmetric operator: rmatvec should match matvec
    assert np.allclose(op.rmatvec(x), op.matvec(x))


def test_normal_equations_diagonal_matches_explicit_diagonal():
    rng = np.random.default_rng(11)
    M = sparse.random(9, 5, density=0.6, format="csr", random_state=rng)

    diag_fast = _normal_equations_diagonal(M)
    diag_explicit = np.asarray((M.T @ M).diagonal(), dtype=float)

    assert np.allclose(diag_fast, diag_explicit)
