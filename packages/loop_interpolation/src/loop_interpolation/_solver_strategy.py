"""Internal solver strategy helpers for discrete interpolation.

This module centralises backend-specific solve behavior (CG, LSMR, ADMM)
so interpolator classes can focus on orchestration and state management.
"""
from __future__ import annotations

import inspect
from typing import Callable, Optional, Union

import numpy as np
from scipy import sparse


def resolve_solver_choice(
    solver: Optional[Union[Callable[[sparse.csr_matrix, np.ndarray], np.ndarray], str]], logger
):
    if callable(solver):
        return solver
    if isinstance(solver, str) or solver is None:
        if solver not in ["cg", "lsmr", "admm", None]:
            logger.warning(
                "Unknown solver %s using cg. Available solvers are cg, lsmr, admm or a custom callable",
                solver,
            )
            return "cg"
        return "cg" if solver is None else solver
    logger.warning("Unsupported solver type %s, using cg", type(solver).__name__)
    return "cg"


def solve_with_callable(
    solver: Callable[[sparse.csr_matrix, np.ndarray], np.ndarray],
    A: sparse.spmatrix,
    b: np.ndarray,
    timing: dict,
    logger,
) -> tuple[np.ndarray, bool]:
    from time import perf_counter

    logger.warning("Using custom solver")
    solve_step_started = perf_counter()
    c = solver(A.tocsr(), b)
    timing["solve_seconds"] = perf_counter() - solve_step_started
    return np.asarray(c), True


def solve_with_cg(
    A: sparse.spmatrix,
    b: np.ndarray,
    tol: Optional[float],
    solver_kwargs: dict,
    timing: dict,
    logger,
) -> tuple[np.ndarray, bool]:
    from time import perf_counter

    logger.info("Solving using cg")
    if "atol" not in solver_kwargs or "rtol" not in solver_kwargs:
        if tol is not None:
            solver_kwargs["atol"] = tol

    logger.info(f"Solver kwargs: {solver_kwargs}")
    solve_step_started = perf_counter()
    res = sparse.linalg.cg(A.T @ A, A.T @ b, **solver_kwargs)
    timing["solve_seconds"] = perf_counter() - solve_step_started
    if res[1] > 0:
        logger.warning(
            "CG reached iteration limit (%s) and did not converge; using last iteration",
            res[1],
        )
    return np.asarray(res[0]), True


def solve_with_cg_normal_equations(
    N,
    rhs: np.ndarray,
    tol: Optional[float],
    solver_kwargs: dict,
    timing: dict,
    logger,
) -> tuple[np.ndarray, bool]:
    """Solve an already-assembled normal-equations system ``N x = rhs`` with CG.

    Used by the matrix-free-regularisation ``cg`` fast path
    (``DiscreteInterpolator._solve_with_cg_fused_regularisation``): unlike
    :func:`solve_with_cg`, ``N``/``rhs`` here are already ``A^T A`` / ``A^T b``
    (assembled with a fused, boundary-corrected regularisation contribution),
    so this must NOT square an already-square system again.
    """
    from time import perf_counter

    logger.info("Solving using cg (matrix-free fused regularisation normal equations)")
    if "atol" not in solver_kwargs or "rtol" not in solver_kwargs:
        if tol is not None:
            solver_kwargs["atol"] = tol

    logger.info(f"Solver kwargs: {solver_kwargs}")
    solve_step_started = perf_counter()
    res = sparse.linalg.cg(N, rhs, **solver_kwargs)
    timing["solve_seconds"] = perf_counter() - solve_step_started
    if res[1] > 0:
        logger.warning(
            "CG reached iteration limit (%s) and did not converge; using last iteration",
            res[1],
        )
    return np.asarray(res[0]), True


def solve_with_lsmr(
    A: sparse.spmatrix,
    b: np.ndarray,
    tol: Optional[float],
    solver_kwargs: dict,
    timing: dict,
    logger,
) -> tuple[np.ndarray, bool]:
    from time import perf_counter

    logger.info("Solving using lsmr")
    if "btol" not in solver_kwargs:
        if tol is not None:
            solver_kwargs["btol"] = tol
            solver_kwargs["atol"] = 0.0
            logger.info(f"Setting lsmr btol to {tol}")
    logger.info(f"Solver kwargs: {solver_kwargs}")
    solve_step_started = perf_counter()
    res = sparse.linalg.lsmr(A, b, **solver_kwargs)
    timing["solve_seconds"] = perf_counter() - solve_step_started

    if res[1] in (1, 2, 4, 5, 7):
        return np.asarray(res[0]), True
    if res[1] == 0:
        logger.warning("Solution to least squares problem is all zeros, check input data")
        return np.asarray(res[0]), True
    if res[1] in (3, 6):
        logger.warning("COND(A) seems to be greater than CONLIM, check input data")
        return np.asarray(res[0]), True
    return np.asarray(res[0]), True


def extract_admm_kwargs(solver_kwargs: dict, admm_solve) -> tuple[str, dict]:
    linsys_solver = solver_kwargs.pop("linsys_solver", "lsmr")
    admm_kwargs = {
        "batch_size": solver_kwargs.pop("batch_size", None),
        "batch_fraction": solver_kwargs.pop("batch_fraction", None),
        "random_seed": solver_kwargs.pop("random_seed", None),
        "adaptive_rho": solver_kwargs.pop("adaptive_rho", False),
        "adaptive_rho_mu": solver_kwargs.pop("adaptive_rho_mu", 10.0),
        "adaptive_rho_tau": solver_kwargs.pop("adaptive_rho_tau", 2.0),
        "adaptive_rho_min": solver_kwargs.pop("adaptive_rho_min", 1e-4),
        "adaptive_rho_max": solver_kwargs.pop("adaptive_rho_max", 1e3),
        "admm_weight_final": solver_kwargs.pop("admm_weight_final", None),
        "admm_weight_schedule": solver_kwargs.pop("admm_weight_schedule", "geometric"),
        "admm_abs_tol": solver_kwargs.pop("admm_abs_tol", 1e-4),
        "admm_rel_tol": solver_kwargs.pop("admm_rel_tol", 1e-3),
        "min_iterations": solver_kwargs.pop("min_iterations", 5),
        "reuse_inner_solve": solver_kwargs.pop("reuse_inner_solve", True),
        "active_set": solver_kwargs.pop("active_set", False),
        "active_set_padding": solver_kwargs.pop("active_set_padding", 1e-3),
        "active_set_min_size": solver_kwargs.pop("active_set_min_size", 0),
        "active_set_max_size": solver_kwargs.pop("active_set_max_size", 0),
        "active_set_hysteresis": solver_kwargs.pop("active_set_hysteresis", True),
        "active_set_refresh_interval": solver_kwargs.pop("active_set_refresh_interval", 1),
        "cg_preconditioner": solver_kwargs.pop("cg_preconditioner", True),
        "cg_preconditioner_shift": solver_kwargs.pop("cg_preconditioner_shift", 1e-12),
        "matrix_free": solver_kwargs.pop("matrix_free", False),
        "inner_rtol_start": solver_kwargs.pop("inner_rtol_start", None),
        "inner_rtol_end": solver_kwargs.pop("inner_rtol_end", None),
        "inner_atol_start": solver_kwargs.pop("inner_atol_start", None),
        "inner_atol_end": solver_kwargs.pop("inner_atol_end", None),
        "inner_maxiter_start": solver_kwargs.pop("inner_maxiter_start", None),
        "inner_maxiter_end": solver_kwargs.pop("inner_maxiter_end", None),
        "inner_maxiter_schedule": solver_kwargs.pop("inner_maxiter_schedule", "linear"),
        "model_update_tol": solver_kwargs.pop("model_update_tol", 0.0),
        "return_history": solver_kwargs.pop("return_history", False),
    }
    supported_optional = set(inspect.signature(admm_solve).parameters.keys())
    admm_kwargs = {k: v for k, v in admm_kwargs.items() if k in supported_optional}
    return linsys_solver, admm_kwargs


def solve_with_admm(
    A: sparse.spmatrix,
    b: np.ndarray,
    Q: sparse.spmatrix,
    bounds: np.ndarray,
    solver_kwargs: dict,
    timing: dict,
    support,
    logger,
) -> tuple[np.ndarray, Optional[list], bool]:
    from time import perf_counter

    from .loopsolver import admm_solve

    if "x0" in solver_kwargs:
        x0 = solver_kwargs["x0"](support)
    else:
        x0 = np.zeros(A.shape[1])
    solver_kwargs.pop("x0", None)

    linsys_solver, admm_kwargs = extract_admm_kwargs(solver_kwargs, admm_solve)
    solve_step_started = perf_counter()
    res = admm_solve(
        A,
        b,
        Q,
        bounds,
        x0=x0,
        admm_weight=solver_kwargs.pop("admm_weight", 0.01),
        nmajor=solver_kwargs.pop("nmajor", 200),
        linsys_solver_kwargs=solver_kwargs,
        linsys_solver=linsys_solver,
        **admm_kwargs,
    )
    timing["solve_seconds"] = perf_counter() - solve_step_started

    if isinstance(res, tuple):
        return np.asarray(res[0]), res[1], True
    return np.asarray(res), None, True
