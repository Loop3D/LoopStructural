"""Internal solver strategy helpers for discrete interpolation.

This module centralises backend-specific solve behavior (CG, LSMR, ADMM)
so interpolator classes can focus on orchestration and state management.
"""
from __future__ import annotations

import inspect
from typing import Callable

import numpy as np
from scipy import sparse

# scipy's default rtol (1e-5) stops P2 solves up to 0.6 from the converged field
DEFAULT_CG_RTOL = 1e-8

# The stop rule of a CG solve with a position tolerance (cg_to_position_tolerance).
# The stop is tested every CG_POSITION_CHECK iterations, with the change of the
# solution over the last CG_POSITION_DELAY iterations.
CG_POSITION_DELAY = 50
CG_POSITION_CHECK = 10
# CG on these systems often starts with a plateau: the solution hardly moves
# for hundreds of iterations and the smallest Ritz value is still far above
# the smallest eigenvalue. The stop is not tested before the residual has
# dropped by this factor.
CG_POSITION_GUARD_RTOL = 1e-6
# Nothing more can be gained below this relative residual.
CG_POSITION_FLOOR_RTOL = 1e-13


def resolve_solver_choice(
    solver: Callable[[sparse.csr_matrix, np.ndarray], np.ndarray] | str | None, logger
):
    if callable(solver):
        return solver
    if isinstance(solver, str) or solver is None:
        if solver not in ["cg", "lsmr", "admm", "native", None]:
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


def _set_cg_tolerances(solver_kwargs: dict, tol: float | None) -> None:
    """Fill in the CG stopping tolerances that the caller did not give.

    scipy's cg stops at ``max(rtol * |b|, atol)``. A ``tol`` sets ``atol``;
    without an ``rtol`` scipy would use its own default 1e-5, which then
    decides the stop and makes ``tol`` meaningless. So ``rtol`` defaults to
    ``DEFAULT_CG_RTOL`` in both cases.
    """
    if tol is not None and "atol" not in solver_kwargs:
        solver_kwargs["atol"] = tol
    if "rtol" not in solver_kwargs and (tol is not None or "atol" not in solver_kwargs):
        solver_kwargs["rtol"] = DEFAULT_CG_RTOL


def _smallest_ritz_value(diagonal: list, off_diagonal: list) -> float:
    """Smallest eigenvalue of the Lanczos matrix of the CG steps so far."""
    from scipy.linalg import eigvalsh_tridiagonal

    if len(diagonal) == 1:
        return float(diagonal[0])
    return float(
        eigvalsh_tridiagonal(
            np.asarray(diagonal), np.asarray(off_diagonal), select="i", select_range=(0, 0)
        )[0]
    )


def cg_to_position_tolerance(
    N,
    rhs: np.ndarray,
    position_tolerance: float,
    position_error: Callable[[np.ndarray, np.ndarray], float],
    x0: np.ndarray | None = None,
    maxiter: int | None = None,
    callback: Callable[[np.ndarray], None] | None = None,
) -> tuple[np.ndarray, int, dict]:
    """Conjugate gradients on ``N x = rhs`` that stop on an estimate of the
    position error of the isosurfaces.

    The iteration is the one of ``scipy.sparse.linalg.cg`` without a
    preconditioner (the same operations in the same order), so the iterates
    are the iterates of scipy. Only the stop differs. Every
    ``CG_POSITION_CHECK`` iterations, once the residual is below
    ``CG_POSITION_GUARD_RTOL * |rhs|``, the solve stops if::

        max(1, bound / |change|) * position_error(x_k, change) <= position_tolerance

        change = x_k - x_(k - CG_POSITION_DELAY)
        bound  = |r_k| / theta_k

    Two estimates of the error ``x - x_k`` are combined:

    * ``change`` is the delay estimate of the error: the vector form of the
      estimate of Strakos and Tichy (2002) for the A-norm of the error, from
      the same CG steps. It is a **lower** estimate and it is far too small
      while CG stagnates (measured: up to 1000 times).
    * ``bound`` is the classical bound ``|x - x_k| <= |r_k| / lambda_min``
      of the 2-norm of the error, with the smallest eigenvalue replaced by
      ``theta_k``, the smallest Ritz value of the Lanczos matrix that the CG
      coefficients define (Golub and Meurant; Meurant and Tichy 2019).
      ``theta_k`` comes down to ``lambda_min`` from above, so ``bound`` is
      an upper bound only once ``theta_k`` has converged; before that it is
      too small by the factor ``theta_k / lambda_min``.

    The bound gives a size but no place. The rule assumes that the error is
    spread over the unknowns like the recent change, and scales the change
    up to the size of the bound. This assumption has no proof. So the result
    is an **estimate, not a guarantee** (the measured cases are given in
    ``DiscreteInterpolator.solve_system``).

    Parameters
    ----------
    N
        Symmetric positive definite matrix or ``LinearOperator``.
    rhs
        Right-hand side.
    position_tolerance
        Requested position error, in the length units of ``position_error``.
    position_error
        ``position_error(x, dx)``: the position error, as a length, that a
        change ``dx`` of the solution ``x`` means. It must be a norm of
        ``dx`` (it is scaled linearly). ``inf`` if it has no meaning (for
        example a field without a gradient).
    x0, maxiter, callback
        As in ``scipy.sparse.linalg.cg``.

    Returns
    -------
    x, info, report
        ``info`` is 0 when a stop rule was met and ``maxiter`` otherwise, as
        in scipy. ``report`` has ``iterations``, ``relative_residual``,
        ``position_error_estimate`` (the last tested estimate; NaN if none
        was tested), ``position_error_change`` (the lower estimate in it),
        ``smallest_ritz_value`` and ``stop``: ``"position"``, ``"residual"``
        (the floor ``CG_POSITION_FLOOR_RTOL`` was reached first),
        ``"zero_rhs"`` or ``"maxiter"``.
    """
    rhs = np.asarray(rhs, dtype=float).reshape(-1)
    n = rhs.shape[0]
    matvec = N.matvec if hasattr(N, "matvec") else (lambda v: N @ v)
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float).reshape(-1)
    report = {
        "iterations": 0,
        "relative_residual": 0.0,
        "position_error_estimate": float("nan"),
        "position_error_change": float("nan"),
        "smallest_ritz_value": float("nan"),
        "stop": "zero_rhs",
    }
    bnrm2 = np.linalg.norm(rhs)
    if bnrm2 == 0:
        return rhs.copy(), 0, report
    if maxiter is None:
        maxiter = n * 10
    r = rhs - np.asarray(matvec(x)).reshape(-1) if x.any() else rhs.copy()
    guard = CG_POSITION_GUARD_RTOL * bnrm2
    floor = CG_POSITION_FLOOR_RTOL * bnrm2
    lag = CG_POSITION_DELAY // CG_POSITION_CHECK
    # x at the last checks: history[0] is CG_POSITION_DELAY iterations old
    # when it holds lag entries
    history: list[np.ndarray] = [x.copy()]
    # the Lanczos matrix of the CG steps (Golub and Van Loan, section 11.3)
    diagonal: list[float] = []
    off_diagonal: list[float] = []
    rho_prev, p, alpha_prev = None, None, None
    info, stop = maxiter, "maxiter"
    iteration = 0
    for iteration in range(maxiter):
        rnorm = np.linalg.norm(r)
        if rnorm < floor:
            info, stop = 0, "residual"
            break
        if iteration > 0 and iteration % CG_POSITION_CHECK == 0:
            if len(history) == lag and rnorm < guard:
                change = x - history[0]
                lower = float(position_error(x, change))
                theta = _smallest_ritz_value(diagonal, off_diagonal)
                change_norm = float(np.linalg.norm(change))
                report["position_error_change"] = lower
                report["smallest_ritz_value"] = theta
                if theta > 0 and change_norm > 0:
                    estimate = max(1.0, float(rnorm) / theta / change_norm) * lower
                    report["position_error_estimate"] = estimate
                    if estimate <= position_tolerance:
                        info, stop = 0, "position"
                        break
            history.append(x.copy())
            del history[:-lag]
        rho_cur = np.dot(r, r)
        if iteration > 0:
            beta = rho_cur / rho_prev
            off_diagonal.append(float(np.sqrt(beta) / alpha_prev))
            p *= beta
            p += r
        else:
            beta = 0.0
            p = r.copy()
        q = np.asarray(matvec(p)).reshape(-1)
        alpha = rho_cur / np.dot(p, q)
        diagonal.append(float(1.0 / alpha + (beta / alpha_prev if iteration > 0 else 0.0)))
        x += alpha * p
        r -= alpha * q
        rho_prev, alpha_prev = rho_cur, alpha
        if callback:
            callback(x)
    else:
        iteration = maxiter
    report["iterations"] = iteration
    report["relative_residual"] = float(np.linalg.norm(r) / bnrm2)
    report["stop"] = stop
    return x, info, report


def _solve_cg_position(N, rhs, position_stop, solver_kwargs, timing, logger, x0):
    """Run cg_to_position_tolerance for solve_with_cg(_normal_equations)."""
    from time import perf_counter

    position_tolerance, position_error = position_stop
    if solver_kwargs.get("M") is not None:
        raise ValueError("position_tolerance does not support a cg preconditioner (M)")
    ignored = sorted(k for k in solver_kwargs if k not in ("maxiter", "callback", "M"))
    if ignored:
        logger.warning(
            "position_tolerance sets the cg stop rule; these solver_kwargs are ignored: %s",
            ignored,
        )
    logger.info(f"cg position tolerance: {position_tolerance}")
    solve_step_started = perf_counter()
    x, info, report = cg_to_position_tolerance(
        N,
        rhs,
        position_tolerance,
        position_error,
        x0=x0,
        maxiter=solver_kwargs.get("maxiter"),
        callback=solver_kwargs.get("callback"),
    )
    timing["solve_seconds"] = perf_counter() - solve_step_started
    timing["cg_iterations"] = report["iterations"]
    timing["cg_stop"] = report["stop"]
    timing["cg_relative_residual"] = report["relative_residual"]
    timing["position_tolerance"] = float(position_tolerance)
    timing["position_error_estimate"] = report["position_error_estimate"]
    timing["position_error_change"] = report["position_error_change"]
    if info > 0:
        logger.warning(
            "CG reached iteration limit (%s) and did not converge; using last iteration",
            info,
        )
    elif report["stop"] == "residual":
        logger.warning(
            "CG reached the relative residual %s before the position tolerance %s "
            "(last estimate %s); using last iteration",
            CG_POSITION_FLOOR_RTOL,
            position_tolerance,
            report["position_error_estimate"],
        )
    return np.asarray(x), info == 0


def solve_with_cg(
    A: sparse.spmatrix,
    b: np.ndarray,
    tol: float | None,
    solver_kwargs: dict,
    timing: dict,
    logger,
    x0: np.ndarray | None = None,
    position_stop: tuple[float, Callable[[np.ndarray, np.ndarray], float]] | None = None,
) -> tuple[np.ndarray, bool]:
    from time import perf_counter

    logger.info("Solving using cg")
    if position_stop is not None:
        # tol is not used: the position tolerance sets the stop
        return _solve_cg_position(
            A.T @ A, A.T @ b, position_stop, solver_kwargs, timing, logger, x0
        )
    _set_cg_tolerances(solver_kwargs, tol)

    logger.info(f"Solver kwargs: {solver_kwargs}")
    solve_step_started = perf_counter()
    res = sparse.linalg.cg(A.T @ A, A.T @ b, x0=x0, **solver_kwargs)
    timing["solve_seconds"] = perf_counter() - solve_step_started
    if res[1] > 0:
        logger.warning(
            "CG reached iteration limit (%s) and did not converge; using last iteration",
            res[1],
        )
    return np.asarray(res[0]), res[1] == 0


def solve_with_cg_normal_equations(
    N,
    rhs: np.ndarray,
    tol: float | None,
    solver_kwargs: dict,
    timing: dict,
    logger,
    x0: np.ndarray | None = None,
    position_stop: tuple[float, Callable[[np.ndarray, np.ndarray], float]] | None = None,
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
    if position_stop is not None:
        return _solve_cg_position(N, rhs, position_stop, solver_kwargs, timing, logger, x0)
    _set_cg_tolerances(solver_kwargs, tol)

    logger.info(f"Solver kwargs: {solver_kwargs}")
    solve_step_started = perf_counter()
    res = sparse.linalg.cg(N, rhs, x0=x0, **solver_kwargs)
    timing["solve_seconds"] = perf_counter() - solve_step_started
    if res[1] > 0:
        logger.warning(
            "CG reached iteration limit (%s) and did not converge; using last iteration",
            res[1],
        )
    return np.asarray(res[0]), res[1] == 0


def solve_with_lsmr(
    A: sparse.spmatrix,
    b: np.ndarray,
    tol: float | None,
    solver_kwargs: dict,
    timing: dict,
    logger,
    x0: np.ndarray | None = None,
) -> tuple[np.ndarray, bool]:
    from time import perf_counter

    logger.info("Solving using lsmr")
    if "btol" not in solver_kwargs and tol is not None:
        solver_kwargs["btol"] = tol
        solver_kwargs["atol"] = 0.0
        logger.info(f"Setting lsmr btol to {tol}")
    if x0 is not None:
        solver_kwargs["x0"] = x0
    logger.info(f"Solver kwargs: {solver_kwargs}")
    solve_step_started = perf_counter()
    res = sparse.linalg.lsmr(A, b, **solver_kwargs)
    timing["solve_seconds"] = perf_counter() - solve_step_started

    if res[1] == 0:
        logger.warning("Solution to least squares problem is all zeros, check input data")
    elif res[1] in (3, 6):
        logger.warning("COND(A) seems to be greater than CONLIM, check input data")
    elif res[1] == 7:
        logger.warning("LSMR reached iteration limit (maxiter) and did not converge")
    return np.asarray(res[0]), res[1] != 7


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
        "warn_on_nonconvergence": solver_kwargs.pop("warn_on_nonconvergence", True),
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
    x0: np.ndarray | None = None,
) -> tuple[np.ndarray, list | None, bool]:
    from time import perf_counter

    from .loopsolver import admm_solve

    legacy_x0 = solver_kwargs.pop("x0", None)
    if x0 is None and legacy_x0 is not None:
        x0 = legacy_x0(support) if callable(legacy_x0) else legacy_x0
    if x0 is None:
        x0 = np.zeros(A.shape[1])
    x0 = np.array(x0, dtype=float)

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
