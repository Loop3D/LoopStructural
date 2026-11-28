"""Internal pipeline helpers for discrete interpolation solve flow.

This module owns solve-stage orchestration primitives that are backend-agnostic:
- timing payload setup/finalisation
- extraction of common solver options
- matrix assembly and preprocessing
- inequality assembly
"""

from time import perf_counter

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator


def extract_constant_norm_options(solver_kwargs: dict, logger):
    """Pop constant-norm options from solver kwargs.

    Returns
    -------
    tuple[int, float, Optional[float]]
        iterations, base weight, optional target norm
    """
    constant_norm_iterations = max(0, int(solver_kwargs.pop("constant_norm_iterations", 0)))
    constant_norm_weight = float(solver_kwargs.pop("constant_norm_weight", 0.0))
    constant_norm_target = solver_kwargs.pop("constant_norm_target", None)
    if constant_norm_target is not None:
        constant_norm_target = float(constant_norm_target)
        if constant_norm_target <= 0.0:
            logger.warning(
                "constant_norm_target must be > 0; disabling explicit target and using auto mode"
            )
            constant_norm_target = None
    return constant_norm_iterations, constant_norm_weight, constant_norm_target


def extract_ridge_options(
    solver_kwargs: dict, add_ridge_default: bool, ridge_factor_default: float, logger
):
    """Pop per-solve ridge options from solver kwargs.

    ``ridge_factor`` scales the identity rows appended to the system, so it
    adds ``ridge_factor**2`` to the diagonal of the normal equations. A
    ``ridge_factor`` of 0 disables the ridge term.

    Returns
    -------
    tuple[bool, float]
        whether to add the ridge term, and the ridge factor
    """
    add_ridge = solver_kwargs.pop("add_ridge_regularisation", None)
    ridge_factor = solver_kwargs.pop("ridge_factor", None)
    if ridge_factor is None:
        ridge_factor = float(ridge_factor_default)
        if add_ridge is None:
            add_ridge = bool(add_ridge_default)
    else:
        ridge_factor = float(ridge_factor)
        if ridge_factor < 0.0:
            logger.warning("ridge_factor must be >= 0; using %s", abs(ridge_factor))
            ridge_factor = abs(ridge_factor)
        if add_ridge is None:
            add_ridge = ridge_factor > 0.0
    return bool(add_ridge) and ridge_factor > 0.0, ridge_factor


def extract_initial_guess(solver_kwargs: dict, previous_c, support, dof: int, logger):
    """Pop ``x0`` and ``warm_start`` from solver kwargs and resolve the
    initial guess for the solver.

    ``x0`` is an array of length ``dof`` or a callable that takes the support
    and returns one. It has priority over ``warm_start``. With
    ``warm_start=True`` the previous solution is the initial guess, if it is
    finite and has the correct length.

    Returns
    -------
    Optional[np.ndarray]
        initial guess, or None for a zero start
    """
    x0 = solver_kwargs.pop("x0", None)
    warm_start = bool(solver_kwargs.pop("warm_start", False))
    if x0 is not None:
        if callable(x0):
            x0 = x0(support)
        x0 = np.array(x0, dtype=float).reshape(-1)
        if x0.shape[0] != dof:
            logger.warning(
                "x0 has %d values but the system has %d degrees of freedom; ignoring x0",
                x0.shape[0],
                dof,
            )
            return None
        return x0
    if not warm_start:
        return None
    if previous_c is None:
        return None
    previous_c = np.array(previous_c, dtype=float).reshape(-1)
    if previous_c.shape[0] != dof or not np.all(np.isfinite(previous_c)):
        logger.info("No usable previous solution for warm start; starting from zero")
        return None
    logger.info("Warm starting solver from previous solution")
    return previous_c


def init_timing(solver_choice) -> dict:
    return {
        "solver": "callable" if callable(solver_choice) else solver_choice,
        "backend": "python",
    }


def assemble_main_system(build_matrix_fn, timing: dict):
    assembly_started = perf_counter()
    A, b = build_matrix_fn()
    timing["assembly_seconds"] = perf_counter() - assembly_started
    timing["matrix_rows"] = int(A.shape[0])
    timing["matrix_cols"] = int(A.shape[1])
    # A matrix-free regularisation system (see FiniteDifferenceInterpolator's
    # `regularisation_matrix_free`) returns a scipy LinearOperator here instead
    # of a sparse matrix, which has no `.nnz`. Record None rather than crashing.
    timing["matrix_nnz"] = int(A.nnz) if hasattr(A, "nnz") else None
    return A, b


def _add_ridge_to_linear_operator(A: LinearOperator, b: np.ndarray, ridge_factor: float):
    """Append ridge-regularisation rows (``ridge_factor * I``) to a matrix-free
    system ``LinearOperator`` without materialising it as a sparse block.

    Equivalent to the explicit ``sparse.vstack([A, sparse.eye(dof) * ridge_factor])``
    used on the concrete-matrix path, but composed via matvec/rmatvec so a
    matrix-free regularisation operator (from
    ``DiscreteInterpolator._combine_explicit_matrix_with_linear_operator``) never
    needs to be converted to a concrete matrix.
    """
    dof = A.shape[1]
    n_rows = A.shape[0]

    def matvec(x):
        x = np.asarray(x).reshape(-1)
        return np.concatenate([A.matvec(x), ridge_factor * x])

    def rmatvec(y):
        y = np.asarray(y).reshape(-1)
        return A.rmatvec(y[:n_rows]) + ridge_factor * y[n_rows:]

    combined = LinearOperator(
        shape=(n_rows + dof, dof), matvec=matvec, rmatvec=rmatvec, dtype=float
    )
    combined_b = np.concatenate([np.asarray(b, dtype=float).reshape(-1), np.zeros(dof)])
    return combined, combined_b


def preprocess_main_system(
    A,
    b,
    add_ridge_regularisation: bool,
    ridge_factor: float,
    apply_scaling_matrix: bool,
    compute_column_scaling_matrix_fn,
    logger,
    timing: dict,
):
    preprocess_started = perf_counter()
    scaling_matrix = None
    is_matrix_free = isinstance(A, LinearOperator)
    if add_ridge_regularisation:
        if is_matrix_free:
            A, b = _add_ridge_to_linear_operator(A, b, ridge_factor)
            logger.info("Adding ridge regularisation to matrix-free interpolation operator")
        else:
            ridge = sparse.eye(A.shape[1]) * ridge_factor
            A = sparse.vstack([A, ridge])
            b = np.hstack([b, np.zeros(A.shape[1])])
            logger.info("Adding ridge regularisation to interpolation matrix")
    if apply_scaling_matrix:
        if is_matrix_free:
            # Should already be prevented upstream (FiniteDifferenceInterpolator.
            # setup_interpolator reverts regularisation_matrix_free to False when
            # apply_scaling_matrix=True), but guard here too: column scaling needs
            # explicit per-column norms and cannot be computed for a LinearOperator.
            logger.warning(
                "apply_scaling_matrix=True requested but the assembled system is a "
                "matrix-free LinearOperator; column scaling requires an explicit sparse "
                "matrix and cannot be applied to a LinearOperator. Skipping column "
                "scaling for this solve."
            )
        else:
            scaling_matrix = compute_column_scaling_matrix_fn(A)
            A = A @ scaling_matrix
    timing["preprocess_seconds"] = perf_counter() - preprocess_started
    return A, b, scaling_matrix


def assemble_inequality_system(build_inequality_matrix_fn, timing: dict):
    inequality_started = perf_counter()
    Q, bounds = build_inequality_matrix_fn()
    timing["inequality_seconds"] = perf_counter() - inequality_started
    timing["inequality_rows"] = int(Q.shape[0])
    return Q, bounds


def finalize_timing(timing: dict, solve_started: float, up_to_date: bool) -> dict:
    timing["total_seconds"] = perf_counter() - solve_started
    timing["up_to_date"] = bool(up_to_date)
    return timing
