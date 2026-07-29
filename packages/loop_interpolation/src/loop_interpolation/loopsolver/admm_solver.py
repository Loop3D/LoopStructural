import numpy as np
import importlib
import inspect
from .admm_method import ADMM
from dataclasses import dataclass
from scipy.sparse.linalg import lsmr, lsqr, cg, LinearOperator
from scipy.sparse import vstack, csr_matrix, diags


@dataclass
class Config:
    verbose: bool = False
    progress: bool = True


progressbar = lambda x: x


def _normal_equations_operator(matrix: csr_matrix) -> LinearOperator:
    """Return a LinearOperator representing ``matrix.T @ matrix`` (the Gram/normal-
    equations matrix) without ever materializing that product.

    CG only ever needs the *action* of the normal-equations matrix on a vector,
    so we compose the two matvecs (``matrix @ x`` then ``matrix.T @ (...)``)
    instead of forming the (denser, fill-in prone) sparse-sparse product.
    """
    n = matrix.shape[1]

    def matvec(x):
        return matrix.T @ (matrix @ x)

    return LinearOperator((n, n), matvec=matvec, rmatvec=matvec, dtype=float)


def _normal_equations_diagonal(matrix: csr_matrix) -> np.ndarray:
    """Compute ``diag(matrix.T @ matrix)`` directly as a column sum of squares,
    i.e. ``diag[j] = sum_i matrix[i, j] ** 2``, without materializing the full
    matrix product.
    """
    return np.asarray(matrix.multiply(matrix).sum(axis=0), dtype=float).ravel()


try:
    tqdm_module = importlib.import_module("tqdm")

    if Config.progress:
        progressbar = tqdm_module.tqdm
    else:
        progressbar = lambda x: x
except ModuleNotFoundError:
    Config.progress = False
    progressbar = lambda x: x


def admm_solve(
    A: csr_matrix,
    b: np.ndarray,
    Q: csr_matrix,
    bounds: np.ndarray,
    x0: np.ndarray,
    admm_weight: float = 0.1,
    nmajor=200,
    linsys_solver_kwargs=None,
    linsys_solver="lsmr",
    reuse_inner_solve: bool = True,
    adaptive_rho: bool = False,
    adaptive_rho_mu: float = 10.0,
    adaptive_rho_tau: float = 2.0,
    adaptive_rho_min: float = 1e-4,
    adaptive_rho_max: float = 1e3,
    admm_weight_final: float | None = None,
    admm_weight_schedule: str = "geometric",
    active_set: bool = False,
    active_set_padding: float = 1e-3,
    active_set_min_size: int = 0,
    active_set_max_size: int = 0,
    active_set_hysteresis: bool = True,
    active_set_refresh_interval: int = 1,
    cg_preconditioner: bool = True,
    cg_preconditioner_shift: float = 1e-12,
    matrix_free: bool = False,
    inner_rtol_start: float | None = None,
    inner_rtol_end: float | None = None,
    inner_atol_start: float | None = None,
    inner_atol_end: float | None = None,
    inner_maxiter_start: int | None = None,
    inner_maxiter_end: int | None = None,
    inner_maxiter_schedule: str = "linear",
    admm_abs_tol: float = 1e-4,
    admm_rel_tol: float = 1e-3,
    model_update_tol: float = 0.0,
    min_iterations: int = 5,
    return_history: bool = False,
):
    if A.shape[1] != x0.shape[0]:
        raise ValueError("Number of columns in interpolation matrix does not match x0")
    if A.shape[1] != Q.shape[1]:
        raise ValueError(
            "Number of columns in interpolation matrix and inequality matrix are different "
        )
    if Q.shape[0] != bounds.shape[0]:
        raise ValueError("Number of rows in inequality matrix and bounds are different")
    if bounds.shape[1] == 2:
        bounds = np.hstack([bounds, np.ones((bounds.shape[0], 1))])
    if bounds.shape[1] != 3:
        raise ValueError("Bounds must have two columns")
    if A.shape[0] != b.shape[0]:
        raise ValueError("Number of rows in interpolation matrix and b are different")
    n_ie = bounds.shape[0]
    qx_val = np.zeros((Q.shape[0], 1))
    model = np.zeros(A.shape[1])
    model[:] = x0[:]
    # initialise the admm method, sets up the u and v matrices as 0s
    admm_method = ADMM(n_ie)
    b0 = np.zeros(b.shape)
    b0[:] = b[:]
    # the b vector used for the lsqr soln is the size of A + Q
    b = np.zeros(A.shape[0] + Q.shape[0])
    A_size = A.shape[0]
    xmin = bounds[:, [0]]
    xmax = bounds[:, [1]]
    x0_ADMM = np.zeros(Q.shape[0])
    Q_base = Q.tocsr().copy()
    if linsys_solver_kwargs is None:
        linsys_solver_kwargs = {}
    solver_name = linsys_solver.lower() if isinstance(linsys_solver, str) else linsys_solver

    nmajor = int(nmajor)
    if nmajor <= 0:
        raise ValueError("nmajor must be positive")
    if model_update_tol < 0:
        raise ValueError("model_update_tol must be >= 0")
    if inner_maxiter_schedule not in ("linear", "geometric"):
        raise ValueError("inner_maxiter_schedule must be 'linear' or 'geometric'")

    if admm_weight_final is not None:
        if admm_weight_final <= 0:
            raise ValueError("admm_weight_final must be > 0")
        if admm_weight <= 0:
            raise ValueError("admm_weight must be > 0")
        if admm_weight_schedule == "linear":
            rho_schedule = np.linspace(float(admm_weight), float(admm_weight_final), nmajor)
        else:
            rho_schedule = np.geomspace(float(admm_weight), float(admm_weight_final), nmajor)
    else:
        rho_schedule = np.full(nmajor, float(admm_weight), dtype=float)

    adaptive_rho_enabled = bool(adaptive_rho) and admm_weight_final is None

    def _schedule(start, end, n, geometric=False):
        if start is None or end is None:
            return None
        if geometric:
            if start <= 0 or end <= 0:
                raise ValueError("Geometric schedule endpoints must be > 0")
            return np.geomspace(float(start), float(end), n)
        return np.linspace(float(start), float(end), n)

    def _schedule_int(start, end, n, geometric=False):
        if start is None and end is None:
            return None
        if start is None:
            start = end
        if end is None:
            end = start
        start = int(start)
        end = int(end)
        if start <= 0 or end <= 0:
            raise ValueError("inner_maxiter_start/end must be positive when provided")
        if geometric:
            vals = np.geomspace(float(start), float(end), n)
        else:
            vals = np.linspace(float(start), float(end), n)
        return np.maximum(np.rint(vals).astype(int), 1)

    rtol_sched = _schedule(inner_rtol_start, inner_rtol_end, nmajor, geometric=True)
    atol_sched = _schedule(inner_atol_start, inner_atol_end, nmajor, geometric=True)
    if rtol_sched is not None:
        if solver_name in ("lsmr", "lsqr"):
            linsys_solver_kwargs.setdefault("btol", rtol_sched.tolist())
        else:
            linsys_solver_kwargs.setdefault("rtol", rtol_sched.tolist())
    if atol_sched is not None:
        linsys_solver_kwargs.setdefault("atol", atol_sched.tolist())
    maxiter_sched = _schedule_int(
        inner_maxiter_start,
        inner_maxiter_end,
        nmajor,
        geometric=(inner_maxiter_schedule == "geometric"),
    )
    if maxiter_sched is not None:
        maxiter_key = "iter_lim" if solver_name == "lsqr" else "maxiter"
        linsys_solver_kwargs.setdefault(maxiter_key, maxiter_sched.tolist())

    lsmr_params = set(inspect.signature(lsmr).parameters.keys())
    lsqr_params = set(inspect.signature(lsqr).parameters.keys())
    cg_params = set(inspect.signature(cg).parameters.keys())

    def _inject_x0_if_supported(kwargs, x0_guess, supported_params):
        if x0_guess is None:
            return kwargs
        if "x0" in supported_params and "x0" not in kwargs:
            kwargs = dict(kwargs)
            kwargs["x0"] = x0_guess
        return kwargs

    def _active_rows(qx: np.ndarray, xmin_arr: np.ndarray, xmax_arr: np.ndarray, prev_active=None):
        lower_gap = xmin_arr[:, 0] - qx
        upper_gap = qx - xmax_arr[:, 0]
        violation = np.maximum(lower_gap, 0.0) + np.maximum(upper_gap, 0.0)
        distance_to_bounds = np.minimum(np.abs(qx - xmin_arr[:, 0]), np.abs(qx - xmax_arr[:, 0]))

        active = np.logical_or(violation > 0.0, distance_to_bounds <= active_set_padding)
        if prev_active is not None and active_set_hysteresis:
            active = np.logical_or(active, prev_active)

        if active_set_min_size > 0 and np.sum(active) < active_set_min_size:
            score = violation + np.maximum(active_set_padding - distance_to_bounds, 0.0)
            k = int(min(active_set_min_size, qx.shape[0]))
            if k > 0:
                topk = np.argpartition(score, -k)[-k:]
                active[topk] = True
        if active_set_max_size > 0 and np.sum(active) > active_set_max_size:
            score = violation + np.maximum(active_set_padding - distance_to_bounds, 0.0)
            k = int(min(active_set_max_size, qx.shape[0]))
            keep = np.argpartition(score, -k)[-k:]
            constrained = np.zeros_like(active, dtype=bool)
            constrained[keep] = True
            active = constrained
        if np.sum(active) == 0 and qx.shape[0] > 0:
            active[np.argmax(violation)] = True
        return active

    matrix = None
    matrix_transpose = None
    precomputed_lhs = None
    precomputed_M = None

    def _update_system(rho_value: float):
        nonlocal matrix, matrix_transpose, precomputed_lhs, precomputed_M
        Q_scaled = (Q_base * rho_value).tocsr()
        matrix = vstack([A, Q_scaled]).tocsr()
        matrix_transpose = matrix.T
        precomputed_lhs = None
        precomputed_M = None
        if solver_name == "cg" and not active_set:
            if matrix_free:
                precomputed_lhs = _normal_equations_operator(matrix)
                if cg_preconditioner:
                    diag_vals = _normal_equations_diagonal(matrix)
                    diag_vals = np.maximum(diag_vals, cg_preconditioner_shift)
                    inv_diag = 1.0 / diag_vals
                    precomputed_M = diags(inv_diag)
            else:
                precomputed_lhs = matrix_transpose @ matrix
                if cg_preconditioner:
                    diag_vals = np.asarray(precomputed_lhs.diagonal(), dtype=float)
                    diag_vals = np.maximum(diag_vals, cg_preconditioner_shift)
                    inv_diag = 1.0 / diag_vals
                    precomputed_M = diags(inv_diag)

    rho = float(rho_schedule[0])
    _update_system(rho)

    def _solve_linsys(system_matrix, rhs, kwargs):
        x0_guess = kwargs.pop("_x0_guess", None)
        if callable(solver_name):
            try:
                return np.asarray(
                    solver_name(system_matrix, rhs, x0=x0_guess, **kwargs), dtype=float
                )
            except TypeError:
                return np.asarray(solver_name(system_matrix, rhs, **kwargs), dtype=float)
        if solver_name == "lsmr":
            kwargs = _inject_x0_if_supported(kwargs, x0_guess, lsmr_params)
            res = lsmr(system_matrix, rhs, **kwargs)
            return res[0]
        if solver_name == "lsqr":
            kwargs = _inject_x0_if_supported(kwargs, x0_guess, lsqr_params)
            res = lsqr(system_matrix, rhs, **kwargs)
            return res[0]
        if solver_name == "cg":
            cg_kwargs = dict(kwargs)
            if "btol" in cg_kwargs and "rtol" not in cg_kwargs:
                cg_kwargs["rtol"] = cg_kwargs.pop("btol")
            cg_kwargs = _inject_x0_if_supported(cg_kwargs, x0_guess, cg_params)
            lhs = cg_kwargs.pop("_lhs", None)
            rhs_normal = cg_kwargs.pop("_rhs_normal", None)
            M = cg_kwargs.pop("_M", None)
            if lhs is None:
                if matrix_free and not active_set:
                    lhs = _normal_equations_operator(system_matrix)
                else:
                    lhs = system_matrix.T @ system_matrix
            if rhs_normal is None:
                rhs_normal = system_matrix.T @ rhs
            if M is None and cg_preconditioner:
                if isinstance(lhs, LinearOperator):
                    diag_vals = _normal_equations_diagonal(system_matrix)
                else:
                    diag_vals = np.asarray(lhs.diagonal(), dtype=float)
                diag_vals = np.maximum(diag_vals, cg_preconditioner_shift)
                inv_diag = 1.0 / diag_vals
                M = diags(inv_diag)
            if M is not None and "M" not in cg_kwargs:
                if isinstance(M, LinearOperator):
                    cg_kwargs["M"] = M
                else:
                    cg_kwargs["M"] = M
            dx, info = cg(lhs, rhs_normal, **cg_kwargs)
            if info < 0:
                raise ValueError("CG inner solve failed")
            return dx
        raise ValueError(
            f"Unknown linsys_solver '{linsys_solver}'. Use one of: 'lsmr', 'lsqr', 'cg', or a callable."
        )

    history = []
    inner_x0 = None
    active_rows_prev = None
    use_cached_normal_eq = solver_name == "cg" and not active_set
    for k in linsys_solver_kwargs:
        if (
            not hasattr(linsys_solver_kwargs[k], "__len__")
            or len(linsys_solver_kwargs[k]) != nmajor
        ):
            linsys_solver_kwargs[k] = [linsys_solver_kwargs[k]] * nmajor
    for _i in progressbar(range(nmajor)):
        target_rho = float(rho_schedule[_i])
        if target_rho != rho:
            rho_old = rho
            rho = target_rho
            if rho_old > 0:
                admm_method.u *= rho_old / rho
            _update_system(rho)

        z_prev = admm_method.z.copy()
        # current model value
        Mx = matrix @ model  # np.dot(A, model)
        b[:A_size] = b0[:A_size] - Mx[:A_size]

        if Q.shape[0] > 0:
            qx_val[:, 0] = Mx[A_size:,] / rho
            x0_ADMM = admm_method.admm_method_iterate_admm_array(xmin, xmax, qx_val)
            # print(x0_ADMM, qx_val.shape)
            # raise Exception
            b[A_size:] = -rho * (qx_val[:, 0] - x0_ADMM)
        # cost_data1 = np.linalg.norm(b[:A_size])
        # cost_data2 = np.linalg.norm(b0[A_size:])
        # model_norm = np.linalg.norm(model)
        if Config.verbose:
            cost_data1 = np.linalg.norm(b[:A_size])
            cost_data2 = np.linalg.norm(b0[:A_size])
            model_norm = np.linalg.norm(model)
            cost_data = -1.0
            cost_data_model = 0.0
            if cost_data2 > 0:
                cost_data = cost_data1 / cost_data2
            if model_norm > 0:
                cost_data_model = cost_data1 / model_norm
            cost_admm1 = np.linalg.norm(qx_val - admm_method.z)
            cost_admm2 = np.linalg.norm(admm_method.z)
            cost_admm = -1.0
            if cost_admm2 > 0:
                cost_admm = cost_admm1 / cost_admm2
            print("----------------------------------------")
            print(f"it = {_i}")
            print("cost_data = ", cost_data)
            print("cost_data_model = ", cost_data_model)
            print("cost_admm = ", cost_admm)
            print("----------------------------------------")
        linsys_kwargs = {k: v[_i] for k, v in linsys_solver_kwargs.items()}
        rhs = b
        solve_matrix = matrix
        if active_set and Q.shape[0] > 0:
            refresh = (
                active_rows_prev is None
                or int(active_set_refresh_interval) <= 1
                or (_i % int(active_set_refresh_interval) == 0)
            )
            if refresh:
                active_rows_prev = _active_rows(
                    qx_val[:, 0], xmin, xmax, prev_active=active_rows_prev
                )
            q_active = (Q_base * rho)[active_rows_prev]
            solve_matrix = vstack([A, q_active]).tocsr()
            rhs = np.hstack([b[:A_size], b[A_size:][active_rows_prev]])

        if use_cached_normal_eq:
            linsys_kwargs["_lhs"] = precomputed_lhs
            linsys_kwargs["_rhs_normal"] = matrix_transpose @ rhs
            if precomputed_M is not None:
                linsys_kwargs["_M"] = precomputed_M
        if reuse_inner_solve and inner_x0 is not None:
            linsys_kwargs["_x0_guess"] = inner_x0
        dx = _solve_linsys(solve_matrix, rhs, linsys_kwargs)
        dx_norm = float(np.linalg.norm(dx))
        model_norm_before = float(np.linalg.norm(model))
        model += dx
        if reuse_inner_solve:
            inner_x0 = np.asarray(dx, dtype=float)

        if Q.shape[0] > 0:
            primal_residual = qx_val[:, 0] - admm_method.z
            primal_norm = float(np.linalg.norm(primal_residual))
            dual_norm = float(rho * np.linalg.norm(admm_method.z - z_prev))
            qx_norm = float(np.linalg.norm(qx_val[:, 0]))
            z_norm = float(np.linalg.norm(admm_method.z))
            u_norm = float(np.linalg.norm(rho * admm_method.u))
            n_ie_sqrt = float(np.sqrt(n_ie))
            eps_pri = float(n_ie_sqrt * admm_abs_tol + admm_rel_tol * max(qx_norm, z_norm))
            eps_dual = float(n_ie_sqrt * admm_abs_tol + admm_rel_tol * u_norm)

            if return_history:
                history.append(
                    {
                        "iteration": int(_i + 1),
                        "primal_norm": primal_norm,
                        "dual_norm": dual_norm,
                        "eps_pri": eps_pri,
                        "eps_dual": eps_dual,
                        "rho": float(rho),
                        "active_count": int(np.sum(active_rows_prev))
                        if active_set and active_rows_prev is not None
                        else int(n_ie),
                    }
                )

            if (_i + 1) >= min_iterations and primal_norm <= eps_pri and dual_norm <= eps_dual:
                break

            if model_update_tol > 0.0 and (_i + 1) >= min_iterations:
                model_scale = max(model_norm_before, 1.0)
                if dx_norm <= model_update_tol * model_scale:
                    if return_history:
                        history[-1]["stopped_on_model_update"] = True
                        history[-1]["dx_norm"] = dx_norm
                        history[-1]["model_norm"] = model_norm_before
                    break

            if adaptive_rho_enabled and (_i + 1) >= min_iterations:
                rho_new = rho
                if primal_norm > adaptive_rho_mu * dual_norm:
                    rho_new = min(rho * adaptive_rho_tau, adaptive_rho_max)
                elif dual_norm > adaptive_rho_mu * primal_norm:
                    rho_new = max(rho / adaptive_rho_tau, adaptive_rho_min)
                if rho_new != rho:
                    if rho > 0:
                        admm_method.u *= rho / rho_new
                    rho = float(rho_new)
                    _update_system(rho)
        elif model_update_tol > 0.0 and (_i + 1) >= min_iterations:
            model_scale = max(model_norm_before, 1.0)
            if dx_norm <= model_update_tol * model_scale:
                break

    if return_history:
        return model, history
    return model
