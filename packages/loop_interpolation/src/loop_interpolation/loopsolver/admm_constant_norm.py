import importlib
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.sparse import csr_matrix, vstack
from scipy.sparse.linalg import lsmr

from .admm_method import ADMM


@dataclass
class Config:
    verbose: bool = False
    progress: bool = True


progressbar = lambda x: x

try:
    tqdm_module = importlib.import_module("tqdm")

    if Config.progress:
        progressbar = tqdm_module.tqdm
    else:
        progressbar = lambda x: x
except ModuleNotFoundError:
    Config.progress = False
    progressbar = lambda x: x


def admm_solve_constant_norm(
    A: csr_matrix,
    b: np.ndarray,
    Q: csr_matrix,
    bounds: np.ndarray,
    t: np.ndarray,
    update_r: Callable,
    x0: np.ndarray,
    admm_weight: float = 0.1,
    nmajor=200,
    linsys_solver_kwargs={"maxiter": 100},
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
    # if R.shape[0] != t.shape[0]:
    #     raise ValueError("Number of rows in R matrix and t are different")
    # if R.shape[1] != x0.shape[0]:
    #     raise ValueError("Number of columns in R matrix does not match x0")
    n_ie = bounds.shape[0]
    qx_val = np.zeros((Q.shape[0], 1))
    model = np.zeros(A.shape[1])
    model[:] = x0[:]
    # initialise the admm method, sets up the u and v matrices as 0s
    admm_method = ADMM(n_ie)
    b0 = np.zeros(b.shape[0] + t.shape[0])
    b0[: b.shape[0]] = b[:]
    b0[b.shape[0] :] = t[:]
    # the b vector used for the lsqr soln is the size of A + Q
    b = np.zeros(A.shape[0] + t.shape[0] + Q.shape[0])
    A_size = A.shape[0] + t.shape[0]
    xmin = bounds[:, [0]]
    xmax = bounds[:, [1]]
    x0_ADMM = np.zeros(Q.shape[0])
    # scale the Q matrix by the admm f
    Q *= admm_weight

    # matrix = vstack([A, Q])
    for _i in progressbar(range(nmajor)):
        # current model value
        R = update_r(model, _i)
        matrix = vstack([A, R, Q])
        Mx = matrix @ model  # np.dot(A, model)

        qx_val[:, 0] = Mx[A_size:,] / admm_weight
        x0_ADMM = admm_method.admm_method_iterate_admm_array(xmin, xmax, qx_val)
        # print(x0_ADMM, qx_val.shape)
        # raise Exception
        b[:A_size] = b0[:A_size] - Mx[:A_size]
        b[A_size:] = -admm_weight * (qx_val[:, 0] - x0_ADMM)
        cost_data1 = np.linalg.norm(b[:A_size])
        cost_data2 = np.linalg.norm(b0[A_size:])
        model_norm = np.linalg.norm(model)
        if Config.verbose:
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
        x = lsmr(matrix, b, **linsys_solver_kwargs)
        model += x[0]
    return model
