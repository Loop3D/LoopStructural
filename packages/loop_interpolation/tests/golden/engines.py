"""Build, set up and solve one interpolator for one case. Deterministic."""
from __future__ import annotations

import time

import numpy as np
from loop_common.math import rng as _package_rng
from loop_interpolation import (
    DiscreteFoldInterpolator,
    FDFoldInterpolator,
    FiniteDifferenceInterpolator,
    FoldEvent,
    P1Interpolator,
    P2Interpolator,
    P2UnstructuredTetMesh,
    StructuredGrid,
    TetMesh,
)

N_NODES = 17  # nodes per axis for FD and P1 (16 cells)
P2_CELLS = 8  # P2 has ~8 nodes per cell: 8 cells gives ~4k unknowns, similar to FD/P1
ENGINES = ("FD", "P1", "P2", "DFI", "FDFold")
# Pinned. The default CG stops early and still returns True (review bug 3), so the golden
# fields use a converged solve. `solve_default` below reproduces what a user gets today.
SOLVER = "cg"
SOLVER_KWARGS = {"maxiter": 20000, "rtol": 1e-10}


class _Feature:
    """Analytic fold-frame coordinate: value = coordinate `axis`, gradient = unit vector."""

    def __init__(self, axis: int):
        self.axis = axis

    def evaluate_value(self, p):
        return np.asarray(p)[:, self.axis].astype(float)

    def evaluate_gradient(self, p):
        g = np.zeros((len(p), 3))
        g[:, self.axis] = 1.0
        return g


class _Frame:
    features = (_Feature(0), _Feature(1))


def make_fold_event(case) -> FoldEvent:
    f = case.fold
    if f.axis_angle_deg is not None:
        return FoldEvent(_Frame(), fold_limb_rotation=f.limb_angle_deg, fold_axis_rotation=f.axis_angle_deg)
    return FoldEvent(_Frame(), fold_limb_rotation=f.limb_angle_deg, fold_axis=f.axis)


def reseed():
    """The fold interpolators shuffle rows with a shared generator. Reset it."""
    _package_rng.bit_generator.state = np.random.default_rng(12345).bit_generator.state


def make_interpolator(engine: str, case):
    n = N_NODES
    origin, step = np.zeros(3), np.full(3, 1.0 / (n - 1))
    if engine in ("FD", "FDFold"):
        grid = StructuredGrid(origin=origin, nsteps=np.full(3, n), step_vector=step)
        if engine == "FD":
            return FiniteDifferenceInterpolator(grid)
        return FDFoldInterpolator(grid, fold=make_fold_event(case))
    if engine in ("P1", "DFI"):
        mesh = TetMesh(origin=origin, nsteps=np.full(3, n), step_vector=step)
        if engine == "P1":
            return P1Interpolator(mesh)
        return DiscreteFoldInterpolator(mesh, fold=make_fold_event(case))
    if engine == "P2":
        mesh = P2UnstructuredTetMesh(origin=origin, step_vector=np.full(3, 1.0 / P2_CELLS), nsteps_cells=np.full(3, P2_CELLS))
        return P2Interpolator(mesh)
    raise ValueError(engine)


def add_data(interp, case):
    interp.set_value_constraints(case.values)
    interp.set_normal_constraints(case.normals)
    if len(case.tangents):
        interp.set_tangent_constraints(case.tangents)


def setup_kwargs(engine: str) -> dict:
    if engine in ("DFI", "FDFold"):
        return {}
    return {}


def run(engine: str, case, solver_kwargs: dict | None = None):
    """Return (interpolator, info). `info` holds timing and the solver flag.

    `solver_kwargs=None` uses the pinned converged solve; pass `{}` for the package default.
    """
    solver_kwargs = SOLVER_KWARGS if solver_kwargs is None else solver_kwargs
    reseed()
    interp = make_interpolator(engine, case)
    add_data(interp, case)
    t0 = time.perf_counter()
    interp.setup_interpolator(**setup_kwargs(engine))
    t1 = time.perf_counter()
    ok = interp.solve_system(SOLVER, solver_kwargs=dict(solver_kwargs))
    t2 = time.perf_counter()
    return interp, {"solved": bool(ok), "t_setup": t1 - t0, "t_solve": t2 - t1, "dof": int(interp.dof)}


def constraint_fingerprint(interp) -> dict:
    """Compact, comparable summary of every constraint family in the assembled system.

    Row order inside a family is part of the contract, so `row_hash` is order sensitive.
    """
    out = {}
    for name, c in interp.constraints.items():
        m = c["matrix"].tocsr()
        m.sort_indices()
        b = np.asarray(c["b"], float)
        w = np.asarray(c["w"], float)
        wm = m.multiply(w[:, None]).tocsr() if w.ndim else m
        out[name] = {
            "shape": list(m.shape),
            "nnz": int(m.nnz),
            "A_sum": float(m.sum()),
            "A_fro": float(np.sqrt(m.multiply(m).sum())),
            "b_sum": float(b.sum()),
            "b_norm": float(np.linalg.norm(b)),
            "w_min": float(np.min(w)),
            "w_max": float(np.max(w)),
            "WA_fro": float(np.sqrt(wm.multiply(wm).sum())),
            # order-sensitive probe: first moment of the row index
            "A_rowmoment": float(m.multiply(np.arange(m.shape[0])[:, None] + 1.0).sum()),
        }
    return out
