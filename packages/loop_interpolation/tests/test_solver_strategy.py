import logging

import numpy as np
from loop_interpolation import _solver_strategy as strategy
from scipy import sparse


def test_resolve_solver_choice_fallbacks_to_cg_for_unknown_name():
    logger = logging.getLogger("test_solver_strategy")
    resolved = strategy.resolve_solver_choice("not-a-solver", logger)
    assert resolved == "cg"


def test_extract_admm_kwargs_filters_to_supported_signature():
    def fake_admm_solve(
        A,
        b,
        Q,
        bounds,
        x0,
        admm_weight,
        nmajor,
        linsys_solver_kwargs,
        linsys_solver,
        adaptive_rho=False,
        return_history=False,
    ):
        return x0

    kwargs = {
        "linsys_solver": "cg",
        "adaptive_rho": True,
        "return_history": True,
        "batch_size": 8,
        "inner_rtol_start": 1e-4,
    }

    linsys_solver, admm_kwargs = strategy.extract_admm_kwargs(kwargs, fake_admm_solve)

    assert linsys_solver == "cg"
    assert "adaptive_rho" in admm_kwargs
    assert "return_history" in admm_kwargs
    assert "batch_size" not in admm_kwargs
    assert "inner_rtol_start" not in admm_kwargs


def test_solve_with_lsmr_applies_tol_defaults(monkeypatch):
    logger = logging.getLogger("test_solver_strategy")
    captured = {}

    def fake_lsmr(A, b, **kwargs):
        captured["kwargs"] = kwargs
        return (np.array([1.0, 2.0]), 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    monkeypatch.setattr(strategy.sparse.linalg, "lsmr", fake_lsmr)

    timing = {}
    A = sparse.eye(2, format="csr")
    b = np.array([1.0, 2.0])
    c, ok = strategy.solve_with_lsmr(A, b, tol=1e-4, solver_kwargs={}, timing=timing, logger=logger)

    assert ok is True
    assert np.allclose(c, np.array([1.0, 2.0]))
    assert captured["kwargs"]["btol"] == 1e-4
    assert captured["kwargs"]["atol"] == 0.0
    assert "solve_seconds" in timing


def test_solve_with_admm_passes_x0_and_returns_history(monkeypatch):
    logger = logging.getLogger("test_solver_strategy")
    captured = {}

    def fake_admm_solve(
        A,
        b,
        Q,
        bounds,
        x0,
        admm_weight,
        nmajor,
        linsys_solver_kwargs,
        linsys_solver,
        return_history=False,
    ):
        captured["x0"] = x0
        captured["admm_weight"] = admm_weight
        captured["nmajor"] = nmajor
        captured["linsys_solver"] = linsys_solver
        captured["return_history"] = return_history
        return np.array([0.25, 0.75]), [{"iteration": 1}]

    from loop_interpolation import loopsolver

    monkeypatch.setattr(loopsolver, "admm_solve", fake_admm_solve)

    timing = {}
    A = sparse.eye(2, format="csr")
    b = np.array([0.0, 0.0])
    Q = sparse.csr_matrix((0, 2), dtype=float)
    bounds = np.zeros((0, 3), dtype=float)
    solver_kwargs = {
        "x0": lambda _support: np.array([1.0, 2.0]),
        "admm_weight": 0.2,
        "nmajor": 5,
        "linsys_solver": "lsmr",
        "return_history": True,
    }

    class _Support:
        pass

    c, history, ok = strategy.solve_with_admm(
        A=A,
        b=b,
        Q=Q,
        bounds=bounds,
        solver_kwargs=solver_kwargs,
        timing=timing,
        support=_Support(),
        logger=logger,
    )

    assert ok is True
    assert np.allclose(c, np.array([0.25, 0.75]))
    assert history == [{"iteration": 1}]
    assert np.allclose(captured["x0"], np.array([1.0, 2.0]))
    assert captured["admm_weight"] == 0.2
    assert captured["nmajor"] == 5
    assert captured["linsys_solver"] == "lsmr"
    assert captured["return_history"] is True
    assert "solve_seconds" in timing
