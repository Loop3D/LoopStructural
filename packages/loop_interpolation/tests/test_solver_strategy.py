import logging

import numpy as np
from loop_interpolation import _solver_strategy as strategy
from scipy import sparse


def test_resolve_solver_choice_fallbacks_to_cg_for_unknown_name():
    logger = logging.getLogger("test_solver_strategy")
    resolved = strategy.resolve_solver_choice("not-a-solver", logger)
    assert resolved == "cg"


def test_resolve_solver_choice_accepts_native_without_warning(caplog):
    logger = logging.getLogger("test_solver_strategy")
    with caplog.at_level(logging.WARNING, logger="test_solver_strategy"):
        resolved = strategy.resolve_solver_choice("native", logger)
    assert resolved == "native"
    assert not [record for record in caplog.records if record.levelno >= logging.WARNING]


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


def test_cg_iteration_limit_reports_not_converged():
    import logging

    import numpy as np
    from loop_interpolation import _solver_strategy
    from scipy import sparse

    rng = np.random.default_rng(0)
    A = sparse.csr_matrix(rng.normal(size=(60, 30)))
    b = rng.normal(size=60)
    _, ok = _solver_strategy.solve_with_cg(A, b, None, {"maxiter": 2}, {}, logging.getLogger())
    assert ok is False
    _, ok = _solver_strategy.solve_with_cg(A, b, None, {"maxiter": 500}, {}, logging.getLogger())
    assert ok is True


def _p1_interpolator():
    from loop_interpolation import P1Interpolator, TetMesh

    mesh = TetMesh(origin=np.zeros(3), nsteps=np.full(3, 4), step_vector=np.ones(3))
    interp = P1Interpolator(mesh)
    interp.set_value_constraints(np.array([[1.0, 1.0, 1.0, 0.5, 1.0], [2.0, 2.0, 2.0, 1.0, 1.0]]))
    interp.setup_interpolator(regularisation=1.0)
    return interp


def test_solve_system_ridge_factor_changes_solution():
    interp = _p1_interpolator()
    interp.solve_system("lsmr", solver_kwargs={"ridge_factor": 0.0})
    no_ridge = interp.c.copy()
    interp.solve_system("lsmr", solver_kwargs={"ridge_factor": 1.0})
    ridge = interp.c.copy()
    assert np.all(np.isfinite(ridge))
    # a strong ridge pulls the solution towards zero
    assert np.linalg.norm(ridge) < np.linalg.norm(no_ridge)


def test_solve_system_warm_start_uses_previous_solution(monkeypatch):
    interp = _p1_interpolator()
    for solver in ("cg", "lsmr"):
        interp.solve_system(solver)
        previous = interp.c.copy()
        captured = {}
        name = "solve_with_cg" if solver == "cg" else "solve_with_lsmr"
        original = getattr(strategy, name)

        def spy(*args, x0=None, _original=original, **kwargs):
            captured["x0"] = x0
            return _original(*args, x0=x0, **kwargs)

        monkeypatch.setattr(strategy, name, spy)
        interp.solve_system(solver, solver_kwargs={"warm_start": True})
        monkeypatch.undo()

        assert captured["x0"] is not None
        assert np.all(np.isfinite(captured["x0"]))
        assert np.allclose(interp.c, previous, atol=1e-6)


def test_solve_system_warm_start_needs_fewer_iterations():
    interp = _p1_interpolator()
    interp.apply_scaling_matrix = False
    interp.solve_system("cg")
    iterations = {"n": 0}

    def count(_xk):
        iterations["n"] += 1

    interp.solve_system("cg", solver_kwargs={"warm_start": True, "callback": count})
    warm = iterations["n"]
    iterations["n"] = 0
    interp.solve_system("cg", solver_kwargs={"callback": count})
    cold = iterations["n"]
    assert warm < cold


def test_update_reuses_setup_kwargs():
    import numpy as np
    from loop_interpolation import P1Interpolator, TetMesh

    mesh = TetMesh(origin=np.zeros(3), nsteps=np.full(3, 4), step_vector=np.ones(3))
    interp = P1Interpolator(mesh)
    interp.set_value_constraints(np.array([[1.0, 1.0, 1.0, 0.5, 1.0], [2.0, 2.0, 2.0, 1.0, 1.0]]))
    interp.setup_interpolator(regularisation=3.0)
    interp.solve_system("cg")
    before = {k: float(np.mean(v["w"])) for k, v in interp.constraints.items()}
    interp.up_to_date = False
    interp.update()
    after = {k: float(np.mean(v["w"])) for k, v in interp.constraints.items()}
    assert before == after
