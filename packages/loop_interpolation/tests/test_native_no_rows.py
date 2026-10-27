"""Tests for the native solver and the sparse rows (FR-019, plan D7.2).

Tests (a) to (c) need the native core and skip when it is not installed.
Test (d) does not need the core.
"""

import logging
import subprocess
import sys

import numpy as np
import pytest
from loop_interpolation import (
    ConstraintDiagnosticsReport,
    FiniteDifferenceInterpolator,
    StructuredGrid,
)

SETUP_KWARGS = {"cpw": 1.0, "gpw": 0.0, "npw": 0.0, "tpw": 0.0, "ipw": 0.0}


def _make_interpolator(solver="cg"):
    grid = StructuredGrid(
        origin=np.array([0.0, 0.0, 0.0]),
        nsteps=np.array([6, 6, 6]),
        step_vector=np.array([1.0, 1.0, 1.0]),
    )
    interpolator = FiniteDifferenceInterpolator(grid)
    # The solver is set before setup. Plain attribute, as in the factory.
    interpolator.solver = solver
    points = np.array(
        [
            [1.5, 2.5, 2.5, 0.0, 1.0],
            [4.5, 2.5, 2.5, 1.0, 1.0],
            [2.5, 3.5, 1.5, 0.5, 1.0],
            [3.5, 1.5, 4.5, 0.7, 1.0],
        ]
    )
    interpolator.set_value_constraints(points)
    return interpolator


def _field(interpolator):
    return np.asarray(interpolator.c, dtype=float).copy()


def test_native_before_setup_builds_no_rows():
    """(a) Native solver set before setup: no sparse rows, report exists."""
    pytest.importorskip("loop_interpolation_core")
    interpolator = _make_interpolator("native")
    report = interpolator.setup(**SETUP_KWARGS)
    assert len(interpolator.constraints) == 0
    assert isinstance(report, ConstraintDiagnosticsReport)
    assert interpolator.latest_diagnostics_report is report


def test_native_after_setup_works_with_one_warning(caplog):
    """(b) Native solver set after setup: rows exist, field is correct, one warning."""
    pytest.importorskip("loop_interpolation_core")
    reference = _make_interpolator("cg")
    reference.setup(**SETUP_KWARGS)
    reference.solve_system("cg", tol=1e-12)
    expected = _field(reference)

    interpolator = _make_interpolator("cg")
    interpolator.setup(**SETUP_KWARGS)
    assert len(interpolator.constraints) > 0
    with caplog.at_level(logging.WARNING):
        interpolator.solve_system("native")
    memory_warnings = [
        r
        for r in caplog.records
        if r.levelno == logging.WARNING and "memory" in r.getMessage().lower()
    ]
    assert len(memory_warnings) == 1
    assert len(interpolator.constraints) > 0
    scale = max(float(np.ptp(expected)), 1e-12)
    assert np.max(np.abs(_field(interpolator) - expected)) <= 1e-3 * scale


def test_update_after_native_setup_keeps_field():
    """(c) update() after (a) does not raise and the field is the same."""
    pytest.importorskip("loop_interpolation_core")
    interpolator = _make_interpolator("native")
    interpolator.setup(**SETUP_KWARGS)
    interpolator._setup_kwargs = dict(SETUP_KWARGS)
    interpolator.solve_system("native")
    first = _field(interpolator)
    interpolator.up_to_date = False
    interpolator.update()
    assert len(interpolator.constraints) == 0
    assert np.allclose(_field(interpolator), first, rtol=1e-9, atol=1e-9)


def test_import_does_not_import_core():
    """(d) Import of loop_interpolation never imports the core."""
    code = (
        "import sys, loop_interpolation; "
        "bad = [m for m in sys.modules if m.startswith('loop_interpolation_core')]; "
        "sys.exit(1 if bad else 0)"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True)
    assert result.returncode == 0, result.stderr.decode()


@pytest.mark.parametrize("solver", ["cg", "lsmr"])
def test_existing_solvers_unchanged_without_core(solver):
    """(d) cg and lsmr give the same field, with no core needed."""
    reference = _make_interpolator("cg")
    reference.setup(**SETUP_KWARGS)
    reference.solve_system("cg", tol=1e-12)
    expected = _field(reference)

    interpolator = _make_interpolator(solver)
    interpolator.setup(**SETUP_KWARGS)
    assert len(interpolator.constraints) > 0
    interpolator.solve_system(solver)
    field = _field(interpolator)
    assert np.all(np.isfinite(field))
    scale = max(float(np.ptp(expected)), 1e-12)
    assert np.max(np.abs(field - expected)) <= 1e-2 * scale


def test_admm_solver_unchanged_without_core():
    """(d) The admm solver runs and gives the same field as cg, with no core.

    The ADMM path needs inequality rows, so the test calls the solver
    on a small box-bounded problem and compares with a direct solve.
    """
    from loop_interpolation.loopsolver import admm_solve
    from scipy import sparse

    rng = np.random.default_rng(0)
    A = sparse.csr_matrix(rng.normal(size=(30, 8)))
    x_true = np.linspace(0.1, 0.8, 8)
    b = A @ x_true
    Q = sparse.identity(8, format="csr")
    bounds = np.array([[-10.0, 10.0]] * 8)
    x = admm_solve(
        A,
        b,
        Q,
        bounds,
        np.zeros(8),
        admm_weight=0.05,
        nmajor=200,
        linsys_solver="cg",
    )
    assert np.all(np.isfinite(x))
    assert np.allclose(x, x_true, atol=1e-2)
