import logging

import numpy as np
from scipy import sparse

from loop_interpolation import _solver_pipeline as pipeline


class _DummyScaling:
    def __init__(self, arr):
        self.arr = arr

    def __rmatmul__(self, other):
        return other @ self.arr


def test_extract_constant_norm_options_invalid_target_disables_target():
    logger = logging.getLogger("test_solver_pipeline")
    kwargs = {
        "constant_norm_iterations": 3,
        "constant_norm_weight": 0.2,
        "constant_norm_target": -1.0,
    }
    iters, weight, target = pipeline.extract_constant_norm_options(kwargs, logger)
    assert iters == 3
    assert weight == 0.2
    assert target is None
    assert kwargs == {}


def test_preprocess_main_system_adds_ridge_and_scaling():
    logger = logging.getLogger("test_solver_pipeline")
    A = sparse.eye(2, format="csr")
    b = np.array([1.0, 2.0])

    def fake_scaling(mat):
        return sparse.diags([2.0, 3.0])

    timing = {}
    A2, b2, S = pipeline.preprocess_main_system(
        A=A,
        b=b,
        add_ridge_regularisation=True,
        ridge_factor=1e-8,
        apply_scaling_matrix=True,
        compute_column_scaling_matrix_fn=fake_scaling,
        logger=logger,
        timing=timing,
    )

    assert A2.shape[0] == 4
    assert A2.shape[1] == 2
    assert b2.shape[0] == 4
    assert S is not None
    assert "preprocess_seconds" in timing


def test_assemble_inequality_system_records_timing():
    def fake_build_ineq():
        return sparse.csr_matrix((3, 2), dtype=float), np.zeros((3, 3), dtype=float)

    timing = {}
    Q, bounds = pipeline.assemble_inequality_system(fake_build_ineq, timing)
    assert Q.shape == (3, 2)
    assert bounds.shape == (3, 3)
    assert timing["inequality_rows"] == 3
    assert "inequality_seconds" in timing


def test_finalize_timing_sets_total_and_status():
    timing = {"solver": "cg"}
    out = pipeline.finalize_timing(timing=timing, solve_started=0.0, up_to_date=True)
    assert out["up_to_date"] is True
    assert "total_seconds" in out
