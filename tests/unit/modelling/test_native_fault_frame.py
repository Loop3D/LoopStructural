"""A fault frame with solver="native": all three coordinates use the native core.

The coordinates 1 and 2 have gradient orthogonal constraints. The native core
gets them as raw data, so it solves the full frame.
The tests skip while the native core (loop_interpolation_core) is not installed.
"""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("loop_interpolation_core")

from LoopStructural import GeologicalModel


def _fault_model(solver):
    data = pd.DataFrame(
        {
            "X": [500.0, 500.0, 500.0, 500.0],
            "Y": [100.0, 500.0, 900.0, 500.0],
            "Z": [-500.0, -500.0, -500.0, -500.0],
            "val": [0.0, 0.0, 0.0, np.nan],
            "gx": [np.nan, np.nan, np.nan, 1.0],
            "gy": [np.nan, np.nan, np.nan, 0.0],
            "gz": [np.nan, np.nan, np.nan, 0.0],
            "feature_name": "f1",
        }
    )
    model = GeologicalModel([0, 0, -1000], [1000, 1000, 0])
    model.data = data
    fault = model.create_and_add_fault("f1", 50, nelements=20_000, solver=solver)
    fault.builder.up_to_date()
    return fault


def test_all_coordinates_of_a_fault_use_the_native_solver():
    fault = _fault_model("native")
    for coordinate in range(3):
        builder = fault[coordinate].builder
        assert builder.build_arguments["solver"] == "native"
        assert builder.interpolator.latest_native_report["converged"]


def test_the_native_fault_frame_is_the_cg_fault_frame():
    native = _fault_model("native")
    matrix = _fault_model("cg")
    points = np.random.default_rng(0).random((200, 3)) * [1000, 1000, 1000] + [0, 0, -1000]
    for coordinate in range(3):
        a = native[coordinate].evaluate_value(points)
        b = matrix[coordinate].evaluate_value(points)
        scale = np.nanmax(np.abs(b))
        assert np.nanmax(np.abs(a - b)) <= 1e-3 * scale, coordinate
