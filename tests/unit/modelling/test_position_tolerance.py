"""The solver tolerance in length units of the model (LS-F1)."""

import numpy as np
import pandas as pd
import pytest

from LoopStructural import GeologicalModel
from tests.fixtures.fault_models import COLUMNS

EXTENT = 2000.0  # metres


def _model(model_tolerance=None, **feature_kwargs):
    """A noisy dipping series in a 2 km box that starts at a large offset."""
    origin = np.array([500000.0, 7000000.0, -1000.0])
    model = GeologicalModel(origin, origin + EXTENT)
    if model_tolerance is not None:
        model.position_tolerance = model_tolerance
    rng = np.random.default_rng(11)
    local = rng.uniform(0.0, EXTENT, size=(300, 3))
    value = local[:, 2] + 0.2 * local[:, 0] + 60.0 * np.sin(local[:, 1] / 300.0)
    points = local + origin
    df = pd.DataFrame(
        {
            "X": points[:, 0], "Y": points[:, 1], "Z": points[:, 2], "val": value,
            "nx": np.nan, "ny": np.nan, "nz": np.nan, "coord": np.nan, "feature_name": "strat",
        }
    )[COLUMNS]
    model.create_and_add_foliation("strat", data=df, nelements=8000, buffer=0.2,
                                   **feature_kwargs)
    model.update()
    return model


def _interpolator(model):
    return model["strat"].builder.interpolator


@pytest.fixture(scope="module")
def reference():
    return _interpolator(_model(solver_kwargs={"rtol": 1e-13, "atol": 0.0}))


def _position_error(interpolator, reference):
    nodes = reference.support.nodes
    gradient = np.linalg.norm(reference.support.evaluate_gradient(nodes, reference.c), axis=1)
    return np.max(np.abs(interpolator.c - reference.c) / gradient)


def test_default_is_none_and_keeps_the_residual_stop():
    model = _model()
    assert model.position_tolerance is None
    assert "cg_stop" not in _interpolator(model).get_last_solve_timing()


@pytest.mark.parametrize("metres", [20.0, 0.2])
def test_model_position_tolerance_is_in_model_units(reference, metres):
    interpolator = _interpolator(_model(metres))
    timing = interpolator.get_last_solve_timing()
    assert timing["cg_stop"] == "position"
    assert timing["position_tolerance"] == metres
    assert _position_error(interpolator, reference) <= metres


def test_a_tighter_tolerance_is_closer_to_the_converged_field(reference):
    loose = _position_error(_interpolator(_model(20.0)), reference)
    tight = _position_error(_interpolator(_model(0.02)), reference)
    assert tight < loose
    assert tight <= 0.02


def test_feature_keyword_has_priority_over_the_model():
    model = _model(20.0)
    assert _interpolator(model).get_last_solve_timing()["position_tolerance"] == 20.0
    model = _model(20.0, position_tolerance=2.0)
    assert _interpolator(model).get_last_solve_timing()["position_tolerance"] == 2.0
    model = _model(position_tolerance=2.0)
    assert _interpolator(model).get_last_solve_timing()["position_tolerance"] == 2.0


def test_two_builds_are_bit_equal():
    a = _interpolator(_model(1.0))
    b = _interpolator(_model(1.0))
    np.testing.assert_array_equal(a.c, b.c)
