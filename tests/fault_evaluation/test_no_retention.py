"""Tests that no state stays between evaluation calls (FR-009, user story 2, scenario 6, T047).

After a change of the model, the next call must give the result of a model that has the
same change from the start. The model keeps no result of an earlier call.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest

from benchmarks.fault_evaluation.instrumentation import EvaluationInstrumentation
from LoopStructural import GeologicalModel
from tests.fault_evaluation.golden import assert_matches, capture
from tests.fault_evaluation.models import build_model, sample_points

DATA = Path(__file__).parent / "data"
MANIFEST = json.loads((DATA / "manifest.json").read_text(encoding="utf8"))
POINT_COUNT = 200


def _scale_displacement(model: GeologicalModel) -> None:
    """Multiply the displacement of the first fault by 1.5."""
    model.get_feature_by_name("fault_0").displacement *= 1.5


def _disable_faults(model: GeologicalModel) -> None:
    """Disable the faults of the oldest group."""
    model.get_feature_by_name("group_0").faults_enabled = False


def _set_steps(model: GeologicalModel) -> None:
    """Set the steps of the first fault to 5."""
    model.get_feature_by_name("fault_0").steps = 5


@pytest.mark.parametrize(
    "name,change",
    [
        ("f1_u0", _scale_displacement),
        ("f1_u1", _disable_faults),
        ("f1_u0", _set_steps),
    ],
    ids=["displacement", "faults_enabled", "steps"],
)
def test_change_after_a_call_gives_the_result_of_a_fresh_model(
    name: str, change: Callable[[GeologicalModel], None]
) -> None:
    """A call, a change and a call equal a model that has the change from the start."""
    points = sample_points(name, POINT_COUNT, seed=1)
    model = build_model(name)
    model.evaluate_model(points.copy())
    change(model)
    fresh = build_model(name)
    change(fresh)
    assert_matches(
        name,
        capture(model, points),
        capture(fresh, points),
        MANIFEST["models"][name]["box_size"],
    )


def test_repeat_call_applies_the_faults_again() -> None:
    """Two equal calls make the same number of fault calls and give equal results."""
    model = build_model("f1_u0")
    points = sample_points("f1_u0", POINT_COUNT, seed=1)
    counts = []
    results = []
    for _ in range(2):
        with EvaluationInstrumentation() as instrumentation:
            results.append(model.evaluate_model(points.copy()))
        counts.append(instrumentation.snapshot()["fault_calls_total"])
    assert counts[0] > 0
    assert counts[0] == counts[1], "the second call must apply the faults again"
    np.testing.assert_array_equal(results[0], results[1])


def test_model_has_no_cache_attributes() -> None:
    """The model has no ``evaluation_cache`` and no ``evaluation_revision``."""
    model = build_model("f1_u0")
    assert not hasattr(model, "evaluation_cache")
    assert not hasattr(model, "evaluation_revision")
