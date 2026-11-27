"""Tests for the restore order that a feature keeps between evaluations (FR-001, FR-002).

``restore_order(faults, reverse=..., owner=...)`` gives the order in which faults restore
points. With an owner, it computes the order once for each state of the faults and keeps it.
"""

from __future__ import annotations

from types import SimpleNamespace
from unittest import mock

import numpy as np
import pytest

from LoopStructural import LoopStructuralConfig
from LoopStructural.utils import _fault_block
from LoopStructural.utils._fault_block import cutters_first, restore_order, with_cutters
from tests.fault_evaluation.golden import capture
from tests.fault_evaluation.models import build_model, sample_points


def _fake_fault(name, cutters=()):
    """Make a fault-like object whose cutters are read from its builder."""
    builder = SimpleNamespace(builders=[SimpleNamespace(faults=list(cutters))])
    return SimpleNamespace(name=name, builder=builder, features=[])


def _faulted_feature(model):
    return next(f for f in model.features if getattr(f, "faults", None))


def test_order_is_computed_once_for_one_state():
    model = build_model("cutting_faults")
    feature = _faulted_feature(model)
    points = sample_points("cutting_faults", 50)
    feature.evaluate_value(points)
    with mock.patch.object(_fault_block, "cutters_first", wraps=cutters_first) as spy:
        feature.evaluate_value(points)
    assert spy.call_count == 0


def test_order_changes_when_a_cutter_is_added():
    other = _fake_fault("other")
    fault = _fake_fault("fault")
    owner = SimpleNamespace()
    first = restore_order([fault], owner=owner)
    assert first == [fault]
    assert restore_order([fault], owner=owner) is first
    fault.builder.builders[0].faults.append(other)
    second = restore_order([fault], owner=owner)
    assert second != first
    assert second.index(other) < second.index(fault)


@pytest.mark.parametrize("name", ["cutting_faults", "f4_u3"])
def test_cache_gives_identical_values(name, monkeypatch):
    # The test checks the Python path bit for bit; the native path gives equal results within the tolerance only.
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", False)
    model = build_model(name)
    points = sample_points(name, 200, seed=1)
    cached = capture(model, points)
    original = restore_order
    with mock.patch.object(
        _fault_block, "restore_order", lambda *a, **k: original(*a, **{**k, "owner": None})
    ):
        plain = capture(build_model(name), points)
    assert cached.keys() == plain.keys()
    for key in cached:
        np.testing.assert_array_equal(cached[key], plain[key])


def test_reverse_order_rule():
    a = _fake_fault("a")
    b = _fake_fault("b", cutters=[a])
    c = _fake_fault("c")
    faults = [b, c]
    expected = cutters_first(with_cutters(faults)[::-1])[::-1]
    assert restore_order(faults, reverse=True) == expected
    assert restore_order(faults, reverse=True, owner=SimpleNamespace()) == expected
    assert restore_order(faults) == cutters_first(with_cutters(faults))
    assert expected != restore_order(faults)


@pytest.mark.parametrize("name", ["cutting_faults", "f4_u3"])
def test_shared_frame_restoration_gives_identical_values(name, monkeypatch):
    # The test checks the Python path bit for bit; the native path gives equal results within the tolerance only.
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", False)
    from LoopStructural.modelling.features.fault import _fault_segment

    points = sample_points(name, 200, seed=1)
    shared = capture(build_model(name), points)
    with mock.patch.object(_fault_segment, "_SHARE_CUT_FRAME_RESTORATION", False):
        legacy = capture(build_model(name), points)
    assert shared.keys() == legacy.keys()
    for key in shared:
        np.testing.assert_array_equal(shared[key], legacy[key])
