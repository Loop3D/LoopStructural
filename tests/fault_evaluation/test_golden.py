"""Replay all frozen pre-performance fault-evaluation references."""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pytest

from tests.fault_evaluation.golden import assert_matches, capture
from tests.fault_evaluation.models import build_model, model_names

DATA = Path(__file__).parent / "data"
MANIFEST = json.loads((DATA / "manifest.json").read_text(encoding="utf8"))


def _reference(name: str) -> dict[str, np.ndarray]:
    path = DATA / f"{name}.npz"
    metadata = MANIFEST["models"][name]
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == metadata["sha256"], f"{name} frozen reference hash changed"
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in archive.files}


def test_frozen_model_catalog_matches_manifest():
    """Keep the complete frozen model catalogue synchronized with its manifest."""
    names = model_names()
    assert len(names) == 21
    assert set(names) == set(MANIFEST["models"])


@pytest.mark.parametrize("name", model_names())
def test_frozen_model_matches_legacy_output(name, record_property):
    """Build and replay one frozen case without invoking a solver or writer."""
    started = time.perf_counter()
    expected = _reference(name)
    model = build_model(name)
    actual = capture(model, expected["points"])
    assert_matches(name, actual, expected, MANIFEST["models"][name]["box_size"])
    record_property(f"legacy_golden_capture_seconds_{name}", time.perf_counter() - started)


def test_one_percent_fault_movement_change_is_rejected(monkeypatch):
    """The comparator detects a small scoped change in actual fault movement."""
    from LoopStructural.modelling.features.fault import FaultSegment

    name = "f1_u0"
    expected = _reference(name)
    model = build_model(name)
    original = FaultSegment.apply_to_points

    def perturbed(self, points, *args, **kwargs):
        points = np.asarray(points)
        moved = original(self, points, *args, **kwargs)
        return points + 1.01 * (moved - points)

    with monkeypatch.context() as scoped:
        scoped.setattr(FaultSegment, "apply_to_points", perturbed)
        actual = capture(model, expected["points"])

    with pytest.raises(AssertionError) as error:
        assert_matches(name, actual, expected, MANIFEST["models"][name]["box_size"])
    message = str(error.value)
    assert "f1_u0" in message
    assert "max difference=" in message
    assert "at point" in message
