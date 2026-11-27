"""Tests that the native fault path is safe to use (story 5, scenarios 2 to 7, edge cases).

The tests need the installed package ``loop_model_core``. They use the model ``f4_u0``.
"""

from __future__ import annotations

import builtins
import threading
from copy import deepcopy

import numpy as np
import pytest

pytest.importorskip("loop_model_core")

from loop_model_core import FaultChain

from LoopStructural import LoopStructuralConfig
from LoopStructural.modelling.features import _native_faults
from tests.fault_evaluation.golden import _faulted_features
from tests.fault_evaluation.models import build_model, sample_points

NAME = "f4_u0"
POINT_COUNT = 2000
BIG_COUNT = 10_000_000


@pytest.fixture(autouse=True)
def settings(monkeypatch):
    """Reset the native settings and the environment variable for each test."""
    monkeypatch.setattr(LoopStructuralConfig, "native_faults", None)
    monkeypatch.setattr(LoopStructuralConfig, "native_threads", None)
    monkeypatch.delenv("LOOPSTRUCTURAL_NATIVE_FAULTS", raising=False)


@pytest.fixture(scope="module")
def shared():
    """Return the model, a faulted feature and sample points. Tests must not change them."""
    model = build_model(NAME)
    feature = _faulted_features(model)[0][0]
    return model, feature, sample_points(NAME, POINT_COUNT)


@pytest.fixture
def fresh():
    """Return a new model, a faulted feature and sample points, for tests that change state."""
    model = build_model(NAME)
    feature = _faulted_features(model)[0][0]
    return model, feature, sample_points(NAME, POINT_COUNT)


def _tolerance(model) -> float:
    # Box size rule of tests/fault_evaluation/regenerate.py: the longest side of the box.
    box = model.bounding_box
    return 1e-8 * float(np.max(np.asarray(box.maximum) - np.asarray(box.origin)))


def _python_values(feature, pts):
    LoopStructuralConfig.native_faults = False
    try:
        return feature.evaluate_value(pts)
    finally:
        LoopStructuralConfig.native_faults = None


def _python_restored(feature, pts):
    LoopStructuralConfig.native_faults = False
    try:
        return feature._apply_faults(pts)
    finally:
        LoopStructuralConfig.native_faults = None


def test_native_path_is_used(shared):
    _, feature, pts = shared
    assert _native_faults.restore(feature, pts) is not None
    feature.evaluate_value(pts)
    assert hasattr(feature, "_native_chain")


def test_plain_python_function_uses_python_path(fresh):
    model, feature, pts = fresh
    # The setter accepts any callable, so the public attribute is used.
    feature.faults[0].faultfunction = lambda gx, gy, gz: np.ones_like(gx)
    assert _native_faults.chain_for(feature) is None
    assert _native_faults.restore(feature, pts) is None
    expected = _python_values(feature, pts)
    np.testing.assert_allclose(
        feature.evaluate_value(pts), expected, atol=_tolerance(model), equal_nan=True
    )


def test_second_call_reads_no_file_and_keeps_chain(fresh, monkeypatch):
    _, feature, pts = fresh
    feature.evaluate_value(pts)
    first_chain = feature._native_chain[2]

    def fail(*args, **kwargs):
        raise AssertionError("a file was opened")

    with monkeypatch.context() as patch:
        patch.setattr(builtins, "open", fail)
        assert _native_faults.restore(feature, pts) is not None
        feature.evaluate_value(pts)
    assert feature._native_chain[2] is first_chain


def test_changed_displacement_is_used(fresh):
    model, feature, pts = fresh
    assert _native_faults.restore(feature, pts) is not None
    first_chain = feature._native_chain[2]
    feature.faults[0].displacement = feature.faults[0].displacement * 1.5 + 1.0
    native = feature.evaluate_value(pts)
    assert feature._native_chain[2] is not first_chain
    np.testing.assert_allclose(
        native, _python_values(feature, pts), atol=_tolerance(model), equal_nan=True
    )


def test_gil_is_released_during_restore(shared):
    _, feature, pts = shared
    big = np.ascontiguousarray(np.tile(pts, (BIG_COUNT // len(pts), 1)), dtype=np.float64)
    counter = [0]
    stop = threading.Event()

    def count():
        while not stop.is_set():
            counter[0] += 1

    worker = threading.Thread(target=count)
    worker.start()
    try:
        before = counter[0]
        result = _native_faults.restore(feature, big)
        after = counter[0]
    finally:
        stop.set()
        worker.join()
    assert result is not None
    assert after - before >= 1000


def test_two_threads_equal_single_thread(fresh):
    _, feature, pts = fresh
    assert _native_faults.restore(feature, pts) is not None
    other = pts[::-1].copy()
    alone = [feature.evaluate_value(pts), feature.evaluate_value(other)]
    results = [None, None]

    def run(index, points):
        for _ in range(5):
            results[index] = feature.evaluate_value(points)

    threads = [threading.Thread(target=run, args=(i, p)) for i, p in enumerate((pts, other))]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    for got, want in zip(results, alone):
        assert np.array_equal(got, want, equal_nan=True)


def test_one_thread_equals_all_threads(shared):
    _, feature, pts = shared
    LoopStructuralConfig.native_threads = 0
    all_threads = _native_faults.restore(feature, pts)
    LoopStructuralConfig.native_threads = 1
    one_thread = _native_faults.restore(feature, pts)
    assert all_threads is not None
    assert np.array_equal(all_threads, one_thread, equal_nan=True)


def test_empty_points(shared):
    _, feature, _ = shared
    result = _native_faults.restore(feature, np.zeros((0, 3)))
    assert result is not None
    assert result.shape == (0, 3)


def test_nan_coordinate(shared):
    model, feature, pts = shared
    pts = pts.copy()
    pts[3, 1] = np.nan
    result = _native_faults.restore(feature, pts)
    assert result is not None
    # The Python path keeps a row with a NaN coordinate as it is. The native path must do the same.
    expected = _python_restored(feature, pts)
    np.testing.assert_array_equal(np.isnan(result), np.isnan(expected))
    np.testing.assert_allclose(result, expected, atol=_tolerance(model))


def test_bad_spec_raises_and_state_is_kept(fresh):
    _, feature, pts = fresh
    before = deepcopy(_native_faults.restore(feature, pts))
    assert before is not None
    first_chain = feature._native_chain[2]
    with pytest.raises(ValueError):
        FaultChain({"version": 99})
    after = _native_faults.restore(feature, pts)
    assert np.array_equal(before, after, equal_nan=True)
    assert feature._native_chain[2] is first_chain
