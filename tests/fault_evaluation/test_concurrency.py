"""Tests for simultaneous evaluation calls on models (D4, D5, FR-010a, T022, T047).

The plugin evaluates in background threads. These tests define the behaviour of
evaluation calls that run at the same time.

* Two threads that evaluate one model give the results of serial calls.
* Two calls on one model run at the same time. No lock of the model serializes them.
* The scope of one thread is not visible in another thread.
* An exception in one thread releases the scope of that thread.

No test uses a sleep to synchronize. Events and barriers do the synchronization. Every
join has a time limit, so a deadlock fails the test.
"""

from __future__ import annotations

import threading
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from LoopStructural import GeologicalModel
from LoopStructural.modelling.features import _evaluation
from LoopStructural.modelling.features.fault._fault_segment import FaultSegment
from tests.fault_evaluation.models import build_model

DATA = Path(__file__).parent / "data"
JOIN_SECONDS = 60


def _points(name: str) -> np.ndarray:
    """Return the stored golden points of a model."""
    with np.load(DATA / f"{name}.npz", allow_pickle=False) as archive:
        return archive["points"].copy()


def _group(model: GeologicalModel, index: int):
    return model.get_feature_by_name(f"group_{index}")


class _Worker(threading.Thread):
    """Thread that stores the result or the exception of its function."""

    def __init__(self, function: Callable[[], Any], name: str = "worker") -> None:
        super().__init__(name=name, daemon=True)
        self._function = function
        self.result: Any = None
        self.error: BaseException | None = None

    def run(self) -> None:
        try:
            self.result = self._function()
        except BaseException as error:  # noqa: BLE001 - the test inspects every error
            self.error = error

    def finish(self) -> Any:
        """Join with a time limit. Fail on a deadlock or an exception. Return the result."""
        self.join(JOIN_SECONDS)
        assert not self.is_alive(), f"thread {self.name!r} did not finish: deadlock"
        assert self.error is None, f"thread {self.name!r} raised {self.error!r}"
        return self.result


def _join_all(*workers: _Worker) -> None:
    for worker in workers:
        worker.finish()


def _hook_fault_calls(monkeypatch, before: Callable[[], None]) -> None:
    """Call ``before()`` at the start of each fault call, then run the real call."""
    original = FaultSegment.apply_to_points

    def apply_to_points(fault, points, *args, **kwargs):
        before()
        return original(fault, points, *args, **kwargs)

    monkeypatch.setattr(FaultSegment, "apply_to_points", apply_to_points)


# ---------------------------------------------------------------------------
# Simultaneous calls on one model
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name,group_index", [("f1_u0", 0), ("f4_u3", 3)])
def test_two_threads_on_one_model_match_serial_references(name, group_index):
    """evaluate_model in one thread and evaluate_value in another give the serial results."""
    model = build_model(name)
    points = _points(name)
    group = _group(model, group_index)
    # Make the serial references first.
    model_reference = model.evaluate_model(points)
    value_reference = group.evaluate_value(points)
    for _ in range(3):
        barrier = threading.Barrier(2, timeout=JOIN_SECONDS)

        def run_model(barrier=barrier):
            barrier.wait()
            return model.evaluate_model(points)

        def run_value(barrier=barrier):
            barrier.wait()
            return group.evaluate_value(points)

        first, second = _Worker(run_model, "model"), _Worker(run_value, "value")
        first.start()
        second.start()
        _join_all(first, second)
        np.testing.assert_array_equal(first.result, model_reference)
        np.testing.assert_array_equal(second.result, value_reference)


def test_two_threads_with_the_same_call_match_the_serial_reference():
    """Two equal calls at the same time both give the serial result."""
    model = build_model("f4_u3")
    points = _points("f4_u3")
    reference = model.evaluate_model(points)
    barrier = threading.Barrier(2, timeout=JOIN_SECONDS)

    def run():
        barrier.wait()
        return model.evaluate_model(points.copy())

    workers = [_Worker(run, f"same-{index}") for index in range(2)]
    for worker in workers:
        worker.start()
    _join_all(*workers)
    for worker in workers:
        np.testing.assert_array_equal(worker.result, reference)


def test_two_calls_on_one_model_run_at_the_same_time(monkeypatch):
    """Two threads are inside a fault call of one model at the same time.

    Each thread waits at a barrier for the other thread, inside its first fault call.
    A lock of the model would keep the second thread out, and the barrier would break.
    """
    points = _points("f1_u0")
    reference = build_model("f1_u0").evaluate_model(points)  # a second model: no warm state
    model = build_model("f1_u0")
    barrier = threading.Barrier(2, timeout=10)
    waited: set[int] = set()
    guard = threading.Lock()
    original = FaultSegment.apply_to_points

    def apply_to_points(fault, pts, *args, **kwargs):
        with guard:
            first = threading.get_ident() not in waited
            waited.add(threading.get_ident())
        if first:
            barrier.wait()  # raises threading.BrokenBarrierError when no second thread comes
        return original(fault, pts, *args, **kwargs)

    monkeypatch.setattr(FaultSegment, "apply_to_points", apply_to_points)
    workers = [_Worker(lambda: model.evaluate_model(points.copy()), f"same-{i}") for i in range(2)]
    for worker in workers:
        worker.start()
    _join_all(*workers)
    assert len(waited) == 2, "both threads must reach the fault call"
    for worker in workers:
        np.testing.assert_array_equal(worker.result, reference)


# ---------------------------------------------------------------------------
# Two models
# ---------------------------------------------------------------------------


def test_evaluation_in_progress_on_one_model_does_not_block_another_model(monkeypatch):
    """Model B evaluates to the end while a thread is inside an evaluation of model A."""
    first_model, second_model = build_model("f1_u0"), build_model("f1_u0")
    points = _points("f1_u0")
    reference = second_model.evaluate_model(points)
    inside, release = threading.Event(), threading.Event()

    def before():
        if threading.current_thread().name == "model-a":
            inside.set()
            assert release.wait(JOIN_SECONDS)

    _hook_fault_calls(monkeypatch, before)
    first = _Worker(lambda: first_model.evaluate_model(points), "model-a")
    first.start()
    assert inside.wait(JOIN_SECONDS)
    second = _Worker(lambda: second_model.evaluate_model(points), "model-b")
    try:
        second.start()
        second.finish()
    finally:
        release.set()
    first.finish()
    np.testing.assert_array_equal(second.result, reference)
    np.testing.assert_array_equal(first.result, reference)


# ---------------------------------------------------------------------------
# Isolation of the scope
# ---------------------------------------------------------------------------


def test_scopes_of_two_threads_are_separate(monkeypatch):
    """Each thread has its own context. No context is seen in two threads."""
    model = build_model("f4_u3")
    points = _points("f4_u3")
    seen: dict[str, list] = {}
    original = _evaluation.EvaluationContext.restore
    guard = threading.Lock()

    def recording(context, feature, pts, reverse=False):
        with guard:
            seen.setdefault(threading.current_thread().name, []).append(context)
        return original(context, feature, pts, reverse)

    monkeypatch.setattr(_evaluation.EvaluationContext, "restore", recording)
    barrier = threading.Barrier(2, timeout=JOIN_SECONDS)

    def run():
        assert _evaluation._ACTIVE.get() is None
        barrier.wait()
        result = model.evaluate_model(points)
        assert _evaluation._ACTIVE.get() is None, "the scope stayed active in the thread"
        return result

    workers = [_Worker(run, f"scope-{index}") for index in range(2)]
    for worker in workers:
        worker.start()
    _join_all(*workers)
    assert _evaluation._ACTIVE.get() is None
    assert set(seen) == {"scope-0", "scope-1"}
    # The lists hold the context objects, so no id can be reused.
    first_ids = {id(context) for context in seen["scope-0"]}
    second_ids = {id(context) for context in seen["scope-1"]}
    assert first_ids.isdisjoint(second_ids), "a context was shared between threads"
    np.testing.assert_array_equal(workers[0].result, workers[1].result)


def test_scope_of_one_thread_is_not_visible_in_another_thread():
    """A scope that is open in one thread is not the active scope of another thread."""
    model = build_model("f1_u0")
    inside, release = threading.Event(), threading.Event()

    def hold():
        with _evaluation.evaluation_scope(model):
            inside.set()
            assert release.wait(JOIN_SECONDS)

    holder = _Worker(hold, "scope-holder")
    holder.start()
    assert inside.wait(JOIN_SECONDS)
    other = _Worker(lambda: _evaluation._ACTIVE.get(), "scope-reader")
    try:
        other.start()
        other.finish()
    finally:
        release.set()
    holder.finish()
    assert other.result is None
    assert _evaluation._ACTIVE.get() is None


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class _Boom(RuntimeError):
    """Error that the test raises on purpose."""


def test_exception_in_one_thread_releases_the_scope(monkeypatch):
    """A failing thread leaves no scope. Other threads and later calls work."""
    model = build_model("f4_u3")
    points = _points("f4_u3")
    reference = build_model("f4_u3").evaluate_model(points)  # a second model: no warm state
    inside, release = threading.Event(), threading.Event()

    def before():
        if threading.current_thread().name == "failing":
            inside.set()
            assert release.wait(JOIN_SECONDS)
            raise _Boom("on purpose")

    _hook_fault_calls(monkeypatch, before)
    scope_after_failure: list[Any] = []

    def failing():
        try:
            return model.evaluate_model(points)
        finally:
            scope_after_failure.append(_evaluation._ACTIVE.get())

    bad = _Worker(failing, "failing")
    bad.start()
    assert inside.wait(JOIN_SECONDS)
    good = _Worker(lambda: model.evaluate_model(points), "good")
    good.start()
    release.set()
    bad.join(JOIN_SECONDS)
    assert not bad.is_alive(), "the failing thread did not finish"
    assert isinstance(bad.error, _Boom)
    assert scope_after_failure == [None], "the scope was not reset after the exception"
    good.finish()
    np.testing.assert_array_equal(good.result, reference)
    assert _evaluation._ACTIVE.get() is None
    again = _Worker(lambda: model.evaluate_model(points), "again")
    again.start()
    np.testing.assert_array_equal(again.finish(), reference)
