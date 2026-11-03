"""Tests for simultaneous evaluation calls on models (D4, D5, FR-010a, T022).

The plugin evaluates in background threads. These tests define the behaviour that T023
must give.

* Two threads that evaluate one model give the results of serial calls.
* A lazy build of a dirty builder happens exactly once, also when two threads start.
* One model serializes its evaluations, its builds and its changes under the lock
  ``model.evaluation_cache.lock``. Different models do not share a lock.
* The scope of one thread is not visible in another thread.
* An exception in one thread releases the lock and the scope, and it publishes nothing.
* A change through a public setter, while another thread evaluates, waits for the
  evaluation. No call gives a mixed result.

No test uses a sleep to synchronize. Events and barriers do the synchronization. Some
tests wait for an event that must NOT be set. That wait is short and it cannot fail a
correct implementation. Every join has a time limit, so a deadlock fails the test.
"""

from __future__ import annotations

import threading
from collections.abc import Callable
from typing import Any

import numpy as np
import pytest

from LoopStructural import GeologicalModel
from LoopStructural.modelling.features import _evaluation
from LoopStructural.modelling.features.builders._geological_feature_builder import (
    GeologicalFeatureBuilder,
)
from LoopStructural.modelling.features.fault._fault_segment import FaultSegment
from tests.fault_evaluation.models import build_model
from tests.fault_evaluation.test_invalidation import (
    _evaluate_disabled,
    _fault,
    _group,
    _install_counters,
    _points,
    _tiny_model,
    _tiny_points,
)

JOIN_SECONDS = 60
# Time to wait for an event that a correct implementation must not set.
ABSENT_SECONDS = 0.5


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


def _lock(model: GeologicalModel):
    """Return the lock of the model (the lock of its evaluation cache)."""
    return model.evaluation_cache.lock


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
    # Make the references first. A cache-disabled call clears all entries of the model.
    model_reference, _ = _evaluate_disabled(model, lambda: model.evaluate_model(points))
    value_reference, _ = _evaluate_disabled(model, lambda: group.evaluate_value(points))
    for _ in range(3):
        model.evaluation_cache.clear()
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
    reference, _ = _evaluate_disabled(model, lambda: model.evaluate_model(points))
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


def test_first_lazy_build_happens_once_with_two_threads(monkeypatch):
    """Two threads on a model with dirty builders build each builder exactly once.

    The first build waits until a second build could enter. Without a lock, the second
    thread enters the same builder while it is still dirty, and the build runs twice.
    """
    builds, solves = _install_counters(monkeypatch)
    counted_build = GeologicalFeatureBuilder.build  # the counting wrapper
    state = {"first": True}
    in_build, release, second_entered = threading.Event(), threading.Event(), threading.Event()
    guard = threading.Lock()

    def gated_build(self, *args, **kwargs):
        with guard:
            is_first = state["first"]
            state["first"] = False
        if is_first:
            in_build.set()
            release.wait(JOIN_SECONDS)
        else:
            second_entered.set()
        return counted_build(self, *args, **kwargs)

    monkeypatch.setattr(GeologicalFeatureBuilder, "build", gated_build)
    model = _tiny_model()
    group, fault, points = _group(model, 0), _fault(model, 0), _tiny_points()
    assert group.faults == [], "the model must start with dirty builders"
    first = _Worker(lambda: group.evaluate_value(points), "builder-1")
    second = _Worker(lambda: group.evaluate_value(points), "builder-2")
    first.start()
    assert in_build.wait(JOIN_SECONDS), "the first build did not start"
    second.start()
    # A second build in the other thread while the first one runs is the failure.
    concurrent = second_entered.wait(ABSENT_SECONDS)
    release.set()
    for worker in (first, second):
        worker.join(JOIN_SECONDS)
        assert not worker.is_alive(), f"thread {worker.name!r} did not finish: deadlock"
    # Check the overlap first. A racing build can also raise in a worker.
    assert not concurrent, "a second thread built a builder that is already being built"
    _join_all(first, second)
    assert set(builds.values()) == {1}, f"a builder was built more than once: {builds}"
    expected = {id(group.builder)} | {id(builder) for builder in fault.builder.builders}
    assert set(builds) == expected
    assert solves[0] == len(expected)
    np.testing.assert_array_equal(first.result, second.result)


# ---------------------------------------------------------------------------
# The lock of one model
# ---------------------------------------------------------------------------


def test_held_lock_blocks_evaluation_also_with_the_cache_disabled(monkeypatch):
    """While one thread holds the lock of a model, an evaluation of the model waits.

    The cache is off, so the evaluation uses no persistent entry. The full outer
    evaluation must still take the lock of the model. This test passes before T023 only
    because reading ``cache.enabled`` takes the same lock. It guards the rule after T023.
    """
    model = build_model("f1_u0")
    points = _points("f1_u0")
    reference, _ = _evaluate_disabled(model, lambda: model.evaluate_model(points))
    model.evaluation_cache.enabled = False
    entered = threading.Event()
    _hook_fault_calls(monkeypatch, entered.set)
    with _lock(model):
        worker = _Worker(lambda: model.evaluate_model(points), "blocked")
        worker.start()
        started = entered.wait(ABSENT_SECONDS)
    assert not started, "the evaluation ran while another thread held the lock of the model"
    worker.finish()
    assert entered.is_set()
    np.testing.assert_array_equal(worker.result, reference)


def test_nested_and_repeated_calls_reuse_the_lock_in_one_thread():
    """A thread that holds the lock can evaluate the model. The lock is re-entrant."""
    model = build_model("f1_u0")
    points = _points("f1_u0")
    group = _group(model, 0)
    reference, _ = _evaluate_disabled(model, lambda: model.evaluate_model(points))

    def run():
        with _lock(model), _lock(model):
            return model.evaluate_model(points), group.evaluate_value(points)

    worker = _Worker(run, "reentrant")
    worker.start()
    worker.finish()
    np.testing.assert_array_equal(worker.result[0], reference)


def test_thread_holding_the_lock_of_one_model_does_not_block_another_model():
    """There is no global lock. A held lock of model A leaves model B free."""
    first_model, second_model = build_model("f1_u0"), build_model("f1_u0")
    assert _lock(first_model) is not _lock(second_model)
    points = _points("f1_u0")
    reference, _ = _evaluate_disabled(second_model, lambda: second_model.evaluate_model(points))
    held, release = threading.Event(), threading.Event()

    def hold():
        with _lock(first_model):
            held.set()
            assert release.wait(JOIN_SECONDS)

    holder = _Worker(hold, "holder")
    holder.start()
    assert held.wait(JOIN_SECONDS)
    other = _Worker(lambda: second_model.evaluate_model(points), "other-model")
    try:
        other.start()
        other.finish()
    finally:
        release.set()
    holder.finish()
    np.testing.assert_array_equal(other.result, reference)
    assert first_model.evaluation_cache.retained_bytes == 0
    assert second_model.evaluation_cache.retained_bytes > 0


def test_evaluation_in_progress_on_one_model_does_not_block_another_model(monkeypatch):
    """Model B evaluates to the end while a thread is inside an evaluation of model A."""
    first_model, second_model = build_model("f1_u0"), build_model("f1_u0")
    points = _points("f1_u0")
    reference, _ = _evaluate_disabled(second_model, lambda: second_model.evaluate_model(points))
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
    # Each model keeps its own entries.
    first_cache, second_cache = first_model.evaluation_cache, second_model.evaluation_cache
    assert first_cache is not second_cache
    assert first_cache.retained_bytes == second_cache.retained_bytes == points.nbytes


# ---------------------------------------------------------------------------
# Isolation of the scope
# ---------------------------------------------------------------------------


def test_scopes_of_two_threads_are_separate(monkeypatch):
    """Each thread has its own context. No context is seen in two threads."""
    model = build_model("f4_u3")
    points = _points("f4_u3")
    seen: dict[str, set[int]] = {}
    original = _evaluation.EvaluationContext.restore
    guard = threading.Lock()

    def recording(context, feature, pts, reverse=False):
        with guard:
            seen.setdefault(threading.current_thread().name, set()).add(id(context))
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
    assert seen["scope-0"].isdisjoint(seen["scope-1"]), "a context was shared between threads"
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


def test_exception_in_one_thread_releases_lock_and_scope(monkeypatch):
    """A failing thread leaves no lock and no scope. Others and later calls work."""
    model = build_model("f4_u3")
    points = _points("f4_u3")
    reference, _ = _evaluate_disabled(model, lambda: model.evaluate_model(points))
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
    # The lock is free for a new thread. A lock that a dead thread still owns would block.
    acquired: list[bool] = []

    def try_lock():
        got = _lock(model).acquire(timeout=JOIN_SECONDS)
        acquired.append(got)
        if got:
            _lock(model).release()

    prober = _Worker(try_lock, "prober")
    prober.start()
    prober.finish()
    assert acquired == [True], "the lock of the model stayed held after the exception"
    assert _evaluation._ACTIVE.get() is None
    again = _Worker(lambda: model.evaluate_model(points), "again")
    again.start()
    np.testing.assert_array_equal(again.finish(), reference)


def test_failed_call_publishes_no_entry(monkeypatch):
    """An exception in the middle of a call leaves no entry of that call in the cache."""
    model = build_model("f4_u3")
    points = _points("f4_u3")
    reference, _ = _evaluate_disabled(model, lambda: model.evaluate_model(points))
    cache = model.evaluation_cache
    cache.clear()
    count = [0]
    guard = threading.Lock()

    def before():
        with guard:
            count[0] += 1
            number = count[0]
        if number == 3:
            raise _Boom("on purpose")

    _hook_fault_calls(monkeypatch, before)
    worker = _Worker(lambda: model.evaluate_model(points), "partial")
    worker.start()
    worker.join(JOIN_SECONDS)
    assert not worker.is_alive()
    assert isinstance(worker.error, _Boom)
    assert cache.retained_bytes == 0, "a failed call published a restoration"
    monkeypatch.undo()
    after = _Worker(lambda: model.evaluate_model(points), "after-failure")
    after.start()
    np.testing.assert_array_equal(after.finish(), reference)


# ---------------------------------------------------------------------------
# Changes while another thread evaluates
# ---------------------------------------------------------------------------


def _set_displacement(model):
    _fault(model, 0).displacement = 20.0


def _set_steps(model):
    _fault(model, 2).steps = 2


@pytest.mark.parametrize("change", [_set_displacement, _set_steps], ids=["displacement", "steps"])
def test_change_during_evaluation_waits_for_the_evaluation(monkeypatch, change):
    """A public setter waits for the evaluation in progress. The result is not mixed.

    The evaluating thread stops inside its first fault call and holds the lock of the
    model. A second thread then changes the model. The change must not finish before the
    evaluation ends. The result of the evaluation is the old reference. A later call
    gives the new reference.
    """
    model = build_model("f4_u3")
    points = _points("f4_u3")
    old_reference, _ = _evaluate_disabled(model, lambda: model.evaluate_model(points))
    inside, release = threading.Event(), threading.Event()
    changed = threading.Event()

    def before():
        if threading.current_thread().name == "evaluating":
            inside.set()
            assert release.wait(JOIN_SECONDS)

    _hook_fault_calls(monkeypatch, before)

    def mutate():
        change(model)
        changed.set()

    evaluating = _Worker(lambda: model.evaluate_model(points), "evaluating")
    evaluating.start()
    assert inside.wait(JOIN_SECONDS)
    mutator = _Worker(mutate, "mutator")
    mutator.start()
    finished_early = changed.wait(ABSENT_SECONDS)
    release.set()
    _join_all(evaluating, mutator)
    assert not finished_early, "the change did not wait for the evaluation in progress"
    # No torn result: the evaluation equals the old state.
    np.testing.assert_array_equal(evaluating.result, old_reference)
    monkeypatch.undo()
    new_reference, _ = _evaluate_disabled(model, lambda: model.evaluate_model(points))
    assert not np.array_equal(new_reference, old_reference, equal_nan=True), "no effect"
    later = _Worker(lambda: model.evaluate_model(points), "later")
    later.start()
    np.testing.assert_array_equal(later.finish(), new_reference)


def test_repeated_changes_and_evaluations_never_give_a_mixed_result():
    """A changing thread and an evaluating thread: every result equals one known state."""
    model = build_model("f4_u3")
    points = _points("f4_u3")
    fault = _fault(model, 0)
    old, new = fault.displacement, 20.0
    old_reference, _ = _evaluate_disabled(model, lambda: model.evaluate_model(points))
    fault.displacement = new
    new_reference, _ = _evaluate_disabled(model, lambda: model.evaluate_model(points))
    fault.displacement = old
    assert not np.array_equal(old_reference, new_reference, equal_nan=True)
    barrier = threading.Barrier(2, timeout=JOIN_SECONDS)

    def flip():
        barrier.wait()
        for index in range(6):
            fault.displacement = new if index % 2 == 0 else old

    def evaluate():
        barrier.wait()
        return [model.evaluate_model(points) for _ in range(6)]

    flipper, evaluator = _Worker(flip, "flipper"), _Worker(evaluate, "evaluator")
    flipper.start()
    evaluator.start()
    _join_all(flipper, evaluator)
    for result in evaluator.result:
        assert np.array_equal(result, old_reference, equal_nan=True) or np.array_equal(
            result, new_reference, equal_nan=True
        ), "a result is neither the old nor the new state"
    # The last change (index 5) set the old value.
    final = _Worker(lambda: model.evaluate_model(points), "final")
    final.start()
    np.testing.assert_array_equal(final.finish(), old_reference)
