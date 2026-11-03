"""Tests for the reuse of fault restorations between evaluation calls (D4, FR-009, SC-010).

The model cache ``GeologicalModel.evaluation_cache`` keeps restored coordinates between
top-level evaluation calls. These tests define the behaviour that T019 must give.

* A second ``evaluate_model`` call on equal points makes zero fault calls.
* The additive controls ``enabled``, ``max_bytes``, ``clear`` and ``retained_bytes`` work.
* The key holds the full content of the points. A copy with equal content hits. A
  reordered or changed array does not hit, and gives the same result as a cache-disabled
  call.
* A change of the model state makes a new restoration. The result equals the result of a
  cache-disabled call on the changed model.
* The scope and the cache share one stored array.
* A model with an opaque dependency (fingerprint None) never uses a persistent entry.

Fault calls are the calls of ``FaultSegment.apply_to_points``. The reference result is
always made with the cache disabled.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pytest

from benchmarks.fault_evaluation.instrumentation import EvaluationInstrumentation
from LoopStructural import GeologicalModel
from LoopStructural.modelling.features import _evaluation
from LoopStructural.modelling.features._evaluation_cache import (
    RestorationCache,
    fault_set_fingerprint,
)
from tests.fault_evaluation.models import build_model
from tests.fault_evaluation.test_invalidation import (
    MUTATION_IDS,
    MUTATIONS,
    _built_tiny_model,
    _evaluate_disabled,
    _fault,
    _group,
    _points,
    _revision,
)

GOLDEN_MODELS = ["f4_u3", "f1_u0"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _calls(call: Callable[[], np.ndarray]):
    """Run ``call``. Return its result and the number of fault calls."""
    with EvaluationInstrumentation() as instrumentation:
        result = call()
    return result, instrumentation.snapshot()["fault_calls_total"]


def _reference(model: GeologicalModel, call: Callable[[], np.ndarray]) -> np.ndarray:
    """Return the result of ``call`` with the cache of ``model`` disabled."""
    result, _ = _evaluate_disabled(model, call)
    return result


def _cache(model: GeologicalModel) -> RestorationCache:
    """Return ``model.evaluation_cache``, or fail with a clear reason when it is absent."""
    cache = getattr(model, "evaluation_cache", None)
    assert cache is not None, "GeologicalModel.evaluation_cache is required by D4"
    return cache


def _warm_model(name: str = "f1_u0"):
    """Return a model, golden points and the first (cache-enabled) result."""
    model = build_model(name)
    points = _points(name)
    first, calls = _calls(lambda: model.evaluate_model(points))
    assert calls > 0, "the first call must restore points"
    return model, points, first


def _single_restoration_model():
    """Return the f1_u0 model, its only faulted group and a function for the group value.

    The group has one fault, so one distinct point batch gives one stored restoration of
    exactly ``points.nbytes`` bytes.
    """
    model = build_model("f1_u0")
    group = _group(model, 0)
    assert len(group.faults) == 1
    return model, group


def _batches(count: int) -> list[np.ndarray]:
    """Return ``count`` different float64 point batches of equal size inside the f1_u0 box."""
    base = _points("f1_u0")[:40]
    return [np.ascontiguousarray(base + 3.0 * index) for index in range(count)]


# ---------------------------------------------------------------------------
# Reuse between calls
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", GOLDEN_MODELS)
def test_second_evaluate_model_makes_no_fault_calls(name):
    """The second call on the same grid restores nothing and gives the same result."""
    model = build_model(name)
    points = _points(name)
    reference = _reference(model, lambda: model.evaluate_model(points))
    first, first_calls = _calls(lambda: model.evaluate_model(points))
    second, second_calls = _calls(lambda: model.evaluate_model(points))
    assert first_calls > 0
    assert second_calls == 0, "the second call must reuse the persistent restorations"
    np.testing.assert_array_equal(first, reference)
    np.testing.assert_array_equal(second, reference)


def test_second_feature_value_call_makes_no_fault_calls():
    """A feature value call on equal points also reuses the entry of the first call."""
    model, group = _single_restoration_model()
    points = _points("f1_u0")
    reference = _reference(model, lambda: group.evaluate_value(points))
    first, first_calls = _calls(lambda: group.evaluate_value(points))
    second, second_calls = _calls(lambda: group.evaluate_value(points))
    assert first_calls == 1
    assert second_calls == 0
    np.testing.assert_array_equal(first, reference)
    np.testing.assert_array_equal(second, reference)


def test_identical_copy_of_the_points_hits():
    """An independent copy with equal content gives a hit. Object identity is not the key."""
    model, points, first = _warm_model()
    copy = points.copy()
    assert copy is not points
    second, calls = _calls(lambda: model.evaluate_model(copy))
    assert calls == 0
    np.testing.assert_array_equal(second, first)


def test_strided_rows_with_equal_content_hit():
    """A strided view and a Fortran-order array with equal content hit the entry."""
    model, group = _single_restoration_model()
    wide = np.repeat(_points("f1_u0"), 2, axis=0)
    strided = wide[::2]
    assert not strided.flags.c_contiguous
    contiguous = np.ascontiguousarray(strided)
    # Make the reference first. A cache-disabled call clears all entries of the model.
    reference = _reference(model, lambda: group.evaluate_value(contiguous))
    group.evaluate_value(contiguous)
    for variant in (strided, np.asfortranarray(contiguous)):
        result, calls = _calls(lambda variant=variant: group.evaluate_value(variant))
        assert calls == 0
        np.testing.assert_array_equal(result, reference)


def test_row_reorder_restores_again_and_matches_the_reference():
    """The same rows in another order are another point batch. They must restore again."""
    model = build_model("f1_u0")
    points = _points("f1_u0")
    order = np.random.default_rng(7).permutation(len(points))
    shuffled = points[order]
    assert not np.array_equal(shuffled, points)
    # Make the reference first. A cache-disabled call clears all entries of the model,
    # so it must not run after the warm call.
    reference = _reference(model, lambda: model.evaluate_model(shuffled))
    first, first_calls = _calls(lambda: model.evaluate_model(points))
    assert first_calls > 0, "the first call must restore points"
    result, calls = _calls(lambda: model.evaluate_model(shuffled))
    assert calls > 0, "a reordered batch must not hit the entry of the original order"
    np.testing.assert_array_equal(result, reference)
    np.testing.assert_array_equal(result, first[order])
    # Both orders are now stored.
    _, again = _calls(lambda: model.evaluate_model(points))
    assert again == 0


def test_input_mutation_after_a_call_restores_again():
    """A change of the content of the input array between calls must not give a stale hit."""
    model, points, first = _warm_model()
    mutated = points.copy()
    model.evaluate_model(mutated)
    mutated[: len(mutated) // 2, 0] += 40.0
    reference = _reference(model, lambda: model.evaluate_model(mutated.copy()))
    result, calls = _calls(lambda: model.evaluate_model(mutated))
    assert calls > 0
    np.testing.assert_array_equal(result, reference)
    # The old content still hits, or at least is correct.
    old, _ = _calls(lambda: model.evaluate_model(points))
    np.testing.assert_array_equal(old, first)


def test_feature_value_input_mutation_restores_again():
    """The same check for the value call of a feature with an array changed in place."""
    model, group = _single_restoration_model()
    points = _points("f1_u0").copy()
    group.evaluate_value(points)
    points[:, 1] += 25.0
    reference = _reference(model, lambda: group.evaluate_value(points.copy()))
    result, calls = _calls(lambda: group.evaluate_value(points))
    assert calls == 1
    np.testing.assert_array_equal(result, reference)


# ---------------------------------------------------------------------------
# Models, fault sets and reverse
# ---------------------------------------------------------------------------


def test_different_models_do_not_share_entries():
    """Two equal models have their own caches. A warm call of one does not help the other."""
    first_model, points, first = _warm_model()
    other = build_model("f1_u0")
    assert _cache(other) is not _cache(first_model)
    stored = _cache(first_model).retained_bytes
    assert stored > 0, "the first model must retain its restoration"
    assert _cache(other).retained_bytes == 0
    result, calls = _calls(lambda: other.evaluate_model(points))
    assert calls > 0, "the other model must restore for itself"
    np.testing.assert_array_equal(result, first)
    assert _cache(first_model).retained_bytes == stored
    assert _cache(other).retained_bytes == stored


def test_different_fault_sets_have_different_entries():
    """Groups with different ORDERED fault sets do not share a stored restoration."""
    model = build_model("f4_u3")
    points = _points("f4_u3")
    full = list(_group(model, 3).faults)
    assert len(full) >= 3
    _group(model, 0).faults = full[:2]
    _group(model, 1).faults = list(reversed(full))
    _group(model, 2).faults = list(full)
    groups = [_group(model, index) for index in range(3)]
    restored = [_evaluation.restore(group, points) for group in groups]
    assert not np.array_equal(restored[0], restored[2]), "subset and full set must differ"
    assert not np.array_equal(restored[1], restored[2]), "the fault order must matter"
    assert _cache(model).retained_bytes == 3 * points.nbytes
    for group, expected in zip(groups, restored):
        result, calls = _calls(lambda group=group: _evaluation.restore(group, points))
        assert calls == 0
        np.testing.assert_array_equal(result, expected)


def test_reverse_flag_has_its_own_entry():
    """``reverse`` is part of the key, so the two directions have two entries.

    The values of the two directions may be equal today, because ``apply_to_points``
    ignores the flag. The entries must stay separate.
    """
    model, group = _single_restoration_model()
    points = _points("f1_u0")
    forward = _evaluation.restore(group, points)
    backward, calls = _calls(lambda: _evaluation.restore(group, points, reverse=True))
    assert calls > 0, "a reverse call must not use the forward entry"
    size = points.nbytes
    assert _cache(model).retained_bytes == 2 * size
    for reverse, expected in ((False, forward), (True, backward)):
        result, calls = _calls(lambda reverse=reverse: _evaluation.restore(group, points, reverse))
        assert calls == 0
        np.testing.assert_array_equal(result, expected)


# ---------------------------------------------------------------------------
# Changes of the model
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("mutation", MUTATIONS, ids=MUTATION_IDS)
def test_model_change_restores_again_and_matches_the_changed_reference(mutation):
    """After a change the next call restores again. Its result equals the changed reference."""
    model = build_model(mutation.model)
    points = _points(mutation.model)
    model.evaluate_model(points)
    _, warm_calls = _calls(lambda: model.evaluate_model(points))
    assert warm_calls == 0, "the call before the change must reuse its restorations"
    undo = mutation.apply(model)
    try:
        reference = _reference(model, lambda: model.evaluate_model(points))
        result, calls = _calls(lambda: model.evaluate_model(points))
        if mutation.touches_faults or mutation.revision:
            assert calls > 0, "a changed model must not use a stale entry"
        np.testing.assert_array_equal(result, reference)
    finally:
        undo()


def test_fault_data_change_with_rebuild_restores_again():
    """A change of the data of a fault and a rebuild gives a new restoration."""
    model, group, points, _ = _built_tiny_model()
    group.evaluate_value(points)
    _, warm_calls = _calls(lambda: group.evaluate_value(points))
    assert warm_calls == 0
    before = _reference(model, lambda: group.evaluate_value(points))
    # The fault builder holds no data. Its member builders hold the data frames.
    builder = _fault(model, 0).builder.builders[0]
    revision = _revision(model)
    builder.data.loc[builder.data["val"] == 0.0, "X"] += 150.0
    builder.set_not_up_to_date(None)
    assert _revision(model) > revision
    reference = _reference(model, lambda: group.evaluate_value(points))
    assert not np.array_equal(reference, before, equal_nan=True), "the change has no effect"
    result, calls = _calls(lambda: group.evaluate_value(points))
    assert calls > 0
    np.testing.assert_array_equal(result, reference)


def test_faults_enabled_change_restores_again():
    """Switching ``faults_enabled`` off and on again follows the current state each time."""
    model, group = _single_restoration_model()
    points = _points("f1_u0")
    enabled = group.evaluate_value(points)
    group.faults_enabled = False
    disabled, calls = _calls(lambda: group.evaluate_value(points))
    assert calls == 0, "no fault is applied when faults are off"
    assert not np.array_equal(enabled, disabled, equal_nan=True)
    np.testing.assert_array_equal(disabled, _reference(model, lambda: group.evaluate_value(points)))
    group.faults_enabled = True
    result, _ = _calls(lambda: group.evaluate_value(points))
    np.testing.assert_array_equal(result, enabled)


# ---------------------------------------------------------------------------
# Controls of the cache
# ---------------------------------------------------------------------------


def test_default_controls():
    """The cache is enabled with 64 MiB and holds nothing before the first call."""
    cache = _cache(build_model("f1_u0"))
    assert isinstance(cache, RestorationCache)
    assert cache.enabled is True
    assert cache.max_bytes == 64 * 1024 * 1024
    assert cache.retained_bytes == 0


def test_a_call_stores_the_restoration():
    """After a call with faults the cache retains bytes. The count is read-only."""
    model, _, _ = _warm_model()
    cache = _cache(model)
    assert cache.retained_bytes > 0
    with pytest.raises(AttributeError):
        cache.retained_bytes = 0


def test_disabled_cache_gives_no_reuse_and_clears_entries():
    """``enabled=False`` removes the entries. Every call restores and nothing is kept."""
    model, points, first = _warm_model()
    cache = _cache(model)
    assert cache.retained_bytes > 0
    cache.enabled = False
    assert cache.retained_bytes == 0
    results, counts = [], []
    for _ in range(2):
        result, calls = _calls(lambda: model.evaluate_model(points))
        results.append(result)
        counts.append(calls)
    assert counts[0] > 0 and counts[1] == counts[0], "a disabled cache gives no reuse"
    assert cache.retained_bytes == 0
    for result in results:
        np.testing.assert_array_equal(result, first)
    cache.enabled = True
    _, calls = _calls(lambda: model.evaluate_model(points))
    assert calls > 0, "the entries were cleared when the cache was disabled"
    _, calls = _calls(lambda: model.evaluate_model(points))
    assert calls == 0


def test_clear_removes_entries_but_keeps_the_settings():
    """``clear()`` forces a new restoration. It does not change the revision or the limits."""
    model, points, first = _warm_model()
    cache = _cache(model)
    limit, revision = cache.max_bytes, _revision(model)
    cache.clear()
    assert cache.retained_bytes == 0
    assert cache.enabled is True and cache.max_bytes == limit
    assert _revision(model) == revision
    result, calls = _calls(lambda: model.evaluate_model(points))
    assert calls > 0
    np.testing.assert_array_equal(result, first)
    _, calls = _calls(lambda: model.evaluate_model(points))
    assert calls == 0


def test_result_larger_than_max_bytes_is_not_stored():
    """An entry larger than ``max_bytes`` is not retained, so the next call restores again."""
    model, group = _single_restoration_model()
    points = _points("f1_u0")
    cache = _cache(model)
    cache.max_bytes = points.nbytes - 1
    reference = _reference(model, lambda: group.evaluate_value(points))
    for _ in range(2):
        result, calls = _calls(lambda: group.evaluate_value(points))
        assert calls == 1
        assert cache.retained_bytes == 0
        np.testing.assert_array_equal(result, reference)
    cache.max_bytes = points.nbytes
    group.evaluate_value(points)
    assert cache.retained_bytes == points.nbytes
    _, calls = _calls(lambda: group.evaluate_value(points))
    assert calls == 0


def test_negative_limit_is_refused():
    """A negative ``max_bytes`` raises ValueError and keeps the old limit."""
    cache = _cache(build_model("f1_u0"))
    limit = cache.max_bytes
    with pytest.raises(ValueError):
        cache.max_bytes = -1
    assert cache.max_bytes == limit


def test_lru_eviction_order():
    """The least recently used restoration leaves first. A hit marks an entry as recent."""
    model, group = _single_restoration_model()
    cache = _cache(model)
    first, second, third = _batches(3)
    size = first.nbytes
    cache.max_bytes = 2 * size

    def count(points: np.ndarray) -> int:
        return _calls(lambda: group.evaluate_value(points))[1]

    assert count(first) == 1
    assert count(second) == 1
    assert cache.retained_bytes == 2 * size
    assert count(first) == 0, "a hit marks the first entry as recent"
    assert count(third) == 1, "the third entry evicts the second one"
    assert cache.retained_bytes == 2 * size
    assert count(first) == 0
    assert count(third) == 0
    assert count(second) == 1, "the second entry was the least recently used"
    # Shrinking the limit evicts at once, oldest first. The newest entry stays.
    cache.max_bytes = size
    assert cache.retained_bytes == size
    assert count(second) == 0


# ---------------------------------------------------------------------------
# Shared storage between the scope and the cache
# ---------------------------------------------------------------------------


def test_scope_and_cache_share_one_stored_array():
    """One restoration is retained once. The array handed out is read-only and shared."""
    model, group = _single_restoration_model()
    points = _points("f1_u0")
    cache = _cache(model)
    with _evaluation.evaluation_scope(model):
        inside = _evaluation.restore(group, points)
        again = _evaluation.restore(group, points)
        assert inside is again
        assert cache.retained_bytes == inside.nbytes == points.nbytes
    assert not inside.flags.writeable
    with pytest.raises(ValueError):
        inside[0, 0] = 0.0
    assert cache.retained_bytes == points.nbytes, "the scope must not keep a second copy"
    later = _evaluation.restore(group, points)
    assert np.shares_memory(later, inside), "a later call must hand out the stored array"
    assert not later.flags.writeable
    assert cache.retained_bytes == points.nbytes


# ---------------------------------------------------------------------------
# Full lookup hashing
# ---------------------------------------------------------------------------


def test_points_are_hashed_in_the_real_call_path_of_a_warm_call(monkeypatch):
    """A warm call computes the digest of the full points. No timing shortcut skips it."""
    model, points, first = _warm_model()
    digested: list[np.ndarray] = []
    original = _evaluation._digest

    def counting_digest(array):
        digested.append(np.array(array, copy=True))
        return original(array)

    monkeypatch.setattr(_evaluation, "_digest", counting_digest)
    second, calls = _calls(lambda: model.evaluate_model(points))
    assert calls == 0, "the warm call must be a hit"
    np.testing.assert_array_equal(second, first)
    assert any(
        array.shape == points.shape and np.array_equal(array, points) for array in digested
    ), "the digest of the full points must be computed on the second call"


def test_mutated_points_are_hashed_again_and_not_trusted(monkeypatch):
    """A mutated array is digested again on each call. Equal identity or time is no proof."""
    model, group = _single_restoration_model()
    points = _points("f1_u0").copy()
    group.evaluate_value(points)
    digested = []
    original = _evaluation._digest

    def counting_digest(array):
        digested.append(len(digested))
        return original(array)

    monkeypatch.setattr(_evaluation, "_digest", counting_digest)
    for shift in (0.0, 10.0, 20.0):
        points[:, 2] += shift
        expected = _reference(model, lambda: group.evaluate_value(points.copy()))
        before = len(digested)
        result = group.evaluate_value(points)
        assert len(digested) > before
        np.testing.assert_array_equal(result, expected)


# ---------------------------------------------------------------------------
# Opaque dependencies
# ---------------------------------------------------------------------------


def test_opaque_dependency_bypasses_persistent_entries():
    """A plain function as a region of a fault has hidden state. The key is None, so no reuse.

    The closure is on a frame member of a fault, because the restoration reads the regions of
    the frame members and not ``FaultSegment.regions``.
    """
    model = build_model("cutting_faults")
    points = _points("cutting_faults")
    group = _group(model, 0)
    state = {"limit": 1.0e9}

    def predicate(xyz):
        return xyz[:, 0] < state["limit"]

    _fault(model, 0)[0].regions.append(predicate)
    assert fault_set_fingerprint(tuple(group.faults)) is None
    cache = _cache(model)
    reference = _reference(model, lambda: group.evaluate_value(points))
    first, first_calls = _calls(lambda: group.evaluate_value(points))
    second, second_calls = _calls(lambda: group.evaluate_value(points))
    assert first_calls > 0
    assert second_calls == first_calls, "an opaque model must restore again in each call"
    assert cache.retained_bytes == 0, "no persistent entry may be kept"
    np.testing.assert_array_equal(first, reference)
    np.testing.assert_array_equal(second, reference)
    # A hidden change of the closure state is seen, because nothing is stored.
    state["limit"] = float(np.median(points[:, 0]))
    changed = _reference(model, lambda: group.evaluate_value(points))
    assert not np.array_equal(changed, reference, equal_nan=True)
    result, _ = _calls(lambda: group.evaluate_value(points))
    np.testing.assert_array_equal(result, changed)


def _evaluation_fingerprint(feature):
    """Return ``dependency_fingerprint(feature)``."""
    from LoopStructural.modelling.features._evaluation_cache import dependency_fingerprint

    return dependency_fingerprint(feature)
