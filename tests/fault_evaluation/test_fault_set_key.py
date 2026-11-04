"""Tests for the fault-set cache key of the restoration cache (D4a, FR-009, SC-003).

D4a puts the ordered fault set in the persistent key, and nothing of the root feature.
The names are private API. This file fixes them for task T041.

``fault_set_fingerprint(faults, build=True)`` in ``LoopStructural.modelling.features._evaluation_cache``
    A 32-byte SHA256 digest of the state that the restoration of the ordered faults
    depends on, or None when a dependency is opaque.

``RestorationCache.sync_revision(revision)``
    Clears all entries when the revision differs from the stored one.

Each test imports the new names inside the test, so a test fails alone while the name is
absent. The tests for ``_update_frame`` and the pickle defect use present names.
"""

from __future__ import annotations

import hashlib
import threading
from collections.abc import Callable

import numpy as np
import pandas as pd
import pytest

from benchmarks.fault_evaluation.instrumentation import EvaluationInstrumentation
from LoopStructural.modelling.features import _evaluation_cache
from LoopStructural.modelling.features._evaluation_cache import (
    RestorationCache,
    _update_frame,
    dependency_fingerprint,
)
from tests.fault_evaluation.models import build_model
from tests.fault_evaluation.test_invalidation import (
    _fault,
    _group,
    _points,
    _revision,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _fault_set_fingerprint(faults, build: bool = True):
    """Call ``fault_set_fingerprint``, or fail with a clear reason when it is absent."""
    from LoopStructural.modelling.features import _evaluation_cache

    function = getattr(_evaluation_cache, "fault_set_fingerprint", None)
    assert function is not None, "fault_set_fingerprint is required by D4a"
    return function(faults, build)


def _sync_revision(cache: RestorationCache, revision: int) -> None:
    """Call ``cache.sync_revision``, or fail with a clear reason when it is absent."""
    method = getattr(cache, "sync_revision", None)
    assert method is not None, "RestorationCache.sync_revision is required by D4a"
    method(revision)


def _fault_calls(call: Callable[[], np.ndarray]) -> int:
    """Run ``call`` and return the number of fault calls."""
    with EvaluationInstrumentation() as instrumentation:
        call()
    return instrumentation.snapshot()["fault_calls_total"]


def _run_with_timeout(call: Callable, seconds: float = 60.0):
    """Run ``call`` in a daemon thread. Fail when it does not end, re-raise its error."""
    outcome: dict = {}

    def target() -> None:
        try:
            outcome["value"] = call()
        except BaseException as error:  # noqa: BLE001 - the error is re-raised below
            outcome["error"] = error

    thread = threading.Thread(target=target, daemon=True)
    thread.start()
    thread.join(timeout=seconds)
    assert not thread.is_alive(), "the call did not end (a lock that cannot be acquired?)"
    if "error" in outcome:
        raise outcome["error"]
    return outcome.get("value")


# ---------------------------------------------------------------------------
# (1) The fingerprint of a fault set
# ---------------------------------------------------------------------------


def test_groups_with_equal_ordered_faults_have_equal_fault_set_bytes():
    """The 4 groups of ``f4_u3`` have the same ordered faults and different root data."""
    model = build_model("f4_u3")
    groups = [_group(model, index) for index in range(4)]
    assert all(list(group.faults) == list(groups[0].faults) for group in groups)
    digests = [_fault_set_fingerprint(group.faults) for group in groups]
    assert digests[0] is not None and len(digests[0]) == 32
    assert digests == [digests[0]] * 4


def test_other_fault_order_or_subset_gives_other_bytes():
    """A reversed list and a subset of the faults give two other digests."""
    model = build_model("f4_u3")
    faults = list(_group(model, 0).faults)
    assert len(faults) >= 3
    full = _fault_set_fingerprint(faults)
    reordered = _fault_set_fingerprint(list(reversed(faults)))
    subset = _fault_set_fingerprint(faults[:2])
    assert None not in (full, reordered, subset)
    assert len({full, reordered, subset}) == 3


def test_opaque_region_on_a_fault_gives_none():
    """A plain function as a region of a fault has unknown state, so there is no digest."""
    model = build_model("f4_u3")
    faults = list(_group(model, 0).faults)
    before = _fault_set_fingerprint(faults)
    assert before is not None
    _fault(model, 0).regions.append(lambda xyz: xyz[:, 0] < 600.0)
    assert _fault_set_fingerprint(faults) is None


# ---------------------------------------------------------------------------
# (2) Sharing between top-level calls (SC-003)
# ---------------------------------------------------------------------------


def test_separate_top_level_calls_share_the_restoration_of_equal_fault_sets():
    """After ``group_0`` a separate call on ``group_1`` makes no fault call."""
    model = build_model("f4_u3")
    points = _points("f4_u3")
    assert model.evaluation_cache.enabled
    assert _fault_calls(lambda: _group(model, 0).evaluate_value(points)) == 4
    assert _fault_calls(lambda: _group(model, 1).evaluate_value(points)) == 0


def test_cold_evaluate_model_makes_exactly_one_call_for_each_fault():
    """A cold ``evaluate_model`` on ``f4_u3`` applies each of the 4 faults once (SC-003)."""
    model = build_model("f4_u3")
    points = _points("f4_u3")
    model.evaluation_cache.clear()
    assert _fault_calls(lambda: model.evaluate_model(points)) == 4


# ---------------------------------------------------------------------------
# (3) The revision term of the cache
# ---------------------------------------------------------------------------


def _array() -> np.ndarray:
    return np.arange(12, dtype=float).reshape(4, 3)


def test_sync_revision_keeps_entries_for_the_first_and_the_same_revision():
    """The first call stores the revision. The same revision keeps the entries."""
    cache = RestorationCache()
    _sync_revision(cache, 5)
    cache.put("a", _array())
    size = cache.retained_bytes
    assert size == _array().nbytes
    _sync_revision(cache, 5)
    assert cache.retained_bytes == size
    assert cache.get("a") is not None


def test_sync_revision_clears_all_entries_for_another_revision():
    """A different revision, higher or lower, removes all entries."""
    cache = RestorationCache()
    _sync_revision(cache, 5)
    cache.put("a", _array())
    cache.put("b", _array() + 1.0)
    assert cache.retained_bytes > 0
    _sync_revision(cache, 6)
    assert cache.retained_bytes == 0
    assert cache.get("a") is None
    cache.put("c", _array())
    _sync_revision(cache, 4)
    assert cache.retained_bytes == 0


# ---------------------------------------------------------------------------
# (4) The publish guard
# ---------------------------------------------------------------------------


def test_revision_change_during_the_faults_blocks_the_store(monkeypatch):
    """When a fault call bumps the revision, the result stays out of the cache.

    The bump is made through the setter ``GeologicalModel.data``. The fault then runs its real code. The
    result must equal the result of a call with the cache off.
    """
    model = build_model("f1_u0")
    group = _group(model, 0)
    fault = group.faults[0]
    points = _points("f1_u0")
    original = fault.apply_to_points

    def bumping(*args, **kwargs):
        model.data = model.data.copy()  # this setter always bumps the revision
        return original(*args, **kwargs)

    monkeypatch.setattr(fault, "apply_to_points", bumping)
    revision = _revision(model)
    cache = model.evaluation_cache
    result = group.evaluate_value(points)
    assert _revision(model) > revision, "the setter must bump the revision"
    assert cache.retained_bytes == 0, "a result made while the revision changed is not stored"
    cache.enabled = False
    reference = group.evaluate_value(points)
    np.testing.assert_array_equal(result, reference)


# ---------------------------------------------------------------------------
# (5) The bounding-box transform
# ---------------------------------------------------------------------------


def _shift_box(model) -> None:
    """Change the transform of the box of the fault interpolator in place."""
    box = _fault(model, 0).features[0].interpolator.bounding_box
    box._world_to_local[0, 3] += 5.0


def test_box_transform_change_changes_the_fault_set_fingerprint():
    """An in-place edit of the transform of the interpolator box changes the digest."""
    model = build_model("f4_u3")
    faults = list(_group(model, 0).faults)
    before = _fault_set_fingerprint(faults)
    _shift_box(model)
    assert _fault_set_fingerprint(faults) != before


def test_box_transform_change_changes_the_dependency_fingerprint():
    """The same edit changes the digest of a feature that reaches the fault."""
    model = build_model("f4_u3")
    group = _group(model, 0)
    before = dependency_fingerprint(group)
    _shift_box(model)
    assert dependency_fingerprint(group) != before


def test_box_transform_change_gives_no_stale_restoration():
    """After a change of the transform, the cache gives no stale hit.

    The first call stores the restoration. After the change, a call with the cache on must
    apply the faults again, and its result must equal the result with the cache off. The
    policy of the coordinate owner (accept and rebuild, or reject before mutation) belongs to
    spec 005. This test covers only the stale-hit rule, which both policies need.
    """
    model = build_model("f4_u3")
    points = _points("f4_u3")
    group = _group(model, 0)
    group.evaluate_value(points)
    assert model.evaluation_cache.retained_bytes > 0
    _shift_box(model)
    with EvaluationInstrumentation() as instrumentation:
        cached = group.evaluate_value(points)
    assert instrumentation.snapshot()["fault_calls_total"] > 0, "stale hit after the change"
    model.evaluation_cache.enabled = False
    reference = group.evaluate_value(points)
    np.testing.assert_array_equal(cached, reference)


# ---------------------------------------------------------------------------
# (6) The digest of a DataFrame
# ---------------------------------------------------------------------------


def _frame_digest(frame: pd.DataFrame) -> bytes:
    hasher = hashlib.sha256()
    _update_frame(hasher, frame)
    return hasher.digest()


def _frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "X": [1.0, 2.0, 3.0],
            "n": np.array([1, 2, 3], dtype=np.int64),
            "name": ["a", "b", "c"],
        }
    )


def _small_frame_cells() -> int:
    """Return the limit of cells for the small-frame path of the digest.

    Returns
    -------
    int
        The value of ``_SMALL_FRAME_CELLS`` in the module ``_evaluation_cache``.
    """
    limit = getattr(_evaluation_cache, "_SMALL_FRAME_CELLS", None)
    assert limit is not None, "D4b needs _SMALL_FRAME_CELLS"
    return limit


def _sized_frame(rows: int) -> pd.DataFrame:
    """Make a frame with the columns of ``_frame()`` and ``rows`` rows.

    Parameters
    ----------
    rows : int
        The number of rows.

    Returns
    -------
    pd.DataFrame
        A frame with one float, one integer and one text column.
    """
    return pd.DataFrame(
        {
            "X": np.arange(rows) * 0.5,
            "n": np.arange(rows, dtype=np.int64),
            "name": [f"s{i}" for i in range(rows)],
        }
    )


def _rows_for(size: str) -> int:
    """Return the number of rows for the size ``"small"`` or ``"large"``."""
    if size == "small":
        return 3
    return _small_frame_cells() // 3 + 1


def _edit_value(frame: pd.DataFrame) -> pd.DataFrame:
    changed = frame.copy()
    changed.iloc[0, 0] = 99.5
    return changed


def _edit_object_value(frame: pd.DataFrame) -> pd.DataFrame:
    changed = frame.copy()
    changed.loc[0, "name"] = "edited"
    return changed


@pytest.mark.parametrize("size", ["small", "large"])
def test_equal_frames_have_equal_digests(size):
    """Two frames with equal content, also an independent copy, give equal bytes."""
    frame = _sized_frame(_rows_for(size))
    assert _frame_digest(frame) == _frame_digest(frame.copy())


@pytest.mark.parametrize("size", ["small", "large"])
@pytest.mark.parametrize(
    "change",
    [
        _edit_value,
        lambda frame: frame.rename(columns={"X": "Y"}),
        lambda frame: frame.astype({"n": np.float64}),
        _edit_object_value,
    ],
    ids=["value", "column_name", "dtype", "object_value"],
)
def test_changed_frame_has_another_digest(size, change):
    """A value, a column name, a dtype or a value of an object column changes the digest."""
    frame = _sized_frame(_rows_for(size))
    assert _frame_digest(change(frame)) != _frame_digest(frame)


@pytest.mark.parametrize("size", ["small", "large"])
@pytest.mark.parametrize(
    "edit",
    [
        lambda frame: frame.iloc.__setitem__((1, 0), 99.25),
        lambda frame: frame.loc.__setitem__((1, "name"), "edited"),
    ],
    ids=["number", "text"],
)
def test_in_place_edit_changes_the_digest(size, edit):
    """An edit of the same frame object changes the digest."""
    frame = _sized_frame(_rows_for(size))
    before = _frame_digest(frame)
    edit(frame)
    assert _frame_digest(frame) != before


def test_object_values_of_other_types_have_other_digests():
    """An integer and a float with an equal number in an object column differ."""
    first = pd.DataFrame({"v": pd.Series([1, "a"], dtype=object)})
    second = pd.DataFrame({"v": pd.Series([1.0, "a"], dtype=object)})
    assert _frame_digest(first) != _frame_digest(second)


def test_small_frame_does_not_hash_each_column(monkeypatch):
    """A small frame is hashed in one step, without a loop over the columns."""

    def _forbidden(*args, **kwargs):
        raise AssertionError("per-column path used for a small frame")

    monkeypatch.setattr(pd.util, "hash_pandas_object", _forbidden)
    monkeypatch.setattr(pd.DataFrame, "items", _forbidden)
    _frame_digest(_frame())


# ---------------------------------------------------------------------------
# (7) Pickle defect: the lock of a cache that is loaded
# ---------------------------------------------------------------------------


def test_cache_keeps_its_settings_and_a_usable_lock_after_a_dill_round_trip():
    """``enabled`` and ``max_bytes`` stay, the entries go, and the lock can be acquired."""
    import dill

    cache = RestorationCache(enabled=False, max_bytes=1234)

    def round_trip() -> dict:
        loaded = dill.loads(dill.dumps(cache))
        # Test the lock first: the properties below wait for it.
        acquired = loaded.lock.acquire(timeout=2)
        if not acquired:
            return {"acquired": False}
        loaded.lock.release()
        return {
            "acquired": True,
            "enabled": loaded.enabled,
            "max_bytes": loaded.max_bytes,
            "retained_bytes": loaded.retained_bytes,
        }

    state = _run_with_timeout(round_trip)
    assert state["acquired"] is True, "the lock of a loaded cache cannot be acquired"
    assert state["enabled"] is False
    assert state["max_bytes"] == 1234
    assert state["retained_bytes"] == 0


def test_loaded_model_evaluates_with_the_same_ids():
    """A model that is saved and loaded again does not hang and gives the same ids."""
    import dill

    model = build_model("f1_u0")
    points = _points("f1_u0")
    before = model.evaluate_model(points)

    def round_trip() -> np.ndarray:
        loaded = dill.loads(dill.dumps(model))
        return loaded.evaluate_model(points)

    after = _run_with_timeout(round_trip)
    np.testing.assert_array_equal(after, before)
