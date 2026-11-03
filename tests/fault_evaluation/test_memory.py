"""Tests for the memory accounting of the fault evaluation (D4, D5, FR-010, T022).

* ``retained_bytes`` is the exact sum of ``nbytes`` of the stored arrays.
* One restoration is stored once. The scope and the cache share one array.
* ``max_bytes`` holds across many different point batches.
* The workspace of the digest of the points and of the dependency fingerprint is
  bounded. It is far smaller than the arrays that are hashed.
* A disabled cache retains zero bytes.

The unit tests use small batches. The large case with 10 million points and the peak RSS
is not part of normal pytest (T029, T030).
"""

from __future__ import annotations

import tracemalloc
from collections.abc import Callable

import numpy as np

from LoopStructural.modelling.features import _evaluation
from LoopStructural.modelling.features._evaluation_cache import (
    RestorationCache,
    dependency_fingerprint,
)
from tests.fault_evaluation.models import build_model
from tests.fault_evaluation.test_invalidation import _fault, _group, _points

MEGABYTE = 1024 * 1024


def _batches(count: int, rows: int = 40) -> list[np.ndarray]:
    """Return ``count`` different float64 point batches of equal size inside the f1_u0 box."""
    base = _points("f1_u0")[:rows]
    return [np.ascontiguousarray(base + 3.0 * index) for index in range(count)]


def _peak(call: Callable[[], object]) -> int:
    """Return the peak of the Python and NumPy allocations of ``call`` in bytes."""
    call()  # warm call: imports and lazy tables are not counted
    tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        call()
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return peak


# ---------------------------------------------------------------------------
# Exact accounting of the storage
# ---------------------------------------------------------------------------


def test_retained_bytes_is_the_exact_sum_of_nbytes():
    """Put, replace, discard and clear keep the count equal to the sum of the stored arrays."""
    cache = RestorationCache()
    small = np.zeros((5, 3))
    medium = np.zeros((7, 3), dtype=np.float32)
    large = np.zeros((11, 3))
    assert cache.retained_bytes == 0
    cache.put("a", small)
    assert cache.retained_bytes == small.nbytes == 120
    cache.put("b", medium)
    assert cache.retained_bytes == small.nbytes + medium.nbytes
    cache.put("a", large)  # the same key replaces the old array
    assert cache.retained_bytes == large.nbytes + medium.nbytes
    cache.put("a", large)  # an equal put does not add a second copy
    assert cache.retained_bytes == large.nbytes + medium.nbytes
    cache.discard("missing")
    assert cache.retained_bytes == large.nbytes + medium.nbytes
    cache.discard("b")
    assert cache.retained_bytes == large.nbytes
    cache.clear()
    assert cache.retained_bytes == 0


def test_retained_bytes_counts_the_stored_array_not_the_view_of_the_caller():
    """A view of a large array is stored as a compact copy. The count is the view size."""
    cache = RestorationCache()
    wide = np.zeros((100, 3))
    view = wide[::2]
    cache.put("view", view)
    stored = cache.get("view")
    assert stored is not None
    assert cache.retained_bytes == stored.nbytes == view.nbytes == 50 * 3 * 8
    assert not np.shares_memory(stored, wide), "the cache must not keep the array of the caller"
    assert not stored.flags.writeable


def test_eviction_and_limit_change_keep_the_count_exact():
    """Eviction by the limit and by a smaller limit leaves the exact sum of the arrays."""
    size = np.zeros((10, 3)).nbytes
    cache = RestorationCache(max_bytes=3 * size)
    for index in range(5):
        cache.put(index, np.full((10, 3), float(index)))
        assert cache.retained_bytes == min(index + 1, 3) * size
    assert cache.get(0) is None and cache.get(1) is None
    for index in (2, 3, 4):
        np.testing.assert_array_equal(cache.get(index), np.full((10, 3), float(index)))
    cache.max_bytes = size
    assert cache.retained_bytes == size
    cache.max_bytes = size - 1
    assert cache.retained_bytes == 0


def test_model_retains_one_array_for_each_distinct_batch_and_fault_set():
    """Small batches through the model give an exact count: batch bytes for each distinct key."""
    model = build_model("f1_u0")
    group = _group(model, 0)
    cache = model.evaluation_cache
    batches = _batches(4)
    for number, batch in enumerate(batches, start=1):
        group.evaluate_value(batch)
        assert cache.retained_bytes == number * batch.nbytes
    # A repeated batch, also an equal copy, adds nothing.
    group.evaluate_value(batches[0])
    group.evaluate_value(batches[1].copy())
    assert cache.retained_bytes == 4 * batches[0].nbytes


def test_nested_fault_sets_are_counted_for_each_set():
    """Three different ordered fault sets on one batch retain three arrays."""
    model = build_model("f4_u3")
    points = _points("f4_u3")
    full = list(_group(model, 3).faults)
    assert len(full) >= 3
    _group(model, 0).faults = full[:2]
    _group(model, 1).faults = list(reversed(full))
    _group(model, 2).faults = list(full)
    for index in range(3):
        _evaluation.restore(_group(model, index), points)
    assert model.evaluation_cache.retained_bytes == 3 * points.nbytes


# ---------------------------------------------------------------------------
# No duplicate arrays
# ---------------------------------------------------------------------------


def test_one_stored_array_for_each_key_shared_with_the_scope():
    """The scope and the cache share one array. A later call hands out the same memory."""
    model = build_model("f1_u0")
    group = _group(model, 0)
    points = _points("f1_u0")
    cache = model.evaluation_cache
    with _evaluation.evaluation_scope(model):
        inside = _evaluation.restore(group, points)
        again = _evaluation.restore(group, points.copy())
        assert inside is again
        assert cache.retained_bytes == points.nbytes
    assert cache.retained_bytes == points.nbytes, "the scope must not keep a second copy"
    later = _evaluation.restore(group, points)
    assert np.shares_memory(later, inside)
    assert not later.flags.writeable
    assert cache.retained_bytes == points.nbytes


def test_two_features_with_one_fault_set_retain_one_array():
    """The group and its model call use the same fault set. One array is retained."""
    model = build_model("f1_u0")
    points = _points("f1_u0")
    group = _group(model, 0)
    model.evaluate_model(points)
    group.evaluate_value(points)
    assert model.evaluation_cache.retained_bytes == points.nbytes
    first = _evaluation.restore(group, points)
    second = _evaluation.restore(group, points.copy())
    assert np.shares_memory(first, second)
    assert model.evaluation_cache.retained_bytes == points.nbytes


def test_scope_keeps_no_second_array_when_the_result_is_too_large():
    """A result that the cache refuses stays in the scope only. Nothing is retained."""
    model = build_model("f1_u0")
    group = _group(model, 0)
    points = _points("f1_u0")
    cache = model.evaluation_cache
    cache.max_bytes = points.nbytes - 1
    with _evaluation.evaluation_scope(model):
        first = _evaluation.restore(group, points)
        second = _evaluation.restore(group, points)
        assert first is second
    assert cache.retained_bytes == 0


# ---------------------------------------------------------------------------
# The byte limit
# ---------------------------------------------------------------------------


def test_max_bytes_holds_across_many_different_batches():
    """The count never exceeds the limit, also for many different batches."""
    model = build_model("f1_u0")
    group = _group(model, 0)
    cache = model.evaluation_cache
    batches = _batches(25)
    size = batches[0].nbytes
    limit = 3 * size + size // 2
    cache.max_bytes = limit
    for number, batch in enumerate(batches, start=1):
        group.evaluate_value(batch)
        assert cache.retained_bytes <= limit
        assert cache.retained_bytes == min(number, 3) * size
    # The newest batches stay. The oldest were evicted and restore again.
    assert cache.retained_bytes == 3 * size


def test_max_bytes_holds_when_batches_have_different_sizes():
    """Batches of growing size never push the count over the limit."""
    model = build_model("f1_u0")
    group = _group(model, 0)
    cache = model.evaluation_cache
    base = _points("f1_u0")
    limit = 40 * 3 * 8 * 4
    cache.max_bytes = limit
    for rows in (10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 20, 10):
        batch = np.ascontiguousarray(base[:rows] + float(rows))
        group.evaluate_value(batch)
        assert cache.retained_bytes <= limit


def test_zero_limit_retains_nothing():
    """``max_bytes = 0`` retains nothing and the evaluation still gives a result."""
    model = build_model("f1_u0")
    group = _group(model, 0)
    points = _points("f1_u0")
    model.evaluation_cache.max_bytes = 0
    assert np.isfinite(group.evaluate_value(points)).any()
    assert model.evaluation_cache.retained_bytes == 0


# ---------------------------------------------------------------------------
# A disabled cache
# ---------------------------------------------------------------------------


def test_disabled_cache_retains_zero_bytes():
    """A disabled cache keeps no array, for direct puts and for model evaluations."""
    model = build_model("f1_u0")
    group = _group(model, 0)
    cache = model.evaluation_cache
    group.evaluate_value(_points("f1_u0"))
    assert cache.retained_bytes > 0
    cache.enabled = False
    assert cache.retained_bytes == 0, "disabling must drop the stored arrays"
    cache.put("direct", np.zeros((10, 3)))
    assert cache.retained_bytes == 0
    assert cache.get("direct") is None
    for batch in _batches(3):
        model.evaluate_model(batch)
        group.evaluate_value(batch)
        assert cache.retained_bytes == 0
    with _evaluation.evaluation_scope(model):
        _evaluation.restore(group, _points("f1_u0"))
        assert cache.retained_bytes == 0
    assert RestorationCache(enabled=False).retained_bytes == 0


def test_disabled_cache_put_retains_zero_bytes_from_the_start():
    """A cache that starts disabled retains nothing after a put."""
    cache = RestorationCache(enabled=False)
    cache.put("a", np.zeros((10, 3)))
    assert cache.retained_bytes == 0
    assert cache.get("a") is None


# ---------------------------------------------------------------------------
# Bounded hashing workspace
# ---------------------------------------------------------------------------


def test_points_digest_uses_a_bounded_workspace():
    """The digest of 24 MB of points allocates far less than the points.

    A copy of the array (for example ``tobytes``) would need 24 MB. The limit of 6 MB
    leaves room for the chunk buffer, the digest state and Python objects.
    """
    points = np.linspace(0.0, 1.0, 3_000_000).reshape(-1, 3)
    size = points.nbytes
    assert size >= 24_000_000
    peak = _peak(lambda: _evaluation._digest(points))
    assert peak < 6 * MEGABYTE < size, f"peak {peak} bytes for an array of {size} bytes"


def test_points_digest_of_a_view_uses_a_bounded_workspace():
    """The digest of a large non-contiguous view must not copy the whole array.

    A Fortran-order array is a plausible input (for example a transposed array). Chunks
    of bounded size must be copied, not the full array.
    """
    points = np.asfortranarray(np.linspace(0.0, 1.0, 3_000_000).reshape(-1, 3))
    assert not points.flags.c_contiguous
    size = points.nbytes
    peak = _peak(lambda: _evaluation._digest(points))
    assert peak < 6 * MEGABYTE < size, f"peak {peak} bytes for an array of {size} bytes"


def test_points_digest_is_the_same_for_the_chunked_view():
    """The chunked digest of a non-contiguous array equals the digest of its C-order copy."""
    contiguous = np.linspace(0.0, 1.0, 600_000).reshape(-1, 3)
    for variant in (np.asfortranarray(contiguous), np.repeat(contiguous, 2, axis=0)[::2]):
        assert _evaluation._digest(variant) == _evaluation._digest(contiguous)


def test_dependency_fingerprint_uses_a_bounded_workspace():
    """The fingerprint of a field with 3 million coefficients allocates far less than 24 MB."""
    model = build_model("f1_u0")
    group = _group(model, 0)
    interpolator = _fault(model, 0).features[0].interpolator
    interpolator.c = np.linspace(0.0, 1.0, 3_000_000)
    size = interpolator.c.nbytes
    peak = _peak(lambda: dependency_fingerprint(group))
    assert peak < 6 * MEGABYTE < size, f"peak {peak} bytes for an array of {size} bytes"


def test_dependency_fingerprint_of_a_non_contiguous_array_uses_a_bounded_workspace():
    """A strided coefficient array is hashed in chunks, without a full contiguous copy."""
    model = build_model("f1_u0")
    group = _group(model, 0)
    interpolator = _fault(model, 0).features[0].interpolator
    wide = np.linspace(0.0, 1.0, 6_000_000)
    interpolator.c = wide[::2]
    size = interpolator.c.nbytes
    assert not interpolator.c.flags.c_contiguous
    peak = _peak(lambda: dependency_fingerprint(group))
    assert peak < 6 * MEGABYTE < size, f"peak {peak} bytes for an array of {size} bytes"
