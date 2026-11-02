"""Tests for the storage and the exact byte limits of the restoration cache (D4).

The tests use small arrays only. No test builds a model.
"""

from __future__ import annotations

import importlib

import numpy as np
import pytest

MODULE_NAME = "LoopStructural.modelling.features._evaluation_cache"
DEFAULT_MAX_BYTES = 64 * 1024 * 1024


def _api(name: str):
    """Return an object of the cache module, or fail with a clear reason when it is absent."""
    try:
        module = importlib.import_module(MODULE_NAME)
    except ImportError as error:
        pytest.fail(f"{MODULE_NAME} is required by D4 (import failed: {error})")
    value = getattr(module, name, None)
    assert value is not None, f"{MODULE_NAME}.{name} is required by D4"
    return value


def _cache(**kwargs):
    return _api("RestorationCache")(**kwargs)


def _array(count: int, fill: float = 1.0) -> np.ndarray:
    """Return a float64 array with exactly 8 * count bytes."""
    return np.full(count, fill, dtype=np.float64)


def test_defaults():
    """A new cache is enabled, empty and has the 64 MiB default limit."""
    cache = _cache()
    assert cache.enabled is True
    assert cache.max_bytes == DEFAULT_MAX_BYTES
    assert cache.retained_bytes == 0
    assert cache.get("missing") is None


def test_constructor_arguments():
    cache = _cache(enabled=False, max_bytes=1024)
    assert cache.enabled is False
    assert cache.max_bytes == 1024


def test_put_get_returns_equal_values():
    cache = _cache()
    values = np.arange(12, dtype=np.float64).reshape(4, 3)
    cache.put("a", values)
    result = cache.get("a")
    assert isinstance(result, np.ndarray)
    np.testing.assert_array_equal(result, values)
    assert result.dtype == values.dtype
    assert result.shape == values.shape


def test_get_miss_returns_none():
    cache = _cache()
    cache.put("a", _array(2))
    assert cache.get("b") is None


def test_get_returns_the_same_object_without_copy():
    """A hit gives the stored array. It does not make a new copy for each call."""
    cache = _cache()
    cache.put("a", _array(4))
    assert cache.get("a") is cache.get("a")


def test_retained_bytes_is_sum_of_nbytes():
    cache = _cache()
    first = np.zeros((5, 3), dtype=np.float64)
    second = np.zeros(7, dtype=np.float32)
    cache.put("a", first)
    cache.put("b", second)
    assert cache.retained_bytes == first.nbytes + second.nbytes == 120 + 28
    assert isinstance(cache.retained_bytes, int)


def test_retained_bytes_is_read_only():
    cache = _cache()
    with pytest.raises(AttributeError):
        cache.retained_bytes = 0
    assert cache.retained_bytes == 0


def test_stored_array_is_read_only():
    cache = _cache()
    cache.put("a", _array(4))
    stored = cache.get("a")
    assert stored.flags.writeable is False
    with pytest.raises(ValueError):
        stored[0] = 5.0


def test_stored_array_is_not_the_writable_array_of_the_caller():
    """Rule: a later write to the array of the caller does not change the stored values."""
    cache = _cache()
    original = _array(4, fill=2.0)
    cache.put("a", original)
    original[:] = -9.0
    np.testing.assert_array_equal(cache.get("a"), np.full(4, 2.0))
    stored = cache.get("a")
    assert stored is not original
    assert not np.shares_memory(stored, original) or stored.flags.writeable is False
    assert original.flags.writeable is True


def test_put_does_not_change_the_flags_of_the_caller_array():
    cache = _cache()
    original = _array(4)
    cache.put("a", original)
    assert original.flags.writeable is True


def test_put_of_a_view_of_a_writable_base_is_isolated_from_the_base():
    """A read-only view of a writable base is not safe. A write to the base must not leak."""
    cache = _cache()
    base = _array(8, fill=3.0)
    view = base[:4]
    view.flags.writeable = False
    cache.put("a", view)
    base[:] = 0.0
    np.testing.assert_array_equal(cache.get("a"), np.full(4, 3.0))


def test_no_second_retained_copy():
    """One put keeps one array. The counted bytes equal the nbytes of one array."""
    cache = _cache()
    values = _array(100)
    cache.put("a", values)
    assert cache.retained_bytes == values.nbytes == 800
    cache.put("a", values)
    assert cache.retained_bytes == 800


def test_non_contiguous_input_is_counted_by_nbytes():
    cache = _cache()
    values = np.arange(20, dtype=np.float64).reshape(5, 4)[:, ::2]
    cache.put("a", values)
    assert cache.retained_bytes == cache.get("a").nbytes == values.nbytes
    np.testing.assert_array_equal(cache.get("a"), values)


def test_overwrite_replaces_entry_and_fixes_count():
    cache = _cache()
    cache.put("a", _array(10, fill=1.0))
    cache.put("b", _array(5))
    cache.put("a", _array(3, fill=7.0))
    np.testing.assert_array_equal(cache.get("a"), np.full(3, 7.0))
    assert cache.retained_bytes == 3 * 8 + 5 * 8


def test_overwrite_with_larger_array_updates_count():
    cache = _cache()
    cache.put("a", _array(2))
    cache.put("a", _array(20))
    assert cache.retained_bytes == 160


def test_lru_evicts_least_recently_put_first():
    cache = _cache(max_bytes=24)
    cache.put("a", _array(1))
    cache.put("b", _array(1))
    cache.put("c", _array(1))
    cache.put("d", _array(1))
    assert cache.get("a") is None
    assert cache.get("b") is not None
    assert cache.get("c") is not None
    assert cache.get("d") is not None
    assert cache.retained_bytes == 24


def test_get_marks_entry_as_most_recently_used():
    cache = _cache(max_bytes=24)
    cache.put("a", _array(1))
    cache.put("b", _array(1))
    cache.put("c", _array(1))
    assert cache.get("a") is not None
    cache.put("d", _array(1))
    assert cache.get("b") is None
    assert cache.get("a") is not None
    assert cache.get("c") is not None
    assert cache.get("d") is not None


def test_overwrite_marks_entry_as_most_recently_used():
    cache = _cache(max_bytes=24)
    cache.put("a", _array(1))
    cache.put("b", _array(1))
    cache.put("c", _array(1))
    cache.put("a", _array(1, fill=5.0))
    cache.put("d", _array(1))
    assert cache.get("b") is None
    assert cache.get("a") is not None


def test_miss_does_not_change_order():
    cache = _cache(max_bytes=16)
    cache.put("a", _array(1))
    cache.put("b", _array(1))
    assert cache.get("zzz") is None
    cache.put("c", _array(1))
    assert cache.get("a") is None
    assert cache.get("b") is not None


def test_one_put_evicts_several_entries_when_needed():
    cache = _cache(max_bytes=40)
    for name in "abcd":
        cache.put(name, _array(1))
    cache.put("big", _array(4))
    assert cache.get("a") is None
    assert cache.get("b") is None
    assert cache.get("c") is None
    assert cache.get("d") is not None
    assert cache.get("big") is not None
    assert cache.retained_bytes == 40


def test_entry_equal_to_limit_is_retained():
    cache = _cache(max_bytes=32)
    cache.put("a", _array(4))
    assert cache.get("a") is not None
    assert cache.retained_bytes == 32


def test_oversize_entry_is_not_retained():
    cache = _cache(max_bytes=31)
    cache.put("a", _array(4))
    assert cache.get("a") is None
    assert cache.retained_bytes == 0


def test_oversize_entry_does_not_evict_other_entries():
    cache = _cache(max_bytes=32)
    cache.put("a", _array(2))
    cache.put("big", _array(100))
    assert cache.get("big") is None
    assert cache.get("a") is not None
    assert cache.retained_bytes == 16


def test_oversize_overwrite_does_not_keep_stale_value():
    """A put with a too large array for a known key must not leave the old value as a hit."""
    cache = _cache(max_bytes=32)
    cache.put("a", _array(2, fill=1.0))
    cache.put("a", _array(100, fill=2.0))
    result = cache.get("a")
    assert result is None
    assert cache.retained_bytes == 0


def test_shrinking_max_bytes_evicts_at_once_least_recently_used_first():
    cache = _cache(max_bytes=64)
    for name in "abcd":
        cache.put(name, _array(2))
    assert cache.get("a") is not None
    cache.max_bytes = 32
    assert cache.max_bytes == 32
    assert cache.retained_bytes == 32
    assert cache.get("b") is None
    assert cache.get("c") is None
    assert cache.get("d") is not None
    assert cache.get("a") is not None


def test_shrinking_to_zero_removes_everything():
    cache = _cache()
    cache.put("a", _array(2))
    cache.max_bytes = 0
    assert cache.retained_bytes == 0
    assert cache.get("a") is None


def test_growing_max_bytes_keeps_entries():
    cache = _cache(max_bytes=16)
    cache.put("a", _array(2))
    cache.max_bytes = 1024
    assert cache.get("a") is not None
    assert cache.retained_bytes == 16


def test_negative_limit_in_constructor_raises_value_error():
    with pytest.raises(ValueError):
        _cache(max_bytes=-1)


def test_negative_limit_by_attribute_raises_value_error_and_keeps_state():
    cache = _cache(max_bytes=64)
    cache.put("a", _array(2))
    with pytest.raises(ValueError):
        cache.max_bytes = -1
    assert cache.max_bytes == 64
    assert cache.get("a") is not None
    assert cache.retained_bytes == 16


def test_zero_limit_stops_retention():
    cache = _cache(max_bytes=0)
    cache.put("a", _array(1))
    assert cache.get("a") is None
    assert cache.retained_bytes == 0


def test_disabled_cache_stores_nothing():
    cache = _cache(enabled=False)
    cache.put("a", _array(1))
    assert cache.get("a") is None
    assert cache.retained_bytes == 0


def test_disabling_stops_retention_and_hides_entries():
    """After disable, put stores nothing and get returns None. Enable again works."""
    cache = _cache()
    cache.put("a", _array(1))
    cache.enabled = False
    assert cache.get("a") is None
    cache.put("b", _array(1))
    assert cache.get("b") is None
    cache.enabled = True
    cache.put("c", _array(1, fill=4.0))
    np.testing.assert_array_equal(cache.get("c"), np.full(1, 4.0))
    assert cache.get("b") is None


def test_disabled_cache_keeps_count_exact():
    cache = _cache()
    cache.put("a", _array(2))
    cache.enabled = False
    cache.put("b", _array(2))
    cache.enabled = True
    expected = 16 * sum(1 for key in ("a", "b") if cache.get(key) is not None)
    assert cache.retained_bytes == expected


def test_zero_limit_can_be_raised_again():
    cache = _cache(max_bytes=0)
    cache.max_bytes = 64
    cache.put("a", _array(2))
    assert cache.get("a") is not None


def test_clear_removes_all_entries_and_keeps_settings():
    cache = _cache(enabled=True, max_bytes=128)
    cache.put("a", _array(2))
    cache.put("b", _array(2))
    cache.clear()
    assert cache.get("a") is None
    assert cache.get("b") is None
    assert cache.retained_bytes == 0
    assert cache.enabled is True
    assert cache.max_bytes == 128
    cache.put("c", _array(2))
    assert cache.get("c") is not None


def test_invalidate_removes_all_entries_and_keeps_settings():
    cache = _cache(enabled=True, max_bytes=128)
    cache.put("a", _array(2))
    cache.put("b", _array(2))
    cache.invalidate()
    assert cache.get("a") is None
    assert cache.get("b") is None
    assert cache.retained_bytes == 0
    assert cache.enabled is True
    assert cache.max_bytes == 128


def test_clear_and_invalidate_on_empty_cache_are_safe():
    cache = _cache()
    cache.clear()
    cache.invalidate()
    assert cache.retained_bytes == 0


def test_clear_keeps_disabled_state():
    cache = _cache(enabled=False)
    cache.clear()
    assert cache.enabled is False


class _FailingArray:
    """Object that fails when the cache converts it to an array."""

    def __array__(self, dtype=None, copy=None):
        raise RuntimeError("conversion failed")

    @property
    def nbytes(self):
        raise RuntimeError("nbytes failed")


def test_failed_put_publishes_no_entry():
    cache = _cache()
    cache.put("keep", _array(2))
    with pytest.raises(Exception):  # noqa: B017
        cache.put("bad", _FailingArray())
    assert cache.get("bad") is None
    assert cache.get("keep") is not None
    assert cache.retained_bytes == 16


def test_failed_put_keeps_old_entry_of_the_same_key():
    cache = _cache()
    cache.put("a", _array(2, fill=6.0))
    with pytest.raises(Exception):  # noqa: B017
        cache.put("a", _FailingArray())
    np.testing.assert_array_equal(cache.get("a"), np.full(2, 6.0))
    assert cache.retained_bytes == 16


class _NoCopyArray(np.ndarray):
    """Real ndarray view where every route to a private copy raises."""

    def copy(self, *args, **kwargs):
        raise RuntimeError("copy failed")

    def __array__(self, *args, **kwargs):
        raise RuntimeError("conversion failed")

    def __copy__(self):
        raise RuntimeError("copy failed")

    def __deepcopy__(self, memo=None):
        raise RuntimeError("copy failed")

    @property
    def nbytes(self):
        raise RuntimeError("nbytes failed")


def test_failed_copy_step_keeps_count_exact():
    """A failure in the copy step must not publish a partial entry or change the count."""
    cache = _cache()
    cache.put("keep", _array(3))
    cache.put("same", _array(2, fill=4.0))
    bad = np.zeros(4).view(_NoCopyArray)
    with pytest.raises(Exception):  # noqa: B017
        cache.put("new", bad)
    with pytest.raises(Exception):  # noqa: B017
        cache.put("same", bad)
    assert cache.get("new") is None
    np.testing.assert_array_equal(cache.get("same"), np.full(2, 4.0))
    assert cache.get("keep") is not None
    assert cache.retained_bytes == 40


def test_failed_put_does_not_evict_other_entries():
    cache = _cache(max_bytes=32)
    cache.put("a", _array(2))
    cache.put("b", _array(2))
    with pytest.raises(Exception):  # noqa: B017
        cache.put("c", _FailingArray())
    assert cache.get("a") is not None
    assert cache.get("b") is not None
    assert cache.retained_bytes == 32


def test_many_operations_keep_count_equal_to_stored_nbytes():
    """The count equals the sum of the nbytes of all entries that a get can still return."""
    cache = _cache(max_bytes=200)
    rng = np.random.default_rng(0)
    keys = [f"k{index}" for index in range(12)]
    for step in range(200):
        key = keys[int(rng.integers(len(keys)))]
        action = int(rng.integers(4))
        if action == 0:
            cache.put(key, _array(int(rng.integers(1, 40))))
        elif action == 1:
            cache.get(key)
        elif action == 2 and step % 25 == 0:
            cache.max_bytes = int(rng.integers(0, 300))
        elif action == 3 and step % 40 == 0:
            cache.clear()
        assert cache.retained_bytes <= cache.max_bytes
    total = sum(cache.get(key).nbytes for key in keys if cache.get(key) is not None)
    assert cache.retained_bytes == total
