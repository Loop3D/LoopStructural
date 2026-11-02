"""Bounded storage for restored fault evaluation results.

The module defines :class:`RestorationCache`. The cache keeps read-only NumPy arrays
under hashable keys. It removes the least recently used entries when the byte limit is
exceeded.

Memory accounting
-----------------
``retained_bytes`` is the exact sum of ``ndarray.nbytes`` of the stored arrays. It does
not include the metadata of the cache. The metadata are the key objects, the
``OrderedDict`` itself, the ndarray headers and the lock. The workspace of a later
hashing step is also outside this count. The metadata and the workspace must be
documented and limited separately from the stored arrays.

The cache keeps no input array, no mask and no final output of the caller. ``put`` stores
one private read-only copy.
"""

from __future__ import annotations

import threading
from collections import OrderedDict
from collections.abc import Hashable

import numpy as np

DEFAULT_MAX_BYTES = 64 * 1024 * 1024


def _check_max_bytes(max_bytes: int) -> int:
    """Return ``max_bytes`` as an int, or raise ValueError when it is negative."""
    value = int(max_bytes)
    if value < 0:
        raise ValueError(f"max_bytes must not be negative, got {max_bytes}")
    return value


class RestorationCache:
    """Thread-safe least-recently-used cache for restored coordinate arrays.

    Parameters
    ----------
    enabled : bool, default True
        When False, the cache retains nothing and ``get`` returns None.
    max_bytes : int, default 64 MiB
        Largest total ``nbytes`` of the stored arrays. Zero stops retention.

    Raises
    ------
    ValueError
        If ``max_bytes`` is negative.
    """

    def __init__(self, enabled: bool = True, max_bytes: int = DEFAULT_MAX_BYTES) -> None:
        self._lock = threading.RLock()
        self._entries: OrderedDict[Hashable, np.ndarray] = OrderedDict()
        self._retained_bytes = 0
        self._max_bytes = _check_max_bytes(max_bytes)
        self._enabled = bool(enabled)

    @property
    def enabled(self) -> bool:
        """bool: Whether the cache retains and returns entries."""
        with self._lock:
            return self._enabled

    @enabled.setter
    def enabled(self, value: bool) -> None:
        with self._lock:
            self._enabled = bool(value)
            if not self._enabled:
                self._clear_locked()

    @property
    def max_bytes(self) -> int:
        """int: Largest total ``nbytes`` of the stored arrays."""
        with self._lock:
            return self._max_bytes

    @max_bytes.setter
    def max_bytes(self, value: int) -> None:
        new_limit = _check_max_bytes(value)
        with self._lock:
            self._max_bytes = new_limit
            self._evict_locked(new_limit)

    @property
    def retained_bytes(self) -> int:
        """int: Exact sum of ``nbytes`` of the stored arrays (read-only)."""
        with self._lock:
            return self._retained_bytes

    def get(self, key: Hashable) -> np.ndarray | None:
        """Return the stored read-only array for ``key``, or None.

        A hit marks the entry as most recently used. The stored array is returned
        without a copy.
        """
        with self._lock:
            if not self._enabled:
                return None
            restored = self._entries.get(key)
            if restored is None:
                return None
            self._entries.move_to_end(key)
            return restored

    def put(self, key: Hashable, restored: np.ndarray) -> None:
        """Store a read-only copy of ``restored`` under ``key``.

        The copy is made before the store changes. If a step fails, the store and the
        byte count stay as they were. An array larger than ``max_bytes`` is not
        retained, and it removes an older entry of the same key so that no stale value
        stays. The array of the caller is not changed and is not kept.
        """
        with self._lock:
            if not self._enabled or self._max_bytes == 0:
                return
            size = int(restored.nbytes)
            if size > self._max_bytes:
                self._discard_locked(key)
                return
            stored = np.array(restored, copy=True, order="C")
            stored.flags.writeable = False
            size = int(stored.nbytes)
            # All steps that can fail are done. The next steps only change the store.
            self._discard_locked(key)
            self._entries[key] = stored
            self._retained_bytes += size
            self._evict_locked(self._max_bytes)

    def clear(self) -> None:
        """Remove all entries. The settings do not change."""
        with self._lock:
            self._clear_locked()

    def invalidate(self) -> None:
        """Remove all entries because their source data changed."""
        self.clear()

    def _discard_locked(self, key: Hashable) -> None:
        old = self._entries.pop(key, None)
        if old is not None:
            self._retained_bytes -= int(old.nbytes)

    def _clear_locked(self) -> None:
        self._entries.clear()
        self._retained_bytes = 0

    def _evict_locked(self, limit: int) -> None:
        while self._retained_bytes > limit and self._entries:
            _, old = self._entries.popitem(last=False)
            self._retained_bytes -= int(old.nbytes)
