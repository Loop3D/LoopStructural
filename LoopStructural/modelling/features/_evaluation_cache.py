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

Invalidation
------------
Two guards protect the stored restorations (plan D4).

* ``GeologicalModel.evaluation_revision`` is a monotonic counter. :func:`bump_revision`
  increases it. Setters and methods that change a fault, a builder, an interpolator or
  the model call it.
* :func:`dependency_fingerprint` is a SHA256 digest of the current state of all objects
  that the restoration of the faults of a feature depends on. It finds the changes that
  no setter can see: in-place edits of lists, dictionaries, DataFrames and arrays.
"""

from __future__ import annotations

import hashlib
import threading
from collections import OrderedDict
from collections.abc import Hashable
from typing import Protocol

import numpy as np

DEFAULT_MAX_BYTES = 64 * 1024 * 1024
_CHUNK_BYTES = 1 << 20


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


# ---------------------------------------------------------------------------
# Revision hooks
# ---------------------------------------------------------------------------


def _model_of(owner: object) -> object | None:
    """Return the model of a feature, builder or interpolator owner, or None.

    A builder can have no model of its own. It then uses the model of its feature. A
    structural frame can have no model of its own. It then uses the model of its members.
    """
    for _ in range(3):
        if owner is None:
            return None
        model = getattr(owner, "_model", None)
        if model is not None:
            return model
        for member in getattr(owner, "features", None) or ():
            model = getattr(member, "_model", None)
            if model is not None:
                return model
        owner = getattr(owner, "_feature", None)
    return None


def bump_revision(owner: object) -> None:
    """Increase the evaluation revision of the model of ``owner``.

    The function does nothing when ``owner`` has no model yet.
    """
    model = _model_of(owner)
    bump = getattr(model, "_bump_evaluation_revision", None)
    if bump is not None:
        bump()


def values_differ(old: object, new: object) -> bool:
    """Return True when ``new`` is not equal to ``old``, also for array-like values."""
    try:
        return bool(old != new)
    except (TypeError, ValueError):
        return True


# ---------------------------------------------------------------------------
# Dependency fingerprint
# ---------------------------------------------------------------------------


class _Updatable(Protocol):
    """The part of a ``hashlib`` object that the fingerprint uses."""

    def update(self, data: bytes | np.ndarray, /) -> None:
        """Add bytes to the digest."""


class _Opaque(Exception):
    """A dependency has a state that the fingerprint cannot read."""


def _update_array(hasher: _Updatable, array: object) -> None:
    """Add the dtype, the shape and the bytes of ``array`` to ``hasher`` in chunks.

    A C-contiguous array is read through a byte view. Any other array is read in blocks
    of rows. No copy of the whole array is made.
    """
    array = np.asarray(array)
    if array.dtype.hasobject:
        raise _Opaque
    hasher.update(f"A{array.dtype.str}{array.shape};".encode())
    if array.size == 0:
        return
    if array.flags.c_contiguous:
        flat = array.reshape(-1).view(np.uint8)
        for start in range(0, flat.size, _CHUNK_BYTES):
            hasher.update(flat[start : start + _CHUNK_BYTES])
        return
    row_bytes = array[0].nbytes if array.ndim > 1 else array.itemsize
    rows = max(1, _CHUNK_BYTES // max(1, row_bytes))
    for start in range(0, array.shape[0], rows):
        block = np.ascontiguousarray(array[start : start + rows])
        hasher.update(block.reshape(-1).view(np.uint8))


def _encode(hasher: _Updatable, value: object) -> None:
    """Add a plain value (scalar, text, array, list, tuple or dict) to ``hasher``."""
    if value is None:
        hasher.update(b"N;")
    elif isinstance(value, (bool, np.bool_)):
        hasher.update(b"B1;" if value else b"B0;")
    elif isinstance(value, (int, np.integer)):
        hasher.update(f"I{int(value)};".encode())
    elif isinstance(value, (float, np.floating)):
        hasher.update(f"F{float(value).hex()};".encode())
    elif isinstance(value, (str, bytes)):
        raw = value.encode() if isinstance(value, str) else value
        hasher.update(f"S{len(raw)}:".encode())
        hasher.update(raw)
    elif isinstance(value, np.ndarray):
        _update_array(hasher, value)
    elif isinstance(value, (list, tuple)):
        hasher.update(f"L{len(value)};".encode())
        for item in value:
            _encode(hasher, item)
    elif isinstance(value, dict):
        hasher.update(f"M{len(value)};".encode())
        for key in sorted(value, key=repr):
            _encode(hasher, key)
            _encode(hasher, value[key])
    else:
        raise _Opaque


def _update_frame(hasher: _Updatable, frame: object) -> None:
    """Add the column names, dtypes and values of a DataFrame to ``hasher``."""
    import pandas as pd

    if not isinstance(frame, pd.DataFrame):
        raise _Opaque
    hasher.update(f"D{frame.shape};".encode())
    for position in range(frame.shape[1]):
        column = frame.iloc[:, position]
        _encode(hasher, repr(frame.columns[position]))
        values = column.to_numpy()
        if values.dtype.hasobject:
            values = pd.util.hash_pandas_object(column, index=False).to_numpy()
        _update_array(hasher, values)


def _update_profile(hasher: _Updatable, function: object) -> None:
    """Add the state of a fault displacement function to ``hasher``.

    Only the known function classes of LoopStructural are read. Any other callable can
    use state that is not visible, so it is opaque.
    """
    from .fault._fault_function import (
        Composite,
        CubicFunction,
        FaultDisplacement,
        Ones,
        Zeros,
        smooth_peak,
    )

    if function is smooth_peak:
        hasher.update(b"Psmooth_peak;")
    elif type(function) is FaultDisplacement:
        hasher.update(b"PFaultDisplacement;")
        _encode(hasher, function.scale)
        for part in (function.gx, function.gy, function.gz):
            _update_profile(hasher, part)
    elif type(function) is Composite:
        hasher.update(b"PComposite;")
        _update_profile(hasher, function.positive)
        _update_profile(hasher, function.negative)
    elif type(function) is CubicFunction:
        hasher.update(b"PCubic;")
        for name in ("A", "B", "lim", "max_v", "min_v"):
            _encode(hasher, getattr(function, name))
        # ``w`` is the cached solution of A and B. It is the only source without A.
        if len(function.A) == 0:
            _encode(hasher, np.asarray(function.w))
    elif type(function) in (Ones, Zeros):
        hasher.update(f"P{type(function).__name__};".encode())
    else:
        raise _Opaque


def _update_support(hasher: _Updatable, support: object) -> None:
    """Add the geometry of an interpolation support to ``hasher``."""
    from loop_common.supports import UnStructuredTetMesh
    from loop_common.supports._3d_base_structured import BaseStructuredSupport

    if isinstance(support, BaseStructuredSupport):
        hasher.update(f"G{type(support).__name__};".encode())
        for name in ("origin", "nsteps", "step_vector", "rotation_xy"):
            _encode(hasher, np.asarray(getattr(support, name)))
    elif isinstance(support, UnStructuredTetMesh):
        hasher.update(b"GUnStructuredTetMesh;")
        for name in ("nodes", "elements", "neighbours"):
            _encode(hasher, np.asarray(getattr(support, name)))
    else:
        raise _Opaque


def _update_interpolator(hasher: _Updatable, interpolator: object) -> None:
    """Add the coefficients and the support of a discrete interpolator to ``hasher``."""
    from loop_interpolation._discrete_interpolator import DiscreteInterpolator

    if not isinstance(interpolator, DiscreteInterpolator):
        raise _Opaque
    coefficients = getattr(interpolator, "c", None)
    if not isinstance(coefficients, np.ndarray):
        raise _Opaque
    hasher.update(f"C{type(interpolator).__name__};".encode())
    _update_array(hasher, coefficients)
    _update_support(hasher, getattr(interpolator, "support", None))


def _known_feature_types() -> tuple:
    """Return the feature classes whose complete state the fingerprint reads."""
    from ._geological_feature import GeologicalFeature
    from ._structural_frame import StructuralFrame
    from ._unconformity_feature import UnconformityFeature
    from .fault import FaultSegment
    from .fold import FoldFrame

    return (GeologicalFeature, UnconformityFeature, StructuralFrame, FaultSegment, FoldFrame)


def _is_feature(value: object) -> bool:
    from ._base_geological_feature import BaseFeature

    return isinstance(value, BaseFeature)


def _region_feature(region: object) -> object | None:
    """Return the feature that a region reads, or None."""
    if _is_feature(region):
        return region
    feature = getattr(region, "feature", None)
    return feature if _is_feature(feature) else None


def _dependencies(node: object) -> list:
    """Return the features that ``node`` reaches directly, in a fixed order."""
    result = list(getattr(node, "faults", None) or [])
    members = getattr(node, "features", None)
    if isinstance(members, (list, tuple)):
        result.extend(member for member in members if member is not None)
    for region in getattr(node, "regions", None) or []:
        feature = _region_feature(region)
        if feature is not None:
            result.append(feature)
        parent = getattr(region, "parent", None)
        if _is_feature(parent):
            result.append(parent)
    parent = getattr(node, "parent", None)
    if _is_feature(parent):
        result.append(parent)
    return [item for item in result if _is_feature(item)]


def _build_reachable(root: object) -> None:
    """Build every dirty builder that ``root`` reaches, once, before any state is read.

    The faults of a lazy feature exist only after its builder has run. A node is
    therefore built before its dependencies are listed.
    """
    seen: set = set()
    builders: set = set()
    stack = [root]
    while stack:
        node = stack.pop()
        if id(node) in seen:
            continue
        seen.add(id(node))
        builder = getattr(node, "builder", None)
        if builder is not None and id(builder) not in builders:
            builders.add(id(builder))
            up_to_date = getattr(builder, "up_to_date", None)
            if callable(up_to_date):
                up_to_date()
        stack.extend(reversed(_dependencies(node)))


class _StateHasher:
    """Write the state of a feature graph into one SHA256 digest.

    A node is written once for each mode (shallow and deep). A second visit in the same
    mode writes only the index of the first visit, so a cycle ends. A deep visit of a
    node that has only a shallow record writes the full state. The digest holds no object id, only the state and the visit order.
    """

    def __init__(self) -> None:
        self.hasher = hashlib.sha256()
        self._visited: dict = {}

    def node(self, node: object, deep: bool) -> None:
        """Write one feature.

        ``deep`` is False for the root. The restoration of the faults of the root does not
        depend on the field of the root, so its coefficients, support and displacement
        are not written. Its data, build arguments and regions are written.
        """
        hasher = self.hasher
        # A deep record holds all of the state, so it also answers a shallow visit.
        # A shallow record leaves out the field, so a deep visit must write it in full.
        for key in ((id(node), True), (id(node), deep)):
            if key in self._visited:
                hasher.update(f"R{self._visited[key]};".encode())
                return
        self._visited[(id(node), deep)] = len(self._visited)
        kind = type(node)
        if deep and kind not in _known_feature_types():
            raise _Opaque
        hasher.update(f"({kind.__module__}.{kind.__qualname__};".encode())
        _encode(hasher, getattr(node, "faults_enabled", True))
        for name in ("value", "sign"):  # unconformity threshold and side
            if hasattr(node, name):
                _encode(hasher, getattr(node, name))
        self._builder(getattr(node, "builder", None))
        faults = list(getattr(node, "faults", None) or [])
        hasher.update(f"F{len(faults)};".encode())
        for fault in faults:
            if not _is_feature(fault):
                raise _Opaque
            self.node(fault, True)
        regions = list(getattr(node, "regions", None) or [])
        hasher.update(f"G{len(regions)};".encode())
        for region in regions:
            self._region(region)
        members = getattr(node, "features", None)
        if members is not None:
            if getattr(node, "fold", None) is not None:
                raise _Opaque
            hasher.update(f"M{len(members)};".encode())
            for member in members:
                if member is None:
                    hasher.update(b"N;")
                else:
                    self.node(member, deep)
        parent = getattr(node, "parent", None)
        if _is_feature(parent):
            self.node(parent, True)
        if deep:
            self._deep_state(node)
        hasher.update(b");")

    def _deep_state(self, node: object) -> None:
        """Write the state that only a dependency of the root contributes."""
        from .fault import FaultSegment

        if isinstance(node, FaultSegment):
            for name in ("displacement", "steps", "fault_offset"):
                _encode(self.hasher, getattr(node, name))
            _update_profile(self.hasher, node.faultfunction)
        if not isinstance(getattr(node, "features", None), (list, tuple)):
            _update_interpolator(self.hasher, getattr(node, "interpolator", None))

    def _builder(self, builder: object) -> None:
        """Write the data and the build arguments of a builder."""
        if builder is None:
            return
        # A frame builder owns member builders. The member features write their own.
        if getattr(builder, "builders", None) is None:
            data = getattr(builder, "data", None)
            if data is not None:
                _update_frame(self.hasher, data)
        arguments = getattr(builder, "build_arguments", None)
        if arguments is not None:
            _encode(self.hasher, arguments)

    def _region(self, region: object) -> None:
        """Write one region. A plain callable is opaque: its closure state is unknown."""
        from ...utils.regions import NegativeRegion, PositiveRegion, RegionEverywhere
        from ._region import Region

        hasher = self.hasher
        if _is_feature(region):
            hasher.update(b"Gfeature;")
            self.node(region, True)
        elif type(region) in (PositiveRegion, NegativeRegion):
            hasher.update(f"G{type(region).__name__};".encode())
            _encode(hasher, region.vector)
            _encode(hasher, region.point)
            self.node(region.feature, True)
            if _is_feature(region.parent):
                self.node(region.parent, True)
        elif type(region) is Region:
            hasher.update(b"GRegion;")
            _encode(hasher, region.value)
            _encode(hasher, region.sign)
            self.node(region.feature, True)
        elif type(region) is RegionEverywhere:
            hasher.update(b"GEverywhere;")
        else:
            raise _Opaque


def dependency_fingerprint(feature: object) -> bytes | None:
    """Return a SHA256 digest of the state that the faults of ``feature`` depend on.

    The digest covers the ordered faults, the frames of the faults, the regions and
    their parents, the data and the build arguments of the builders, the coefficients
    and the geometry of the supports, and the displacement, the steps and the function
    of each fault. It holds no revision and no object id.

    A dirty builder that the feature reaches is built first, once. The function reads
    the model and changes nothing else. Arrays are hashed in chunks, with no copy of
    a whole array. A cycle in the graph ends the walk.

    Parameters
    ----------
    feature : BaseFeature
        The feature whose restoration is stored.

    Returns
    -------
    bytes or None
        A 32-byte digest, or None when a dependency is opaque (for example a plain
        function as a region, a custom interpolator or an unknown support). Persistent
        entries must then be bypassed.
    """
    _build_reachable(feature)
    state = _StateHasher()
    try:
        state.node(feature, False)
    except _Opaque:
        return None
    return state.hasher.digest()
