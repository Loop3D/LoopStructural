"""Evaluation scope that shares fault restorations inside one call.

An evaluation scope is a short-lived :class:`EvaluationContext`. It stores the restored
coordinates of the faults of a feature. The key of an entry holds the ordered faults,
``reverse``, a state token of the faults and the shape of the points. The entry also holds
the exact content of the points. Two features that use the same ordered faults on equal
points therefore share one restoration. Each key keeps only its last point batch, so the
scope holds at most one restoration for each fault set and batch shape (FR-010, OP1).

The scope is stored in a ``contextvars.ContextVar``. It is always reset when the
outermost ``with`` block ends, also after an exception. The scope keeps
its entries for one call only. It does not keep any entry after it ends. A partial restoration is never stored.
"""

from __future__ import annotations

import hashlib
from collections.abc import Hashable, Iterator
from contextlib import contextmanager
from contextvars import ContextVar
from typing import TYPE_CHECKING

import numpy as np

from LoopStructural.utils._fault_block import restore_through

if TYPE_CHECKING:
    from LoopStructural import GeologicalModel
    from LoopStructural.modelling.features._base_geological_feature import BaseFeature

_CHUNK_BYTES = 1 << 20
_ACTIVE: ContextVar = ContextVar("loopstructural_evaluation_context", default=None)


def _digest(points: np.ndarray) -> bytes:
    """Return a digest of the exact content of ``points``.

    The digest is made in blocks of bounded size. A contiguous array is not copied.
    A non-contiguous array is copied one block of rows at a time, never as a whole.
    The digest is the same for equal content in any memory layout.
    """
    array = np.asarray(points)
    hasher = hashlib.blake2b(digest_size=16)
    if array.size == 0:
        return hasher.digest()
    if array.ndim == 0:
        array = array.reshape(1)
    if array.flags.c_contiguous:
        blocks = [array]
    else:
        row_bytes = max(1, array.itemsize * int(np.prod(array.shape[1:], dtype=np.int64)))
        rows = max(1, _CHUNK_BYTES // row_bytes)
        blocks = (array[start : start + rows] for start in range(0, array.shape[0], rows))
    for block in blocks:
        flat = np.ascontiguousarray(block).reshape(-1).view(np.uint8)
        for start in range(0, flat.size, _CHUNK_BYTES):
            hasher.update(flat[start : start + _CHUNK_BYTES])
    return hasher.digest()


def _hashable(value: object) -> Hashable:
    """Return the value, or a token of its bytes if it cannot be a dict key.

    A displacement can be a NumPy array, for example from ``df[...].values``.
    """
    try:
        hash(value)
    except TypeError:
        array = np.asarray(value)
        return (array.dtype.str, array.shape, array.tobytes())
    return value


def _fault_state(fault: BaseFeature) -> Hashable:
    """Return a token for the identity and the state of one fault."""
    return (
        id(fault),
        _hashable(getattr(fault, "displacement", None)),
        _hashable(getattr(fault, "steps", None)),
    )


class EvaluationContext:
    """Storage for the restorations of one evaluation scope.

    Parameters
    ----------
    model : GeologicalModel, optional
        The model of the scope. None for a standalone feature.
    """

    def __init__(self, model: GeologicalModel | None = None) -> None:
        self.model = model
        self._restorations: dict = {}
        self._batches: dict = {}

    def restore(
        self, feature: BaseFeature, points: np.ndarray, reverse: bool = False
    ) -> np.ndarray:
        """Return the points restored by the ordered faults of ``feature``.

        Parameters
        ----------
        feature : BaseFeature
            The feature whose ``faults`` are applied.
        points : np.ndarray
            Nx3 array of points.
        reverse : bool, default False
            Apply the faults in the reverse direction.

        Returns
        -------
        np.ndarray
            Read-only array of the restored points.
        """
        points = np.asarray(points)
        faults = tuple(feature.faults) if feature.faults_enabled else ()
        if not faults:
            unchanged = points.view()
            unchanged.flags.writeable = False
            return unchanged
        # Build every reachable builder before the first lookup, so a lazy build
        # happens once and not between two lookups.
        for fault in faults:
            builder = getattr(fault, "builder", None)
            if builder is not None:
                builder.up_to_date()
        if points.dtype.hasobject:
            # An array of Python objects has no bytes to compare. Restore it without
            # sharing, as the legacy path does (apply_to_points makes it float).
            return restore_through(points, faults, reverse=reverse)
        digest = _digest(points)
        # The scope key holds the faults only. Features with the same ordered faults share
        # one restoration inside the scope. The model does not change inside one scope.
        # One slot for each fault set and batch shape, so the N points and the 4N
        # tetrahedron corners of a gradient do not push each other out.
        fault_set = (tuple(_fault_state(fault) for fault in faults), bool(reverse), points.shape)
        batch = (digest, points.dtype.str)
        if self._batches.get(fault_set) == batch:
            return self._restorations[fault_set]
        # Free the old batch of this slot before the new one is made.
        self._restorations.pop(fault_set, None)
        self._batches.pop(fault_set, None)
        restored = np.asarray(restore_through(points, faults, reverse=reverse))
        if restored is points or np.shares_memory(restored, points):
            restored = restored.copy()
        restored.flags.writeable = False
        # The entry is stored only after all faults are applied. Each slot keeps only
        # its last point batch (FR-010). The batches of a fault frame move with each
        # step, so an older batch is not asked for again.
        self._restorations[fault_set] = restored
        self._batches[fault_set] = batch
        return restored


@contextmanager
def evaluation_scope(model: GeologicalModel | None = None) -> Iterator[EvaluationContext]:
    """Open an evaluation scope and give its context.

    A nested scope uses the active context when ``model`` is the model of the active
    context, or when ``model`` is None. A scope for another model is a separate context.
    The ContextVar is reset when the scope that opened the context ends.

    Parameters
    ----------
    model : GeologicalModel, optional
        The model of the scope, by default None.

    Yields
    ------
    EvaluationContext
        The active context.
    """
    active = _ACTIVE.get()
    if active is not None and (model is None or active.model is model):
        yield active
        return
    context = EvaluationContext(model)
    token = _ACTIVE.set(context)
    try:
        yield context
    finally:
        _ACTIVE.reset(token)


def restore(feature: BaseFeature, points: np.ndarray, reverse: bool = False) -> np.ndarray:
    """Restore ``points`` with the ordered faults of ``feature``, shared in the scope.

    Parameters
    ----------
    feature : BaseFeature
        The feature whose faults are applied.
    points : np.ndarray
        Nx3 array of points.
    reverse : bool, default False
        Apply the faults in the reverse direction.

    Returns
    -------
    np.ndarray
        Read-only array of the restored points.
    """
    with evaluation_scope(getattr(feature, "model", None)) as context:
        return context.restore(feature, points, reverse)
