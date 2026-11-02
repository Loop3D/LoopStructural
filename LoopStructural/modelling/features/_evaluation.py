"""Evaluation scope that shares fault restorations inside one top-level call.

An evaluation scope is a short-lived :class:`EvaluationContext`. It stores the restored
coordinates of the faults of a feature. The key of an entry holds the exact content of the
points, the ordered faults, ``reverse``, ``faults_enabled`` and a state token of the faults.
Two features that use the same ordered faults on equal points therefore share one
restoration.

The scope is stored in a ``contextvars.ContextVar``. It is always reset when the
outermost ``with`` block ends, also after an exception. The scope does not keep any
entry after it ends. A partial restoration is never stored.
"""

from __future__ import annotations

import hashlib
from collections.abc import Hashable, Iterator
from contextlib import contextmanager
from contextvars import ContextVar
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from LoopStructural import GeologicalModel
    from LoopStructural.modelling.features._base_geological_feature import BaseFeature

_CHUNK_BYTES = 1 << 20
_ACTIVE: ContextVar = ContextVar("loopstructural_evaluation_context", default=None)


def _digest(points: np.ndarray) -> bytes:
    """Return a digest of the exact content of ``points``.

    The digest is made in chunks of bounded size, so no large copy is made for a
    contiguous array.
    """
    array = np.ascontiguousarray(points)
    hasher = hashlib.blake2b(digest_size=16)
    flat = array.reshape(-1).view(np.uint8) if array.size else array.reshape(-1)
    for start in range(0, flat.size, _CHUNK_BYTES):
        hasher.update(flat[start : start + _CHUNK_BYTES])
    return hasher.digest()


def _fault_state(fault: BaseFeature) -> Hashable:
    """Return a token for the identity and the state of one fault."""
    return (id(fault), getattr(fault, "displacement", None), getattr(fault, "steps", None))


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
        key = (
            _digest(points),
            points.shape,
            points.dtype.str,
            tuple(_fault_state(fault) for fault in faults),
            bool(reverse),
            True,
        )
        restored = self._restorations.get(key)
        if restored is not None:
            return restored
        restored = points
        for fault in faults:
            restored = fault.apply_to_points(restored, reverse=reverse)
        restored = np.asarray(restored)
        if restored is points or np.shares_memory(restored, points):
            restored = restored.copy()
        restored.flags.writeable = False
        # The entry is stored only after all faults are applied.
        self._restorations[key] = restored
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
