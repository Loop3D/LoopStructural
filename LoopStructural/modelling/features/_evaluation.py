"""Evaluation scope that shares fault restorations inside and between calls.

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

from ._evaluation_cache import dependency_fingerprint

if TYPE_CHECKING:
    from LoopStructural import GeologicalModel
    from LoopStructural.modelling.features._base_geological_feature import BaseFeature
    from LoopStructural.modelling.features._evaluation_cache import RestorationCache

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
        self._fingerprints: dict = {}
        self._published: list = []
        self._opaque = False
        self._model_state: tuple | None = None

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
        digest = _digest(points)
        # The scope key holds the faults only. Features with the same ordered faults share
        # one restoration inside the scope. The model does not change inside one scope.
        key = (
            digest,
            points.shape,
            points.dtype.str,
            tuple(_fault_state(fault) for fault in faults),
            bool(reverse),
            True,
        )
        restored = self._restorations.get(key)
        if restored is not None:
            return restored
        # The persistent key holds the full dependency state of the feature. An opaque
        # state or a cache that is off gives no persistent key.
        cache = self._persistent_cache()
        fingerprint = self._fingerprint(feature, cache)
        persistent_key = None
        if fingerprint is not None:
            persistent_key = (
                "persistent",
                digest,
                points.shape,
                points.dtype.str,
                bool(reverse),
                fingerprint,
                self._state_of_model(cache),
            )
            restored = cache.get(persistent_key)
            if restored is not None:
                self._restorations[key] = restored
                return restored
        restored = points
        for fault in faults:
            restored = fault.apply_to_points(restored, reverse=reverse)
        restored = np.asarray(restored)
        if restored is points or np.shares_memory(restored, points):
            restored = restored.copy()
        restored.flags.writeable = False
        if persistent_key is not None:
            # The cache keeps one private read-only copy. The scope shares that copy,
            # so the result is retained once. A result that is too large stays in the
            # scope only.
            cache.put(persistent_key, restored)
            stored = cache.get(persistent_key)
            if stored is not None:
                restored = stored
                self._published.append(persistent_key)
        # The entry is stored only after all faults are applied.
        self._restorations[key] = restored
        return restored

    def _persistent_cache(self) -> RestorationCache | None:
        """Return the cache of the model when it can keep entries, else None.

        After the scope meets an opaque dependency, it uses no persistent entry.
        """
        if self._opaque:
            return None
        cache = getattr(self.model, "evaluation_cache", None)
        if cache is None or not cache.enabled or cache.max_bytes == 0:
            return None
        return cache

    def _fingerprint(self, feature: BaseFeature, cache: RestorationCache | None) -> bytes | None:
        """Return the dependency fingerprint of ``feature``, once for each scope.

        The function returns None when the cache is off or a dependency is opaque.
        The lazy build of dirty builders runs under the lock of the cache, so it
        happens once.
        """
        if cache is None:
            return None
        known = self._fingerprints.get(id(feature))
        if known is not None:
            return known[1]
        with cache.lock:
            fingerprint = dependency_fingerprint(feature)
        # The scope holds the feature, so its id stays valid until the scope ends.
        self._fingerprints[id(feature)] = (feature, fingerprint)
        if fingerprint is None:
            self._opaque = True
        return fingerprint

    def _state_of_model(self, cache: RestorationCache) -> bytes:
        """Return a digest of the state of all features of the model, once for each scope.

        An entry is valid only while the whole model has the state in which it was made.
        A change in one feature, also an edit in place that no setter sees, therefore
        makes all older entries unreachable. The function builds nothing, so it does not
        change which builders are built. A feature that cannot be read makes the digest
        None, and the scope then uses no persistent entry.
        """
        if self._model_state is None:
            hasher = hashlib.sha256()
            with cache.lock:
                for item in getattr(self.model, "features", None) or ():
                    fingerprint = dependency_fingerprint(item, build=False)
                    hasher.update(fingerprint if fingerprint is not None else b"opaque")
            self._model_state = (hasher.digest(),)
        return self._model_state[0]

    def _withdraw(self) -> None:
        """Remove the entries that this scope published, after an exception.

        A call that fails must not leave the result of an unfinished evaluation in
        the cache of the model.
        """
        cache = getattr(self.model, "evaluation_cache", None)
        if cache is not None:
            for key in self._published:
                cache.discard(key)
        self._published.clear()


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
    except BaseException:
        context._withdraw()
        raise
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
