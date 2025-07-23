from __future__ import annotations

from collections.abc import Callable
from contextlib import contextmanager
from typing import Any, Generic, Protocol, TypeAlias, TypeVar, runtime_checkable
import threading
import weakref

__all__ = ["Observer", "Observable", "Disposable"]


@runtime_checkable
class Observer(Protocol):
    """Objects implementing an *update* method can subscribe."""

    def update(self, observable: "Observable", event: str, *args: Any, **kwargs: Any) -> None:  # noqa: D401,E501
        """Receive a notification."""


Callback: TypeAlias = Callable[["Observable", str, Any], None]
T = TypeVar("T", bound="Observable")


class Disposable:
    """A small helper that detaches an observer when disposed."""

    __slots__ = ("_detach",)

    def __init__(self, detach: Callable[[], None]):
        self._detach = detach

    def dispose(self) -> None:  # noqa: D401
        """Detach the associated observer immediately."""

        self._detach()

    # Allow use as a context‑manager for temporary subscriptions
    def __enter__(self) -> "Disposable":
        return self

    def __exit__(self, exc_type, exc, tb):  # noqa: D401
        self.dispose()
        return False  # do not swallow exceptions


class Observable(Generic[T]):
    """Base‑class that provides Observer pattern plumbing."""

    #: Internal storage: mapping *event* → WeakSet[Callback]
    _observers: dict[str, weakref.WeakSet[Callback]]
    _any_observers: weakref.WeakSet[Callback]
    
    def __init__(self) -> None:
        self._lock = threading.RLock()
        self._observers = {}
        self._any_observers = weakref.WeakSet()
        self._frozen = 0
        self._pending: list[tuple[str, tuple[Any, ...], dict[str, Any]]] = []

    # ‑‑‑ subscription api --------------------------------------------------
    def attach(self, listener: Observer | Callback, event: str | None = None) -> Disposable:  # noqa: D401,E501
        """Register *listener* for *event* (all events if *event* is None).

        Returns a :class:`Disposable` so the caller can easily detach again.
        """
        callback: Callback = (
            listener.update  # type: ignore[attr‑defined]
            if isinstance(listener, Observer)  # type: ignore[misc]
            else listener  # already a callable
        )

        with self._lock:
            if event is None:
                self._any_observers.add(callback)
            else:
                self._observers.setdefault(event, weakref.WeakSet()).add(callback)

        return Disposable(lambda: self.detach(listener, event))

    def detach(self, listener: Observer | Callback, event: str | None = None) -> None:  # noqa: D401,E501
        """Unregister a previously attached *listener*."""

        callback: Callback = (
            listener.update  # type: ignore[attr‑defined]
            if isinstance(listener, Observer)  # type: ignore[misc]
            else listener
        )

        with self._lock:
            if event is None:
                self._any_observers.discard(callback)
                for s in self._observers.values():
                    s.discard(callback)
            else:
                self._observers.get(event, weakref.WeakSet()).discard(callback)

    # ‑‑‑ notification api --------------------------------------------------
    def notify(self: T, event: str, *args: Any, **kwargs: Any) -> None:  # noqa: D401,E501
        """Notify observers that *event* happened."""

        with self._lock:
            if self._frozen:
                # defer until freeze_notifications() exits
                self._pending.append((event, args, kwargs))
                return

            observers = list(self._any_observers)
            observers.extend(self._observers.get(event, ()))

        # Call outside lock — prevent deadlocks if observers trigger other
        # notifications.
        for cb in observers:
            try:
                cb(self, event, *args, **kwargs)
            except Exception:  # pragma: no cover
                # Optionally log; never allow an observer error to break flow.
                import logging

                logging.getLogger(__name__).exception(
                    "Unhandled error in observer %s for event %s", cb, event
                )

    # ‑‑‑ batching ----------------------------------------------------------
    @contextmanager
    def freeze_notifications(self):  # noqa: D401
        """Context manager that batches notifications until exit."""

        with self._lock:
            self._frozen += 1
        try:
            yield self
        finally:
            with self._lock:
                self._frozen -= 1
                if self._frozen == 0 and self._pending:
                    pending = self._pending[:]
                    self._pending.clear()
            for event, args, kw in pending:  # type: ignore[has‑type]
                self.notify(event, *args, **kw)
