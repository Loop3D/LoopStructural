from __future__ import annotations

from collections.abc import Callable
from contextlib import contextmanager
from typing import Any, Generic, Protocol, TypeVar, runtime_checkable
import threading
import weakref

__all__ = ["Observer", "Observable", "Disposable"]


@runtime_checkable
class Observer(Protocol):
    """Protocol for objects that can observe events from Observable objects.

    Classes implementing this protocol must provide an update method that
    will be called when observed events occur.
    """

    def update(self, observable: "Observable", event: str, *args: Any, **kwargs: Any) -> None:
        """Receive a notification from an observable object.

        Parameters
        ----------
        observable : Observable
            The observable object that triggered the event
        event : str
            The name of the event that occurred
        *args : Any
            Positional arguments associated with the event
        **kwargs : Any
            Keyword arguments associated with the event
        """


Callback = Callable[["Observable", str, Any], None]
T = TypeVar("T", bound="Observable")


class Disposable:
    """A helper class that manages detachment of observers.

    This class provides a convenient way to detach observers from observables.
    It can be used as a context manager for temporary subscriptions.

    Parameters
    ----------
    detach : Callable[[], None]
        Function to call when disposing of the observer
    """

    __slots__ = ("_detach",)

    def __init__(self, detach: Callable[[], None]):
        self._detach = detach

    def dispose(self) -> None:  
        """Detach the associated observer immediately."""

        self._detach()

    # Allow use as a context‑manager for temporary subscriptions
    def __enter__(self) -> "Disposable":
        return self

    def __exit__(self, exc_type, exc, tb):  
        self.dispose()
        return False  # do not swallow exceptions


class Observable(Generic[T]):
    """Base class that implements the Observer pattern.

    This class provides the infrastructure for managing observers and 
    notifying them of events. Observers can be attached to specific events
    or to all events.

    Attributes
    ----------
    _observers : dict[str, weakref.WeakSet[Callback]]
        Internal storage mapping event names to sets of callbacks
    _any_observers : weakref.WeakSet[Callback]
        Set of callbacks that listen to all events
    """

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
    def attach(self, listener: Observer | Callback, event: str | None = None) -> Disposable:  
        """Register a listener for specific event or all events.

        Parameters
        ----------
        listener : Observer | Callback
            The observer object or callback function to attach
        event : str | None, optional
            The specific event to listen for. If None, listens to all events, by default None

        Returns
        -------
        Disposable
            A disposable object that can be used to detach the listener
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

    def detach(self, listener: Observer | Callback, event: str | None = None) -> None:  
        """Unregister a previously attached listener.

        Parameters
        ----------
        listener : Observer | Callback
            The observer object or callback function to detach
        event : str | None, optional
            The specific event to stop listening for. If None, detaches from all events, by default None
        """

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
    def __getstate__(self):
        """Prepare object state for pickling by removing unpicklable attributes.

        Returns
        -------
        dict
            Object state dictionary with thread locks and weak references removed
        """
        state = self.__dict__.copy()
        state.pop('_lock', None)  # RLock cannot be pickled
        state.pop('_observers', None)  # WeakSet cannot be pickled
        state.pop('_any_observers', None)
        return state
    
    def __setstate__(self, state):
        """Restore object state after unpickling and reinitialize locks and observers.

        Parameters
        ----------
        state : dict
            The restored object state dictionary
        """
        self.__dict__.update(state)
        self._lock = threading.RLock()
        self._observers = {}
        self._any_observers = weakref.WeakSet()
        self._frozen = 0
    # ‑‑‑ notification api --------------------------------------------------
    def notify(self: T, event: str, *args: Any, **kwargs: Any) -> None:  
        """Notify all observers that an event has occurred.

        Parameters
        ----------
        event : str
            The name of the event that occurred
        *args : Any
            Positional arguments to pass to the observers
        **kwargs : Any
            Keyword arguments to pass to the observers
        """

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
    def freeze_notifications(self):  
        """Context manager that batches notifications until exit.

        While in this context, notifications are queued rather than sent
        immediately. When the context exits, all queued notifications are
        sent in order.

        Yields
        ------
        Observable
            Self reference for method chaining
        """

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
