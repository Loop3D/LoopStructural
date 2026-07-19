import gc
import pickle

import pytest

from LoopStructural.utils.observer import Observable, Disposable


class Recorder:
    """Simple Observer implementation used across tests.

    Implements the `update` method required by the Observer protocol and
    just records every call it receives so tests can assert on them.
    """

    def __init__(self):
        self.calls = []

    def update(self, observable, event, *args, **kwargs):
        self.calls.append((observable, event, args, kwargs))


def test_attach_callback_and_notify():
    obs = Observable()
    received = []

    def callback(observable, event, *args, **kwargs):
        received.append((observable, event, args, kwargs))

    obs.attach(callback)
    obs.notify("changed", 1, 2, key="value")

    assert len(received) == 1
    assert received[0][0] is obs
    assert received[0][1] == "changed"
    assert received[0][2] == (1, 2)
    assert received[0][3] == {"key": "value"}


def test_attach_observer_object_receives_notifications():
    """`attach()` tracks bound-method listeners (e.g. `listener.update`) via
    `weakref.WeakMethod`, which is keyed on the owning instance rather than
    the transient bound-method wrapper object. As long as the observer
    object itself (`recorder`) is kept alive, it continues to receive
    notifications through the "Observer protocol" pattern (attaching an
    object that implements `update`).
    """
    obs = Observable()
    recorder = Recorder()

    obs.attach(recorder)
    gc.collect()
    obs.notify("event_a")

    assert len(recorder.calls) == 1
    assert recorder.calls[0][1] == "event_a"


def test_attach_specific_event_only_triggers_for_that_event():
    obs = Observable()
    calls = []

    def callback(observable, event, *args, **kwargs):
        calls.append(event)

    obs.attach(callback, event="specific")
    obs.notify("other")
    obs.notify("specific")

    assert calls == ["specific"]


def test_detach_removes_listener():
    obs = Observable()
    calls = []

    def callback(observable, event, *args, **kwargs):
        calls.append(event)

    obs.attach(callback)
    obs.notify("first")
    obs.detach(callback)
    obs.notify("second")

    assert calls == ["first"]


def test_detach_event_specific_listener():
    obs = Observable()
    recorder = Recorder()

    obs.attach(recorder, event="my_event")
    obs.detach(recorder, event="my_event")
    obs.notify("my_event")

    assert recorder.calls == []


def test_attach_returns_disposable_that_detaches():
    obs = Observable()
    calls = []

    def callback(observable, event, *args, **kwargs):
        calls.append(event)

    disposable = obs.attach(callback)
    assert isinstance(disposable, Disposable)
    obs.notify("first")
    disposable.dispose()
    obs.notify("second")

    assert calls == ["first"]


def test_disposable_as_context_manager_detaches_on_exit():
    obs = Observable()
    calls = []

    def callback(observable, event, *args, **kwargs):
        calls.append(event)

    with obs.attach(callback) as disposable:
        assert isinstance(disposable, Disposable)
        obs.notify("inside")

    obs.notify("outside")

    assert calls == ["inside"]


def test_disposable_context_manager_does_not_swallow_exceptions():
    obs = Observable()
    recorder = Recorder()

    with pytest.raises(ValueError):
        with obs.attach(recorder):
            raise ValueError("boom")


def test_multiple_observers_all_notified():
    obs = Observable()
    calls1 = []
    calls2 = []

    def cb1(observable, event, *args, **kwargs):
        calls1.append(event)

    def cb2(observable, event, *args, **kwargs):
        calls2.append(event)

    obs.attach(cb1)
    obs.attach(cb2)
    obs.notify("broadcast")

    assert calls1 == ["broadcast"]
    assert calls2 == ["broadcast"]


def test_observer_exception_does_not_break_notification_of_others():
    obs = Observable()
    calls = []

    def bad(observable, event, *args, **kwargs):
        raise RuntimeError("observer failed")

    def good(observable, event, *args, **kwargs):
        calls.append(event)

    obs.attach(bad)
    obs.attach(good)

    # should not raise even though `bad` raises internally
    obs.notify("event")

    assert calls == ["event"]


def test_freeze_notifications_batches_and_replays_in_order():
    obs = Observable()
    calls = []

    def callback(observable, event, *args, **kwargs):
        calls.append(event)

    obs.attach(callback)

    with obs.freeze_notifications():
        obs.notify("first")
        obs.notify("second")
        # nothing delivered yet while frozen
        assert calls == []

    assert calls == ["first", "second"]


def test_freeze_notifications_yields_self():
    obs = Observable()
    # at least one notification must occur inside the block, otherwise
    # exiting freeze_notifications() hits the UnboundLocalError bug
    # documented in test_freeze_notifications_with_no_pending_events_bug
    with obs.freeze_notifications() as ctx:
        assert ctx is obs
        obs.notify("noop")


def test_freeze_notifications_with_no_pending_events_bug():
    """Documents a real bug in Observable.freeze_notifications().

    On exit, the generator only assigns the local variable `pending` inside
    `if self._frozen == 0 and self._pending:`, but then unconditionally
    iterates over `pending` afterwards. If nothing was notified while frozen
    (`self._pending` is empty), `pending` is never assigned and exiting the
    context manager raises UnboundLocalError - even though nothing else
    about the usage was incorrect.
    """
    obs = Observable()
    with pytest.raises(UnboundLocalError):
        with obs.freeze_notifications():
            pass


def test_nested_freeze_notifications_bug():
    """Documents the same freeze_notifications bug as above, triggered by
    nesting: the inner context manager exits while the outer one is still
    active, so `self._frozen` is not yet back to 0 and `pending` is never
    assigned, even though a notification did occur.
    """
    obs = Observable()
    with pytest.raises(UnboundLocalError):
        with obs.freeze_notifications():
            with obs.freeze_notifications():
                obs.notify("nested")


def test_weakref_callback_stops_receiving_after_garbage_collection():
    """Plain callback functions (as opposed to bound `.update` methods, see
    test_attach_observer_object_is_dropped_immediately_bug) are held
    correctly by the internal WeakSet: they keep receiving notifications as
    long as something else keeps them alive, and are silently dropped (no
    error) once they are garbage collected.
    """
    obs = Observable()
    calls = []

    def make_callback():
        def callback(observable, event, *args, **kwargs):
            calls.append(event)

        return callback

    callback = make_callback()
    obs.attach(callback)
    obs.notify("first")
    assert calls == ["first"]

    del callback
    gc.collect()

    # should not raise even though the callback has been garbage collected
    obs.notify("second")
    assert calls == ["first"]


def test_pickling_drops_lock_and_observers_but_restores_state():
    obs = Observable()
    recorder = Recorder()
    obs.attach(recorder)

    data = pickle.dumps(obs)
    restored = pickle.loads(data)

    # the restored object should have fresh internal bookkeeping
    assert restored._observers == {}
    assert list(restored._any_observers) == []
    assert restored._frozen == 0

    # notify should work fine on the restored object even though the
    # original observer registration was not (and cannot be) preserved
    restored.notify("event")
