"""Tests for the plugin push API."""

from __future__ import annotations

import itertools
import threading
import time
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

import numpy as np
from loop_link.push import PushResult, push_objects

from loop_link.errors import Refusal
from loop_link.objects import Points
from loop_link.testing import FakeViewer
from loop_link.viewer import Viewer


def _points(offset: float) -> Points:
    """Return a small point set that depends on the offset."""
    vertices = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]) + offset
    return Points(vertices=vertices, properties={"value": np.array([1.0, 2.0, 3.0]) + offset})


def _wait_for(condition: Callable[[], bool], timeout: float = 5.0) -> bool:
    """Wait until the condition is true or the timeout ends."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if condition():
            return True
        time.sleep(0.01)
    return condition()


def _header_without_version(header: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in header.items() if key != "version"}


def test_push_frames_equal_to_send(
    fake_viewer: Any, viewer_factory: Callable[..., Any], tmp_path: Path
) -> None:
    objects = [(f"obj/{index}", _points(float(index))) for index in range(3)]
    second_fake = FakeViewer()
    try:
        second_settings = tmp_path / "second"
        second_fake.write_instance_file(second_settings)
        pushed = viewer_factory()
        result = push_objects(pushed, objects)
        assert result.sent == 3
        with Viewer(settings_dir=second_settings, heartbeat_interval=0.1) as direct:
            for object_id, obj in objects:
                direct.send(obj, id=object_id)
            direct.flush()
        assert _wait_for(lambda: len(second_fake.objects) == 3)
        assert _wait_for(lambda: len(fake_viewer.viewer.objects) == 3)
        assert set(fake_viewer.viewer.objects) == set(second_fake.objects)
        for object_id, _ in objects:
            first = fake_viewer.viewer.objects[object_id]
            second = second_fake.objects[object_id]
            assert _header_without_version(first.header) == _header_without_version(second.header)
            assert set(first.buffers) == set(second.buffers)
            for name, array in first.buffers.items():
                assert np.array_equal(array, second.buffers[name])
    finally:
        second_fake.close()


def test_replace_prefix_removes_old_objects(
    fake_viewer: Any, viewer_factory: Callable[..., Any]
) -> None:
    viewer = viewer_factory()
    viewer.send(_points(0.0), id="model/old")
    viewer.send(_points(1.0), id="other/keep")
    viewer.flush()
    assert _wait_for(lambda: len(fake_viewer.viewer.objects) == 2)
    viewer.release_data()  # the next connect reads an inventory with both old objects
    result = push_objects(viewer, [("model/new", _points(2.0))], replace_prefix="model/")
    assert result.removed == 1
    assert result.sent == 1
    assert _wait_for(lambda: "model/old" not in fake_viewer.viewer.objects)
    ids = {entry["id"] for entry in fake_viewer.viewer.inventory()}
    assert ids == {"other/keep", "model/new"}


def test_invalid_object_is_skipped_and_push_goes_on(
    fake_viewer: Any, viewer_factory: Callable[..., Any]
) -> None:
    viewer = viewer_factory()
    objects = [("first", _points(0.0)), ("bad", 42), ("last", _points(1.0))]
    result = push_objects(viewer, objects)
    assert result.sent == 2
    assert len(result.skipped) == 1
    skipped_id, reason = result.skipped[0]
    assert skipped_id == "bad"
    assert isinstance(reason, str)
    assert reason != ""
    assert _wait_for(lambda: {"first", "last"} <= set(fake_viewer.viewer.objects))
    assert "bad" not in fake_viewer.viewer.objects


def test_cancel_stops_after_first_object(
    fake_viewer: Any, viewer_factory: Callable[..., Any]
) -> None:
    viewer = viewer_factory()
    cancel = threading.Event()

    def objects() -> Iterator[tuple[str, Points]]:
        yield "one", _points(0.0)
        cancel.set()  # the push asks for the next object after it sent the first one
        yield "two", _points(1.0)
        yield "three", _points(2.0)

    result = push_objects(viewer, objects(), cancel=cancel)
    assert result.cancelled is True
    assert result.sent == 1
    assert _wait_for(lambda: "one" in fake_viewer.viewer.objects)
    assert "two" not in fake_viewer.viewer.objects


def test_progress_text_before_each_object_on_calling_thread(
    fake_viewer: Any, viewer_factory: Callable[..., Any]
) -> None:
    viewer = viewer_factory()
    texts: list[str] = []
    threads: list[int] = []

    def progress(text: str) -> None:
        texts.append(text)
        threads.append(threading.get_ident())

    objects = [("a", _points(0.0)), ("b", _points(1.0)), ("c", _points(2.0))]
    push_objects(viewer, objects, progress=progress)
    assert texts == ["Sending a...", "Sending b...", "Sending c..."]
    assert set(threads) == {threading.get_ident()}


def test_main_thread_is_not_blocked_by_slow_push(tmp_path: Path, monkeypatch: Any) -> None:
    monkeypatch.setenv("GEOVIEWER_TOKEN", "secret")
    slow = FakeViewer(read_delay=0.05)
    try:
        slow_settings = tmp_path / "slow"
        slow.write_instance_file(slow_settings)
        with Viewer(settings_dir=slow_settings, heartbeat_interval=0.1) as viewer:
            objects = [(f"slow/{index}", _points(float(index))) for index in range(5)]
            results: list[PushResult] = []
            worker = threading.Thread(target=lambda: results.append(push_objects(viewer, objects)))
            start = time.monotonic()
            worker.start()
            ticks = 0
            tick_times = [start]
            for _ in range(20):
                time.sleep(0.01)
                ticks += 1
                tick_times.append(time.monotonic())
            elapsed = time.monotonic() - start
            push_running = worker.is_alive()
            worker.join(timeout=10.0)
            assert not worker.is_alive()
            assert ticks == 20
            largest_gap = max(b - a for a, b in itertools.pairwise(tick_times))
            assert largest_gap < 0.2
            assert push_running or elapsed < 0.5
            assert results[0].sent == 5
    finally:
        slow.close()


def test_summary_texts() -> None:
    refusal = Refusal(object_id="x", code="invalid_object", message="No.")
    skipped = [("a", "reason")]
    cases = [
        (PushResult(sent=3), "3 objects sent to GeoViewer."),
        (PushResult(sent=3, removed=1), "3 objects sent to GeoViewer, 1 old objects removed."),
        (PushResult(sent=3, skipped=skipped), "3 objects sent to GeoViewer, 1 objects skipped."),
        (
            PushResult(sent=3, removed=1, skipped=skipped),
            "3 objects sent to GeoViewer, 1 old objects removed, 1 objects skipped.",
        ),
        (
            PushResult(sent=3, refused=[refusal, refusal]),
            "3 objects sent to GeoViewer. The viewer refused 2 messages.",
        ),
        (
            PushResult(sent=3, removed=1, refused=[refusal, refusal]),
            "3 objects sent to GeoViewer, 1 old objects removed. The viewer refused 2 messages.",
        ),
    ]
    for result, expected in cases:
        assert result.summary() == expected


def test_push_result_defaults() -> None:
    result = PushResult()
    assert result.sent == 0
    assert result.removed == 0
    assert result.refused == []
    assert result.skipped == []
    assert result.cancelled is False
