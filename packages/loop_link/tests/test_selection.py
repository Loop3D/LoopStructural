"""Tests for the selection in the viewer (``loop_link.selection`` and the ``Viewer`` methods)."""

from __future__ import annotations

import gc
import threading
import time
from collections.abc import Callable, Iterator
from typing import Any

import numpy as np
import pytest
from loop_link.selection import Selection, SelectionWatcher, parse_selection_result

from loop_link.encode import encode_object, hello_frame
from loop_link.errors import (
    InvalidInputError,
    MissingFeatureError,
    StaleRevisionError,
    ViewerClosedError,
)
from loop_link.objects import BlockModel
from loop_link.protocol import Frame, ProtocolError
from loop_link.testing import FakeViewer

WAIT = 3.0
THREAD_NAME = "loop-link-selection-poll"


def _grid(cells: tuple[int, int, int] = (4, 4, 4)) -> BlockModel:
    """Return a block model with ``cells`` cells on each axis."""
    return BlockModel((0, 0, 0), (1, 1, 1), cells)


def _wait_for(check: Callable[[], bool], timeout: float = WAIT) -> bool:
    """Return True when ``check`` is true. Return False after ``timeout`` seconds."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if check():
            return True
        time.sleep(0.01)
    return check()


def _poll_threads() -> list[threading.Thread]:
    """Return the running selection poll threads."""
    return [thread for thread in threading.enumerate() if thread.name == THREAD_NAME]


def _put_grid(fake: FakeViewer, cells: tuple[int, int, int] = (4, 4, 4)) -> None:
    """Put a block model with id ``grid`` into the fake viewer. Open no data connection."""
    fake.accept_frame(hello_frame("", (0.0, 0.0, 0.0), "test"))
    data = encode_object(
        _grid(cells), id="grid", name=None, parent=None, revision=1, origin=(0.0, 0.0, 0.0)
    )
    fake.accept_frame(data)


@pytest.fixture
def viewer(viewer_factory: Callable[..., Any]) -> Any:
    """Return a viewer that has sent a 4x4x4 block model with id ``grid``."""
    made = viewer_factory(timeout=WAIT)
    made.send(_grid(), id="grid")
    made.flush()
    return made


@pytest.fixture
def selector(fake_viewer: Any, viewer_factory: Callable[..., Any]) -> Any:
    """Return a viewer that only selects. The scene has a block model ``grid``."""
    _put_grid(fake_viewer.viewer)
    return viewer_factory(timeout=WAIT)


@pytest.fixture
def watchers() -> Iterator[list[SelectionWatcher]]:
    """Return a list. Each watcher in it is stopped at teardown."""
    made: list[SelectionWatcher] = []
    yield made
    for watcher in made:
        watcher.stop()


# -- Selection -------------------------------------------------------------------------


def _selection(
    object_id: str | None = "grid",
    cells: list[int] | None = None,
    vertices: list[int] | None = None,
    revision: int = 1,
) -> Selection:
    return Selection(
        object_id,
        revision,
        np.array(cells or [], dtype=np.uint32),
        np.array(vertices or [], dtype=np.uint32),
    )


def test_selection_dataclass_empty() -> None:
    assert _selection(None).empty
    assert not _selection(cells=[1]).empty
    assert not _selection(vertices=[1]).empty


def test_selection_is_frozen() -> None:
    selection = _selection(cells=[1])
    with pytest.raises(AttributeError):
        selection.object_id = "other"  # type: ignore[misc]


def test_same_as_ignores_revision_and_order() -> None:
    first = _selection(cells=[3, 1, 2], vertices=[9, 8], revision=1)
    second = _selection(cells=[1, 2, 3], vertices=[8, 9], revision=7)
    assert first.same_as(second)
    assert second.same_as(first)


def test_same_as_false_for_other_object_or_indices() -> None:
    base = _selection(cells=[1, 2], vertices=[5])
    assert not base.same_as(_selection("other", cells=[1, 2], vertices=[5]))
    assert not base.same_as(_selection(cells=[1, 3], vertices=[5]))
    assert not base.same_as(_selection(cells=[1, 2], vertices=[6]))
    assert not base.same_as(_selection(cells=[1, 2]))


# -- parse_selection_result ------------------------------------------------------------


def _result(entries: list[dict[str, Any]], buffers: dict[str, np.ndarray]) -> Frame:
    header = {
        "type": "selection.result",
        "request": 1,
        "object": "grid",
        "revision": 2,
        "entries": entries,
    }
    return Frame(header, buffers)


def test_parse_selection_result_maps_kinds() -> None:
    frame = _result(
        [
            {"kind": "cells", "count": 2, "indices": "c"},
            {"kind": "nodes", "count": 3, "indices": "n"},
        ],
        {"c": np.array([4, 5], dtype=np.uint32), "n": np.array([1, 2, 3], dtype=np.uint32)},
    )
    selection = parse_selection_result(frame)
    assert selection.object_id == "grid"
    assert selection.revision == 2
    assert selection.cells.tolist() == [4, 5]
    assert selection.vertices.tolist() == [1, 2, 3]
    assert selection.cells.dtype == np.uint32
    assert selection.vertices.dtype == np.uint32


def test_parse_selection_result_count_without_buffer() -> None:
    with pytest.raises(ProtocolError):
        parse_selection_result(_result([{"kind": "cells", "count": 2}], {}))


def test_parse_selection_result_buffer_not_uint32() -> None:
    frame = _result(
        [{"kind": "cells", "count": 2, "indices": "c"}],
        {"c": np.array([4, 5], dtype=np.int32)},
    )
    with pytest.raises(ProtocolError):
        parse_selection_result(frame)


def test_parse_selection_result_length_not_count() -> None:
    frame = _result(
        [{"kind": "cells", "count": 3, "indices": "c"}],
        {"c": np.array([4, 5], dtype=np.uint32)},
    )
    with pytest.raises(ProtocolError):
        parse_selection_result(frame)


# -- story 2: select -------------------------------------------------------------------


def test_select_cells(fake_viewer: Any, viewer: Any) -> None:
    result = viewer.select("grid", cells=[12, 13])
    assert result == {"cell": 2, "vertex": 0}
    selection = fake_viewer.viewer.selection
    assert selection is not None
    assert selection["object"] == "grid"
    assert selection["cells"].tolist() == [12, 13]


def test_select_vertices(fake_viewer: Any, viewer: Any) -> None:
    result = viewer.select("grid", vertices=[0, 5, 124])
    assert result == {"cell": 0, "vertex": 3}
    selection = fake_viewer.viewer.selection
    assert selection is not None
    assert selection["nodes"].tolist() == [0, 5, 124]


def test_select_cells_and_vertices_together(fake_viewer: Any, viewer: Any) -> None:
    result = viewer.select("grid", cells=[1, 2], vertices=[3])
    assert result == {"cell": 2, "vertex": 1}
    selection = fake_viewer.viewer.selection
    assert selection is not None
    assert selection["cells"].tolist() == [1, 2]
    assert selection["nodes"].tolist() == [3]


def test_select_duplicates_count_one_time(fake_viewer: Any, viewer: Any) -> None:
    result = viewer.select("grid", cells=[5, 5, 6, 5])
    assert result["cell"] == 2
    assert fake_viewer.viewer.selection["cells"].tolist() == [5, 6]


def test_select_no_indices_clears(fake_viewer: Any, viewer: Any) -> None:
    viewer.select("grid", cells=[1, 2])
    assert fake_viewer.viewer.selection is not None
    result = viewer.select("grid", cells=[])
    assert result == {"cell": 0, "vertex": 0}
    assert fake_viewer.viewer.selection is None


def test_select_accepts_numpy_array(fake_viewer: Any, viewer: Any) -> None:
    viewer.select("grid", cells=np.array([7, 8], dtype=np.int64))
    assert fake_viewer.viewer.selection["cells"].tolist() == [7, 8]


def test_select_index_out_of_range_keeps_selection(fake_viewer: Any, viewer: Any) -> None:
    viewer.select("grid", cells=[12, 13])
    with pytest.raises(InvalidInputError):
        viewer.select("grid", cells=[64])
    selection = fake_viewer.viewer.selection
    assert selection is not None
    assert selection["cells"].tolist() == [12, 13]


def test_select_with_current_revision(fake_viewer: Any, viewer: Any) -> None:
    revision = viewer.revisions["grid"]
    viewer.select("grid", cells=[1], revision=revision)
    assert fake_viewer.viewer.selection["cells"].tolist() == [1]


def test_select_with_old_revision_gives_stale_revision_error(fake_viewer: Any, viewer: Any) -> None:
    viewer.select("grid", cells=[12, 13])
    viewer.send(_grid(), id="grid")
    viewer.flush()
    with pytest.raises(StaleRevisionError):
        viewer.select("grid", cells=[1], revision=1)
    assert fake_viewer.viewer.selection["cells"].tolist() == [12, 13]


def test_select_without_indices_names_the_argument(viewer: Any) -> None:
    with pytest.raises(InvalidInputError) as info:
        viewer.select("grid")
    assert "cells" in str(info.value)
    assert "vertices" in str(info.value)


@pytest.mark.parametrize(
    ("kwargs", "name"),
    [
        ({"cells": [1.5]}, "cells"),
        ({"vertices": [2.5]}, "vertices"),
        ({"cells": [True]}, "cells"),
        ({"vertices": [False, True]}, "vertices"),
        ({"cells": [-1]}, "cells"),
        ({"vertices": [-1]}, "vertices"),
        ({"cells": [2**32]}, "cells"),
        ({"vertices": [2**32]}, "vertices"),
    ],
)
def test_select_bad_index_names_the_argument(
    fake_viewer: Any, viewer: Any, kwargs: dict[str, Any], name: str
) -> None:
    viewer.select("grid", cells=[3])
    with pytest.raises(InvalidInputError) as info:
        viewer.select("grid", **kwargs)
    assert name in str(info.value)
    # The check is local, so the selection in the viewer stays.
    assert fake_viewer.viewer.selection["cells"].tolist() == [3]


def test_select_largest_index_is_accepted(fake_viewer: Any, viewer: Any) -> None:
    # 2**32-1 is a valid value for the client. The viewer refuses it as out of range.
    with pytest.raises(InvalidInputError):
        viewer.select("grid", cells=[2**32 - 1])
    assert fake_viewer.viewer.selection is None


def test_select_100000_cells_is_one_frame(fake_viewer: Any, viewer_factory: Any) -> None:
    big = viewer_factory(timeout=10.0)
    big.send(_grid((50, 50, 50)), id="grid")
    big.flush()
    before = len(fake_viewer.viewer.tool_requests)
    result = big.select("grid", cells=np.arange(100_000))
    assert result == {"cell": 100_000, "vertex": 0}
    frames = fake_viewer.viewer.tool_requests[before:]
    sets = [frame for frame in frames if frame.header["type"] == "selection.set"]
    assert len(sets) == 1
    assert len(frames) == 1
    buffer = sets[0].buffers[sets[0].header["entries"][0]["indices"]]
    assert buffer.dtype == np.uint32
    assert buffer.shape == (100_000,)


def test_select_headless_viewer_gives_missing_feature(fake_viewer: Any) -> None:
    from loop_link.viewer import Viewer

    with FakeViewer(kind="headless") as fake:
        _put_grid(fake)
        selector = Viewer(port=fake.port, token="secret", timeout=WAIT)
        try:
            with pytest.raises(MissingFeatureError):
                selector.select("grid", cells=[1])
        finally:
            selector.close()


def test_select_old_protocol_gives_missing_feature(fake_viewer: Any) -> None:
    from loop_link.viewer import Viewer

    with FakeViewer(protocol_minor=14) as fake:
        _put_grid(fake)
        selector = Viewer(port=fake.port, token="secret", timeout=WAIT)
        try:
            with pytest.raises(MissingFeatureError):
                selector.select("grid", cells=[1])
            with pytest.raises(MissingFeatureError):
                selector.selection()
        finally:
            selector.close()


def test_selecting_viewer_does_not_open_data_connection(fake_viewer: Any, selector: Any) -> None:
    assert not selector.tool_connected
    result = selector.select("grid", cells=[2, 3])
    assert result == {"cell": 2, "vertex": 0}
    assert selector.tool_connected
    assert not selector.connected
    assert fake_viewer.viewer.selection["cells"].tolist() == [2, 3]
    assert selector.selection().cells.tolist() == [2, 3]
    assert not selector.connected


def test_send_viewer_has_tool_connection_only_after_select(viewer: Any) -> None:
    assert viewer.connected
    assert not viewer.tool_connected
    viewer.select("grid", cells=[1])
    assert viewer.tool_connected


# -- story 3: read the selection -------------------------------------------------------


def test_selection_returns_user_selection(fake_viewer: Any, selector: Any) -> None:
    fake_viewer.viewer.set_user_selection("grid", cells=[3, 1, 2], nodes=[10, 11])
    selection = selector.selection()
    assert isinstance(selection, Selection)
    assert selection.object_id == "grid"
    assert selection.revision == 1
    assert selection.cells.tolist() == [1, 2, 3]
    assert selection.vertices.tolist() == [10, 11]
    assert selection.cells.dtype == np.uint32
    assert selection.vertices.dtype == np.uint32
    assert not selection.empty


def test_selection_of_100000_cells(fake_viewer: Any, viewer_factory: Any) -> None:
    _put_grid(fake_viewer.viewer, (50, 50, 50))
    reader = viewer_factory(timeout=10.0)
    fake_viewer.viewer.set_user_selection("grid", cells=range(100_000))
    selection = reader.selection()
    assert selection.cells.dtype == np.uint32
    assert selection.cells.shape == (100_000,)
    assert selection.cells[-1] == 99_999


def test_selection_none_has_no_object(selector: Any) -> None:
    selection = selector.selection()
    assert selection.object_id is None
    assert selection.empty
    assert selection.cells.size == 0
    assert selection.vertices.size == 0


def test_selection_after_select_matches(selector: Any) -> None:
    selector.select("grid", cells=[9, 8], vertices=[1])
    assert selector.selection().same_as(_selection(cells=[8, 9], vertices=[1]))


def test_selection_after_close_gives_viewer_closed_error(selector: Any) -> None:
    selector.selection()
    selector.close()
    with pytest.raises(ViewerClosedError):
        selector.selection()


# -- story 3: watch_selection ----------------------------------------------------------


def test_watch_calls_once_for_each_change(
    fake_viewer: Any, selector: Any, watchers: list[SelectionWatcher]
) -> None:
    calls: list[Selection] = []
    assert _poll_threads() == []
    watcher = selector.watch_selection(calls.append, interval=0.05)
    watchers.append(watcher)
    assert isinstance(watcher, SelectionWatcher)
    assert watcher.running
    assert len(_poll_threads()) == 1
    fake_viewer.viewer.set_user_selection("grid", cells=[1, 2])
    assert _wait_for(lambda: len(calls) == 1)
    assert calls[0].object_id == "grid"
    assert calls[0].cells.tolist() == [1, 2]
    fake_viewer.viewer.set_user_selection("grid", cells=[5])
    assert _wait_for(lambda: len(calls) == 2)
    assert calls[1].cells.tolist() == [5]
    fake_viewer.viewer.set_user_selection("grid", cells=[1], nodes=[2])
    assert _wait_for(lambda: len(calls) == 3)
    assert calls[2].vertices.tolist() == [2]


def test_watch_does_not_call_without_change(
    fake_viewer: Any, selector: Any, watchers: list[SelectionWatcher]
) -> None:
    fake_viewer.viewer.set_user_selection("grid", cells=[1, 2])
    calls: list[Selection] = []
    watchers.append(selector.watch_selection(calls.append, interval=0.05))
    # Many polls pass. Nothing changes, so the callback is not called.
    time.sleep(0.4)
    assert calls == []
    # The same selection again is not a change.
    fake_viewer.viewer.set_user_selection("grid", cells=[2, 1])
    time.sleep(0.3)
    assert calls == []
    fake_viewer.viewer.set_user_selection("grid", cells=[3])
    assert _wait_for(lambda: len(calls) == 1)


def test_watch_initial_calls_once_with_first_selection(
    fake_viewer: Any, selector: Any, watchers: list[SelectionWatcher]
) -> None:
    fake_viewer.viewer.set_user_selection("grid", cells=[4, 5])
    calls: list[Selection] = []
    watchers.append(selector.watch_selection(calls.append, interval=0.05, initial=True))
    assert _wait_for(lambda: len(calls) >= 1)
    time.sleep(0.3)
    assert len(calls) == 1
    assert calls[0].cells.tolist() == [4, 5]


@pytest.mark.parametrize("interval", [0, -1.0])
def test_watch_interval_not_above_zero(selector: Any, interval: float) -> None:
    with pytest.raises(InvalidInputError) as info:
        selector.watch_selection(lambda selection: None, interval=interval)
    assert "interval" in str(info.value)
    assert _poll_threads() == []


def test_watch_stop_ends_thread(selector: Any) -> None:
    assert _poll_threads() == []
    watcher = selector.watch_selection(lambda selection: None, interval=0.05)
    try:
        assert len(_poll_threads()) == 1
    finally:
        watcher.stop()
    assert not watcher.running
    assert _wait_for(lambda: _poll_threads() == [])
    watcher.stop()  # A second stop does no harm.


def test_watcher_is_a_context_manager(selector: Any) -> None:
    with selector.watch_selection(lambda selection: None, interval=0.05) as watcher:
        assert watcher.running
    assert not watcher.running
    assert _wait_for(lambda: _poll_threads() == [])


def test_watch_callback_error_sets_error_and_ends_thread(
    fake_viewer: Any, selector: Any, watchers: list[SelectionWatcher]
) -> None:
    def callback(selection: Selection) -> None:
        raise RuntimeError("callback failed")

    watcher = selector.watch_selection(callback, interval=0.05)
    watchers.append(watcher)
    assert watcher.error is None
    fake_viewer.viewer.set_user_selection("grid", cells=[1])
    assert _wait_for(lambda: watcher.error is not None)
    assert isinstance(watcher.error, RuntimeError)
    assert _wait_for(lambda: not watcher.running)
    assert _wait_for(lambda: _poll_threads() == [])


def test_watch_ends_when_viewer_is_dropped(fake_viewer: Any) -> None:
    from loop_link.viewer import Viewer

    _put_grid(fake_viewer.viewer)
    dropped = Viewer(settings_dir=fake_viewer.settings, timeout=WAIT, heartbeat_interval=0.1)
    watcher = dropped.watch_selection(lambda selection: None, interval=0.05)
    try:
        assert watcher.running
        del dropped
        gc.collect()
        assert _wait_for(lambda: not watcher.running)
        assert _wait_for(lambda: _poll_threads() == [])
    finally:
        watcher.stop()


def test_close_stops_watchers(selector: Any, watchers: list[SelectionWatcher]) -> None:
    first = selector.watch_selection(lambda selection: None, interval=0.05)
    second = selector.watch_selection(lambda selection: None, interval=0.05)
    watchers.extend([first, second])
    assert first.running
    assert second.running
    selector.close()
    assert _wait_for(lambda: not first.running and not second.running)
    assert _wait_for(lambda: _poll_threads() == [])
    with pytest.raises(ViewerClosedError):
        selector.selection()
