"""Tests for the tool connection to the viewer (``loop_link.tool``)."""

from __future__ import annotations

import socket
import threading
import time
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from loop_link.tool import ToolConnection, raise_for_request_error

from loop_link import encode, protocol
from loop_link.errors import (
    BusyError,
    CancelledError,
    InvalidInputError,
    MissingFeatureError,
    NoViewerError,
    Refusal,
    StaleRevisionError,
    TokenError,
    ViewerClosedError,
    ViewerRefusedError,
)
from loop_link.protocol import Frame
from loop_link.testing import FakeViewer
from loop_link.testing.fake_viewer import StoredObject, _send

TIMEOUT = 5.0


class _Tap:
    """Wrap a reader. Keep each message that the reader returns."""

    def __init__(self, reader: Any, sink: list[bytes]) -> None:
        self._reader = reader
        self._sink = sink

    def message(self) -> bytes | None:
        """Return the next message and keep it."""
        data = self._reader.message()
        if data is not None:
            self._sink.append(data)
        return data


class _ScriptedViewer(FakeViewer):
    """A fake viewer that breaks the protocol in one way, set by ``mode``.

    Parameters
    ----------
    mode : str
        ``"normal"`` does not break anything. ``"stale"`` sends a frame with request 99
        before each real reply. ``"unsupported"`` answers each tool message with an
        ``unsupported_version`` error. ``"silent"`` never answers a tool message.
        ``"drop_first"`` closes the first tool connection after it reads one message.
        ``"refuse"`` answers each tool message with an ``error`` frame of an unknown code.
        ``"ignore_cancel"`` sends one ``command.progress`` frame, waits for a
        ``command.cancel`` and then sends ``selection.result`` as if it never saw the cancel.
        ``"count_cancel"`` sends ``command.progress`` frames for 0.5 s, counts the
        ``command.cancel`` messages that arrive and then sends a ``cancelled`` error.
        ``"slow_progress"`` sends ``command.progress`` every 0.1 s for 0.5 s and then
        ``selection.result``.

    Attributes
    ----------
    cancel_frames : list of bytes
        The messages that the viewer read between the parts of a fetch.
    tool_connections : int
        The number of tool connections.
    cancel_count : int
        The number of ``command.cancel`` messages that mode ``"count_cancel"`` read.
    """

    def __init__(self, mode: str = "normal", **kwargs: Any) -> None:
        self.mode = mode
        self.cancel_frames: list[bytes] = []
        self.tool_connections = 0
        self.cancel_count = 0
        super().__init__(**kwargs)

    def _tool_reply(self, frame: Frame) -> list[Frame]:
        if self.mode == "stale":
            stale = Frame(
                header={"type": "selection.result", "request": 99, "entries": []}, buffers={}
            )
            return [stale, *super()._tool_reply(frame)]
        if self.mode == "unsupported":
            header = {
                "type": "error",
                "code": "unsupported_version",
                "message": "Not supported.",
                "supported": {"major": 1, "minor": 14},
            }
            return [Frame(header=header, buffers={})]
        if self.mode == "refuse":
            header = {"type": "error", "code": "something_new", "message": "x"}
            return [Frame(header=header, buffers={})]
        return super()._tool_reply(frame)

    def _cancelled(self, conn: socket.socket, reader: Any, request: object) -> bool:
        return FakeViewer._cancelled(conn, _Tap(reader, self.cancel_frames), request)

    @staticmethod
    def _push(conn: socket.socket, header: dict[str, Any]) -> None:
        """Send one tool frame without buffers."""
        _send(conn, protocol.encode_frame(header, [], version=protocol.TOOL_PROTOCOL_VERSION))

    @staticmethod
    def _poll(conn: socket.socket, reader: Any, wait: float) -> dict[str, Any] | None:
        """Return the header of the next message, or None if none arrives in ``wait`` s."""
        conn.settimeout(wait)
        try:
            data = reader.message()
        except (TimeoutError, BlockingIOError):
            return None
        finally:
            conn.settimeout(None)
        return None if data is None else protocol.decode_frame(data).header

    def _scripted_request(self, conn: socket.socket, reader: Any) -> None:
        """Answer one request in a mode that sends progress frames."""
        data = reader.message()
        if data is None:
            return
        request = protocol.decode_frame(data).header["request"]
        progress = {"type": "command.progress", "request": request, "fraction": 0.5, "text": "x"}
        result = {"type": "selection.result", "request": request, "entries": []}
        if self.mode == "slow_progress":
            for _ in range(5):
                time.sleep(0.1)
                self._push(conn, progress)
            self._push(conn, result)
        elif self.mode == "ignore_cancel":
            self._push(conn, progress)
            deadline = time.monotonic() + 2.0
            while time.monotonic() < deadline:
                header = self._poll(conn, reader, 0.05)
                if header is not None and header.get("type") == "command.cancel":
                    break
            self._push(conn, result)
        else:
            deadline = time.monotonic() + 0.5
            while time.monotonic() < deadline:
                self._push(conn, progress)
                header = self._poll(conn, reader, 0.05)
                if header is not None and header.get("type") == "command.cancel":
                    self.cancel_count += 1
            error = {
                "type": "request.error",
                "request": request,
                "error": {"code": "cancelled", "message": "The client cancelled the request."},
            }
            self._push(conn, error)
        # Wait for the client to close, so a reset does not destroy the last frames.
        while self._poll(conn, reader, 0.05) is not None:
            pass

    def _run(self, conn: socket.socket, reader: Any, tool: bool) -> None:
        if tool:
            self.tool_connections += 1
            if self.mode in ("ignore_cancel", "count_cancel", "slow_progress"):
                self._scripted_request(conn, reader)
                return
            if self.mode == "silent":
                while reader.message() is not None:
                    pass
                return
            if self.mode == "drop_first" and self.tool_connections == 1:
                reader.message()
                conn.shutdown(socket.SHUT_RDWR)
                return
        super()._run(conn, reader, tool)


def _points(count: int = 10) -> StoredObject:
    """Return a stored Points object with ``count`` vertices."""
    header = {
        "id": "pts",
        "object_type": "Points",
        "revision": 1,
        "geometry": {"vertices": "geometry/vertices"},
    }
    vertices = np.arange(count * 3, dtype=np.float32).reshape(count, 3)
    return StoredObject(header=header, buffers={"geometry/vertices": vertices})


@pytest.fixture
def scripted(tmp_path: Path) -> Iterator[Callable[..., _ScriptedViewer]]:
    """Return a function that starts scripted viewers. Each viewer is closed at teardown."""
    made: list[_ScriptedViewer] = []

    def make(mode: str = "normal", **kwargs: Any) -> _ScriptedViewer:
        viewer = _ScriptedViewer(mode, **kwargs)
        made.append(viewer)
        return viewer

    yield make
    for viewer in made:
        viewer.close()


@pytest.fixture
def connection(fake_viewer: Any) -> Iterator[ToolConnection]:
    """Return a tool connection that finds the fake viewer through its instance file."""
    tool = ToolConnection(settings_dir=fake_viewer.settings, timeout=TIMEOUT)
    try:
        yield tool
    finally:
        tool.close()


@pytest.fixture
def direct() -> Iterator[Callable[..., ToolConnection]]:
    """Return a function that makes tool connections with a port. Each is closed at teardown."""
    made: list[ToolConnection] = []

    def make(viewer: FakeViewer, **kwargs: Any) -> ToolConnection:
        options: dict[str, Any] = {"port": viewer.port, "token": "secret", "timeout": TIMEOUT}
        options.update(kwargs)
        tool = ToolConnection(**options)
        made.append(tool)
        return tool

    yield make
    for tool in made:
        tool.close()


# -- connect and request numbers -------------------------------------------------------


def test_connect_finds_viewer_through_instance_file(connection: ToolConnection) -> None:
    assert not connection.connected
    connection.connect()
    assert connection.connected
    connection.close()
    assert not connection.connected


def test_request_connects_when_needed(connection: ToolConnection) -> None:
    connection.request(encode.selection_get_frame, until=("selection.result",), feature="selection")
    assert connection.connected


def test_request_numbers_start_at_one_and_rise(
    fake_viewer: Any, connection: ToolConnection
) -> None:
    for _ in range(3):
        connection.request(
            encode.selection_get_frame, until=("selection.result",), feature="selection"
        )
    numbers = [frame.header["request"] for frame in fake_viewer.viewer.tool_requests]
    assert numbers == [1, 2, 3]


def test_request_returns_selection_result(fake_viewer: Any, connection: ToolConnection) -> None:
    fake_viewer.viewer.objects["pts"] = _points()
    fake_viewer.viewer.set_user_selection("pts", nodes=[1, 3])
    frame = connection.request(
        encode.selection_get_frame, until=("selection.result",), feature="selection"
    )
    assert frame.header["type"] == "selection.result"
    assert frame.header["object"] == "pts"
    assert frame.buffers["nodes"].tolist() == [1, 3]


def test_request_passes_request_number_to_build(
    fake_viewer: Any, connection: ToolConnection
) -> None:
    fake_viewer.viewer.objects["pts"] = _points()
    seen: list[int] = []

    def build(request: int) -> bytes:
        seen.append(request)
        return encode.selection_set_frame(request, "pts", None, np.array([0, 2]), 1)

    frame = connection.request(
        build, until=("selection.result",), feature="selection", object_id="pts"
    )
    assert seen == [1]
    assert frame.header["request"] == 1
    assert fake_viewer.viewer.selection["nodes"].tolist() == [0, 2]


def test_frame_with_other_request_number_is_dropped(
    scripted: Callable[..., _ScriptedViewer], direct: Callable[..., ToolConnection]
) -> None:
    viewer = scripted("stale")
    viewer.objects["pts"] = _points()
    viewer.set_user_selection("pts", nodes=[4])
    seen: list[Frame] = []
    frame = direct(viewer).request(
        encode.selection_get_frame,
        until=("selection.result",),
        feature="selection",
        on_frame=seen.append,
    )
    assert frame.header["request"] == 1
    assert frame.buffers["nodes"].tolist() == [4]
    assert all(item.header["request"] != 99 for item in seen)


# -- progress, on_frame and cancel -----------------------------------------------------


def test_on_frame_gets_each_frame_of_a_fetch(
    scripted: Callable[..., _ScriptedViewer], direct: Callable[..., ToolConnection]
) -> None:
    viewer = scripted(part_bytes=64)
    viewer.objects["pts"] = _points()
    seen: list[str] = []
    frame = direct(viewer).request(
        lambda request: encode.fetch_frame(request, "pts"),
        until=("object.end",),
        feature="fetch",
        object_id="pts",
        on_frame=lambda item: seen.append(item.header["type"]),
    )
    assert frame.header["type"] == "object.end"
    assert seen[0] == "object.header"
    assert seen.count("object.part") == -(-120 // 64)


def test_progress_gets_rising_floats(
    scripted: Callable[..., _ScriptedViewer], direct: Callable[..., ToolConnection]
) -> None:
    viewer = scripted(part_bytes=64)
    viewer.objects["pts"] = _points()
    values: list[float] = []
    direct(viewer).request(
        lambda request: encode.fetch_frame(request, "pts"),
        until=("object.end",),
        feature="fetch",
        object_id="pts",
        progress=values.append,
    )
    assert len(values) >= 2
    assert all(isinstance(value, float) for value in values)
    assert values == sorted(set(values))
    assert 0.0 < values[0] and values[-1] <= 1.0


def test_cancel_gives_cancelled_error_and_one_cancel_frame(
    scripted: Callable[..., _ScriptedViewer], direct: Callable[..., ToolConnection]
) -> None:
    viewer = scripted(part_bytes=64, fetch_delay=0.05)
    viewer.objects["pts"] = _points(40)
    cancel = threading.Event()

    def on_frame(frame: Frame) -> None:
        if frame.header["type"] == "object.part":
            cancel.set()

    with pytest.raises(CancelledError):
        direct(viewer).request(
            lambda request: encode.fetch_frame(request, "pts"),
            until=("object.end",),
            feature="fetch",
            object_id="pts",
            on_frame=on_frame,
            cancel=cancel,
        )
    headers = [protocol.decode_frame(data).header for data in viewer.cancel_frames]
    cancels = [header for header in headers if header.get("type") == "command.cancel"]
    assert len(cancels) == 1
    assert cancels[0]["request"] == 1


def test_cancel_frame_is_sent_once(
    scripted: Callable[..., _ScriptedViewer], direct: Callable[..., ToolConnection]
) -> None:
    viewer = scripted("count_cancel")
    cancel = threading.Event()
    cancel.set()
    with pytest.raises(CancelledError):
        direct(viewer).request(
            encode.selection_get_frame,
            until=("selection.result",),
            feature="selection",
            cancel=cancel,
        )
    assert viewer.cancel_count == 1


def test_until_frame_after_cancel_is_returned(
    scripted: Callable[..., _ScriptedViewer], direct: Callable[..., ToolConnection]
) -> None:
    viewer = scripted("ignore_cancel")
    cancel = threading.Event()
    cancel.set()
    frame = direct(viewer).request(
        encode.selection_get_frame,
        until=("selection.result",),
        feature="selection",
        cancel=cancel,
    )
    assert frame.header["type"] == "selection.result"


def test_connection_is_usable_after_cancel(
    scripted: Callable[..., _ScriptedViewer], direct: Callable[..., ToolConnection]
) -> None:
    viewer = scripted(part_bytes=64, fetch_delay=0.05)
    viewer.objects["pts"] = _points(40)
    cancel = threading.Event()
    cancel.set()
    tool = direct(viewer)
    with pytest.raises(CancelledError):
        tool.request(
            lambda request: encode.fetch_frame(request, "pts"),
            until=("object.end",),
            feature="fetch",
            object_id="pts",
            cancel=cancel,
        )
    frame = tool.request(
        encode.selection_get_frame, until=("selection.result",), feature="selection"
    )
    assert frame.header["type"] == "selection.result"


# -- connection failures ---------------------------------------------------------------


def test_socket_closed_in_the_middle_of_a_request(
    scripted: Callable[..., _ScriptedViewer], direct: Callable[..., ToolConnection]
) -> None:
    viewer = scripted("drop_first")
    tool = direct(viewer)
    with pytest.raises(ViewerClosedError):
        tool.request(encode.selection_get_frame, until=("selection.result",), feature="selection")
    assert not tool.connected
    frame = tool.request(
        encode.selection_get_frame, until=("selection.result",), feature="selection"
    )
    assert frame.header["type"] == "selection.result"
    assert viewer.tool_connections == 2


def test_wrong_token_gives_token_error(direct: Callable[..., ToolConnection]) -> None:
    with FakeViewer() as viewer:
        tool = direct(viewer, token="wrong")
        with pytest.raises(TokenError):
            tool.connect()
        assert not tool.connected


def test_no_viewer_gives_no_viewer_error(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GEOVIEWER_TOKEN", "secret")
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    tool = ToolConnection(port=port, token="secret", timeout=TIMEOUT)
    try:
        with pytest.raises(NoViewerError):
            tool.connect()
    finally:
        tool.close()
    empty = ToolConnection(settings_dir=tmp_path / "empty", timeout=TIMEOUT)
    try:
        with pytest.raises(NoViewerError):
            empty.connect()
    finally:
        empty.close()


def test_no_frame_within_timeout_gives_viewer_closed_error(
    scripted: Callable[..., _ScriptedViewer], direct: Callable[..., ToolConnection]
) -> None:
    viewer = scripted("silent")
    tool = direct(viewer)
    start = time.monotonic()
    with pytest.raises(ViewerClosedError):
        tool.request(
            encode.selection_get_frame,
            until=("selection.result",),
            feature="selection",
            timeout=0.2,
        )
    # The default timeout is 5 s, so a short time shows that the argument was used.
    assert time.monotonic() - start < 2.0
    assert not tool.connected


def test_timeout_restarts_at_each_frame(
    scripted: Callable[..., _ScriptedViewer], direct: Callable[..., ToolConnection]
) -> None:
    viewer = scripted("slow_progress")
    frame = direct(viewer).request(
        encode.selection_get_frame,
        until=("selection.result",),
        feature="selection",
        timeout=0.3,
    )
    assert frame.header["type"] == "selection.result"


# -- errors from the viewer ------------------------------------------------------------


def test_request_error_stale_revision(fake_viewer: Any, connection: ToolConnection) -> None:
    fake_viewer.viewer.objects["pts"] = _points()
    with pytest.raises(StaleRevisionError) as info:
        connection.request(
            lambda request: encode.selection_set_frame(request, "pts", None, np.array([0]), 7),
            until=("selection.result",),
            feature="selection",
            object_id="pts",
        )
    assert info.value.object_id == "pts"


def test_request_error_invalid_input(fake_viewer: Any, connection: ToolConnection) -> None:
    fake_viewer.viewer.objects["pts"] = _points()
    with pytest.raises(InvalidInputError) as info:
        connection.request(
            lambda request: encode.selection_set_frame(request, "pts", None, np.array([999]), 1),
            until=("selection.result",),
            feature="selection",
            object_id="pts",
        )
    assert str(info.value).startswith("The viewer refused the request:")


def test_headless_viewer_gives_missing_feature(
    scripted: Callable[..., _ScriptedViewer], direct: Callable[..., ToolConnection]
) -> None:
    viewer = scripted(kind="headless")
    with pytest.raises(MissingFeatureError):
        direct(viewer).request(
            encode.selection_get_frame, until=("selection.result",), feature="selection"
        )


def test_unknown_message_gives_missing_feature_with_feature_name(
    scripted: Callable[..., _ScriptedViewer], direct: Callable[..., ToolConnection]
) -> None:
    viewer = scripted(protocol_minor=14)
    with pytest.raises(MissingFeatureError) as info:
        direct(viewer).request(
            encode.selection_get_frame, until=("selection.result",), feature="selection"
        )
    assert info.value.feature == "selection"


def test_unsupported_version_gives_missing_feature(
    scripted: Callable[..., _ScriptedViewer], direct: Callable[..., ToolConnection]
) -> None:
    viewer = scripted("unsupported")
    with pytest.raises(MissingFeatureError) as info:
        direct(viewer).request(
            encode.selection_get_frame, until=("selection.result",), feature="fetch"
        )
    assert info.value.feature.startswith("protocol")
    assert info.value.feature != "fetch"
    assert "14" in info.value.feature


def test_other_error_frame_gives_refusal(
    scripted: Callable[..., _ScriptedViewer], direct: Callable[..., ToolConnection]
) -> None:
    viewer = scripted("refuse")
    with pytest.raises(ViewerRefusedError) as info:
        direct(viewer).request(
            encode.selection_get_frame, until=("selection.result",), feature="selection"
        )
    assert info.value.refusals == [Refusal("", "something_new", "x")]


# -- raise_for_request_error -----------------------------------------------------------


def test_error_table_invalid_input_revision() -> None:
    error = {"code": "invalid_input", "message": "Old revision.", "field": "revision"}
    with pytest.raises(StaleRevisionError) as info:
        raise_for_request_error(error, object_id="pts")
    assert info.value.object_id == "pts"


def test_error_table_invalid_input_other_field() -> None:
    error = {"code": "invalid_input", "message": "Bad index.", "field": "indices"}
    with pytest.raises(InvalidInputError) as info:
        raise_for_request_error(error, object_id="pts")
    assert str(info.value).startswith("The viewer refused the request:")
    assert "Bad index." in str(info.value)


def test_error_table_needs_window() -> None:
    error = {"code": "needs_window", "message": "No window."}
    with pytest.raises(MissingFeatureError):
        raise_for_request_error(error, object_id="pts")


def test_error_table_cancelled() -> None:
    error = {"code": "cancelled", "message": "Cancelled."}
    with pytest.raises(CancelledError):
        raise_for_request_error(error, object_id="pts")


def test_error_table_busy() -> None:
    error = {"code": "busy", "message": "Busy."}
    with pytest.raises(BusyError):
        raise_for_request_error(error, object_id="pts")


def test_error_table_other_code_gives_refusal() -> None:
    error = {"code": "something_new", "message": "Not allowed."}
    with pytest.raises(ViewerRefusedError) as info:
        raise_for_request_error(error, object_id="pts")
    assert info.value.refusals == [Refusal("pts", "something_new", "Not allowed.")]


def test_error_table_other_code_without_object_id() -> None:
    error = {"code": "something_new", "message": "Not allowed."}
    with pytest.raises(ViewerRefusedError) as info:
        raise_for_request_error(error, object_id=None)
    assert info.value.refusals == [Refusal("", "something_new", "Not allowed.")]
