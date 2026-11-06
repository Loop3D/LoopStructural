"""Tests for the Viewer class (story 1)."""

from __future__ import annotations

import base64
import contextlib
import gc
import hashlib
import socket
import struct
import subprocess
import threading
import time
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from loop_link import discovery
from loop_link.errors import (
    BusyError,
    MissingFeatureError,
    NoViewerError,
    Refusal,
    TokenError,
    ViewerClosedError,
    ViewerRefusedError,
)
from loop_link.objects import (
    BlockModel,
    Categorical,
    Drillhole,
    Lines,
    Points,
    Surface,
    TetraMesh,
    Vectors,
)
from loop_link.protocol import Frame, encode_frame
from loop_link.testing import FakeViewer
from loop_link.viewer import InventoryEntry, Viewer

FACTORY = Callable[..., Any]


def wait_for(condition: Callable[[], bool], timeout: float = 5.0) -> bool:
    """Poll a condition until it is true or the deadline passes."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if condition():
            return True
        time.sleep(0.01)
    return condition()


def stored_array(stored: Any, name: str) -> np.ndarray:
    """Return the buffer of the data description with this name."""
    entry = next(item for item in stored.header["data"] if item["name"] == name)
    return stored.buffers[entry["buffer"]]


@contextlib.contextmanager
def own_viewer(tmp_path: Path, **kwargs: Any) -> Iterator[tuple[FakeViewer, Path]]:
    """Start a fake viewer with its own settings folder."""
    fake = FakeViewer(**kwargs)
    settings = tmp_path / "own_settings"
    try:
        fake.write_instance_file(settings)
        yield fake, settings
    finally:
        fake.close()


def small_points() -> Points:
    """Return a small Points object."""
    return Points(np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))


# -- scenario 1: discovery ----------------------------------------------------


def test_discovery_without_port_or_token(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """The viewer finds the port and token alone. The object arrives in the scene."""
    viewer = viewer_factory()
    viewer.send(small_points(), id="pts")
    viewer.flush()
    assert wait_for(lambda: "pts" in fake_viewer.viewer.objects)
    assert viewer.connected


# -- scenario 2: revisions ----------------------------------------------------


def test_same_id_twice_gives_revisions(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """A second send with the same id gives revision 2. The stored header has revision 2."""
    viewer = viewer_factory()
    first = viewer.send(small_points(), id="pts")
    second = viewer.send(small_points(), id="pts")
    viewer.flush()
    assert (first, second) == (1, 2)
    assert viewer.revisions["pts"] == 2
    assert wait_for(
        lambda: "pts" in fake_viewer.viewer.objects
        and fake_viewer.viewer.objects["pts"].header["revision"] == 2
    )


# -- scenario 3: no viewer ----------------------------------------------------


def test_no_viewer_gives_error_and_starts_no_process(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Without a viewer, connect gives NoViewerError. No process starts."""
    started: list[Any] = []

    def refuse(*args: Any, **kwargs: Any) -> None:
        started.append(args)
        raise AssertionError("A process must not start.")

    monkeypatch.setattr(subprocess, "Popen", refuse)
    # Use a free port as the default port, so a real viewer on this machine does not matter.
    probe = socket.socket()
    probe.bind(("127.0.0.1", 0))
    free_port = probe.getsockname()[1]
    probe.close()
    monkeypatch.setattr(discovery, "DEFAULT_PORT", free_port)
    viewer = Viewer(settings_dir=tmp_path / "empty", token="secret")
    try:
        with pytest.raises(NoViewerError):
            viewer.connect()
    finally:
        viewer.close()
    assert started == []


# -- scenario 4: properties ---------------------------------------------------


def test_properties_keep_names_and_values(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """The stored buffers have the same property names and values."""
    grade = np.array([0.5, 1.5, 2.5])
    count = np.array([1.0, 2.0, 3.0])
    points = Points(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
        properties={"grade": grade, "count": count},
    )
    viewer = viewer_factory()
    viewer.send(points, id="pts")
    viewer.flush()
    assert wait_for(lambda: "pts" in fake_viewer.viewer.objects)
    stored = fake_viewer.viewer.objects["pts"]
    names = {item["name"] for item in stored.header["data"]}
    assert {"grade", "count"} <= names
    np.testing.assert_allclose(stored_array(stored, "grade"), grade)
    np.testing.assert_allclose(stored_array(stored, "count"), count)


# -- scenario 5: seven object types -------------------------------------------


def test_seven_object_types(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """One object of each type gives 7 stored objects with the right object_type."""
    square = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    objects: dict[str, tuple[Any, str]] = {
        "points": (Points(square), "Points"),
        "vectors": (
            Vectors.from_directions(square, np.tile([0.0, 0.0, 1.0], (4, 1)), kind="normal"),
            "Points",
        ),
        "lines": (Lines(square, [[0, 1], [1, 2]]), "Curve"),
        "surface": (Surface(square, [[0, 1, 2], [0, 1, 3]]), "Surface"),
        "tetra": (TetraMesh(square, [[0, 1, 2, 3]]), "TetraMesh"),
        "block": (BlockModel((0, 0, 0), (1, 1, 1), (2, 2, 2)), "BlockModel"),
        "hole": (
            Drillhole(
                (0.0, 0.0, 0.0),
                [[0.0, 0.0, 90.0], [10.0, 0.0, 90.0]],
                [[0.0, 5.0], [5.0, 10.0]],
            ),
            "Drillhole",
        ),
    }
    viewer = viewer_factory()
    for object_id, (obj, _) in objects.items():
        viewer.send(obj, id=object_id)
    viewer.flush()
    assert wait_for(lambda: len(fake_viewer.viewer.objects) == 7)
    for object_id, (_, expected) in objects.items():
        assert fake_viewer.viewer.objects[object_id].header["object_type"] == expected


# -- scenario 6: orientation values -------------------------------------------


def test_vectors_from_normals_store_dip_direction_and_dip(
    fake_viewer: Any, viewer_factory: FACTORY
) -> None:
    """Normals give the stored dip_direction and dip that geology rules say."""
    dip_direction = np.array([90.0, 180.0])
    dip = np.array([30.0, 60.0])
    dd, dp = np.radians(dip_direction), np.radians(dip)
    # An upward plane normal points to the dip direction, tilted by the dip from vertical.
    normals = np.column_stack([np.sin(dp) * np.sin(dd), np.sin(dp) * np.cos(dd), np.cos(dp)])
    vectors = Vectors.from_directions(
        np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]), normals, kind="normal"
    )
    viewer = viewer_factory()
    viewer.send(vectors, id="vec")
    viewer.flush()
    assert wait_for(lambda: "vec" in fake_viewer.viewer.objects)
    stored = fake_viewer.viewer.objects["vec"]
    np.testing.assert_allclose(stored_array(stored, "orientation/dip_direction"), dip_direction)
    np.testing.assert_allclose(stored_array(stored, "orientation/dip"), dip)
    assert stored_array(stored, "orientation/dip").dtype == np.float64


# -- scenario 7: second viewer is busy ----------------------------------------


def test_second_viewer_is_busy_until_release(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """A second Viewer gets BusyError. After release_data the second one can send."""
    first = viewer_factory()
    first.send(small_points(), id="one")
    first.flush()
    assert wait_for(lambda: "one" in fake_viewer.viewer.objects)
    before = fake_viewer.viewer.objects["one"]

    second = viewer_factory()
    with pytest.raises(BusyError):
        second.connect()
    assert set(fake_viewer.viewer.objects) == {"one"}
    assert fake_viewer.viewer.objects["one"] is before

    first.release_data()

    def second_sends() -> bool:
        try:
            second.send(small_points(), id="two")
            second.flush()
        except BusyError:
            return False
        return True

    assert wait_for(second_sends)
    assert wait_for(lambda: "two" in fake_viewer.viewer.objects)
    assert "one" in fake_viewer.viewer.objects


# -- other cases --------------------------------------------------------------


def test_wrong_token_gives_token_error(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """A wrong token gives TokenError."""
    viewer = viewer_factory(token="wrong")
    with pytest.raises(TokenError):
        viewer.connect()


def test_reply_errors_make_flush_raise(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """A refused object makes flush raise ViewerRefusedError with a Refusal that has the id."""
    fake_viewer.viewer.reply_errors["bad"] = "The object is not good."
    viewer = viewer_factory()
    viewer.send(small_points(), id="bad")
    with pytest.raises(ViewerRefusedError) as info:
        viewer.flush()
    refusals = info.value.refusals
    assert [r.object_id for r in refusals] == ["bad"]
    assert isinstance(refusals[0], Refusal)
    assert "not good" in refusals[0].message


def test_flush_can_return_refusals(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """With raise_on_refusal false, flush returns the refusals."""
    fake_viewer.viewer.reply_errors["bad"] = "No."
    viewer = viewer_factory()
    viewer.send(small_points(), id="bad")
    refusals = viewer.flush(raise_on_refusal=False)
    assert [r.object_id for r in refusals] == ["bad"]


def test_scene_origin_of_viewer_is_used(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The scene origin of the viewer wins over the proposed origin."""
    monkeypatch.setenv("GEOVIEWER_TOKEN", "secret")
    with own_viewer(tmp_path, origin=(100.0, 200.0, 10.0)) as (fake, settings):
        viewer = Viewer(settings_dir=settings, origin=(5.0, 5.0, 5.0), heartbeat_interval=0.1)
        try:
            viewer.connect()
            viewer.send(Points(np.array([[110.0, 210.0, 20.0]])), id="pts")
            viewer.flush()
            assert tuple(viewer.scene_origin) == (100.0, 200.0, 10.0)
            assert viewer.scene_crs == "EPSG:28350"
            assert wait_for(lambda: "pts" in fake.objects)
            stored = fake.objects["pts"]
            buffer = stored.buffers[stored.header["geometry"]["vertices"]]
            np.testing.assert_allclose(buffer, [[10.0, 10.0, 10.0]], atol=1e-4)
        finally:
            viewer.close()


def test_viewer_closes_during_send(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """A viewer that closes in the middle of a send gives ViewerClosedError."""
    viewer = viewer_factory()
    viewer.connect()
    big = Points(np.random.default_rng(1).random((200_000, 3)))
    timer = threading.Timer(0.05, fake_viewer.viewer.close)
    timer.start()
    deadline = time.monotonic() + 5.0
    try:
        with pytest.raises(ViewerClosedError):
            index = 0
            while time.monotonic() < deadline:
                viewer.send(big, id=f"big{index}")
                index += 1
                time.sleep(0.01)
            viewer.flush()
    finally:
        timer.join()


def test_keepalive_holds_connection(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Heartbeats keep the connection open for longer than the viewer timeout."""
    monkeypatch.setenv("GEOVIEWER_TOKEN", "secret")
    with own_viewer(tmp_path, heartbeat_timeout=0.6) as (fake, settings):
        viewer = Viewer(settings_dir=settings, heartbeat_interval=0.1)
        try:
            viewer.connect()
            time.sleep(2.0)  # Real time is needed: the viewer closes idle connections.
            assert viewer.connected
            viewer.send(small_points(), id="late")
            viewer.flush()
            assert wait_for(lambda: "late" in fake.objects)
        finally:
            viewer.close()


def test_unsupported_version_gives_missing_feature(
    fake_viewer: Any, viewer_factory: FACTORY
) -> None:
    """An unsupported_version error from the viewer gives MissingFeatureError."""

    def reply(data: bytes) -> list[Frame]:
        header = {"type": "error", "code": "unsupported_version", "message": "Version."}
        return [Frame(header=header, buffers={})]

    fake_viewer.viewer.accept_frame = reply
    viewer = viewer_factory()
    with pytest.raises(MissingFeatureError):
        viewer.connect()


def test_inventory_lists_objects(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """A new viewer gets the inventory that the viewer sent at hello time."""
    first = viewer_factory()
    first.send(small_points(), id="pts")
    first.flush()
    assert wait_for(lambda: "pts" in fake_viewer.viewer.objects)
    first.release_data()
    assert first.revisions["pts"] == 1  # Revisions are kept after release_data.

    second = viewer_factory()
    deadline = time.monotonic() + 5.0
    while True:
        try:
            second.connect()
            break
        except BusyError:
            assert time.monotonic() < deadline
            time.sleep(0.01)
    entries = second.inventory()
    assert all(isinstance(entry, InventoryEntry) for entry in entries)
    entry = next(entry for entry in entries if entry.id == "pts")
    assert entry.object_type == "Points"
    assert entry.revision == 1


def test_name_defaults_to_id(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """Without a name, the stored header has the id as name."""
    viewer = viewer_factory()
    viewer.send(small_points(), id="x")
    viewer.flush()
    assert wait_for(lambda: "x" in fake_viewer.viewer.objects)
    assert fake_viewer.viewer.objects["x"].header["name"] == "x"


def test_send_after_release_data_reconnects(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """After release_data, a send opens a new connection and the revision goes on."""
    viewer = viewer_factory()
    assert viewer.send(small_points(), id="a") == 1
    viewer.flush()
    viewer.release_data()
    assert wait_for(lambda: "a" in fake_viewer.viewer.objects)
    assert viewer.send(small_points(), id="a") == 2
    viewer.flush()
    assert wait_for(
        lambda: "a" in fake_viewer.viewer.objects
        and fake_viewer.viewer.objects["a"].header["revision"] == 2
    )


def test_send_after_close_raises(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """After close, a send gives ViewerClosedError."""
    viewer = viewer_factory()
    viewer.connect()
    viewer.close()
    with pytest.raises(ViewerClosedError):
        viewer.send(small_points(), id="a")


def test_inventory_follows_sends_removes_and_clear(
    fake_viewer: Any, viewer_factory: FACTORY
) -> None:
    """The inventory shows the known scene: send, send again, remove and clear."""
    viewer = viewer_factory()
    viewer.send(small_points(), id="a")
    viewer.send(small_points(), id="b")
    entries = {entry.id: entry for entry in viewer.inventory()}
    assert set(entries) == {"a", "b"}
    assert entries["a"].revision == 1 and entries["b"].revision == 1
    assert entries["a"].object_type == "Points"
    viewer.send(small_points(), id="a")
    assert {e.id: e.revision for e in viewer.inventory()} == {"a": 2, "b": 1}
    viewer.send_group("g")
    assert {e.id: e.object_type for e in viewer.inventory()}["g"] == "Group"
    viewer.remove("a")
    assert {e.id for e in viewer.inventory()} == {"b", "g"}
    assert "a" not in viewer.revisions
    viewer.clear()
    assert viewer.inventory() == []


def test_inventory_of_new_connection_is_hello_list(
    fake_viewer: Any, viewer_factory: FACTORY
) -> None:
    """A new viewer sees the hello list after the first one released the scene."""
    first = viewer_factory()
    first.send(small_points(), id="a")
    first.flush()
    assert wait_for(lambda: "a" in fake_viewer.viewer.objects)
    first.release_data()
    second = viewer_factory()
    deadline = time.monotonic() + 5.0
    while True:
        try:
            second.connect()
            break
        except BusyError:
            assert time.monotonic() < deadline
            time.sleep(0.01)
    assert [(e.id, e.revision) for e in second.inventory()] == [("a", 1)]


def test_text_and_categorical_properties_are_stored(
    fake_viewer: Any, viewer_factory: FACTORY
) -> None:
    """A text property and a Categorical property reach the scene without refusal."""
    points = Points(
        np.zeros((3, 3)),
        properties={
            "lith": np.array(["a", "bb", "ccc"]),
            "unit": Categorical([0, 1, 1], ["x", "y"]),
        },
    )
    viewer = viewer_factory()
    viewer.send(points, id="pts")
    assert viewer.flush() == []
    assert wait_for(lambda: "pts" in fake_viewer.viewer.objects)


def test_text_offsets_must_match_vertex_count(fake_viewer: Any) -> None:
    """A text array with 2 values on 3 vertices is refused with invalid_object."""
    header = {
        "type": "scene.upsert",
        "id": "bad",
        "object_type": "Points",
        "geometry": {"vertices": "v"},
        "data": [
            {
                "id": "t",
                "name": "t",
                "buffer": "t",
                "association": "vertex",
                "data_type": "text",
                "offsets": "t/offsets",
            }
        ],
    }
    buffers = [
        ("v", np.zeros((3, 3), dtype=np.float32)),
        ("t", np.frombuffer(b"abc", dtype=np.uint8)),
        ("t/offsets", np.array([0, 1, 3], dtype=np.uint32)),
    ]
    replies = fake_viewer.viewer.accept_frame(encode_frame(header, buffers))
    assert [r.header["code"] for r in replies] == ["invalid_object"]
    assert "bad" not in fake_viewer.viewer.objects


def test_connect_returns_viewer(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """connect returns the viewer, so calls can chain."""
    viewer = viewer_factory()
    assert viewer.connect() is viewer
    assert viewer.connect().send(small_points(), id="a") == 1


def test_remove_and_clear_drop_revisions(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """After remove or clear, the next send of the id has revision 1."""
    viewer = viewer_factory()
    viewer.send(small_points(), id="a")
    viewer.send(small_points(), id="b")
    viewer.remove("a")
    assert viewer.send(small_points(), id="a") == 1
    viewer.clear()
    assert viewer.send(small_points(), id="b") == 1


def _serve_busy(server: socket.socket) -> None:
    """Accept one client, send busy and a close frame, and close without reading."""
    try:
        conn, _ = server.accept()
    except OSError:
        return
    crlf = chr(13) + chr(10)
    data = b""
    while (crlf * 2).encode() not in data:
        data += conn.recv(4096)
    key = next(
        line.split(":", 1)[1].strip()
        for line in data.decode().split(crlf)
        if line.lower().startswith("sec-websocket-key")
    )
    guid = b"258EAFA5-E914-47DA-95CA-C5AB0DC85B11"
    accept = base64.b64encode(hashlib.sha1(key.encode() + guid).digest()).decode()
    lines = [
        "HTTP/1.1 101 Switching Protocols",
        "Upgrade: websocket",
        "Connection: Upgrade",
        f"Sec-WebSocket-Accept: {accept}",
        "",
        "",
    ]
    conn.sendall(crlf.join(lines).encode())
    from loop_link.protocol import encode_frame

    body = encode_frame({"type": "error", "code": "busy", "message": "Busy."})
    conn.sendall(bytes([0x82, 126]) + struct.pack(">H", len(body)) + body)
    conn.sendall(bytes([0x88, 0]))
    conn.close()


def test_busy_server_that_closes_at_once_gives_busy_error() -> None:
    """A viewer that sends busy and closes without reading gives BusyError each time."""
    for _ in range(20):
        server = socket.socket()
        server.bind(("127.0.0.1", 0))
        server.listen(1)
        thread = threading.Thread(target=_serve_busy, args=(server,), daemon=True)
        thread.start()
        try:
            viewer = Viewer(port=server.getsockname()[1], token="t", timeout=2.0)
            with pytest.raises(BusyError):
                viewer.connect()
        finally:
            server.close()
            thread.join(timeout=2.0)


def test_dropped_viewer_frees_the_data_slot(fake_viewer: Any, viewer_factory: FACTORY) -> None:
    """A Viewer that is deleted closes its socket, so another Viewer can send."""
    first = Viewer(settings_dir=fake_viewer.settings, heartbeat_interval=0.1)
    first.send(small_points(), id="one")
    first.flush()
    del first
    gc.collect()
    second = viewer_factory()
    deadline = time.monotonic() + 2.0
    while True:
        try:
            second.send(small_points(), id="two")
            break
        except BusyError:
            assert time.monotonic() < deadline
            time.sleep(0.02)


def test_bad_incoming_frame_closes_the_viewer(
    fake_viewer: Any, viewer_factory: FACTORY, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A frame that cannot be decoded gives ViewerClosedError and marks the viewer closed."""
    from loop_link import viewer as viewer_module
    from loop_link.protocol import ProtocolError

    viewer = viewer_factory()
    viewer.connect()
    fake_viewer.viewer.reply_errors["bad"] = "No."
    real = viewer_module.protocol.decode_frame
    main = threading.current_thread()

    def broken(data: bytes) -> Any:
        # The fake viewer decodes frames in other threads, so break only the test thread.
        if threading.current_thread() is not main:
            return real(data)
        monkeypatch.setattr(viewer_module.protocol, "decode_frame", real)
        raise ProtocolError("Garbage.")

    monkeypatch.setattr(viewer_module.protocol, "decode_frame", broken)
    with pytest.raises(ViewerClosedError):
        viewer.send(small_points(), id="bad")
        viewer.flush()
    assert not viewer.connected
    with pytest.raises(ViewerClosedError):
        viewer.send(small_points(), id="again")
