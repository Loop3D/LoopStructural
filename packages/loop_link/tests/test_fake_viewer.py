"""Tests for the validation rules of the fake viewer."""

from __future__ import annotations

import contextlib
import copy
from collections.abc import Iterator
from typing import Any

import numpy as np
import pytest

from loop_link import protocol
from loop_link._websocket import WebSocketClient
from loop_link.encode import (
    cancel_frame,
    encode_object,
    fetch_frame,
    group_frame,
    hello_frame,
    selection_get_frame,
    selection_set_frame,
)
from loop_link.objects import BlockModel, Drillhole, Points, Surface, Vectors
from loop_link.protocol import Frame, decode_frame, encode_frame
from loop_link.testing import FakeViewer


@pytest.fixture
def fake() -> Any:
    """Return a fake viewer that has had a hello frame."""
    viewer = FakeViewer()
    viewer.accept_frame(encode_frame({"type": "hello"}))
    yield viewer
    viewer.close()


def codes(replies: list[Frame]) -> list[str]:
    """Return the error codes of the replies."""
    return [r.header["code"] for r in replies if r.header["type"] == "error"]


def entry(name: str, association: str, **extra: Any) -> dict[str, Any]:
    """Return a float data array description."""
    return {
        "id": name,
        "name": name,
        "buffer": name,
        "association": association,
        "data_type": "float",
        **extra,
    }


def points_frame(data: list[dict[str, Any]], groups: list[dict[str, Any]] | None = None) -> bytes:
    """Return a hand-made Points frame with 3 vertices."""
    header = {
        "type": "scene.upsert",
        "id": "p",
        "object_type": "Points",
        "revision": 1,
        "geometry": {"vertices": "v"},
        "data": data,
        "property_groups": groups or [],
    }
    buffers = [("v", np.zeros((3, 3), dtype=np.float32))]
    buffers += [(d["buffer"], np.zeros(3, dtype=np.float32)) for d in data]
    return encode_frame(header, buffers)


def drillhole_frame(
    association: str,
    log_sets: list[dict[str, Any]] | None = None,
    log_set: str | None = None,
    length: int = 2,
) -> bytes:
    """Return a hand-made Drillhole frame with a data array."""
    header = {
        "type": "scene.upsert",
        "id": "h",
        "object_type": "Drillhole",
        "revision": 1,
        "geometry": {"collar": "c", "survey": "s", "intervals": "i"},
        "data": [entry("x", association, **({"log_set": log_set} if log_set else {}))],
        "log_sets": log_sets or [],
    }
    buffers = [
        ("c", np.zeros(3, dtype=np.float32)),
        ("s", np.zeros((3, 3), dtype=np.float32)),
        ("i", np.zeros((2, 2), dtype=np.float32)),
        ("x", np.zeros(length, dtype=np.float32)),
    ]
    return encode_frame(header, buffers)


def test_drillhole_with_interval_properties_is_accepted(fake: Any) -> None:
    """The encoder sends interval properties as interval arrays. The viewer accepts them."""
    hole = Drillhole(
        collar=[0.0, 0.0, 0.0],
        survey=[[0.0, 0.0, 0.0], [10.0, 0.0, -90.0]],
        intervals=[[0.0, 5.0], [5.0, 10.0]],
        properties={"s": np.array([1.0, 2.0])},
        interval_properties={"grade": np.array([0.5, 1.5]), "lith": np.array(["a", "bb"])},
    )
    data = encode_object(hole, id="h", name=None, parent=None, revision=1, origin=(0.0, 0.0, 0.0))
    assert fake.accept_frame(data) == []
    assert "h" in fake.objects
    header = fake.objects["h"].header
    assert header["log_sets"] == [{"name": "intervals", "kind": "interval"}]
    assert {e["association"] for e in header["data"] if e["name"] == "grade"} == {"interval"}


def test_interval_array_without_log_sets_is_refused(fake: Any) -> None:
    """An interval association needs a log set of kind interval."""
    assert codes(fake.accept_frame(drillhole_frame("interval"))) == ["invalid_object"]
    no_name = drillhole_frame("interval", [{"name": "main", "kind": "interval"}])
    assert codes(fake.accept_frame(no_name)) == ["invalid_object"]
    with_set = drillhole_frame("interval", [{"name": "main", "kind": "interval"}], "main")
    assert fake.accept_frame(with_set) == []


def test_interval_array_with_unknown_or_wrong_kind_log_set_is_refused(fake: Any) -> None:
    """The log set of an interval array must exist and must have the kind interval."""
    sets = [{"name": "main", "kind": "interval"}]
    assert codes(fake.accept_frame(drillhole_frame("interval", sets, "other"))) == [
        "invalid_object"
    ]
    point = [{"name": "main", "kind": "point"}]
    assert codes(fake.accept_frame(drillhole_frame("interval", point, "main"))) == [
        "invalid_object"
    ]


def test_interval_array_with_wrong_length_is_refused(fake: Any) -> None:
    """An interval array needs one value for each row of intervals."""
    sets = [{"name": "main", "kind": "interval"}]
    assert codes(fake.accept_frame(drillhole_frame("interval", sets, "main", 3))) == [
        "invalid_object"
    ]


def test_group_naming_missing_array_is_refused(fake: Any) -> None:
    """A property group member must name an array of the message."""
    group = {
        "id": "orientation",
        "name": "orientation",
        "role": "orientation",
        "members": [
            {"role": "dip_direction", "data": "dd"},
            {"role": "dip", "data": "missing"},
        ],
    }
    frame = points_frame([entry("dd", "vertex")], [group])
    assert codes(fake.accept_frame(frame)) == ["invalid_object"]


def test_duplicate_data_id_is_refused(fake: Any) -> None:
    """Two data arrays with one id are refused."""
    first = entry("a", "vertex")
    second = entry("b", "vertex", id="a")
    assert codes(fake.accept_frame(points_frame([first, second]))) == ["invalid_object"]


def test_revision_must_be_positive_integer(fake: Any) -> None:
    """A revision of 0 is refused."""
    header = {
        "type": "scene.upsert",
        "id": "p",
        "object_type": "Points",
        "revision": 0,
        "geometry": {"vertices": "v"},
    }
    frame = encode_frame(header, [("v", np.zeros((1, 3), dtype=np.float32))])
    assert codes(fake.accept_frame(frame)) == ["invalid_object"]


def test_valid_vectors_are_accepted(fake: Any) -> None:
    """A Vectors object with orientation and lineation is stored."""
    vectors = Vectors(
        np.zeros((2, 3)),
        orientation=([10.0, 20.0], [30.0, 40.0]),
        lineation=([1.0, 2.0], [3.0, 4.0]),
    )
    data = encode_object(
        vectors, id="v", name=None, parent=None, revision=1, origin=(0.0, 0.0, 0.0)
    )
    assert fake.accept_frame(data) == []
    assert "v" in fake.objects


# -- the tool endpoint -------------------------------------------------------------------


# The protocol names of the dtypes, by numpy name.
PROTOCOL_DTYPES = {
    "float32": "f32",
    "float64": "f64",
    "uint32": "u32",
    "int32": "i32",
    "uint8": "u8",
}


def _upsert(viewer: FakeViewer, obj: Any, object_id: str, revision: int) -> None:
    """Send one object to the data side of the fake viewer."""
    data = encode_object(
        obj, id=object_id, name=None, parent=None, revision=revision, origin=(0.0, 0.0, 0.0)
    )
    assert viewer.accept_frame(data) == []


@contextlib.contextmanager
def tool_session(**options: Any) -> Iterator[tuple[FakeViewer, WebSocketClient]]:
    """Start a fake viewer with objects and open a client on its tool endpoint.

    The scene has ``pts`` (4 vertices), ``srf`` (revision 2) and ``srf2`` (4 vertices,
    2 triangles), ``blk`` (2 x 3 x 4 cells), ``dh`` (Drillhole) and ``grp`` (Group).
    """
    viewer = FakeViewer(**options)
    client = None
    try:
        viewer.accept_frame(hello_frame("EPSG:28350", (0.0, 0.0, 0.0), "test"))
        vertices = np.arange(12, dtype=np.float64).reshape(4, 3)
        _upsert(viewer, Points(vertices), "pts", 1)
        _upsert(viewer, Surface(vertices, [[0, 1, 2], [1, 2, 3]]), "srf", 2)
        _upsert(viewer, Surface(vertices, [[0, 1, 2], [1, 2, 3]]), "srf2", 1)
        _upsert(viewer, BlockModel([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2, 3, 4]), "blk", 1)
        hole = Drillhole([0.0, 0.0, 0.0], [[0.0, 0.0, 0.0], [10.0, 0.0, -90.0]], [[0.0, 5.0]])
        _upsert(viewer, hole, "dh", 1)
        assert viewer.accept_frame(group_frame("grp", None, None, 1)) == []
        client = WebSocketClient.connect("127.0.0.1", viewer.port, token="secret", path="/tool")
        yield viewer, client
    finally:
        if client is not None:
            client.close()
        viewer.close()


def reply(client: WebSocketClient) -> Frame:
    """Return the next frame from the tool endpoint. Fail if none arrives in 2 s."""
    data = client.receive(2.0)
    assert data is not None, "The tool endpoint sent no reply."
    return decode_frame(data)


def replies_until_done(client: WebSocketClient) -> list[Frame]:
    """Return the frames up to and including ``object.end`` or an error."""
    frames: list[Frame] = []
    while True:
        frame = reply(client)
        frames.append(frame)
        if frame.header["type"] in ("object.end", "request.error", "error"):
            return frames


def raw_set(entries: list[dict[str, str]], buffers: list[tuple[str, np.ndarray]]) -> bytes:
    """Return a hand-made ``selection.set`` frame for object ``pts``."""
    header = {
        "type": "selection.set",
        "request": 7,
        "object": "pts",
        "revision": 1,
        "entries": entries,
    }
    return encode_frame(header, buffers, version=protocol.TOOL_PROTOCOL_VERSION)


def set_selection(
    client: WebSocketClient,
    object_id: str,
    revision: int,
    *,
    cells: list[int] | None = None,
    nodes: list[int] | None = None,
    request: int = 1,
) -> Frame:
    """Send ``selection.set`` and return the reply."""
    cell_array = None if cells is None else np.array(cells, dtype=np.uint32)
    node_array = None if nodes is None else np.array(nodes, dtype=np.uint32)
    client.send_binary(selection_set_frame(request, object_id, cell_array, node_array, revision))
    return reply(client)


def indices(selection: dict[str, Any] | None, kind: str) -> list[int]:
    """Return the sorted unique indices of one kind of a stored selection."""
    assert selection is not None
    values = selection.get(kind)
    return [] if values is None else sorted({int(v) for v in values})


def error_of(frame: Frame) -> dict[str, Any]:
    """Return the ``error`` object of a ``request.error`` frame."""
    assert frame.header["type"] == "request.error", frame.header
    return frame.header["error"]


def test_protocol_minor_17_answers_unknown_message() -> None:
    """With protocol_minor below 18, each selection or fetch message gets an error frame."""
    with tool_session(protocol_minor=17) as (_, client):
        frames = [
            selection_set_frame(1, "pts", None, np.array([0], dtype=np.uint32), 1),
            selection_get_frame(2),
            fetch_frame(3, "pts"),
        ]
        for data in frames:
            client.send_binary(data)
            answer = reply(client)
            assert answer.header["type"] == "error"
            assert answer.header["code"] == "unknown_message"


def test_headless_viewer_needs_a_window_for_selection() -> None:
    """A headless viewer answers selection.set and selection.get with needs_window."""
    with tool_session(kind="headless") as (_, client):
        client.send_binary(selection_set_frame(1, "pts", None, np.array([0], dtype=np.uint32), 1))
        first = reply(client)
        assert error_of(first)["code"] == "needs_window"
        assert first.header["request"] == 1
        client.send_binary(selection_get_frame(2))
        second = reply(client)
        assert error_of(second)["code"] == "needs_window"
        assert second.header["request"] == 2


def test_selection_set_stores_selection_and_counts_unique_indices() -> None:
    """selection.set replies with the counts. A repeated index counts one time."""
    with tool_session() as (viewer, client):
        answer = set_selection(client, "srf", 2, cells=[0, 1, 1], nodes=[3, 3, 0], request=5)
        assert answer.header["type"] == "selection.result"
        assert answer.header["request"] == 5
        assert answer.header["object"] == "srf"
        assert answer.header["revision"] == 2
        assert answer.header["version"] == protocol.TOOL_PROTOCOL_VERSION
        counts = {e["kind"]: e["count"] for e in answer.header["entries"]}
        assert counts == {"cells": 2, "nodes": 2}
        assert viewer.selection is not None
        assert viewer.selection["object"] == "srf"
        assert viewer.selection["revision"] == 2
        assert indices(viewer.selection, "cells") == [0, 1]
        assert indices(viewer.selection, "nodes") == [0, 3]


def test_selection_set_counts_block_model_nodes_and_cells() -> None:
    """A BlockModel has (nx+1)(ny+1)(nz+1) nodes and nx*ny*nz cells."""
    with tool_session() as (_, client):
        ok = set_selection(client, "blk", 1, cells=[23], nodes=[59])
        assert ok.header["type"] == "selection.result"
        bad_cell = set_selection(client, "blk", 1, cells=[24])
        assert error_of(bad_cell)["field"] == "indices"
        bad_node = set_selection(client, "blk", 1, nodes=[60])
        assert error_of(bad_node)["field"] == "indices"


def test_selection_set_unknown_object_is_invalid_input() -> None:
    """An object id that the viewer lacks gives invalid_input with the field object."""
    with tool_session() as (viewer, client):
        error = error_of(set_selection(client, "missing", 1, nodes=[0], request=9))
        assert error["code"] == "invalid_input"
        assert error["field"] == "object"
        assert isinstance(error["message"], str)
        assert viewer.selection is None


def test_selection_set_wrong_revision_is_invalid_input() -> None:
    """A revision that is not the revision of the object gives the field revision."""
    with tool_session() as (_, client):
        error = error_of(set_selection(client, "srf", 1, nodes=[0]))
        assert error["code"] == "invalid_input"
        assert error["field"] == "revision"


def test_selection_set_kind_the_object_lacks_is_invalid_input() -> None:
    """Points have no cells. Drillhole and Group have no table at all."""
    with tool_session() as (_, client):
        for object_id, kwargs in (
            ("pts", {"cells": [0]}),
            ("dh", {"nodes": [0]}),
            ("dh", {"cells": [0]}),
            ("grp", {"nodes": [0]}),
        ):
            error = error_of(set_selection(client, object_id, 1, **kwargs))
            assert error["code"] == "invalid_input", (object_id, kwargs)
            assert error["field"] == "kind", (object_id, kwargs)


def test_selection_set_index_at_or_above_count_is_invalid_input() -> None:
    """An index equal to the count is refused. The last valid index is accepted."""
    with tool_session() as (viewer, client):
        error = error_of(set_selection(client, "pts", 1, nodes=[0, 4]))
        assert error["code"] == "invalid_input"
        assert error["field"] == "indices"
        assert viewer.selection is None
        error = error_of(set_selection(client, "srf", 2, cells=[2]))
        assert error["field"] == "indices"
        ok = set_selection(client, "srf", 2, cells=[1], nodes=[3])
        assert ok.header["type"] == "selection.result"


def test_selection_set_buffer_not_u32_or_missing_is_invalid_input() -> None:
    """A buffer of another type, or an entry without its buffer, gives the field indices."""
    entries = [{"kind": "nodes", "indices": "nodes"}]
    with tool_session() as (_, client):
        client.send_binary(raw_set(entries, [("nodes", np.array([0, 1], dtype=np.int32))]))
        wrong_type = reply(client)
        assert error_of(wrong_type)["code"] == "invalid_input"
        assert error_of(wrong_type)["field"] == "indices"
        assert wrong_type.header["request"] == 7
        client.send_binary(raw_set(entries, []))
        missing = reply(client)
        assert error_of(missing)["code"] == "invalid_input"
        assert error_of(missing)["field"] == "indices"


def test_new_selection_replaces_the_old_one_also_on_another_object() -> None:
    """A second selection.set replaces the first, on the same object or on another."""
    with tool_session() as (viewer, client):
        set_selection(client, "srf", 2, nodes=[0, 1])
        set_selection(client, "srf", 2, nodes=[2])
        assert indices(viewer.selection, "nodes") == [2]
        set_selection(client, "srf2", 1, cells=[0])
        assert viewer.selection is not None
        assert viewer.selection["object"] == "srf2"
        assert indices(viewer.selection, "nodes") == []
        assert indices(viewer.selection, "cells") == [0]


def test_empty_entries_clear_the_selection() -> None:
    """A selection.set with empty entries clears the selection and replies with no entries."""
    with tool_session() as (viewer, client):
        set_selection(client, "srf", 2, nodes=[0])
        answer = set_selection(client, "srf", 2)
        assert answer.header["type"] == "selection.result"
        assert answer.header["entries"] == []
        assert viewer.selection is None


def test_selection_get_returns_one_u32_entry_for_each_kind() -> None:
    """selection.get replies with object, revision and one u32 buffer for each kind."""
    with tool_session() as (_, client):
        set_selection(client, "srf", 2, cells=[1], nodes=[0, 2])
        client.send_binary(selection_get_frame(4))
        answer = reply(client)
        assert answer.header["type"] == "selection.result"
        assert answer.header["request"] == 4
        assert answer.header["object"] == "srf"
        assert answer.header["revision"] == 2
        assert answer.header["version"] == protocol.TOOL_PROTOCOL_VERSION
        entries = {e["kind"]: e for e in answer.header["entries"]}
        assert set(entries) == {"cells", "nodes"}
        for kind, expected in (("cells", [1]), ("nodes", [0, 2])):
            assert entries[kind]["indices"] == kind
            buffer = answer.buffers[kind]
            assert buffer.dtype == np.dtype("uint32")
            assert sorted(int(v) for v in buffer) == expected


def test_selection_get_without_selection_has_no_entries_and_no_object() -> None:
    """With no selection, selection.get replies with empty entries and no object."""
    with tool_session() as (viewer, client):
        client.send_binary(selection_get_frame(1))
        answer = reply(client)
        assert answer.header["type"] == "selection.result"
        assert answer.header["entries"] == []
        assert "object" not in answer.header
        viewer.set_user_selection("pts", nodes=[1, 2])
        client.send_binary(selection_get_frame(2))
        after = reply(client)
        assert after.header["object"] == "pts"
        assert [e["kind"] for e in after.header["entries"]] == ["nodes"]


def test_fetch_unknown_id_is_invalid_input() -> None:
    """object.fetch for an id that the viewer lacks gives invalid_input with the field id."""
    with tool_session() as (_, client):
        client.send_binary(fetch_frame(3, "missing"))
        answer = reply(client)
        assert answer.header["request"] == 3
        assert error_of(answer)["code"] == "invalid_input"
        assert error_of(answer)["field"] == "id"


def test_fetch_sends_header_parts_and_end_in_order() -> None:
    """object.fetch sends the header, the parts of each array in order, then the end."""
    with tool_session(part_bytes=16) as (viewer, client):
        client.send_binary(fetch_frame(8, "srf"))
        frames = replies_until_done(client)
        assert all(f.header["request"] == 8 for f in frames)
        assert all(f.header["version"] == protocol.TOOL_PROTOCOL_VERSION for f in frames)
        head = frames[0]
        assert head.header["type"] == "object.header"
        assert head.header["crs"] == viewer.crs
        assert head.header["origin"] == list(viewer.origin)
        stored = viewer.objects["srf"]
        assert head.header["object"]["buffers"] == []
        assert head.header["object"]["id"] == "srf"
        arrays = head.header["arrays"]
        assert [a["name"] for a in arrays] == list(stored.buffers)
        for described in arrays:
            source = stored.buffers[described["name"]]
            assert described["dtype"] == PROTOCOL_DTYPES[source.dtype.name]
            assert described["shape"] == [int(n) for n in source.shape]
            assert all(type(n) is int for n in described["shape"])
        parts = [f for f in frames if f.header["type"] == "object.part"]
        assert frames[-1].header["type"] == "object.end"
        assert frames[-1].header["parts"] == len(parts)
        assert len(parts) > 1
        position = 0
        for described in arrays:
            raw = stored.buffers[described["name"]].tobytes()
            joined = b""
            while len(joined) < len(raw):
                part = parts[position]
                position += 1
                assert part.header["array"] == described["name"]
                assert part.header["offset"] == len(joined)
                assert part.header["offset"] % 8 == 0
                data = part.buffers["part"]
                assert data.dtype == np.dtype("uint8")
                assert 0 < part.header["bytes"] <= 16
                assert part.header["bytes"] == data.nbytes
                joined += data.tobytes()
            assert joined == raw
        assert position == len(parts)


def test_fetch_with_more_than_one_part_sends_progress_after_each_part() -> None:
    """command.progress follows each object.part when there is more than one part."""
    with tool_session(part_bytes=16) as (_, client):
        client.send_binary(fetch_frame(2, "pts"))
        frames = replies_until_done(client)
        types = [f.header["type"] for f in frames]
        part_count = types.count("object.part")
        assert part_count > 1
        for index, kind in enumerate(types):
            if kind == "object.part":
                progress = frames[index + 1]
                assert progress.header["type"] == "command.progress"
                assert 0.0 <= progress.header["fraction"] <= 1.0
                assert isinstance(progress.header["text"], str)
        assert types.count("command.progress") == part_count


def test_fetch_cancel_gives_cancelled_error() -> None:
    """A command.cancel with the request number stops a slow fetch with cancelled."""
    with tool_session(part_bytes=8, fetch_delay=0.05) as (_, client):
        client.send_binary(fetch_frame(11, "blk"))
        client.send_binary(cancel_frame(11))
        frames = replies_until_done(client)
        last = frames[-1]
        assert last.header["request"] == 11
        assert error_of(last)["code"] == "cancelled"
        assert "object.end" not in [f.header["type"] for f in frames]


def test_fetch_leaves_scene_and_cursor_unchanged() -> None:
    """object.fetch does not change the objects and does not change the selection."""
    with tool_session(part_bytes=16) as (viewer, client):
        set_selection(client, "srf", 2, nodes=[1])

        def snapshot() -> dict[str, Any]:
            return {
                k: (copy.deepcopy(v.header), {n: a.copy() for n, a in v.buffers.items()})
                for k, v in viewer.objects.items()
            }

        before = snapshot()
        client.send_binary(fetch_frame(1, "srf"))
        replies_until_done(client)
        after = snapshot()
        assert set(after) == set(before)
        for key, (header, arrays) in after.items():
            assert header == before[key][0]
            assert list(arrays) == list(before[key][1])
            for name, array in arrays.items():
                assert np.array_equal(array, before[key][1][name])
        assert viewer.selection is not None
        assert viewer.selection["object"] == "srf"
        assert indices(viewer.selection, "nodes") == [1]


def test_fetch_with_one_part_sends_no_progress() -> None:
    """When all data fits in one part, object.fetch sends no command.progress frame."""
    with tool_session() as (_, client):
        client.send_binary(fetch_frame(4, "pts"))
        frames = replies_until_done(client)
        assert [f.header["type"] for f in frames] == ["object.header", "object.part", "object.end"]
        assert frames[-1].header["parts"] == 1


def test_tool_requests_keep_the_frames_of_the_tool_endpoint() -> None:
    """The fake keeps each tool message in tool_requests."""
    with tool_session() as (viewer, client):
        client.send_binary(selection_get_frame(1))
        reply(client)
        assert [f.header["type"] for f in viewer.tool_requests] == ["selection.get"]
