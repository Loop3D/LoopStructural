"""Tests for the validation rules of the fake viewer."""

from __future__ import annotations

import contextlib
import copy
import time
from collections.abc import Iterator
from typing import Any

import numpy as np
import pytest

from loop_link import protocol
from loop_link._websocket import WebSocketClient
from loop_link.encode import (
    cancel_frame,
    command_frame,
    commands_list_frame,
    encode_object,
    fetch_frame,
    group_frame,
    hello_frame,
    selection_get_frame,
    selection_set_frame,
)
from loop_link.objects import BlockModel, Categorical, Drillhole, Points, Surface, Vectors
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


# -- the commands endpoint ---------------------------------------------------------------

COMMAND_NAMES = [
    "filter.edit",
    "colour.set_active",
    "colour.set_map",
    "colour.set_range",
    "colour.clear_range",
    "colour.set_scale",
    "object.set_colour",
]


def command_points() -> Points:
    """Return 3 points with a float, an integer, a referenced, a boolean and a text array."""
    return Points(
        np.arange(9, dtype=np.float64).reshape(3, 3),
        {
            "val": np.array([1.0, 2.0, 3.0]),
            "num": np.array([1, 2, 3]),
            "cat": Categorical([0, 1, 0], ["a", "b"]),
            "flag": np.array([True, False, True]),
            "txt": np.array(["a", "b", "c"]),
        },
    )


@contextlib.contextmanager
def command_session(**options: Any) -> Iterator[FakeViewer]:
    """Start a fake viewer that has the object ``obj`` (see ``command_points``)."""
    viewer = FakeViewer(**options)
    try:
        viewer.accept_frame(hello_frame("EPSG:28350", (0.0, 0.0, 0.0), "test"))
        _upsert(viewer, command_points(), "obj", 1)
        yield viewer
    finally:
        viewer.close()


def run_command(viewer: FakeViewer, command: str, **fields: Any) -> Frame:
    """Send ``command.run`` to the fake viewer and return the one reply frame."""
    replies = viewer.tool_frame(command_frame(3, command, fields))
    assert len(replies) == 1, replies
    return replies[0]


def result_of(frame: Frame) -> dict[str, Any]:
    """Return the ``result`` of a ``command.result`` frame."""
    assert frame.header["type"] == "command.result", frame.header
    assert frame.header["request"] == 3
    return frame.header["result"]


def edit(viewer: FakeViewer, kind: str, **fields: Any) -> Frame:
    """Send a ``filter.edit`` command for the object ``obj``."""
    return run_command(viewer, "filter.edit", object="obj", edit={"kind": kind, **fields})


def test_fake_command_list_answers_the_seven_commands() -> None:
    """commands.list names the kind of the viewer and the 7 commands, with all keys."""
    with command_session(kind="headless") as viewer:
        (answer,) = viewer.tool_frame(commands_list_frame(4))
        assert answer.header["type"] == "commands"
        assert answer.header["request"] == 4
        assert answer.header["viewer"] == "headless"
        entries = answer.header["commands"]
        assert sorted(e["name"] for e in entries) == sorted(COMMAND_NAMES)
        keys = {
            "name",
            "description",
            "needs_window",
            "tool_available",
            "unavailable_reason",
            "input_schema",
            "result_schema",
        }
        for item in entries:
            assert keys <= set(item)
            assert item["needs_window"] == (item["name"] == "colour.set_scale")


def test_fake_command_list_colour_map_enum_depends_on_protocol_minor() -> None:
    """The names of colour.set_map are listed from minor 18. Below 18 the property is text."""

    def name_property(minor: int) -> dict[str, Any]:
        with command_session(protocol_minor=minor) as viewer:
            (answer,) = viewer.tool_frame(commands_list_frame(1))
        entry = next(e for e in answer.header["commands"] if e["name"] == "colour.set_map")
        return entry["input_schema"]["properties"]["name"]

    assert name_property(18)["enum"] == ["viridis", "grey", "coolwarm"]
    assert name_property(17) == {"type": "string"}


def test_fake_command_below_minor_11_is_unknown_message() -> None:
    """With protocol_minor below 11, both command messages give unknown_message."""
    with command_session(protocol_minor=10) as viewer:
        for data in (commands_list_frame(1), command_frame(2, "colour.clear_range", {"id": "obj"})):
            (answer,) = viewer.tool_frame(data)
            assert answer.header["type"] == "error"
            assert answer.header["code"] == "unknown_message"


def test_fake_command_run_checks_command_and_input() -> None:
    """An unknown command, an unknown or missing field or a wrong type is invalid_input."""
    with command_session() as viewer:
        assert error_of(run_command(viewer, "no.such", id="obj"))["code"] == "invalid_input"
        for fields, field_name in (
            ({"id": "obj", "extra": 1}, "extra"),
            ({}, "id"),
            ({"id": 5}, "id"),
        ):
            error = error_of(run_command(viewer, "colour.clear_range", **fields))
            assert error["code"] == "invalid_input", fields
            assert error["field"] == field_name, fields
        wrong = error_of(run_command(viewer, "colour.set_range", id="obj", min=1, max="2"))
        assert (wrong["code"], wrong["field"]) == ("invalid_input", "min")


def test_fake_command_unknown_object_names_the_id_field() -> None:
    """An unknown object gives the field id, or object for colour.set_scale and filter.edit."""
    with command_session() as viewer:
        error = error_of(run_command(viewer, "colour.clear_range", id="missing"))
        assert (error["code"], error["field"]) == ("invalid_input", "id")
        error = error_of(run_command(viewer, "colour.set_scale", object="missing", scale="log"))
        assert (error["code"], error["field"]) == ("invalid_input", "object")
        error = error_of(
            run_command(viewer, "filter.edit", object="missing", edit={"kind": "cancel"})
        )
        assert (error["code"], error["field"]) == ("invalid_input", "object")


def test_fake_command_headless_needs_a_window_for_set_scale_only() -> None:
    """A headless viewer refuses colour.set_scale with needs_window. Other commands run."""
    with command_session(kind="headless") as viewer:
        error = error_of(run_command(viewer, "colour.set_scale", object="obj", scale="log"))
        assert error["code"] == "needs_window"
        assert result_of(run_command(viewer, "colour.clear_range", id="obj")) == {"status": "done"}


def test_fake_filter_add_and_set_commands_change_the_filter() -> None:
    """add makes an off filter. The set commands change it. remove deletes it."""
    with command_session() as viewer:
        assert result_of(edit(viewer, "add", array="val")) == {"status": "done"}
        assert viewer.filters["obj"]["val"]["kind"] == "range"
        assert viewer.filters["obj"]["val"]["on"] is False
        assert result_of(edit(viewer, "set_on", array="val", on=True)) == {"status": "done"}
        edit(viewer, "set_range", array="val", min="1.5", max="2.5", keep_nan=True)
        state = viewer.filters["obj"]["val"]
        assert (state["on"], state["minimum"], state["maximum"], state["keep_nan"]) == (
            True,
            1.5,
            2.5,
            True,
        )
        edit(viewer, "add", array="cat")
        assert viewer.filters["obj"]["cat"]["kind"] == "values"
        edit(viewer, "set_values", array="cat", values=[1])
        edit(viewer, "set_exclude", array="cat", exclude=True)
        state = viewer.filters["obj"]["cat"]
        assert (list(state["values"]), state["exclude"]) == ([1], True)
        edit(viewer, "remove", array="val")
        assert list(viewer.filters["obj"]) == ["cat"]


def test_fake_filter_property_returns_a_copy() -> None:
    """Changing the dict from .filters does not change the viewer."""
    with command_session() as viewer:
        edit(viewer, "add", array="val")
        viewer.filters["obj"]["val"]["on"] = True
        viewer.filters["obj"].clear()
        assert viewer.filters["obj"]["val"]["on"] is False


def test_fake_filter_set_without_filter_is_refused() -> None:
    """A set command for an array that has no filter gives refused."""
    with command_session() as viewer:
        result = result_of(edit(viewer, "set_on", array="val", on=True))
        assert result["status"] == "refused"
        assert "val" in result["text"]
        assert viewer.filters.get("obj", {}) == {}


def test_fake_filter_refuses_arrays_and_edits_that_do_not_fit() -> None:
    """A filter on a missing or text array, or an edit of the wrong kind, gives refused."""
    with command_session() as viewer:
        for array in ("missing", "txt"):
            result = result_of(edit(viewer, "add", array=array))
            assert result["status"] == "refused", array
            assert array in result["text"]
        for array in ("cat", "flag", "val", "num"):
            edit(viewer, "add", array=array)
        for array in ("cat", "flag"):
            result = result_of(
                edit(viewer, "set_range", array=array, min="0", max="1", keep_nan=False)
            )
            assert result["status"] == "refused", array
            assert array in result["text"]
        for array in ("val", "num"):
            result = result_of(edit(viewer, "set_values", array=array, values=[0]))
            assert result["status"] == "refused", array
            assert array in result["text"]


def test_fake_filter_busy_time_makes_commands_wait_and_cancel_ends_it() -> None:
    """After a filter change, commands wait busy_wait s, then fail with busy."""
    with command_session(filter_delay=30.0, busy_wait=0.05) as viewer:
        edit(viewer, "add", array="val")
        started = time.monotonic()
        error = error_of(run_command(viewer, "colour.set_active", id="obj", array="val"))
        assert error["code"] == "busy"
        assert time.monotonic() - started >= 0.04
        started = time.monotonic()
        error = error_of(run_command(viewer, "colour.set_scale", object="obj", scale="log"))
        assert error["code"] == "busy"
        assert time.monotonic() - started < 0.04
        assert result_of(edit(viewer, "cancel")) == {"status": "done"}
        assert result_of(run_command(viewer, "colour.clear_range", id="obj")) == {"status": "done"}


def test_fake_filter_command_runs_when_busy_time_ends_inside_the_wait() -> None:
    """A waiting command runs when the busy time ends before busy_wait."""
    with command_session(filter_delay=0.1, busy_wait=5.0) as viewer:
        edit(viewer, "add", array="val")
        assert result_of(run_command(viewer, "colour.clear_range", id="obj")) == {"status": "done"}


def test_fake_colour_commands_fill_the_colouring() -> None:
    """The colour commands store their values in .colouring."""
    with command_session() as viewer:
        run_command(viewer, "colour.set_active", id="obj", array="val")
        run_command(viewer, "colour.set_map", id="obj", name="grey")
        run_command(viewer, "colour.set_range", id="obj", min="1", max="3")
        run_command(viewer, "colour.set_scale", object="obj", scale="log")
        state = viewer.colouring["obj"]
        assert state["array"] == "val"
        assert state["colour_map"] == "grey"
        assert tuple(float(v) for v in state["range"]) == (1.0, 3.0)
        assert state["scale"] == "log"
        run_command(viewer, "colour.clear_range", id="obj")
        assert viewer.colouring["obj"]["range"] is None
        run_command(viewer, "object.set_colour", id="obj", colour=[1, 2, 3])
        assert list(viewer.colouring["obj"]["plain"]) == [1, 2, 3]
        run_command(viewer, "object.set_colour", id="obj", colour=None)
        assert viewer.colouring["obj"]["plain"] is None
        viewer.colouring["obj"]["scale"] = "changed"
        assert viewer.colouring["obj"]["scale"] == "log"


def test_fake_colour_set_active_unknown_array_is_refused() -> None:
    """colour.set_active with an array that the object lacks gives refused."""
    with command_session() as viewer:
        result = result_of(run_command(viewer, "colour.set_active", id="obj", array="missing"))
        assert result["status"] == "refused"
        assert "missing" in result["text"]


def test_fake_colour_map_range_and_scale_are_refused_on_referenced_or_boolean_array() -> None:
    """With a referenced or boolean active array, three colour commands give refused."""
    with command_session() as viewer:
        for array in ("cat", "flag"):
            run_command(viewer, "colour.set_active", id="obj", array=array)
            for command, fields in (
                ("colour.set_map", {"id": "obj", "name": "grey"}),
                ("colour.set_range", {"id": "obj", "min": "0", "max": "1"}),
                ("colour.set_scale", {"object": "obj", "scale": "log"}),
            ):
                result = result_of(run_command(viewer, command, **fields))
                assert result["status"] == "refused", (array, command)
                assert array in result["text"]


def test_fake_colour_and_filters_stay_after_a_new_revision() -> None:
    """A new scene.upsert keeps .filters and .colouring and drops a filter of a lost array."""
    with command_session() as viewer:
        edit(viewer, "add", array="val")
        edit(viewer, "add", array="num")
        run_command(viewer, "colour.set_active", id="obj", array="val")
        points = Points(
            np.arange(9, dtype=np.float64).reshape(3, 3), {"val": np.array([1.0, 2.0, 3.0])}
        )
        _upsert(viewer, points, "obj", 2)
        assert list(viewer.filters["obj"]) == ["val"]
        assert viewer.colouring["obj"]["array"] == "val"


# -- the long filter.edit (protocol 1.20) ---------------------------------------------------


@contextlib.contextmanager
def filter_long_session(**options: Any) -> Iterator[tuple[FakeViewer, WebSocketClient]]:
    """Start a fake viewer with the object ``blk`` and open a client on its tool endpoint.

    ``blk`` is a block model of 1,000 cells. The cell property ``grade`` holds the
    values 0 to 999. The viewer uses protocol 1.20 and ``filter_delay`` 0.4 s, unless
    the options say otherwise.
    """
    options = {"protocol_minor": 20, "filter_delay": 0.4, **options}
    viewer = FakeViewer(**options)
    client = None
    try:
        viewer.accept_frame(hello_frame("EPSG:28350", (0.0, 0.0, 0.0), "test"))
        grade = np.arange(1000, dtype=np.float32)
        model = BlockModel(
            [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [10, 10, 10], cell_properties={"grade": grade}
        )
        _upsert(viewer, model, "blk", 1)
        client = WebSocketClient.connect("127.0.0.1", viewer.port, token="secret", path="/tool")
        yield viewer, client
    finally:
        if client is not None:
            client.close()
        viewer.close()


def send_filter_edit(client: WebSocketClient, request: int, kind: str, **fields: Any) -> None:
    """Send a ``filter.edit`` request for ``blk`` on the tool endpoint."""
    edit_fields = {"kind": kind, "array": "grade", **fields}
    client.send_binary(
        command_frame(request, "filter.edit", {"object": "blk", "edit": edit_fields})
    )


def frames_until_end(client: WebSocketClient) -> list[Frame]:
    """Return the frames up to and including a ``command.result`` or an error."""
    frames: list[Frame] = []
    while True:
        frame = reply(client)
        frames.append(frame)
        if frame.header["type"] in ("command.result", "request.error", "error"):
            return frames


def filter_edit_frames(
    client: WebSocketClient, request: int, kind: str, **fields: Any
) -> list[Frame]:
    """Send one ``filter.edit`` request and return all frames of its reply."""
    send_filter_edit(client, request, kind, **fields)
    return frames_until_end(client)


def prepare_filter(client: WebSocketClient) -> None:
    """Add a filter for ``grade`` and set its range 0 to 99. The filter stays off."""
    assert filter_edit_frames(client, 1, "add")[-1].header["result"] == {"status": "done"}
    frames = filter_edit_frames(client, 2, "set_range", min="0.0", max="99.0", keep_nan=False)
    assert frames[-1].header["result"] == {"status": "done"}


def test_filter_long_set_on_sends_progress_and_the_count() -> None:
    """set_on sends command.progress frames, then a result with passing and total."""
    with filter_long_session() as (_, client):
        prepare_filter(client)
        frames = filter_edit_frames(client, 3, "set_on", on=True)
        progress = [f for f in frames if f.header["type"] == "command.progress"]
        assert progress
        for frame in progress:
            assert frame.header["request"] == 3
            assert 0.0 <= frame.header["fraction"] <= 1.0
            assert frame.header["text"] == "filter.edit"
        assert progress[0].header["fraction"] == 0.0
        assert progress[-1].header["fraction"] == 1.0
        last = frames[-1]
        assert last.header["type"] == "command.result"
        assert last.header["request"] == 3
        assert last.header["result"] == {"status": "done", "passing": 100, "total": 1000}


def test_filter_long_set_on_false_passes_all_cells() -> None:
    """set_on with on false gives passing 1000 and total 1000."""
    with filter_long_session() as (_, client):
        prepare_filter(client)
        filter_edit_frames(client, 3, "set_on", on=True)
        frames = filter_edit_frames(client, 4, "set_on", on=False)
        assert frames[-1].header["type"] == "command.result"
        assert frames[-1].header["result"] == {"status": "done", "passing": 1000, "total": 1000}


def test_filter_long_remove_and_add_send_no_progress_and_no_count() -> None:
    """remove and add give only a result with status done."""
    with filter_long_session() as (_, client):
        prepare_filter(client)
        for request, kind in ((3, "remove"), (4, "add")):
            frames = filter_edit_frames(client, request, kind)
            assert [f.header["type"] for f in frames] == ["command.result"], kind
            assert frames[0].header["result"] == {"status": "done"}, kind


def test_filter_long_cancel_gives_cancelled_and_keeps_the_old_filter() -> None:
    """command.cancel of a running set_on gives cancelled. The filter stays off."""
    with filter_long_session() as (viewer, client):
        prepare_filter(client)
        send_filter_edit(client, 3, "set_on", on=True)
        client.send_binary(cancel_frame(3))
        frames = frames_until_end(client)
        last = frames[-1]
        assert last.header["request"] == 3
        assert error_of(last)["code"] == "cancelled"
        assert "command.result" not in [f.header["type"] for f in frames]
        assert viewer.filters["blk"]["grade"]["on"] is False


def test_filter_long_stop_after_gives_state_building() -> None:
    """With filter_stop_after the set_on result is state building after about 0.1 s."""
    with filter_long_session(filter_stop_after=0.1) as (_, client):
        prepare_filter(client)
        started = time.monotonic()
        frames = filter_edit_frames(client, 3, "set_on", on=True)
        elapsed = time.monotonic() - started
        assert frames[-1].header["type"] == "command.result"
        assert frames[-1].header["result"] == {"status": "done", "state": "building"}
        assert 0.05 <= elapsed < 0.35


def test_filter_long_commands_list_marks_only_filter_edit_long() -> None:
    """With minor 20, the commands frame has long true for filter.edit and no other long key."""
    with filter_long_session() as (_, client):
        client.send_binary(commands_list_frame(5))
        answer = reply(client)
        assert answer.header["type"] == "commands"
        for item in answer.header["commands"]:
            if item["name"] == "filter.edit":
                assert item["long"] is True
            else:
                assert "long" not in item, item["name"]


def test_filter_long_below_minor_20_has_no_progress_and_no_long_key() -> None:
    """With minor 18, set_on gives only status done, and the commands have no long key."""
    with filter_long_session(protocol_minor=18) as (_, client):
        prepare_filter(client)
        frames = filter_edit_frames(client, 3, "set_on", on=True)
        assert [f.header["type"] for f in frames] == ["command.result"]
        assert frames[0].header["result"] == {"status": "done"}
        client.send_binary(commands_list_frame(6))
        answer = reply(client)
        assert all("long" not in item for item in answer.header["commands"])
