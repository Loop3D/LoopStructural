"""Check the frames of the geoviewer repository against loop_link (SC-004)."""

from __future__ import annotations

import copy
import json
from typing import Any

import numpy as np
import pytest
from conftest import geoviewer_fixtures

from loop_link.encode import (
    cancel_frame,
    command_frame,
    commands_list_frame,
    encode_object,
    fetch_frame,
    hello_frame,
    selection_get_frame,
    selection_set_frame,
)
from loop_link.objects import Drillhole
from loop_link.protocol import decode_frame
from loop_link.testing import FakeViewer

_SKIP = "set LOOP_LINK_GEOVIEWER_FIXTURES to the folder fixtures/protocol/python of the geoviewer repository"
_DRILLHOLE_FILE = "scene_upsert_drillhole_end_of_hole.bin"


def _folder():
    """Return the fixture folder, or skip the test when it is not set."""
    folder = geoviewer_fixtures()
    if folder is None:
        pytest.skip(_SKIP)
    return folder


def _by_role(header: dict[str, Any], buffers: dict[str, np.ndarray]):
    """Return a header and buffers with each buffer name replaced by its role."""
    name_to_role = {name: role for role, name in header["geometry"].items()}
    result = copy.deepcopy(header)
    result.pop("version", None)
    result["geometry"] = {role: role for role in header["geometry"]}
    for entry in result["buffers"]:
        entry["name"] = name_to_role[entry["name"]]
    return result, {name_to_role[name]: array for name, array in buffers.items()}


def test_viewer_accepts_every_frame() -> None:
    """The fake viewer gives no error reply for any frame of the folder."""
    files = sorted(_folder().glob("*.bin"))
    assert files
    for path in files:
        viewer = FakeViewer()
        try:
            viewer.accept_frame(hello_frame("EPSG:28350", (0.0, 0.0, 0.0), "0"))
            replies = viewer.accept_frame(path.read_bytes())
        finally:
            viewer.close()
        errors = [r.header for r in replies if r.header.get("type") == "error"]
        assert not errors, f"{path.name}: {errors}"


def test_frame_keys_match_encoder() -> None:
    """Each plugin frame has the keys type, object_type and geometry that a frame needs."""
    for path in sorted(_folder().glob("*.bin")):
        if path.name == _DRILLHOLE_FILE:
            continue
        header = decode_frame(path.read_bytes()).header
        if header.get("type") != "scene.upsert":
            continue
        assert isinstance(header["object_type"], str), path.name
        assert isinstance(header["geometry"], dict), path.name


def test_drillhole_frame_matches_encoder() -> None:
    """The encoder builds the drillhole frame of the geoviewer repository."""
    path = _folder() / _DRILLHOLE_FILE
    if not path.exists():
        pytest.skip(f"{_DRILLHOLE_FILE} is not in the folder")
    real = decode_frame(path.read_bytes())
    header = real.header
    buffers = real.buffers
    geometry = header["geometry"]
    obj = Drillhole(
        collar=buffers[geometry["collar"]],
        survey=buffers[geometry["survey"]],
        intervals=buffers[geometry["intervals"]],
        end_of_hole=header["params"]["end_of_hole"][0],
    )
    data = encode_object(
        obj,
        id=header["id"],
        name=header["name"],
        parent=header["parent"],
        revision=header["revision"],
        origin=(0.0, 0.0, 0.0),
    )
    made = decode_frame(data)
    real_header, real_buffers = _by_role(header, buffers)
    made_header, made_buffers = _by_role(made.header, made.buffers)
    assert made_header == real_header
    assert made_buffers.keys() == real_buffers.keys()
    for role, array in real_buffers.items():
        assert np.array_equal(made_buffers[role], array), role


def _tool_header(name: str) -> dict[str, Any]:
    """Return the header of the tool fixture ``name``, or skip when the folder is missing."""
    folder = _folder().parent / "tool"
    if not folder.is_dir():
        pytest.skip(f"the folder {folder} does not exist")
    path = folder / name
    if not path.exists():
        pytest.fail(f"the tool fixture {name} is not in {folder}")
    return json.loads(path.read_text(encoding="utf-8"))


def _check_header(made: bytes, fixture: dict[str, Any]) -> None:
    """Assert that the header of the decoded frame equals the fixture header."""
    assert decode_frame(made).header == fixture


def test_selection_set_matches_fixture() -> None:
    """The encoder builds the selection.set frame of the geoviewer repository."""
    fixture = _tool_header("selection.set.json")
    try:
        counts = {entry["kind"]: 0 for entry in fixture["entries"]}
        shapes = {entry["name"]: entry["shape"] for entry in fixture["buffers"]}
        for entry in fixture["entries"]:
            counts[entry["kind"]] = shapes[entry["indices"]][0]
        unknown = set(counts) - {"cells", "nodes"}
        request, object_id = fixture["request"], fixture["object"]
        revision = fixture.get("revision")
    except (KeyError, IndexError, TypeError) as error:
        pytest.fail(f"selection.set.json has an unexpected form: {error!r}")
    assert not unknown, f"selection.set.json has unknown entry kinds: {unknown}"
    cells = np.arange(counts["cells"]) if "cells" in counts else None
    nodes = np.arange(counts["nodes"]) if "nodes" in counts else None
    _check_header(selection_set_frame(request, object_id, cells, nodes, revision), fixture)


def test_selection_get_matches_fixture() -> None:
    """The encoder builds the selection.get frame of the geoviewer repository."""
    fixture = _tool_header("selection.get.json")
    assert "request" in fixture, "selection.get.json has no request"
    _check_header(selection_get_frame(fixture["request"]), fixture)


def test_object_fetch_matches_fixture() -> None:
    """The encoder builds the object.fetch frame of the geoviewer repository."""
    fixture = _tool_header("object.fetch.json")
    assert "request" in fixture and "id" in fixture, "object.fetch.json has no request or id"
    _check_header(fetch_frame(fixture["request"], fixture["id"]), fixture)


def test_command_cancel_matches_fixture() -> None:
    """The encoder builds the command.cancel frame of the geoviewer repository."""
    fixture = _tool_header("command.cancel.json")
    assert "request" in fixture, "command.cancel.json has no request"
    _check_header(cancel_frame(fixture["request"]), fixture)


def test_commands_list_matches_fixture() -> None:
    """The encoder builds the commands.list frame of the geoviewer repository."""
    fixture = _tool_header("commands.list.json")
    assert "request" in fixture, "commands.list.json has no request"
    _check_header(commands_list_frame(fixture["request"]), fixture)


def test_command_run_matches_fixture() -> None:
    """The encoder builds the command.run frame of the geoviewer repository."""
    fixture = _tool_header("command.run.json")
    assert (
        "request" in fixture and "command" in fixture and "input" in fixture
    ), "command.run.json has no request, command or input"
    _check_header(command_frame(fixture["request"], fixture["command"], fixture["input"]), fixture)
