"""Check the frames of the geoviewer repository against loop_link (SC-004)."""

from __future__ import annotations

import copy
from typing import Any

import numpy as np
import pytest
from conftest import geoviewer_fixtures

from loop_link.encode import encode_object, hello_frame
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
