"""Tests for the validation rules of the fake viewer."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from loop_link.encode import encode_object
from loop_link.objects import Drillhole, Vectors
from loop_link.protocol import Frame, encode_frame
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


def drillhole_frame(association: str, log_sets: list[dict[str, Any]] | None = None) -> bytes:
    """Return a hand-made Drillhole frame with a data array."""
    header = {
        "type": "scene.upsert",
        "id": "h",
        "object_type": "Drillhole",
        "revision": 1,
        "geometry": {"collar": "c", "survey": "s", "intervals": "i"},
        "data": [entry("x", association)],
        "log_sets": log_sets or [],
    }
    buffers = [
        ("c", np.zeros(3, dtype=np.float32)),
        ("s", np.zeros((3, 3), dtype=np.float32)),
        ("i", np.zeros((2, 2), dtype=np.float32)),
        ("x", np.zeros(2, dtype=np.float32)),
    ]
    return encode_frame(header, buffers)


def test_drillhole_with_interval_properties_is_accepted(fake: Any) -> None:
    """The encoder sends interval properties as cell arrays. The viewer accepts them."""
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


def test_interval_array_without_log_sets_is_refused(fake: Any) -> None:
    """An interval association needs a log set of kind interval."""
    assert codes(fake.accept_frame(drillhole_frame("interval"))) == ["invalid_object"]
    with_set = drillhole_frame("interval", [{"name": "main", "kind": "interval"}])
    assert fake.accept_frame(with_set) == []


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
