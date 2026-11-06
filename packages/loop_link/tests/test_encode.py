"""Tests for the encoding of objects into frames."""

from __future__ import annotations

import numpy as np
import pytest

from loop_link import encode
from loop_link.errors import InvalidInputError, ObjectTooLargeError
from loop_link.objects import (
    DEFAULT_CATEGORY_COLOURS,
    BlockModel,
    Categorical,
    Drillhole,
    Lines,
    Points,
    Surface,
    TetraMesh,
    Vectors,
)
from loop_link.protocol import decode_frame

ORIGIN = (0.0, 0.0, 0.0)
GIB2 = 2 * 1024**3


def _encode(obj, origin=ORIGIN, **kwargs):
    """Encode an object with default arguments and decode the frame."""
    args = {"id": "obj-1", "name": "Object 1", "parent": None, "revision": 1}
    args.update(kwargs)
    return decode_frame(encode.encode_object(obj, origin=origin, **args))


def _buffer_dtype(frame, name):
    """Return the dtype name of a buffer from the header."""
    for entry in frame.header["buffers"]:
        if entry["name"] == name:
            return entry["dtype"]
    raise AssertionError(f"no buffer {name}")


def _geometry(frame, role):
    """Return the buffer that the header geometry gives for a role."""
    name = frame.header["geometry"][role]
    assert name in frame.buffers
    return frame.buffers[name]


def _data_entry(frame, name):
    """Return the data description for a property."""
    entries = [e for e in frame.header["data"] if e["name"] == name]
    assert len(entries) == 1
    return entries[0]


def _unchecked_points(vertices):
    """Make a Points object without the input check, to test the encoder checks."""
    obj = object.__new__(Points)
    object.__setattr__(obj, "vertices", vertices)
    object.__setattr__(obj, "properties", {})
    return obj


def test_header_fields_of_upsert():
    frame = _encode(Points(np.zeros((2, 3))), id="p1", name="P", parent="g", revision=4)
    header = frame.header
    assert header["type"] == "scene.upsert"
    assert header["id"] == "p1"
    assert header["object_type"] == "Points"
    assert header["name"] == "P"
    assert header["parent"] == "g"
    assert header["revision"] == 4
    assert header["edit"] == []
    assert header["replaces"] is None
    assert header["params"] == {}
    assert header["data"] == []
    assert header["property_groups"] == []
    assert header["drape"] is None
    assert header["style"] == {}


def test_points_geometry():
    verts = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    frame = _encode(Points(verts))
    got = _geometry(frame, "vertices")
    assert got.shape == (2, 3)
    assert np.allclose(got, verts)


def test_lines_map_to_curve():
    verts = np.arange(9, dtype=float).reshape(3, 3)
    segs = np.array([[0, 1], [1, 2]])
    frame = _encode(Lines(verts, segs))
    assert frame.header["object_type"] == "Curve"
    got = _geometry(frame, "segments")
    assert got.dtype == np.dtype("<u4")
    assert got.tolist() == segs.tolist()
    assert np.allclose(_geometry(frame, "vertices"), verts)


def test_surface_triangles_u32():
    verts = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]])
    tris = np.array([[0, 1, 2]])
    frame = _encode(Surface(verts, tris))
    assert frame.header["object_type"] == "Surface"
    got = _geometry(frame, "triangles")
    assert got.dtype == np.dtype("<u4")
    assert got.tolist() == [[0, 1, 2]]


def test_tetramesh_tetrahedra_u32():
    verts = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    frame = _encode(TetraMesh(verts, np.array([[0, 1, 2, 3]])))
    assert frame.header["object_type"] == "TetraMesh"
    got = _geometry(frame, "tetrahedra")
    assert got.dtype == np.dtype("<u4")
    assert got.tolist() == [[0, 1, 2, 3]]


def test_block_model_params_and_origin_shift():
    block = BlockModel([10.0, 20.0, 30.0], [1.0, 2.0, 3.0], [2, 2, 2])
    frame = _encode(block, origin=(4.0, 5.0, 6.0))
    header = frame.header
    assert header["object_type"] == "BlockModel"
    assert header["geometry"] == {}
    assert header["params"]["origin"] == [6.0, 15.0, 24.0]
    assert header["params"]["cell_sizes"] == [1.0, 2.0, 3.0]
    assert header["params"]["cells"] == [2, 2, 2]


def test_drillhole_shift_and_end_of_hole():
    collar = np.array([1250.0, 3400.0, 312.5])
    survey = np.array([[0.0, 45.0, 60.0], [120.0, 50.0, 65.0]])
    intervals = np.array([[10.0, 240.0]])
    frame = _encode(
        Drillhole(collar, survey, intervals, end_of_hole=250.0),
        origin=(1000.0, 3000.0, 300.0),
    )
    assert frame.header["object_type"] == "Drillhole"
    got = _geometry(frame, "collar")
    assert got.dtype == np.dtype("<f8")
    assert np.allclose(got, [250.0, 400.0, 12.5])
    surv = _geometry(frame, "survey")
    assert surv.dtype == np.dtype("<f8")
    assert surv.shape == (2, 3)
    assert np.allclose(surv, survey)
    inter = _geometry(frame, "intervals")
    assert inter.dtype == np.dtype("<f8")
    assert inter.shape == (1, 2)
    assert np.allclose(inter, intervals)
    assert frame.header["params"] == {"end_of_hole": [250.0]}


def test_drillhole_without_end_of_hole():
    frame = _encode(Drillhole([0, 0, 0], [[0, 0, 0]], [[0, 1]]))
    assert "end_of_hole" not in frame.header["params"]


def test_drillhole_interval_property_association():
    hole = Drillhole(
        [0, 0, 0],
        [[0, 0, 0], [10, 0, 0]],
        [[0, 1], [1, 2]],
        properties={"az": np.array([1.0, 2.0])},
        interval_properties={"lith": np.array([3.0, 4.0])},
    )
    frame = _encode(hole)
    assert _data_entry(frame, "az")["association"] == "vertex"
    assert _data_entry(frame, "lith")["association"] == "interval"
    assert frame.buffers["data/lith"].tolist() == [3.0, 4.0]


def test_vertices_are_f32_after_shift_in_f64():
    # The shift must happen in f64. A shift in f32 would lose the small digits.
    origin = (500000.0, 7000000.0, 0.0)
    verts = np.array([[500000.123, 7000000.456, 12.789], [500100.5, 7000200.25, -30.0]])
    frame = _encode(Points(verts), origin=origin)
    name = frame.header["geometry"]["vertices"]
    assert _buffer_dtype(frame, name) == "f32"
    got = frame.buffers[name]
    assert got.dtype == np.dtype("<f4")
    assert np.allclose(got, verts - np.array(origin), atol=1e-3, rtol=0)


def test_vertices_are_f64_when_far_from_origin():
    verts = np.array([[5000.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    frame = _encode(Points(verts))
    name = frame.header["geometry"]["vertices"]
    assert _buffer_dtype(frame, name) == "f64"
    assert np.allclose(frame.buffers[name], verts)


def test_vertices_at_limit_are_f32():
    verts = np.array([[4000.0, -4000.0, 0.0]])
    frame = _encode(Points(verts))
    name = frame.header["geometry"]["vertices"]
    assert _buffer_dtype(frame, name) == "f32"


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_non_finite_coordinates_are_rejected(bad):
    obj = _unchecked_points(np.array([[0.0, bad, 0.0]]))
    with pytest.raises(InvalidInputError):
        _encode(obj)


def test_float32_property():
    frame = _encode(Points(np.zeros((2, 3)), {"a": np.array([1.5, 2.5], dtype=np.float32)}))
    entry = _data_entry(frame, "a")
    assert entry["id"] == "a"
    assert entry["buffer"] == "data/a"
    assert entry["association"] == "vertex"
    assert entry["data_type"] == "float"
    assert _buffer_dtype(frame, "data/a") == "f32"
    assert frame.buffers["data/a"].tolist() == [1.5, 2.5]


def test_float64_property():
    frame = _encode(Points(np.zeros((2, 3)), {"a": np.array([1.5, 2.5])}))
    assert _data_entry(frame, "a")["data_type"] == "float"
    assert _buffer_dtype(frame, "data/a") == "f64"


def test_bool_property():
    frame = _encode(Points(np.zeros((2, 3)), {"a": np.array([True, False])}))
    assert _data_entry(frame, "a")["data_type"] == "boolean"
    assert _buffer_dtype(frame, "data/a") == "u8"
    assert frame.buffers["data/a"].tolist() == [1, 0]


def test_int_property_in_i32():
    frame = _encode(Points(np.zeros((2, 3)), {"a": np.array([-5, 2**31 - 1], dtype=np.int64)}))
    assert _data_entry(frame, "a")["data_type"] == "integer"
    assert _buffer_dtype(frame, "data/a") == "i32"
    assert frame.buffers["data/a"].tolist() == [-5, 2**31 - 1]


@pytest.mark.parametrize("value", [2**31, -(2**31) - 1])
def test_int_property_outside_i32_is_rejected(value):
    obj = Points(np.zeros((2, 3)), {"a": np.array([0, value], dtype=np.int64)})
    with pytest.raises(InvalidInputError):
        _encode(obj)


def test_text_property():
    frame = _encode(Points(np.zeros((3, 3)), {"t": np.array(["ab", "", "cde"])}))
    entry = _data_entry(frame, "t")
    assert entry["data_type"] == "text"
    assert entry["buffer"] == "data/t"
    assert entry["offsets"] == "data/t/offsets"
    assert _buffer_dtype(frame, "data/t") == "u8"
    assert _buffer_dtype(frame, "data/t/offsets") == "u32"
    assert bytes(frame.buffers["data/t"]) == b"abcde"
    assert frame.buffers["data/t/offsets"].tolist() == [0, 2, 2, 5]


def test_text_property_uses_utf8_bytes():
    frame = _encode(Points(np.zeros((2, 3)), {"t": np.array(["é", "a"])}))
    assert bytes(frame.buffers["data/t"]) == "éa".encode()
    assert frame.buffers["data/t/offsets"].tolist() == [0, 2, 3]


def test_categorical_property_with_default_colours():
    cat = Categorical([0, 2, 1], ["a", "b", "c"])
    frame = _encode(Points(np.zeros((3, 3)), {"c": cat}))
    entry = _data_entry(frame, "c")
    assert entry["data_type"] == "referenced"
    assert _buffer_dtype(frame, "data/c") == "u32"
    assert frame.buffers["data/c"].tolist() == [0, 2, 1]
    expected = [
        {"id": i, "name": n, "colour": list(DEFAULT_CATEGORY_COLOURS[i % 10])}
        for i, n in enumerate(["a", "b", "c"])
    ]
    assert entry["value_map"] == expected


def test_categorical_property_with_given_colours():
    colours = [(1, 2, 3, 4), (5, 6, 7, 8)]
    cat = Categorical([1, 0], ["a", "b"], colours)
    frame = _encode(Points(np.zeros((2, 3)), {"c": cat}))
    assert _data_entry(frame, "c")["value_map"] == [
        {"id": 0, "name": "a", "colour": [1, 2, 3, 4]},
        {"id": 1, "name": "b", "colour": [5, 6, 7, 8]},
    ]


def test_default_colours_wrap_after_ten():
    names = [f"n{i}" for i in range(12)]
    frame = _encode(Points(np.zeros((1, 3)), {"c": Categorical([11], names)}))
    value_map = _data_entry(frame, "c")["value_map"]
    assert value_map[11]["id"] == 11
    assert value_map[11]["colour"] == list(DEFAULT_CATEGORY_COLOURS[1])


def test_object_dtype_property_is_rejected():
    obj = Points(np.zeros((2, 3)), {"a": np.array([None, "x"], dtype=object)})
    with pytest.raises(InvalidInputError):
        _encode(obj)


def test_duplicate_buffer_name_is_rejected():
    # A vertex property and a cell property with the same name give one buffer name.
    obj = Lines(
        np.zeros((3, 3)),
        [[0, 1], [1, 2]],
        properties={"a": np.zeros(3)},
        cell_properties={"a": np.zeros(2)},
    )
    with pytest.raises(InvalidInputError):
        _encode(obj)


def test_cell_property_association():
    obj = Surface(
        np.zeros((3, 3)),
        [[0, 1, 2]],
        properties={"v": np.zeros(3)},
        cell_properties={"c": np.zeros(1)},
    )
    frame = _encode(obj)
    assert _data_entry(frame, "v")["association"] == "vertex"
    assert _data_entry(frame, "c")["association"] == "cell"


def test_vectors_orientation_property_group():
    verts = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    dd = np.array([90.0, 180.0])
    dip = np.array([30.0, 45.0])
    frame = _encode(Vectors(verts, orientation=(dd, dip)))
    assert frame.header["object_type"] == "Points"
    assert frame.header["property_groups"] == [
        {
            "id": "orientation",
            "name": "orientation",
            "role": "orientation",
            "members": [
                {"role": "dip_direction", "data": "orientation/dip_direction"},
                {"role": "dip", "data": "orientation/dip"},
            ],
        }
    ]
    for member in frame.header["property_groups"][0]["members"]:
        assert _data_entry(frame, member["data"])["id"] == member["data"]
    for name, values in (("orientation/dip_direction", dd), ("orientation/dip", dip)):
        entry = _data_entry(frame, name)
        assert entry["association"] == "vertex"
        assert entry["data_type"] == "float"
        assert _buffer_dtype(frame, entry["buffer"]) == "f64"
        assert frame.buffers[entry["buffer"]].tolist() == values.tolist()


def test_vectors_lineation_property_group():
    trend = np.array([10.0])
    plunge = np.array([20.0])
    frame = _encode(Vectors([[0.0, 0, 0]], lineation=(trend, plunge)))
    groups = frame.header["property_groups"]
    assert len(groups) == 1
    assert groups[0]["role"] == "lineation"
    assert [m["role"] for m in groups[0]["members"]] == ["trend", "plunge"]
    for member in groups[0]["members"]:
        assert _data_entry(frame, member["data"])["id"] == member["data"]
    members = {m["role"]: m["data"] for m in groups[0]["members"]}
    assert frame.buffers[_data_entry(frame, members["trend"])["buffer"]].tolist() == [10.0]
    assert frame.buffers[_data_entry(frame, members["plunge"])["buffer"]].tolist() == [20.0]


def test_object_too_large():
    # A broadcast view declares a large size without a real allocation.
    # The encoder must find the size before it makes any copy of the vertices.
    big = np.broadcast_to(np.zeros(3), (2**28, 3))
    obj = _unchecked_points(big)
    with pytest.raises(ObjectTooLargeError) as info:
        encode.encode_object(obj, id="big", name="Big", parent=None, revision=1, origin=ORIGIN)
    assert info.value.object_id == "big"
    assert info.value.limit == GIB2
    assert info.value.size > GIB2


def test_hello_frame():
    frame = decode_frame(encode.hello_frame("EPSG:28350", (1.0, 2, 3), "0.1.0"))
    assert frame.header["type"] == "hello"
    assert frame.header["application"] == {"name": "loop_link", "version": "0.1.0"}
    assert frame.header["crs"] == "EPSG:28350"
    assert frame.header["origin"] == [1.0, 2.0, 3.0]
    assert frame.buffers == {}


def test_heartbeat_frame():
    header = decode_frame(encode.heartbeat_frame(7)).header
    assert header["type"] == "heartbeat"
    assert header["seq"] == 7


def test_remove_frame():
    header = decode_frame(encode.remove_frame("x1")).header
    assert header["type"] == "scene.remove"
    assert header["id"] == "x1"


def test_clear_frame():
    assert decode_frame(encode.clear_frame()).header["type"] == "scene.clear"


def test_group_frame():
    frame = decode_frame(encode.group_frame("g1", "Group", "root", 3))
    header = frame.header
    assert header["type"] == "scene.upsert"
    assert header["object_type"] == "Group"
    assert header["id"] == "g1"
    assert header["name"] == "Group"
    assert header["parent"] == "root"
    assert header["revision"] == 3
    assert frame.buffers == {}


def test_vertices_just_above_limit_are_f64():
    frame = _encode(Points(np.array([[4000.01, 0.0, 0.0]])))
    name = frame.header["geometry"]["vertices"]
    assert _buffer_dtype(frame, name) == "f64"


def test_f64_when_shifted_values_are_large_and_raw_values_small():
    # The rule uses the values after the shift, not the raw values.
    verts = np.array([[1.0, 2.0, 3.0], [-10.0, 0.0, 0.0]])
    frame = _encode(Points(verts), origin=(5000.0, 0.0, 0.0))
    name = frame.header["geometry"]["vertices"]
    assert _buffer_dtype(frame, name) == "f64"
    assert np.allclose(frame.buffers[name], verts - np.array([5000.0, 0.0, 0.0]))


def test_f32_when_raw_values_are_large_and_shifted_values_small():
    verts = np.array([[5001.0, 2.0, 3.0], [5010.0, 0.0, 0.0]])
    frame = _encode(Points(verts), origin=(5000.0, 0.0, 0.0))
    name = frame.header["geometry"]["vertices"]
    assert _buffer_dtype(frame, name) == "f32"
    assert np.allclose(frame.buffers[name], verts - np.array([5000.0, 0.0, 0.0]))


def test_frame_over_frame_limit_is_object_too_large(monkeypatch):
    # The frame limit is small here, so a small object is too large.
    monkeypatch.setattr("loop_link.protocol.MAX_FRAME_BYTES", 1000)
    obj = Points(np.zeros((1000, 3)))
    with pytest.raises(ObjectTooLargeError) as info:
        encode.encode_object(obj, id="frame-big", name="B", parent=None, revision=1, origin=ORIGIN)
    assert info.value.object_id == "frame-big"


def test_block_model_properties_and_cell_properties():
    node_values = np.arange(12, dtype=np.float64)  # (1+1)*(2+1)*(1+1) nodes
    cell_values = np.array([1, 2], dtype=np.int64)  # 1*2*1 cells
    block = BlockModel(
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0],
        [1, 2, 1],
        properties={"n": node_values},
        cell_properties={"c": cell_values},
    )
    frame = _encode(block)
    assert _data_entry(frame, "n")["association"] == "vertex"
    assert _data_entry(frame, "c")["association"] == "cell"
    assert frame.buffers["data/n"].tolist() == node_values.tolist()
    assert frame.buffers["data/c"].tolist() == [1, 2]


def test_float32_nan_property_is_accepted():
    values = np.array([1.0, np.nan], dtype=np.float32)
    frame = _encode(Points(np.zeros((2, 3)), {"a": values}))
    got = frame.buffers["data/a"]
    assert got[0] == 1.0
    assert np.isnan(got[1])


def test_user_property_collides_with_orientation_buffer():
    obj = Vectors(
        [[0.0, 0.0, 0.0]],
        orientation=([90.0], [30.0]),
        properties={"orientation/dip": np.array([1.0])},
    )
    with pytest.raises(InvalidInputError):
        _encode(obj)
