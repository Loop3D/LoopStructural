"""Tests for the decoding of viewer objects into scene objects."""

from __future__ import annotations

import contextlib
import copy
import threading
import time
import warnings
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from loop_link import encode
from loop_link import fetch as fetch_module
from loop_link.errors import (
    BusyError,
    CancelledError,
    InvalidInputError,
    ObjectTooLargeError,
    UnsupportedObjectError,
    ViewerClosedError,
)
from loop_link.fetch import decode_object
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
from loop_link.protocol import Frame, ProtocolError, decode_frame
from loop_link.testing import FakeViewer
from loop_link.tool import ToolConnection
from loop_link.viewer import Viewer

ORIGIN = (1000.0, 2000.0, 3000.0)
WAIT = 3.0
VERTS = np.array([[1001.5, 2002.5, 3003.5], [1010.0, 2020.0, 3030.0], [999.0, 1999.0, 2999.0]])


def _parts(obj, origin=ORIGIN):
    """Encode an object as the viewer holds it. Return a copy of the header and the buffers."""
    frame = decode_frame(
        encode.encode_object(
            obj, id="obj-1", name="Object 1", parent=None, revision=1, origin=origin
        )
    )
    return copy.deepcopy(frame.header), dict(frame.buffers)


def _decode(obj, origin=ORIGIN, decode_origin=ORIGIN):
    """Round-trip an object through the encoder and `decode_object`."""
    header, arrays = _parts(obj, origin)
    return decode_object(header, arrays, decode_origin)


def _datum(header, name):
    """Return the data description with the given name."""
    for desc in header["data"]:
        if desc["name"] == name:
            return desc
    raise AssertionError(f"no data {name}")


def _drillhole(**kwargs):
    """Return a small Drillhole."""
    args = {
        "collar": [1005.0, 2005.0, 3005.0],
        "survey": [[1005.0, 2005.0, 3005.0], [1005.0, 2005.0, 2995.0], [1005.0, 2005.0, 2985.0]],
        "intervals": [[0.0, 5.0], [5.0, 20.0]],
    }
    args.update(kwargs)
    return Drillhole(**args)


# --- geometry ---


def test_points_round_trip():
    result = _decode(Points(VERTS, {"a": np.array([1.0, 2.0, 3.0])}))
    assert isinstance(result, Points)
    assert result.vertices.dtype == np.float64
    np.testing.assert_allclose(result.vertices, VERTS)
    np.testing.assert_allclose(result.properties["a"], [1.0, 2.0, 3.0])


def test_vertices_are_shifted_by_the_origin():
    header, arrays = _parts(Points(VERTS))
    stored = arrays["geometry/vertices"].astype(np.float64)
    result = decode_object(header, arrays, (10.0, 20.0, 30.0))
    np.testing.assert_allclose(result.vertices, stored + np.array([10.0, 20.0, 30.0]))


def test_origin_none_is_zero():
    header, arrays = _parts(Points(VERTS), origin=(0.0, 0.0, 0.0))
    result = decode_object(header, arrays, None)
    np.testing.assert_allclose(result.vertices, VERTS)
    assert result.vertices.dtype == np.float64


def test_f64_vertices_keep_precision():
    far = np.array([[1.0e9 + 0.001, 0.0, 0.0], [1.0e9, 1.0, 1.0]])
    header, arrays = _parts(Points(far), origin=(0.0, 0.0, 0.0))
    assert arrays["geometry/vertices"].dtype == np.float64
    result = decode_object(header, arrays, None)
    np.testing.assert_array_equal(result.vertices, far)


def test_curve_becomes_lines():
    result = _decode(Lines(VERTS, [[0, 1], [1, 2]]))
    assert isinstance(result, Lines)
    np.testing.assert_array_equal(result.segments, [[0, 1], [1, 2]])
    np.testing.assert_allclose(result.vertices, VERTS)


def test_surface():
    result = _decode(Surface(VERTS, [[0, 1, 2]]))
    assert isinstance(result, Surface)
    np.testing.assert_array_equal(result.triangles, [[0, 1, 2]])


def test_tetra_mesh():
    verts = np.vstack([VERTS, [[1003.0, 2003.0, 3010.0]]])
    result = _decode(TetraMesh(verts, [[0, 1, 2, 3]], cell_properties={"c": [4.0]}))
    assert isinstance(result, TetraMesh)
    np.testing.assert_array_equal(result.tetrahedra, [[0, 1, 2, 3]])
    np.testing.assert_allclose(result.cell_properties["c"], [4.0])


def test_block_model_origin_is_params_origin_plus_origin():
    block = BlockModel(
        [1010.0, 2020.0, 3030.0], [1.0, 2.0, 3.0], [1, 1, 1], cell_properties={"c": [5.0]}
    )
    header, arrays = _parts(block)
    assert header["params"]["origin"] == [10.0, 20.0, 30.0]
    result = decode_object(header, arrays, ORIGIN)
    assert isinstance(result, BlockModel)
    np.testing.assert_allclose(result.origin, [1010.0, 2020.0, 3030.0])
    np.testing.assert_allclose(result.cell_sizes, [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(result.cells, [1, 1, 1])
    np.testing.assert_allclose(result.cell_properties["c"], [5.0])
    assert result.properties == {}


def test_block_model_with_origin_none():
    block = BlockModel([1.0, 2.0, 3.0], [1.0, 1.0, 1.0], [1, 1, 1])
    header, arrays = _parts(block, origin=(0.0, 0.0, 0.0))
    result = decode_object(header, arrays, None)
    np.testing.assert_allclose(result.origin, [1.0, 2.0, 3.0])


def test_drillhole_collar_and_end_of_hole():
    hole = _drillhole(end_of_hole=25.0)
    header, arrays = _parts(hole)
    assert header["params"]["end_of_hole"] == [25.0]
    result = decode_object(header, arrays, ORIGIN)
    assert isinstance(result, Drillhole)
    assert result.collar.dtype == np.float64
    np.testing.assert_allclose(result.collar, [1005.0, 2005.0, 3005.0])
    np.testing.assert_allclose(result.survey, hole.survey)
    np.testing.assert_allclose(result.intervals, hole.intervals)
    assert result.end_of_hole == 25.0


def test_drillhole_without_end_of_hole():
    result = _decode(_drillhole())
    assert result.end_of_hole is None


# --- Vectors and groups ---


def test_points_with_orientation_group_become_vectors():
    obj = Vectors(VERTS, orientation=([10.0, 20.0, 30.0], [40.0, 50.0, 60.0]))
    result = _decode(obj)
    assert isinstance(result, Vectors)
    np.testing.assert_allclose(result.orientation[0], [10.0, 20.0, 30.0])
    np.testing.assert_allclose(result.orientation[1], [40.0, 50.0, 60.0])
    assert result.orientation[0].dtype == np.float64
    assert result.lineation is None
    # The group members are not properties.
    assert result.properties == {}


def test_points_with_lineation_group_become_vectors():
    obj = Vectors(VERTS, lineation=([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]))
    result = _decode(obj)
    assert isinstance(result, Vectors)
    np.testing.assert_allclose(result.lineation[0], [1.0, 2.0, 3.0])
    np.testing.assert_allclose(result.lineation[1], [4.0, 5.0, 6.0])
    assert result.orientation is None
    assert result.properties == {}


def test_vectors_with_both_groups_and_a_property():
    obj = Vectors(
        VERTS,
        orientation=([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]),
        lineation=([7.0, 8.0, 9.0], [10.0, 11.0, 12.0]),
        properties={"p": np.array([0.5, 1.5, 2.5])},
    )
    result = _decode(obj)
    assert isinstance(result, Vectors)
    assert result.orientation is not None and result.lineation is not None
    assert list(result.properties) == ["p"]


def test_points_with_other_group_stay_points():
    header, arrays = _parts(Vectors(VERTS, orientation=([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])))
    header["property_groups"][0]["role"] = "other"
    result = decode_object(header, arrays, ORIGIN)
    assert type(result) is Points


# --- unsupported types ---


@pytest.mark.parametrize("object_type", ["Group", "Grid2D", "GeoImage"])
def test_unsupported_types(object_type):
    header, arrays = _parts(Points(VERTS))
    header["object_type"] = object_type
    with pytest.raises(UnsupportedObjectError) as info:
        decode_object(header, arrays, ORIGIN)
    assert "obj-1" in str(info.value)
    assert object_type in str(info.value)


# --- property associations ---


def test_cell_property_on_points_is_unsupported():
    header, arrays = _parts(Points(VERTS, {"a": np.array([1.0, 2.0, 3.0])}))
    _datum(header, "a")["association"] = "cell"
    with pytest.raises(UnsupportedObjectError):
        decode_object(header, arrays, ORIGIN)


def test_cell_property_on_vectors_is_unsupported():
    header, arrays = _parts(
        Vectors(VERTS, orientation=([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]), properties={"a": [1, 2, 3]})
    )
    _datum(header, "a")["association"] = "cell"
    with pytest.raises(UnsupportedObjectError):
        decode_object(header, arrays, ORIGIN)


def test_surface_cell_property_goes_to_cell_properties():
    obj = Surface(VERTS, [[0, 1, 2]], cell_properties={"c": [9.0]}, properties={"v": [1, 2, 3]})
    result = _decode(obj)
    assert list(result.properties) == ["v"]
    assert list(result.cell_properties) == ["c"]


def test_drillhole_properties():
    hole = _drillhole(
        properties={"s": [1.0, 2.0, 3.0]},
        interval_properties={"i": [10.0, 20.0]},
    )
    result = _decode(hole)
    np.testing.assert_allclose(result.properties["s"], [1.0, 2.0, 3.0])
    np.testing.assert_allclose(result.interval_properties["i"], [10.0, 20.0])


def test_drillhole_cell_property_goes_to_interval_properties():
    header, arrays = _parts(_drillhole(interval_properties={"i": [10.0, 20.0]}))
    desc = _datum(header, "i")
    desc["association"] = "cell"
    desc.pop("log_set", None)
    result = decode_object(header, arrays, ORIGIN)
    assert list(result.interval_properties) == ["i"]


def test_drillhole_interval_without_log_set_is_kept():
    header, arrays = _parts(_drillhole(interval_properties={"i": [10.0, 20.0]}))
    _datum(header, "i").pop("log_set")
    result = decode_object(header, arrays, ORIGIN)
    assert list(result.interval_properties) == ["i"]


def _two_log_sets():
    """Return a Drillhole header with a second log set named 'logs'."""
    header, arrays = _parts(_drillhole(interval_properties={"i": [10.0, 20.0]}))
    other = copy.deepcopy(_datum(header, "i"))
    other.update({"id": "j", "name": "j", "buffer": "data/j", "log_set": "logs"})
    header["data"].append(other)
    header["log_sets"].append({"name": "logs", "kind": "interval"})
    arrays["data/j"] = np.array([7.0, 8.0])
    return header, arrays


def test_drillhole_keeps_the_first_interval_set_and_warns_for_the_other():
    header, arrays = _two_log_sets()
    with pytest.warns(UserWarning, match="dropped the property 'j' of object 'obj-1'") as record:
        result = decode_object(header, arrays, ORIGIN)
    assert list(result.interval_properties) == ["i"]
    assert len([w for w in record if "'j'" in str(w.message)]) == 1


@pytest.mark.parametrize("association", ["sample", "object"])
def test_other_associations_are_dropped_with_one_warning(association):
    header, arrays = _parts(Points(VERTS, {"a": np.array([1.0, 2.0, 3.0]), "b": [1, 2, 3]}))
    _datum(header, "a")["association"] = association
    with pytest.warns(UserWarning, match="dropped the property 'a' of object 'obj-1'") as record:
        result = decode_object(header, arrays, ORIGIN)
    assert len(record) == 1
    assert list(result.properties) == ["b"]


def test_no_warning_without_dropped_properties():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _decode(Points(VERTS, {"a": np.array([1.0, 2.0, 3.0])}))


# --- property types ---


def test_float_keeps_f32_and_f64():
    obj = Points(
        VERTS,
        {
            "f32": np.array([1.0, 2.0, 3.0], dtype=np.float32),
            "f64": np.array([1.0, 2.0, 3.0], dtype=np.float64),
        },
    )
    result = _decode(obj)
    assert result.properties["f32"].dtype == np.float32
    assert result.properties["f64"].dtype == np.float64


def test_integer_is_int32():
    result = _decode(Points(VERTS, {"n": np.array([1, -2, 3], dtype=np.int64)}))
    assert result.properties["n"].dtype == np.int32
    np.testing.assert_array_equal(result.properties["n"], [1, -2, 3])


def test_boolean_is_numpy_bool():
    result = _decode(Points(VERTS, {"b": np.array([True, False, True])}))
    assert result.properties["b"].dtype == np.bool_
    np.testing.assert_array_equal(result.properties["b"], [True, False, True])


def test_text_is_unicode_array():
    values = np.array(["alpha", "", "béta"])
    result = _decode(Points(VERTS, {"t": values}))
    assert result.properties["t"].dtype.kind == "U"
    assert result.properties["t"].tolist() == ["alpha", "", "béta"]


def test_referenced_becomes_categorical():
    cat = Categorical([0, 2, 1], ["a", "b", "c"], [(1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12)])
    result = _decode(Points(VERTS, {"c": cat}))
    value = result.properties["c"]
    assert isinstance(value, Categorical)
    assert list(value.names) == ["a", "b", "c"]
    np.testing.assert_array_equal(value.codes, [0, 2, 1])
    assert value.colours == ((1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12))


def _sparse_referenced():
    """Return a header and buffers with the value map ids 5, 2 and 9."""
    header, arrays = _parts(Points(VERTS, {"c": Categorical([0, 1, 2], ["x", "y", "z"])}))
    _datum(header, "c")["value_map"] = [
        {"id": 5, "name": "five", "colour": [10, 20, 30]},
        {"id": 2, "name": "two", "colour": [40, 50, 60, 70]},
        {"id": 9, "name": "nine", "colour": [80, 90, 100]},
    ]
    arrays["data/c"] = np.array([9, 5, 2], dtype=np.uint32)
    return header, arrays


def test_referenced_names_are_ordered_by_id_and_codes_renumbered():
    header, arrays = _sparse_referenced()
    value = decode_object(header, arrays, ORIGIN).properties["c"]
    assert list(value.names) == ["two", "five", "nine"]
    np.testing.assert_array_equal(value.codes, [2, 1, 0])


def test_referenced_three_value_colours_get_alpha_255():
    header, arrays = _sparse_referenced()
    value = decode_object(header, arrays, ORIGIN).properties["c"]
    assert value.colours == ((40, 50, 60, 70), (10, 20, 30, 255), (80, 90, 100, 255))


def test_referenced_code_not_in_value_map_raises():
    header, arrays = _sparse_referenced()
    arrays["data/c"] = np.array([9, 5, 7], dtype=np.uint32)
    with pytest.raises(InvalidInputError):
        decode_object(header, arrays, ORIGIN)


def test_referenced_on_drillhole_interval():
    cat = Categorical([1, 0], ["rock", "soil"])
    result = _decode(_drillhole(interval_properties={"lith": cat}))
    value = result.interval_properties["lith"]
    assert list(value.names) == ["rock", "soil"]
    np.testing.assert_array_equal(value.codes, [1, 0])


# --- the assembler ---

_PROTOCOL_DTYPES = {"float32": "f32", "float64": "f64", "uint32": "u32", "int32": "i32"}


def _points_object() -> dict[str, Any]:
    """Return the object header of a small Points object that has one array, ``v``."""
    return {
        "type": "scene.upsert",
        "id": "obj-1",
        "object_type": "Points",
        "name": "Object 1",
        "parent": None,
        "revision": 4,
        "geometry": {"vertices": "v"},
        "data": [],
        "buffers": [],
    }


def _vertex_bytes() -> bytes:
    """Return the 48 bytes of a (4, 3) f32 array."""
    return np.arange(12, dtype="<f4").tobytes()


def _header_frame(
    obj: dict[str, Any] | None = None,
    arrays: list[dict[str, Any]] | None = None,
    origin: list[float] | None = None,
    crs: str | None = "EPSG:28350",
) -> Frame:
    """Return an ``object.header`` frame. The defaults give one f32 array, ``v``."""
    header: dict[str, Any] = {
        "type": "object.header",
        "request": 1,
        "origin": origin,
        "crs": crs,
        "object": _points_object() if obj is None else obj,
        "arrays": [{"name": "v", "dtype": "f32", "shape": [4, 3]}] if arrays is None else arrays,
    }
    return Frame(header=header, buffers={})


def _part_frame(
    array: str, offset: int, data: bytes, buffer: np.ndarray | None = None, size: int | None = None
) -> Frame:
    """Return an ``object.part`` frame for ``data``."""
    header = {
        "type": "object.part",
        "request": 1,
        "array": array,
        "offset": offset,
        "bytes": len(data) if size is None else size,
    }
    part = np.frombuffer(data, dtype=np.uint8) if buffer is None else buffer
    return Frame(header=header, buffers={"part": part})


def _end_frame(parts: int) -> Frame:
    """Return an ``object.end`` frame."""
    return Frame(header={"type": "object.end", "request": 1, "parts": parts}, buffers={})


def _vertex_parts(sizes: tuple[int, ...] = (16, 16, 16)) -> list[Frame]:
    """Return the ``object.part`` frames of ``v``, cut at ``sizes``."""
    data = _vertex_bytes()
    frames = []
    offset = 0
    for size in sizes:
        frames.append(_part_frame("v", offset, data[offset : offset + size]))
        offset += size
    return frames


def _assemble(frames: list[Frame]) -> Any:
    """Feed all frames to a new assembler and return the result of ``finish``."""
    assembler = fetch_module._Assembler("obj-1")
    for frame in frames:
        assembler.feed(frame)
    return assembler.finish()


def _stream(obj: Any, origin: tuple[float, float, float] = ORIGIN, part_size: int = 24):
    """Return the frames of a good reply to ``object.fetch`` for ``obj``."""
    frame = decode_frame(
        encode.encode_object(
            obj, id="obj-1", name="Object 1", parent="grp", revision=7, origin=origin
        )
    )
    arrays = []
    parts = []
    for name, array in frame.buffers.items():
        raw = np.ascontiguousarray(array).astype(array.dtype.newbyteorder("<")).tobytes()
        dtype = _PROTOCOL_DTYPES[array.dtype.name]
        arrays.append({"name": name, "dtype": dtype, "shape": [int(n) for n in array.shape]})
        for offset in range(0, len(raw), part_size):
            parts.append(_part_frame(name, offset, raw[offset : offset + part_size]))
    header = _header_frame(
        copy.deepcopy(dict(frame.header, buffers=[])), arrays, list(origin), "EPSG:28350"
    )
    return [header, *parts, _end_frame(len(parts))]


def _assert_surface(actual: Any, expected: Surface) -> None:
    """Check that two surfaces have the same vertices, triangles and properties."""
    assert isinstance(actual, Surface)
    np.testing.assert_allclose(actual.vertices, expected.vertices)
    np.testing.assert_array_equal(actual.triangles, expected.triangles)
    assert list(actual.properties) == list(expected.properties)
    for key, values in expected.properties.items():
        np.testing.assert_allclose(actual.properties[key], values)


def _surface(count: int = 40) -> Surface:
    """Return a surface with ``count`` vertices, a few triangles and one property."""
    vertices = np.arange(count * 3, dtype=np.float64).reshape(count, 3) * 0.5
    triangles = np.array([[i, i + 1, i + 2] for i in range(count - 2)])
    return Surface(vertices, triangles, {"p": np.arange(count, dtype=np.float64)})


def test_assembler_good_stream_gives_the_result():
    obj = Surface(VERTS, [[0, 1, 2]], {"v": np.array([1.0, 2.0, 3.0])})
    result = _assemble(_stream(obj))
    assert result.id == "obj-1"
    assert result.object_type == "Surface"
    assert result.revision == 7
    assert result.name == "Object 1"
    assert result.parent == "grp"
    assert result.origin == ORIGIN
    assert isinstance(result.origin, tuple)
    assert result.crs == "EPSG:28350"
    _assert_surface(result.object, obj)


def test_assembler_object_without_arrays():
    block = BlockModel([1.0, 2.0, 3.0], [1.0, 1.0, 1.0], [1, 1, 1])
    result = _assemble(_stream(block, origin=(0.0, 0.0, 0.0)))
    assert isinstance(result.object, BlockModel)
    np.testing.assert_allclose(result.object.origin, [1.0, 2.0, 3.0])


def test_assembler_first_frame_must_be_the_header():
    with pytest.raises(ProtocolError):
        _assemble([*_vertex_parts(), _end_frame(3)])


def test_assembler_header_too_large():
    arrays = [{"name": "v", "dtype": "f32", "shape": [3, 2**29]}]
    with pytest.raises(ObjectTooLargeError):
        _assemble([_header_frame(arrays=arrays)])


def test_assembler_sum_of_arrays_too_large():
    arrays = [
        {"name": "v", "dtype": "f32", "shape": [2**28, 2]},
        {"name": "w", "dtype": "f64", "shape": [2**27, 2]},
    ]
    with pytest.raises(ObjectTooLargeError):
        _assemble([_header_frame(arrays=arrays)])


def test_assembler_small_stream_with_the_revision_of_the_header():
    result = _assemble([_header_frame(), *_vertex_parts(), _end_frame(3)])
    assert result.revision == 4


def test_assembler_accepts_the_same_stream_with_u8_buffer_and_in_order_parts():
    """Positive controls for the cases of the bad streams."""
    assert _assemble([_header_frame(), *_vertex_parts((48,)), _end_frame(1)]).revision == 4
    assert _assemble([_header_frame(), *_vertex_parts((16, 16, 16)), _end_frame(3)]).revision == 4


def _bad_streams() -> dict[str, list[Frame]]:
    """Return streams that each break one rule, by a short name."""
    data = _vertex_bytes()
    wrong_dtype = np.arange(12, dtype=np.float32)
    first = _vertex_parts((16, 16, 16))[:1]
    return {
        "unknown_array": [_header_frame(), _part_frame("w", 0, data[:16]), _end_frame(1)],
        "offset_not_multiple_of_8": [_header_frame(), *_vertex_parts((4, 44)), _end_frame(2)],
        "past_the_end": [
            _header_frame(),
            *_vertex_parts((16, 16)),
            _part_frame("v", 32, data[:32]),
            _end_frame(3),
        ],
        "gap": [_header_frame(), *first, _part_frame("v", 24, data[24:]), _end_frame(2)],
        "overlap": [
            _header_frame(),
            _part_frame("v", 0, data[:24]),
            _part_frame("v", 16, data[16:]),
            _end_frame(2),
        ],
        "ranges_do_not_tile": [
            _header_frame(),
            _part_frame("v", 0, data[:16]),
            _part_frame("v", 8, data[8:24]),
            _part_frame("v", 32, data[32:48]),
            _end_frame(3),
        ],
        "missing_part": [_header_frame(), *_vertex_parts()[:2], _end_frame(2)],
        "too_many_parts_claimed": [_header_frame(), *_vertex_parts(), _end_frame(4)],
        "too_few_parts_claimed": [_header_frame(), *_vertex_parts(), _end_frame(2)],
        "part_not_u8": [
            _header_frame(),
            _part_frame("v", 0, b"", buffer=wrong_dtype, size=48),
            _end_frame(1),
        ],
    }


@pytest.mark.parametrize("case", sorted(_bad_streams()))
def test_assembler_rejects_a_bad_stream(case):
    with pytest.raises(ProtocolError):
        _assemble(_bad_streams()[case])


def test_assembler_geometry_buffer_not_in_arrays():
    obj = _points_object()
    obj["geometry"] = {"vertices": "missing"}
    with pytest.raises(ProtocolError):
        _assemble([_header_frame(obj), *_vertex_parts(), _end_frame(3)])


def test_assembler_data_buffer_not_in_arrays():
    obj = _points_object()
    obj["data"] = [{"id": "a", "name": "a", "buffer": "missing", "association": "vertex"}]
    with pytest.raises(ProtocolError):
        _assemble([_header_frame(obj), *_vertex_parts(), _end_frame(3)])


def test_assembler_offsets_buffer_not_in_arrays():
    obj = _points_object()
    obj["data"] = [
        {"id": "t", "name": "t", "buffer": "v", "offsets": "missing", "data_type": "text"}
    ]
    with pytest.raises(ProtocolError):
        _assemble([_header_frame(obj), *_vertex_parts(), _end_frame(3)])


# --- fetch_object and Viewer.fetch with the fake viewer ---


@pytest.fixture
def fake_factory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> Iterator[Callable[..., tuple[FakeViewer, Path]]]:
    """Return a function that starts a fake viewer with options. It gives the fake and its folder."""
    made: list[FakeViewer] = []

    def make(cls: type[FakeViewer] = FakeViewer, **options: Any) -> tuple[FakeViewer, Path]:
        fake = cls(**options)
        made.append(fake)
        settings = tmp_path / "geoviewer"
        fake.write_instance_file(settings)
        monkeypatch.setenv("GEOVIEWER_TOKEN", "secret")
        return fake, settings

    yield make
    for fake in made:
        fake.close()


@pytest.fixture
def tool_factory() -> Iterator[Callable[[Path], ToolConnection]]:
    """Return a function that makes tool connections. Each one is closed at teardown."""
    made: list[ToolConnection] = []

    def make(settings: Path) -> ToolConnection:
        tool = ToolConnection(settings_dir=settings, token="secret", timeout=WAIT)
        made.append(tool)
        return tool

    yield make
    for tool in made:
        tool.close()


@pytest.fixture
def own_viewers() -> Iterator[Callable[[Path], Viewer]]:
    """Return a function that makes viewers for a settings folder. Each is closed at teardown."""
    made: list[Viewer] = []

    def make(settings: Path) -> Viewer:
        viewer = Viewer(settings_dir=settings, heartbeat_interval=0.1, timeout=WAIT)
        made.append(viewer)
        return viewer

    yield make
    for viewer in made:
        with contextlib.suppress(Exception):
            viewer.close()


def _put(fake: FakeViewer, obj: Any, object_id: str = "s1", revision: int = 3) -> None:
    """Put an object into the fake viewer without a data connection."""
    fake.accept_frame(encode.hello_frame("", (0.0, 0.0, 0.0), "test"))
    fake.accept_frame(
        encode.encode_object(
            obj, id=object_id, name="Name 1", parent=None, revision=revision, origin=fake.origin
        )
    )


class _BadPartFake(FakeViewer):
    """A fake viewer that sends an ``object.part`` with an offset that is not valid."""

    def _fetch_frames(self, frame: Frame) -> Iterator[Frame]:
        for reply in super()._fetch_frames(frame):
            if reply.header["type"] == "object.part":
                reply = Frame(header=dict(reply.header, offset=3), buffers=reply.buffers)
            yield reply


def _same_objects(before: dict[str, Any], fake: FakeViewer) -> bool:
    """Return whether the fake holds the same stored objects as in ``before``."""
    return list(fake.objects) == list(before) and all(
        fake.objects[key] is before[key] for key in before
    )


def test_fetch_object_returns_the_sent_surface(fake_viewer, tool_factory):
    surface = _surface()
    _put(fake_viewer.viewer, surface)
    result = fetch_module.fetch_object(tool_factory(fake_viewer.settings), "s1", timeout=WAIT)
    _assert_surface(result.object, surface)
    assert result.id == "s1"
    assert result.object_type == "Surface"
    assert result.revision == 3
    assert result.name == "Name 1"
    assert result.parent is None
    assert result.origin == (0.0, 0.0, 0.0)
    assert result.crs == "EPSG:28350"


def test_fetch_object_with_the_origin_and_crs_of_the_viewer(fake_factory, tool_factory):
    fake, settings = fake_factory(origin=(500.0, 600.0, 700.0), crs="EPSG:32750")
    _put(fake, Surface(VERTS, [[0, 1, 2]]))
    result = fetch_module.fetch_object(tool_factory(settings), "s1", timeout=WAIT)
    assert result.origin == (500.0, 600.0, 700.0)
    assert result.crs == "EPSG:32750"
    np.testing.assert_allclose(result.object.vertices, VERTS)


def test_fetch_object_block_model_has_no_parts(fake_viewer, tool_factory):
    _put(fake_viewer.viewer, BlockModel([1.0, 2.0, 3.0], [1.0, 1.0, 1.0], [2, 2, 2]))
    result = fetch_module.fetch_object(tool_factory(fake_viewer.settings), "s1", timeout=WAIT)
    assert isinstance(result.object, BlockModel)
    np.testing.assert_array_equal(result.object.cells, [2, 2, 2])


def test_fetch_object_progress_rises_and_ends_at_one(fake_factory, tool_factory):
    fake, settings = fake_factory(part_bytes=64)
    surface = _surface()
    _put(fake, surface)
    calls: list[float] = []
    result = fetch_module.fetch_object(
        tool_factory(settings), "s1", progress=calls.append, timeout=WAIT
    )
    _assert_surface(result.object, surface)
    assert len(calls) > 1
    assert calls == sorted(calls)
    assert calls[-1] == pytest.approx(1.0)
    assert all(0.0 < value <= 1.0 for value in calls)


def test_fetch_object_cancel_set_before_the_call(fake_factory, tool_factory):
    fake, settings = fake_factory(part_bytes=64, fetch_delay=0.05)
    _put(fake, _surface())
    before = dict(fake.objects)
    cancel = threading.Event()
    cancel.set()
    with pytest.raises(CancelledError):
        fetch_module.fetch_object(tool_factory(settings), "s1", cancel=cancel, timeout=WAIT)
    assert _same_objects(before, fake)


def test_fetch_object_cancel_set_from_progress(fake_factory, tool_factory):
    fake, settings = fake_factory(part_bytes=64, fetch_delay=0.05)
    _put(fake, _surface())
    before = dict(fake.objects)
    cancel = threading.Event()
    with pytest.raises(CancelledError):
        fetch_module.fetch_object(
            tool_factory(settings),
            "s1",
            progress=lambda _fraction: cancel.set(),
            cancel=cancel,
            timeout=WAIT,
        )
    assert cancel.is_set()
    assert _same_objects(before, fake)


def test_fetch_object_unknown_id(fake_viewer, tool_factory):
    with pytest.raises(InvalidInputError):
        fetch_module.fetch_object(tool_factory(fake_viewer.settings), "nope", timeout=WAIT)


def test_fetch_object_headless_viewer(fake_factory, tool_factory):
    fake, settings = fake_factory(kind="headless")
    surface = _surface()
    _put(fake, surface)
    result = fetch_module.fetch_object(tool_factory(settings), "s1", timeout=WAIT)
    _assert_surface(result.object, surface)


def test_fetch_object_stream_error_closes_the_connection(fake_factory, tool_factory):
    fake, settings = fake_factory(_BadPartFake)
    _put(fake, _surface())
    tool = tool_factory(settings)
    with pytest.raises(ViewerClosedError) as info:
        fetch_module.fetch_object(tool, "s1", timeout=WAIT)
    assert str(info.value).startswith("The viewer sent an object that is not valid")
    assert not tool.connected


def test_fetch_object_can_run_twice(fake_viewer, tool_factory):
    _put(fake_viewer.viewer, _surface())
    tool = tool_factory(fake_viewer.settings)
    first = fetch_module.fetch_object(tool, "s1", timeout=WAIT)
    second = fetch_module.fetch_object(tool, "s1", timeout=WAIT)
    assert first.revision == second.revision


def test_viewer_fetch_returns_the_sent_surface(viewer_factory):
    viewer = viewer_factory(timeout=WAIT)
    surface = _surface()
    revision = viewer.send(surface, id="s1", name="Name 1")
    viewer.flush()
    result = viewer.fetch("s1")
    _assert_surface(result.object, surface)
    assert result.id == "s1"
    assert result.object_type == "Surface"
    assert result.revision == revision
    assert result.name == "Name 1"


def test_viewer_fetch_progress_and_parts(fake_factory, own_viewers):
    _fake, settings = fake_factory(part_bytes=64)
    viewer = own_viewers(settings)
    surface = _surface()
    viewer.send(surface, id="s1")
    viewer.flush()
    calls: list[float] = []
    result = viewer.fetch("s1", progress=calls.append)
    _assert_surface(result.object, surface)
    assert len(calls) > 1
    assert calls == sorted(calls)
    assert calls[-1] == pytest.approx(1.0)


def test_viewer_fetch_cancel(fake_factory, own_viewers):
    fake, settings = fake_factory(part_bytes=64, fetch_delay=0.05)
    viewer = own_viewers(settings)
    viewer.send(_surface(), id="s1")
    viewer.flush()
    before = dict(fake.objects)
    cancel = threading.Event()
    with pytest.raises(CancelledError):
        viewer.fetch("s1", progress=lambda _fraction: cancel.set(), cancel=cancel)
    assert _same_objects(before, fake)


def test_viewer_fetch_unknown_id(viewer_factory):
    viewer = viewer_factory(timeout=WAIT)
    with pytest.raises(InvalidInputError):
        viewer.fetch("nope")


def test_viewer_fetch_headless_viewer(fake_factory, own_viewers):
    _fake, settings = fake_factory(kind="headless")
    viewer = own_viewers(settings)
    surface = _surface()
    viewer.send(surface, id="s1")
    viewer.flush()
    _assert_surface(viewer.fetch("s1").object, surface)


def test_viewer_fetch_stream_error_closes_the_tool_connection(fake_factory, own_viewers):
    _fake, settings = fake_factory(_BadPartFake)
    viewer = own_viewers(settings)
    viewer.send(_surface(), id="s1")
    viewer.flush()
    with pytest.raises(ViewerClosedError) as info:
        viewer.fetch("s1")
    assert str(info.value).startswith("The viewer sent an object that is not valid")
    assert not viewer.tool_connected


def test_viewer_fetch_gives_the_revision(viewer_factory):
    viewer = viewer_factory(timeout=WAIT)
    viewer.send(_surface(), id="s1")
    viewer.send(_surface(), id="s1")
    revision = viewer.send(_surface(), id="s1")
    viewer.flush()
    assert viewer.fetch("s1").revision == revision == 3


def test_viewer_fetch_while_another_program_holds_the_data_connection(fake_viewer, viewer_factory):
    program_a = viewer_factory(timeout=WAIT)
    program_b = viewer_factory(timeout=WAIT)
    surface = _surface()
    program_a.send(surface, id="a1")
    program_a.flush()
    with pytest.raises(BusyError):
        program_b.send(surface, id="b1")
    assert program_b.select("a1", cells=[0, 1]) == {"cell": 2, "vertex": 0}
    assert program_b.selection().object_id == "a1"
    _assert_surface(program_b.fetch("a1").object, surface)
    program_a.release_data()
    deadline = time.monotonic() + WAIT
    while True:
        try:
            program_b.send(surface, id="b1")
            break
        except BusyError:
            if time.monotonic() > deadline:
                raise
            time.sleep(0.05)
    program_b.flush()
    assert {"a1", "b1"} <= set(fake_viewer.viewer.objects)
