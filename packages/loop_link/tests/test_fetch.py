"""Tests for the decoding of viewer objects into scene objects."""

from __future__ import annotations

import copy
import warnings

import numpy as np
import pytest

from loop_link import encode
from loop_link.errors import InvalidInputError, UnsupportedObjectError
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
from loop_link.protocol import decode_frame

ORIGIN = (1000.0, 2000.0, 3000.0)
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
