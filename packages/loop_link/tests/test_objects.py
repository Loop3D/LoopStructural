"""Tests for the loop_link object model."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest
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

from loop_link.errors import InvalidInputError

VERTS = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]


def make_points(**kw):
    return Points(vertices=VERTS, **kw)


def make_lines(**kw):
    return Lines(vertices=VERTS, segments=[[0, 1], [1, 2]], **kw)


def make_surface(**kw):
    return Surface(vertices=VERTS, triangles=[[0, 1, 2], [0, 1, 3]], **kw)


def make_block(**kw):
    args = dict(origin=[0, 0, 0], cell_sizes=[1, 1, 1], cells=[2, 2, 2])
    args.update(kw)
    return BlockModel(**args)


def make_tetra(**kw):
    return TetraMesh(vertices=VERTS, tetrahedra=[[0, 1, 2, 3]], **kw)


def make_drillhole(**kw):
    args = dict(
        collar=[0, 0, 0],
        survey=[[0, 0, 90], [10, 0, 80], [20, 0, 70]],
        intervals=[[0, 5], [5, 20]],
    )
    args.update(kw)
    return Drillhole(**args)


# --- Valid input ---


def test_valid_objects_are_accepted_and_use_numpy():
    objs = [
        make_points(properties={"a": [1, 2, 3, 4]}),
        Vectors(vertices=VERTS, orientation=([0] * 4, [0] * 4)),
        make_lines(properties={"a": [1, 2, 3, 4]}, cell_properties={"b": [1, 2]}),
        make_surface(properties={"a": [1, 2, 3, 4]}, cell_properties={"b": [1, 2]}),
        make_block(properties={"a": np.arange(27.0)}, cell_properties={"b": np.arange(8.0)}),
        make_tetra(properties={"a": [1, 2, 3, 4]}, cell_properties={"b": [1]}),
        make_drillhole(
            properties={"a": [1, 2, 3]},
            interval_properties={"b": [1, 2]},
            end_of_hole=20.0,
        ),
    ]
    assert len(objs) == 7
    assert isinstance(objs[0].vertices, np.ndarray)
    assert isinstance(objs[0].properties["a"], np.ndarray)
    assert isinstance(objs[2].segments, np.ndarray)
    assert isinstance(objs[3].triangles, np.ndarray)
    assert isinstance(objs[4].origin, np.ndarray)
    assert isinstance(objs[4].cell_sizes, np.ndarray)
    assert isinstance(objs[5].tetrahedra, np.ndarray)
    assert isinstance(objs[6].collar, np.ndarray)
    assert isinstance(objs[6].survey, np.ndarray)
    assert isinstance(objs[6].intervals, np.ndarray)


@pytest.mark.parametrize(
    ("factory", "field"),
    [
        (make_points, "vertices"),
        (lambda: Vectors(vertices=VERTS), "vertices"),
        (make_lines, "vertices"),
        (make_surface, "vertices"),
        (make_block, "cells"),
        (make_tetra, "vertices"),
        (make_drillhole, "collar"),
    ],
)
def test_objects_are_frozen(factory, field):
    obj = factory()
    with pytest.raises(dataclasses.FrozenInstanceError):
        setattr(obj, field, None)


def test_default_properties_are_not_shared():
    first = make_points()
    second = make_points()
    assert first.properties is not second.properties
    first_vec = Vectors(vertices=VERTS)
    second_vec = Vectors(vertices=VERTS)
    assert first_vec.properties is not second_vec.properties


def test_categorical_valid():
    cat = Categorical(codes=[0, 1, 1], names=["a", "b"])
    assert isinstance(cat.codes, np.ndarray)
    assert cat.colours is None


def test_default_category_colours_are_tab10():
    expected = [
        (31, 119, 180),
        (255, 127, 14),
        (44, 160, 44),
        (214, 39, 40),
        (148, 103, 189),
        (140, 86, 75),
        (227, 119, 194),
        (127, 127, 127),
        (188, 189, 34),
        (23, 190, 207),
    ]
    assert [tuple(c) for c in DEFAULT_CATEGORY_COLOURS] == [(*rgb, 255) for rgb in expected]


def test_default_category_colours():
    assert len(DEFAULT_CATEGORY_COLOURS) == 10
    for colour in DEFAULT_CATEGORY_COLOURS:
        assert len(colour) == 4


# --- Vertices ---


@pytest.mark.parametrize(
    "factory",
    [make_points, make_lines, make_surface, make_tetra],
)
def test_vertices_wrong_shape(factory):
    with pytest.raises(InvalidInputError, match="vertices"):
        # Two columns are not valid.
        cls = type(factory())
        kw = {"vertices": [[0, 0], [1, 1]]}
        if cls is Lines:
            kw["segments"] = [[0, 1]]
        elif cls is Surface:
            kw["triangles"] = [[0, 1, 0]]
        elif cls is TetraMesh:
            kw["tetrahedra"] = [[0, 1, 0, 1]]
        cls(**kw)


def test_vertices_one_dimensional():
    with pytest.raises(InvalidInputError, match="vertices"):
        Points(vertices=[0, 0, 0])


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_vertices_not_finite(bad):
    with pytest.raises(InvalidInputError, match="vertices"):
        Points(vertices=[[0, 0, bad]])


# --- Index arrays ---


@pytest.mark.parametrize("bad", [-1, 4])
def test_segments_index_out_of_range(bad):
    with pytest.raises(InvalidInputError, match="segments"):
        Lines(vertices=VERTS, segments=[[0, bad]])


@pytest.mark.parametrize("bad", [-1, 4])
def test_triangles_index_out_of_range(bad):
    with pytest.raises(InvalidInputError, match="triangles"):
        Surface(vertices=VERTS, triangles=[[0, 1, bad]])


@pytest.mark.parametrize("bad", [-1, 4])
def test_tetrahedra_index_out_of_range(bad):
    with pytest.raises(InvalidInputError, match="tetrahedra"):
        TetraMesh(vertices=VERTS, tetrahedra=[[0, 1, 2, bad]])


# --- Property length ---


def test_points_property_wrong_length():
    with pytest.raises(InvalidInputError, match="properties"):
        make_points(properties={"bad": [1, 2, 3]})


def test_categorical_property_valid():
    cat = Categorical(codes=[0, 1, 1, 0], names=["a", "b"])
    pts = make_points(properties={"c": cat})
    assert pts.properties["c"] is cat


def test_categorical_property_wrong_length():
    cat = Categorical(codes=[0, 1, 1], names=["a", "b"])
    with pytest.raises(InvalidInputError, match="properties"):
        make_points(properties={"c": cat})


def test_float_index_array_refused():
    with pytest.raises(InvalidInputError, match="segments"):
        Lines(vertices=VERTS, segments=[[0.5, 1]])


def test_vectors_property_wrong_length():
    with pytest.raises(InvalidInputError, match="properties"):
        Vectors(vertices=VERTS, orientation=([0] * 4, [0] * 4), properties={"bad": [1]})


@pytest.mark.parametrize("factory", [make_lines, make_surface, make_tetra])
def test_vertex_property_wrong_length(factory):
    with pytest.raises(InvalidInputError, match="properties"):
        factory(properties={"bad": [1, 2, 3]})


def test_lines_cell_property_wrong_length():
    with pytest.raises(InvalidInputError, match="cell_properties"):
        make_lines(cell_properties={"bad": [1, 2, 3]})


def test_surface_cell_property_wrong_length():
    with pytest.raises(InvalidInputError, match="cell_properties"):
        make_surface(cell_properties={"bad": [1, 2, 3]})


def test_tetra_cell_property_wrong_length():
    with pytest.raises(InvalidInputError, match="cell_properties"):
        make_tetra(cell_properties={"bad": [1, 2]})


def test_block_node_property_wrong_length():
    # Nodes: (2 + 1) * (2 + 1) * (2 + 1) = 27. Eight values are the cell count.
    with pytest.raises(InvalidInputError, match="properties"):
        make_block(properties={"bad": np.arange(8.0)})


def test_block_cell_property_wrong_length():
    with pytest.raises(InvalidInputError, match="cell_properties"):
        make_block(cell_properties={"bad": np.arange(27.0)})


def test_block_property_lengths_with_unequal_cells():
    # Nodes: 3 * 4 * 5 = 60. Cells: 2 * 3 * 4 = 24.
    block = make_block(
        cells=[2, 3, 4],
        properties={"n": np.arange(60.0)},
        cell_properties={"c": np.arange(24.0)},
    )
    assert block.properties["n"].shape == (60,)
    with pytest.raises(InvalidInputError, match="cell_properties"):
        make_block(cells=[2, 3, 4], cell_properties={"c": np.arange(60.0)})


def test_drillhole_property_wrong_length():
    # Three survey stations are in the survey.
    with pytest.raises(InvalidInputError, match="properties"):
        make_drillhole(properties={"bad": [1, 2]})


def test_drillhole_interval_property_wrong_length():
    # Two intervals are in the data.
    with pytest.raises(InvalidInputError, match="interval_properties"):
        make_drillhole(interval_properties={"bad": [1, 2, 3]})


def test_two_dimensional_property():
    with pytest.raises(InvalidInputError, match="bad"):
        make_points(properties={"bad": np.zeros((4, 2))})


# --- Block model sizes ---


@pytest.mark.parametrize("sizes", [[0, 1, 1], [1, -1, 1], [1, 1, 0]])
def test_block_cell_sizes_not_positive(sizes):
    with pytest.raises(InvalidInputError, match="cell_sizes"):
        make_block(cell_sizes=sizes)


@pytest.mark.parametrize("cells", [[0, 1, 1], [1, 1, -2]])
def test_block_cells_less_than_one(cells):
    with pytest.raises(InvalidInputError, match="cells"):
        make_block(cells=cells)


def test_block_cells_not_integer():
    with pytest.raises(InvalidInputError, match="cells"):
        make_block(cells=[1.5, 1, 1])


# --- Fields that do not exist ---


def test_points_has_no_cell_properties():
    with pytest.raises(TypeError):
        Points(vertices=VERTS, cell_properties={})


def test_vectors_has_no_cell_properties():
    with pytest.raises(TypeError):
        Vectors(vertices=VERTS, cell_properties={})


# --- Categorical ---


def test_categorical_code_too_large():
    with pytest.raises(InvalidInputError, match="codes"):
        Categorical(codes=[0, 2], names=["a", "b"])


@pytest.mark.parametrize("bad", [-1, 2**32])
def test_categorical_code_out_of_range(bad):
    with pytest.raises(InvalidInputError, match="codes"):
        Categorical(codes=[0, bad], names=["a", "b"])


# --- Vectors.from_directions ---


def test_from_directions_bad_kind():
    with pytest.raises(InvalidInputError, match="kind"):
        Vectors.from_directions(VERTS, [[0, 0, 1]] * 4, kind="other")


def test_from_directions_valid_kind():
    vec = Vectors.from_directions(VERTS, [[0, 0, 1]] * 4, kind="normal")
    # Orientation holds two arrays, one value for each vertex.
    assert len(vec.orientation) == 2
    assert all(len(arr) == 4 for arr in vec.orientation)


def test_from_directions_vertical_normal():
    # A normal that points up gives dip direction 0 and dip 0.
    vec = Vectors.from_directions(VERTS, [[0, 0, 1]] * 4, kind="normal")
    np.testing.assert_allclose(vec.orientation[0], [0] * 4, atol=1e-9)
    np.testing.assert_allclose(vec.orientation[1], [0] * 4, atol=1e-9)
    assert vec.lineation is None


def test_from_directions_tilted_normal():
    # Normal formula: n = (sin d sin a, sin d cos a, cos d).
    # Here a is the dip direction and d is the dip, both in degrees.
    a, d = 40.0, 30.0
    ar, dr = np.radians(a), np.radians(d)
    n = [np.sin(dr) * np.sin(ar), np.sin(dr) * np.cos(ar), np.cos(dr)]
    vec = Vectors.from_directions(VERTS, [n] * 4, kind="normal")
    np.testing.assert_allclose(vec.orientation[0], [a] * 4, atol=1e-6)
    np.testing.assert_allclose(vec.orientation[1], [d] * 4, atol=1e-6)


def test_from_directions_lineation():
    # Lineation formula: direction = (cos p sin t, cos p cos t, -sin p).
    # Here t is the trend and p is the plunge, both in degrees.
    t, p = 70.0, 25.0
    tr, pr = np.radians(t), np.radians(p)
    direction = [np.cos(pr) * np.sin(tr), np.cos(pr) * np.cos(tr), -np.sin(pr)]
    vec = Vectors.from_directions(VERTS, [direction] * 4, kind="lineation")
    assert vec.orientation is None
    assert len(vec.lineation) == 2
    np.testing.assert_allclose(vec.lineation[0], [t] * 4, atol=1e-6)
    np.testing.assert_allclose(vec.lineation[1], [p] * 4, atol=1e-6)


def test_from_directions_keeps_properties():
    props = {"a": [1, 2, 3, 4]}
    vec = Vectors.from_directions(VERTS, [[0, 0, 1]] * 4, kind="normal", properties=props)
    np.testing.assert_array_equal(vec.properties["a"], [1, 2, 3, 4])


# --- Vectors orientation and lineation ---


def test_vectors_with_orientation_and_lineation():
    vec = Vectors(vertices=VERTS, orientation=([0] * 4, [10] * 4), lineation=([5] * 4, [20] * 4))
    assert vec.orientation is not None
    assert vec.lineation is not None


@pytest.mark.parametrize("field", ["orientation", "lineation"])
def test_vectors_wrong_length_arrays(field):
    with pytest.raises(InvalidInputError, match=field):
        Vectors(vertices=VERTS, **{field: ([0] * 3, [0] * 3)})
