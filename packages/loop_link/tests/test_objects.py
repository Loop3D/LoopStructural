"""Tests for the loop_link object model."""

from __future__ import annotations

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
        Vectors(vertices=VERTS, orientation=[[1, 0, 0]] * 4),
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


def test_objects_are_frozen():
    pts = make_points()
    with pytest.raises(AttributeError):
        pts.vertices = np.zeros((1, 3))


def test_categorical_valid():
    cat = Categorical(codes=[0, 1, 1], names=["a", "b"])
    assert isinstance(cat.codes, np.ndarray)
    assert cat.colours is None


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
    with pytest.raises(InvalidInputError, match="bad"):
        make_points(properties={"bad": [1, 2, 3]})


def test_vectors_property_wrong_length():
    with pytest.raises(InvalidInputError, match="bad"):
        Vectors(vertices=VERTS, orientation=[[1, 0, 0]] * 4, properties={"bad": [1]})


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


# --- Vectors.from_directions ---


def test_from_directions_bad_kind():
    with pytest.raises(InvalidInputError, match="kind"):
        Vectors.from_directions(VERTS, [[0, 0, 1]] * 4, kind="other")


def test_from_directions_valid_kind():
    vec = Vectors.from_directions(VERTS, [[0, 0, 1]] * 4, kind="normal")
    # Orientation holds two arrays, one value for each vertex.
    assert len(vec.orientation) == 2
    assert all(len(arr) == 4 for arr in vec.orientation)
