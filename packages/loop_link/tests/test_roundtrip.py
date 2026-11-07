"""Tests that an object sent to the fake viewer and fetched again stays the same (SC-002)."""

from __future__ import annotations

import contextlib
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import numpy as np
import pytest

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
from loop_link.testing import FakeViewer
from loop_link.viewer import Viewer

WAIT = 3.0
ORIGIN = (500000.0, 7000000.0, 0.0)
NAMES = ["granite", "basalt", "shale"]
COLOURS = [(10, 20, 30, 255), (40, 50, 60, 128), (70, 80, 90, 0)]
# The dtype that each property has after the round trip.
EXPECTED_DTYPES = {
    "f32": np.float32,
    "f64": np.float64,
    "int": np.int32,
    "flag": np.bool_,
}


@pytest.fixture(params=[(0.0, 0.0, 0.0), ORIGIN], ids=["origin_zero", "origin_far"])
def scene(
    request: pytest.FixtureRequest, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> Iterator[tuple[Viewer, np.ndarray]]:
    """Return a viewer that is connected to a fake viewer, and the origin of the scene."""
    origin = request.param
    fake = FakeViewer(origin=origin)
    settings = tmp_path / "geoviewer"
    viewer = None
    try:
        fake.write_instance_file(settings)
        monkeypatch.setenv("GEOVIEWER_TOKEN", "secret")
        viewer = Viewer(settings_dir=settings, heartbeat_interval=0.1, timeout=WAIT)
        yield viewer, np.array(origin)
    finally:
        if viewer is not None:
            with contextlib.suppress(Exception):  # errors on close do not matter at teardown
                viewer.close()
        fake.close()


def _properties(count: int, prefix: str = "") -> dict[str, Any]:
    """Return one property of each kind, with ``count`` values each.

    A property name is the kind with ``prefix`` before it. Vertex and cell properties need
    different names.
    """
    index = np.arange(count)
    props = {
        "f32": (np.linspace(0.1, 1.7, count) * 1.1).astype(np.float32),
        "f64": np.linspace(-3.3, 9.9, count) / 7.0,
        "int": (index * 1000003 - 5).astype(np.int64),
        "flag": index % 2 == 0,
        "text": np.array([f"item-{i}-é" for i in index]),
        "cat": Categorical(index % 3, NAMES, COLOURS),
    }
    return {prefix + key: value for key, value in props.items()}


def _assert_properties(actual: dict[str, Any], expected: dict[str, Any]) -> None:
    """Check that two mappings of properties have the same names, values and dtypes."""
    assert list(actual) == list(expected)
    for key, want in expected.items():
        got = actual[key]
        if isinstance(want, Categorical):
            assert isinstance(got, Categorical)
            np.testing.assert_array_equal(got.codes, want.codes)
            assert list(got.names) == list(want.names)
            assert got.colours == want.colours
        elif want.dtype.kind == "f":
            assert got.dtype == want.dtype
            np.testing.assert_allclose(got, want, rtol=1e-6)
        elif want.dtype.kind == "U":
            assert got.dtype.kind == "U"
            assert got.tolist() == want.tolist()
        else:
            assert got.dtype == EXPECTED_DTYPES[key.removeprefix("c_")]
            np.testing.assert_array_equal(got, want)


def _assert_coordinates(actual: np.ndarray, expected: np.ndarray) -> None:
    """Check float64 world coordinates against the input."""
    assert actual.dtype == np.float64
    assert actual.shape == expected.shape
    np.testing.assert_allclose(actual, expected, rtol=1e-6)
    # The offset to the scene origin is small, so f32 keeps the error well below one centimetre.
    assert np.max(np.abs(actual - expected)) < 1e-2


def _vertices(origin: np.ndarray, count: int = 6) -> np.ndarray:
    """Return world coordinates near the scene origin."""
    rel = np.arange(count * 3, dtype=np.float64).reshape(count, 3) * 1.37 + 0.123
    rel[:, 2] += np.arange(count) ** 2 * 0.31
    return origin + rel


def _round_trip(viewer: Viewer, obj: Any, object_id: str) -> Any:
    """Send ``obj`` to the viewer, fetch it again and return the fetched object."""
    viewer.send(obj, id=object_id)
    viewer.flush()
    return viewer.fetch(object_id).object


def test_points(scene):
    viewer, origin = scene
    verts = _vertices(origin)
    far = origin + np.array([10000.123456789, -10000.987654321, 10000.5])
    verts = np.vstack([verts, far])
    obj = Points(verts, _properties(len(verts)))
    result = _round_trip(viewer, obj, "points")
    assert isinstance(result, Points)
    _assert_coordinates(result.vertices, verts)
    np.testing.assert_allclose(result.vertices[-1], far, rtol=1e-12)
    _assert_properties(result.properties, obj.properties)


def test_vectors(scene):
    viewer, origin = scene
    verts = _vertices(origin)
    normals = np.array(
        [
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.3, -0.4, 0.8],
            [-1.0, 2.0, 3.0],
            [0, 1, 0.1],
            [1, 1, 1],
        ]
    )
    directions = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, -0.5],
            [1.0, 1.0, 0.0],
            [-1.0, 0.2, -0.3],
            [0, 0, 1],
            [2, 1, 1],
        ]
    )
    orientation = Vectors.from_directions(verts, normals, kind="normal").orientation
    lineation = Vectors.from_directions(verts, directions, kind="lineation").lineation
    obj = Vectors(
        verts, orientation=orientation, lineation=lineation, properties=_properties(len(verts))
    )
    result = _round_trip(viewer, obj, "vectors")
    assert isinstance(result, Vectors)
    _assert_coordinates(result.vertices, verts)
    for name in ("orientation", "lineation"):
        for got, want in zip(getattr(result, name), getattr(obj, name)):
            assert got.dtype == np.float64
            np.testing.assert_allclose(got, want, rtol=0, atol=1e-6)
    _assert_properties(result.properties, obj.properties)


def test_lines(scene):
    viewer, origin = scene
    verts = _vertices(origin)
    segments = np.array([[0, 1], [1, 2], [2, 3], [3, 4], [4, 5]])
    obj = Lines(
        verts, segments, properties=_properties(6), cell_properties=_properties(len(segments), "c_")
    )
    result = _round_trip(viewer, obj, "lines")
    assert isinstance(result, Lines)
    _assert_coordinates(result.vertices, verts)
    np.testing.assert_array_equal(result.segments, segments)
    _assert_properties(result.properties, obj.properties)
    _assert_properties(result.cell_properties, obj.cell_properties)


def test_surface(scene):
    viewer, origin = scene
    verts = _vertices(origin)
    triangles = np.array([[0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 5]])
    obj = Surface(
        verts,
        triangles,
        properties=_properties(6),
        cell_properties=_properties(len(triangles), "c_"),
    )
    result = _round_trip(viewer, obj, "surface")
    assert isinstance(result, Surface)
    _assert_coordinates(result.vertices, verts)
    np.testing.assert_array_equal(result.triangles, triangles)
    _assert_properties(result.properties, obj.properties)
    _assert_properties(result.cell_properties, obj.cell_properties)


def test_block_model(scene):
    viewer, origin = scene
    corner = origin + np.array([100.5, 200.25, -50.125])
    cells = np.array([2, 3, 2])
    obj = BlockModel(
        corner,
        [10.0, 20.5, 5.25],
        cells,
        properties=_properties(int(np.prod(cells + 1))),
        cell_properties=_properties(int(np.prod(cells)), "c_"),
    )
    result = _round_trip(viewer, obj, "block")
    assert isinstance(result, BlockModel)
    assert result.origin.dtype == np.float64
    np.testing.assert_allclose(result.origin, corner, rtol=1e-6)
    np.testing.assert_allclose(result.cell_sizes, obj.cell_sizes, rtol=1e-6)
    np.testing.assert_array_equal(result.cells, cells)
    _assert_properties(result.properties, obj.properties)
    _assert_properties(result.cell_properties, obj.cell_properties)


def test_tetra_mesh(scene):
    viewer, origin = scene
    verts = _vertices(origin)
    tetrahedra = np.array([[0, 1, 2, 3], [1, 2, 3, 4], [2, 3, 4, 5]])
    obj = TetraMesh(
        verts,
        tetrahedra,
        properties=_properties(6),
        cell_properties=_properties(len(tetrahedra), "c_"),
    )
    result = _round_trip(viewer, obj, "tetra")
    assert isinstance(result, TetraMesh)
    _assert_coordinates(result.vertices, verts)
    np.testing.assert_array_equal(result.tetrahedra, tetrahedra)
    _assert_properties(result.properties, obj.properties)
    _assert_properties(result.cell_properties, obj.cell_properties)


def test_drillhole(scene):
    viewer, origin = scene
    collar = origin + np.array([12.5, 30.25, 400.0])
    survey = collar + np.array([[0.0, 0.0, 0.0], [0.5, 0.2, -10.0], [1.1, 0.4, -20.0], [2, 1, -30]])
    intervals = np.array([[0.0, 8.5], [8.5, 17.25], [17.25, 30.0]])
    obj = Drillhole(
        collar,
        survey,
        intervals,
        end_of_hole=30.0,
        interval_properties=_properties(len(intervals)),
    )
    result = _round_trip(viewer, obj, "hole")
    assert isinstance(result, Drillhole)
    _assert_coordinates(result.collar, obj.collar)
    _assert_coordinates(result.survey, survey)
    np.testing.assert_allclose(result.intervals, intervals, rtol=1e-6)
    assert result.end_of_hole == pytest.approx(30.0, rel=1e-6)
    _assert_properties(result.interval_properties, obj.interval_properties)
