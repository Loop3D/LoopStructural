"""Tests for the pyvista adapter."""

from __future__ import annotations

import time
from typing import Any

import numpy as np
import pytest

pv = pytest.importorskip("pyvista")

from loop_link.adapters.pyvista import from_pyvista

from loop_link.errors import InvalidInputError
from loop_link.objects import (
    BlockModel,
    Lines,
    Points,
    Surface,
    TetraMesh,
    Vectors,
)

POINTS = np.array(
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64
)


def wait_for_object(fake_viewer: Any, viewer: Any, object_id: str) -> Any:
    """Flush the viewer and return the stored object with this id."""
    viewer.flush()
    deadline = time.monotonic() + 5.0
    while time.monotonic() < deadline and object_id not in fake_viewer.viewer.objects:
        time.sleep(0.01)
    return fake_viewer.viewer.objects[object_id]


def stored_array(stored: Any, name: str) -> np.ndarray:
    """Return the buffer of the data description with this name."""
    entry = next(item for item in stored.header["data"] if item["name"] == name)
    return stored.buffers[entry["buffer"]]


def test_points_with_point_data() -> None:
    """A vertex-only PolyData gives Points, and point data gives properties."""
    mesh = pv.PolyData(POINTS.copy())
    mesh.point_data["grade"] = np.array([1.0, 2.0, 3.0, 4.0])
    result = from_pyvista(mesh)
    assert isinstance(result, Points)
    np.testing.assert_allclose(result.vertices, POINTS)
    np.testing.assert_allclose(result.properties["grade"], [1.0, 2.0, 3.0, 4.0])


def test_vectors_from_normals() -> None:
    """A point array used as vectors gives Vectors with dip from arccos(n_z)."""
    normals = np.array([[0.0, 0.0, 1.0], [0.6, 0.0, 0.8], [0.0, 0.5, np.sqrt(0.75)]])
    mesh = pv.PolyData(POINTS[:3].copy())
    mesh.point_data["n"] = normals
    mesh.point_data["grade"] = np.array([1.0, 2.0, 3.0])
    result = from_pyvista(mesh, vectors="n")
    assert isinstance(result, Vectors)
    assert result.orientation is not None
    dip_direction, dip = result.orientation
    np.testing.assert_allclose(dip, np.degrees(np.arccos(normals[:, 2])), atol=1e-6)
    # The normal of point 1 leans to +x (east) and of point 2 to +y (north).
    np.testing.assert_allclose(dip_direction[1], 90.0, atol=1e-6)
    np.testing.assert_allclose(dip_direction[2], 0.0, atol=1e-6)
    assert "n" not in result.properties
    np.testing.assert_allclose(result.properties["grade"], [1.0, 2.0, 3.0])


def test_lines() -> None:
    """A polyline gives Lines with segments of shape (m, 2)."""
    mesh = pv.lines_from_points(POINTS.copy())
    result = from_pyvista(mesh)
    assert isinstance(result, Lines)
    np.testing.assert_allclose(result.vertices, POINTS)
    assert np.asarray(result.segments).shape == (3, 2)
    np.testing.assert_array_equal(result.segments, [[0, 1], [1, 2], [2, 3]])


def test_surface_from_triangles() -> None:
    """A triangulated sphere gives Surface, and cell data gives cell_properties."""
    mesh = pv.Sphere().triangulate()
    mesh.clear_data()
    cell_values = np.arange(mesh.n_cells, dtype=np.float64)
    mesh.cell_data["area"] = cell_values
    result = from_pyvista(mesh)
    assert isinstance(result, Surface)
    assert np.asarray(result.vertices).shape == (mesh.n_points, 3)
    assert np.asarray(result.triangles).shape == (mesh.n_cells, 3)
    np.testing.assert_allclose(result.vertices, np.asarray(mesh.points))
    np.testing.assert_allclose(result.cell_properties["area"], cell_values)


def test_surface_from_quads_is_triangulated() -> None:
    """A quad plane is triangulated: 2 triangles for each quad."""
    mesh = pv.Plane(i_resolution=2, j_resolution=2)
    mesh.clear_data()
    assert not mesh.is_all_triangles
    result = from_pyvista(mesh)
    assert isinstance(result, Surface)
    assert np.asarray(result.triangles).shape == (2 * mesh.n_cells, 3)


def test_image_data_gives_block_model() -> None:
    """An ImageData gives BlockModel with node data and cell data."""
    grid = pv.ImageData(dimensions=(4, 5, 6), spacing=(1, 2, 3), origin=(10, 20, 30))
    node_values = np.arange(grid.n_points, dtype=np.float64)
    cell_values = np.arange(grid.n_cells, dtype=np.float64) * 2.0
    grid.point_data["node"] = node_values
    grid.cell_data["cell"] = cell_values
    result = from_pyvista(grid)
    assert isinstance(result, BlockModel)
    np.testing.assert_allclose(result.origin, [10.0, 20.0, 30.0])
    np.testing.assert_allclose(result.cell_sizes, [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(result.cells, [3, 4, 5])
    np.testing.assert_allclose(result.properties["node"], node_values)
    np.testing.assert_allclose(result.cell_properties["cell"], cell_values)


def test_rectilinear_grid_uneven_is_refused() -> None:
    """A RectilinearGrid with uneven spacing gives InvalidInputError."""
    grid = pv.RectilinearGrid([0.0, 1.0, 3.0], [0.0, 1.0, 2.0], [0.0, 1.0, 2.0])
    with pytest.raises(InvalidInputError):
        from_pyvista(grid)


def test_structured_grid_is_refused() -> None:
    """A StructuredGrid gives InvalidInputError."""
    xs, ys, zs = np.meshgrid(np.arange(3.0), np.arange(3.0), np.arange(3.0), indexing="ij")
    grid = pv.StructuredGrid(xs, ys, zs)
    with pytest.raises(InvalidInputError):
        from_pyvista(grid)


def test_tetra_unstructured_grid() -> None:
    """An UnstructuredGrid of tetrahedra gives TetraMesh with tetrahedra (m, 4)."""
    points = np.vstack([POINTS, [[1.0, 1.0, 1.0]]])
    cells = np.array([4, 0, 1, 2, 3, 4, 1, 2, 3, 4])
    grid = pv.UnstructuredGrid(cells, np.array([pv.CellType.TETRA] * 2), points)
    grid.point_data["node"] = np.arange(5, dtype=np.float64)
    grid.cell_data["cell"] = np.array([7.0, 8.0])
    result = from_pyvista(grid)
    assert isinstance(result, TetraMesh)
    np.testing.assert_allclose(result.vertices, points)
    np.testing.assert_array_equal(result.tetrahedra, [[0, 1, 2, 3], [1, 2, 3, 4]])
    np.testing.assert_allclose(result.properties["node"], np.arange(5.0))
    np.testing.assert_allclose(result.cell_properties["cell"], [7.0, 8.0])


def test_mixed_unstructured_grid_is_refused() -> None:
    """An UnstructuredGrid with tetrahedra and triangles gives InvalidInputError."""
    cells = np.array([4, 0, 1, 2, 3, 3, 0, 1, 2])
    types = np.array([pv.CellType.TETRA, pv.CellType.TRIANGLE])
    grid = pv.UnstructuredGrid(cells, types, POINTS.copy())
    with pytest.raises(InvalidInputError):
        from_pyvista(grid)


def test_drillhole_kind_is_refused_with_hint() -> None:
    """The kind Drillhole gives InvalidInputError that names the right function."""
    mesh = pv.PolyData(POINTS.copy())
    with pytest.raises(InvalidInputError) as info:
        from_pyvista(mesh, kind="Drillhole")
    assert "Drillhole" in str(info.value)
    assert "drillholes_from_pandas" in str(info.value)


def test_two_dimensional_property_is_refused() -> None:
    """A point array of shape (n, 2) gives InvalidInputError."""
    mesh = pv.PolyData(POINTS.copy())
    mesh.point_data["pair"] = np.ones((4, 2))
    with pytest.raises(InvalidInputError):
        from_pyvista(mesh)


def test_send_grid_function_in_node_order(fake_viewer: Any, viewer_factory: Any) -> None:
    """The stored buffer of f = x + 10 y + 100 z has the node order x fastest."""
    grid = pv.ImageData(dimensions=(4, 5, 6), spacing=(1, 2, 3), origin=(10, 20, 30))
    nodes = np.asarray(grid.points)
    grid.point_data["f"] = nodes[:, 0] + 10.0 * nodes[:, 1] + 100.0 * nodes[:, 2]
    # Expected values come from explicit loops: x changes fastest, z slowest.
    expected = np.array(
        [
            (10.0 + i) + 10.0 * (20.0 + 2.0 * j) + 100.0 * (30.0 + 3.0 * k)
            for k in range(6)
            for j in range(5)
            for i in range(4)
        ]
    )
    viewer = viewer_factory()
    viewer.send(grid, id="grid")
    stored = wait_for_object(fake_viewer, viewer, "grid")
    assert stored.header["object_type"] == "BlockModel"
    np.testing.assert_allclose(stored.header["params"]["origin"], [10.0, 20.0, 30.0])
    np.testing.assert_allclose(stored.header["params"]["cell_sizes"], [1.0, 2.0, 3.0])
    assert list(stored.header["params"]["cells"]) == [3, 4, 5]
    np.testing.assert_allclose(stored_array(stored, "f"), expected)
