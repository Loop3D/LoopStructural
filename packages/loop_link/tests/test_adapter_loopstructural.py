"""Tests for the LoopStructural adapter (story 1)."""

from __future__ import annotations

import time
from collections.abc import Callable
from typing import Any

import numpy as np
import pytest

geometry = pytest.importorskip("loop_common.geometry")

from loop_link.adapters.loopstructural import from_loopstructural  # noqa: I001
from loop_link.errors import InvalidInputError
from loop_link.objects import BlockModel
from loop_link.objects import Surface as LinkSurface


def wait_for(condition: Callable[[], bool], timeout: float = 5.0) -> bool:
    """Poll a condition until it is true or the deadline passes."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if condition():
            return True
        time.sleep(0.01)
    return condition()


def stored_array(stored: Any, name: str) -> np.ndarray:
    """Return the buffer of the data description with this name."""
    entry = next(item for item in stored.header["data"] if item["name"] == name)
    return stored.buffers[entry["buffer"]]


VERTICES = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 1.0]])
TRIANGLES = np.array([[0, 1, 2], [0, 2, 3]])


def test_surface_keeps_arrays_and_properties() -> None:
    """A Surface gives a loop_link Surface with the same arrays and properties."""
    a = np.array([1.0, 2.0, 3.0, 4.0])
    b = np.array([10.0, 20.0])
    source = geometry.Surface(
        vertices=VERTICES, triangles=TRIANGLES, properties={"a": a}, cell_properties={"b": b}
    )
    result = from_loopstructural(source)
    assert isinstance(result, LinkSurface)
    np.testing.assert_array_equal(result.vertices, VERTICES)
    np.testing.assert_array_equal(result.triangles, TRIANGLES)
    assert list(result.properties) == ["a"]
    assert list(result.cell_properties) == ["b"]
    np.testing.assert_array_equal(result.properties["a"], a)
    np.testing.assert_array_equal(result.cell_properties["b"], b)


def test_surface_without_properties_gives_empty_mappings() -> None:
    """A Surface with properties=None gives an empty properties mapping."""
    source = geometry.Surface(vertices=VERTICES, triangles=TRIANGLES, properties=None)
    result = from_loopstructural(source)
    assert isinstance(result, LinkSurface)
    assert len(result.properties) == 0
    assert len(result.cell_properties) == 0


def test_structured_grid_gives_block_model() -> None:
    """Node counts become cell counts. Origin and step become origin and cell sizes."""
    grid = geometry.StructuredGrid(origin=[1, 2, 3], step_vector=[1, 2, 3], nsteps=[4, 5, 6])
    result = from_loopstructural(grid)
    assert isinstance(result, BlockModel)
    np.testing.assert_array_equal(result.origin, (1, 2, 3))
    np.testing.assert_array_equal(result.cell_sizes, (1, 2, 3))
    np.testing.assert_array_equal(result.cells, (3, 4, 5))


def test_surface_without_cell_properties_gives_empty_mapping() -> None:
    """A Surface with cell_properties=None gives an empty cell_properties mapping."""
    source = geometry.Surface(
        vertices=VERTICES,
        triangles=TRIANGLES,
        properties={"a": np.arange(4.0)},
        cell_properties=None,
    )
    result = from_loopstructural(source)
    assert isinstance(result, LinkSurface)
    assert len(result.cell_properties) == 0


def test_structured_grid_keeps_cell_properties() -> None:
    """A cell property keeps its values and its x-fastest order."""
    nsteps = [4, 5, 6]
    count = 3 * 4 * 5
    values = np.arange(count, dtype=float) * 0.5
    grid = geometry.StructuredGrid(
        origin=[1, 2, 3], step_vector=[1, 2, 3], nsteps=nsteps, cell_properties={"c": values}
    )
    result = from_loopstructural(grid)
    assert isinstance(result, BlockModel)
    assert list(result.cell_properties) == ["c"]
    np.testing.assert_array_equal(result.cell_properties["c"], values)


def test_structured_grid_without_properties_gives_empty_mappings() -> None:
    """A grid with properties=None and cell_properties=None gives empty mappings."""
    grid = geometry.StructuredGrid(
        origin=[1, 2, 3],
        step_vector=[1, 2, 3],
        nsteps=[4, 5, 6],
        properties=None,
        cell_properties=None,
    )
    result = from_loopstructural(grid)
    assert isinstance(result, BlockModel)
    assert len(result.properties) == 0
    assert len(result.cell_properties) == 0


def test_grid_node_property_is_sent_in_x_fastest_order(
    fake_viewer: Any, viewer_factory: Callable[..., Any]
) -> None:
    """A property f = x + 10 y + 100 z reaches the viewer without any reorder."""
    origin = np.array([1.0, 2.0, 3.0])
    step = np.array([1.0, 2.0, 3.0])
    nsteps = np.array([4, 5, 6])
    grid = geometry.StructuredGrid(origin=origin, step_vector=step, nsteps=nsteps)
    # Build the node coordinates here, x fastest, without the class under test.
    z, y, x = np.meshgrid(
        origin[2] + step[2] * np.arange(nsteps[2]),
        origin[1] + step[1] * np.arange(nsteps[1]),
        origin[0] + step[0] * np.arange(nsteps[0]),
        indexing="ij",
    )
    expected = x.ravel() * 1.0 + y.ravel() * 10.0 + z.ravel() * 100.0
    np.testing.assert_array_equal(grid.nodes @ [1, 10, 100], expected)
    grid.properties["f"] = expected
    viewer = viewer_factory()
    viewer.send(grid, id="grid")
    assert wait_for(lambda: "grid" in fake_viewer.viewer.objects)
    stored = fake_viewer.viewer.objects["grid"]
    assert stored.header["object_type"] == "BlockModel"
    np.testing.assert_allclose(stored_array(stored, "f"), expected)


def test_other_loop_common_class_is_refused() -> None:
    """A ValuePoints object is not accepted. The error names the accepted inputs."""
    points = geometry.ValuePoints(locations=np.zeros((2, 3)), values=np.zeros(2))
    with pytest.raises(InvalidInputError, match="Accepted inputs:"):
        from_loopstructural(points)


def test_send_other_loop_common_class_is_refused(
    fake_viewer: Any, viewer_factory: Callable[..., Any]
) -> None:
    """Viewer.send gives the same error for a ValuePoints object."""
    points = geometry.ValuePoints(locations=np.zeros((2, 3)), values=np.zeros(2))
    viewer = viewer_factory()
    with pytest.raises(InvalidInputError, match="Accepted inputs:"):
        viewer.send(points, id="pts")
