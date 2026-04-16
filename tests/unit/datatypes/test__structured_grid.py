import numpy as np
import pytest
from LoopStructural.datatypes._structured_grid import StructuredGrid
from LoopStructural.utils import rng


def test_structured_grid_to_dict():
    origin = np.array([0, 0, 0])
    nsteps = np.array([10, 10, 10])
    step_vector = np.array([1, 1, 1])
    data = {"rng": rng.random(size=(10, 10, 10))}
    name = "grid_data"

    grid = StructuredGrid(
        origin=origin,
        step_vector=step_vector,
        nsteps=nsteps,
        properties=data,
        cell_properties={},
        name=name,
    )
    grid_dict = grid.to_dict()

    assert np.all(grid_dict["origin"] == origin)
    assert np.all(grid_dict["nsteps"] == nsteps)
    assert np.all(grid_dict["step_vector"] == step_vector)
    assert np.all(grid_dict["properties"] == data)
    assert grid_dict["name"] == name


def test_structured_grid_maximum():
    origin = np.array([0, 0, 0])
    nsteps = np.array([10, 10, 10])
    step_vector = np.array([1, 1, 1])

    grid = StructuredGrid(
        origin=origin,
        step_vector=step_vector,
        nsteps=nsteps,
        cell_properties={},
        properties={},
        name=None,
    )
    maximum = grid.maximum

    expected_maximum = origin + nsteps * step_vector
    assert np.array_equal(maximum, expected_maximum)


def test_structured_grid_vtk():
    try:
        import pyvista as pv
    except ImportError:
        pytest.skip("pyvista is required for vtk support")
    origin = np.array([0, 0, 0])
    nsteps = np.array([10, 10, 10])
    step_vector = np.array([1, 1, 1])
    data = {"rng": rng.random(size=(10, 10, 10))}
    name = "grid_data"

    grid = StructuredGrid(
        origin=origin,
        step_vector=step_vector,
        nsteps=nsteps,
        properties=data,
        cell_properties={},
        name=name,
    )
    vtk_grid = grid.vtk()
    grid_points = vtk_grid.points
    grid_origin = np.min(grid_points, axis=0)

    # Add assertions to validate the generated vtk_grid
    assert isinstance(vtk_grid, pv.RectilinearGrid)
    assert np.array_equal(vtk_grid.dimensions, nsteps)
    assert np.array_equal(grid_origin, origin)
    assert np.array_equal(vtk_grid["rng"], data["rng"].flatten(order="F"))


def test_structured_grid_save_gocad_voxet(tmp_path):
    origin = np.array([10.0, 20.0, 30.0])
    nsteps = np.array([2, 3, 4])
    step_vector = np.array([1.0, 2.0, 3.0])
    values = np.arange(np.prod(nsteps), dtype=np.float32).reshape(tuple(nsteps), order="F")
    grid = StructuredGrid(
        origin=origin,
        step_vector=step_vector,
        nsteps=nsteps,
        properties={"scalar": values},
        cell_properties={},
        name="grid_data",
    )

    output = tmp_path / "structured_grid.vo"
    grid.save(output)

    header = output.read_text()
    binary_path = tmp_path / "structured_grid_scalar@@"

    assert "GOCAD Voxet 1" in header
    assert "PROPERTY 1 scalar" in header
    assert f"AXIS_N {nsteps[0]} {nsteps[1]} {nsteps[2]}" in header
    assert "PROP_FILE 1 structured_grid_scalar@@" in header
    assert binary_path.exists()

    exported = np.fromfile(binary_path, dtype=">f4")
    np.testing.assert_array_equal(exported, values.flatten(order="F"))


def test_structured_grid_save_gocad_voxet_uses_cell_properties_when_needed(tmp_path):
    origin = np.array([0.0, 0.0, 0.0])
    nsteps = np.array([3, 3, 3])
    step_vector = np.array([1.0, 1.0, 1.0])
    values = np.arange(np.prod(nsteps - 1), dtype=np.float32).reshape(tuple(nsteps - 1), order="F")
    grid = StructuredGrid(
        origin=origin,
        step_vector=step_vector,
        nsteps=nsteps,
        properties={},
        cell_properties={"cells": values},
        name="cell_grid",
    )

    output = tmp_path / "cell_grid.vo"
    grid.save(output)

    header = output.read_text()
    binary_path = tmp_path / "cell_grid_cells@@"

    assert "PROPERTY 1 cells" in header
    assert "AXIS_N 2 2 2" in header
    assert binary_path.exists()

    exported = np.fromfile(binary_path, dtype=">f4")
    np.testing.assert_array_equal(exported, values.flatten(order="F"))


if __name__ == "__main__":
    test_structured_grid_to_dict()
    test_structured_grid_maximum()
    test_structured_grid_vtk()
