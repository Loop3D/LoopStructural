import numpy as np
import pytest
from LoopStructural.datatypes._structured_grid import StructuredGrid


def test_structured_grid_to_dict():
    origin = np.array([0, 0, 0])
    nsteps = np.array([10, 10, 10])
    step_vector = np.array([1, 1, 1])
    data = np.random.rand(10, 10, 10)
    name = "grid_data"

    grid = StructuredGrid(origin, nsteps, step_vector, data, name)
    grid_dict = grid.to_dict()

    assert np.all(grid_dict["origin"] == origin)
    assert np.all(grid_dict["nsteps"] == nsteps)
    assert np.all(grid_dict["step_vector"] == step_vector)
    assert np.all(grid_dict["data"] == data)
    assert grid_dict["name"] == name


def test_structured_grid_maximum():
    origin = np.array([0, 0, 0])
    nsteps = np.array([10, 10, 10])
    step_vector = np.array([1, 1, 1])

    grid = StructuredGrid(origin, nsteps, step_vector, None, None)
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
    data = np.random.rand(10, 10, 10)
    name = "grid_data"

    grid = StructuredGrid(origin, nsteps, step_vector, data, name)
    vtk_grid = grid.vtk

    # Add assertions to validate the generated vtk_grid
    assert isinstance(vtk_grid, pv.RectilinearGrid)
    assert np.array_equal(vtk_grid.dimensions, nsteps)
    assert np.array_equal(vtk_grid.origin, origin)
    assert np.array_equal(vtk_grid.spacing, step_vector)
    assert np.array_equal(vtk_grid[name], data.flatten(order="F"))
