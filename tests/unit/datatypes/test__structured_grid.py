import numpy as np
import pytest
from LoopStructural.datatypes._structured_grid import StructuredGrid
from LoopStructural.utils import rng


def test_structured_grid_to_dict():
    origin = np.array([0, 0, 0])
    nsteps = np.array([10, 10, 10])
    step_vector = np.array([1, 1, 1])
    data = {'rng': rng.random(size=(10, 10, 10))}
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
    data = {'rng': rng.random(size=(10, 10, 10))}
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
    assert np.array_equal(vtk_grid['rng'], data['rng'].flatten(order="F"))


if __name__ == "__main__":
    test_structured_grid_to_dict()
    test_structured_grid_maximum()
    test_structured_grid_vtk()
