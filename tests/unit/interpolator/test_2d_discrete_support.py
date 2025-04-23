import numpy as np

from LoopStructural.interpolators import StructuredGrid2D


## structured grid 2d tests
def test_create_structured_grid2d():
    grid = StructuredGrid2D()
    assert isinstance(grid, StructuredGrid2D)


def test_create_structured_grid2d_origin_nsteps():
    grid = StructuredGrid2D(origin=np.zeros(2), nsteps=np.array([5, 5]))
    assert grid.n_nodes == 5 * 5
    assert np.sum(grid.maximum - np.ones(2) * 5) == 0


def test_create_structured_grid2d_origin_nsteps_sv():
    grid = StructuredGrid2D(
        origin=np.zeros(2), nsteps=np.array([10, 10]), step_vector=np.array([0.1, 0.1])
    )
    assert np.sum(grid.step_vector - np.array([0.1, 0.1])) == 0
    assert np.sum(grid.maximum - np.ones(2)) == 0


def test_evaluate_value_2d():
    grid = StructuredGrid2D()
    # grid.update_property("X", grid.nodes[:, 0])
    assert (
        np.sum(grid.barycentre[:, 0] - grid.evaluate_value(grid.barycentre, grid.nodes[:, 0])) == 0
    )


def test_evaluate_gradient_2d():
    grid = StructuredGrid2D()
    # grid.update_property("Y", )
    vector = np.mean(grid.evaluate_gradient(grid.barycentre, grid.nodes[:, 1]), axis=0)
    # vector/=np.linalg.norm(vector)
    assert np.sum(vector - np.array([0, grid.step_vector[1]])) == 0


def test_get_element_2d():
    grid = StructuredGrid2D()
    point = grid.barycentre[[0], :]
    idc, inside = grid.position_to_cell_corners(point)
    bary = np.mean(grid.nodes[idc, :], axis=0)
    assert np.sum(point - bary) == 0


def test_global_to_local_coordinates2d():
    grid = StructuredGrid2D()
    point = np.array([[1.2, 1.5, 1.7]])
    local_coords = grid.position_to_local_coordinates(point)
    assert np.isclose(local_coords[0, 0], 0.2)
    assert np.isclose(local_coords[0, 1], 0.5)


def test_get_element_outside2d():
    grid = StructuredGrid2D()
    point = np.array([grid.origin - np.ones(2)])
    idc, inside = grid.position_to_cell_corners(point)
    assert not inside[0]
