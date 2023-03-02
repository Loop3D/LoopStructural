from LoopStructural.interpolators import TetMesh
from LoopStructural.interpolators import StructuredGrid, StructuredGrid2D
import numpy as np





## structured grid 2d tests
def test_create_structured_grid2d():
    grid = StructuredGrid2D()


def test_create_structured_grid2d_origin_nsteps():
    grid = StructuredGrid2D(origin=np.zeros(2), nsteps=np.array([5, 5]))
    assert grid.n_nodes == 5 * 5
    assert np.sum(grid.maximum - np.ones(2) * 5) == 0


def test_create_structured_grid2d_origin_nsteps():
    grid = StructuredGrid2D(
        origin=np.zeros(2), nsteps=np.array([10, 10]), step_vector=np.array([0.1, 0.1])
    )
    assert np.sum(grid.step_vector - np.array([0.1, 0.1])) == 0
    assert np.sum(grid.maximum - np.ones(2)) == 0


def test_evaluate_value_2d():
    grid = StructuredGrid2D()
    grid.update_property("X", grid.nodes[:, 0])
    assert (
        np.sum(grid.barycentre[:, 0] - grid.evaluate_value(grid.barycentre, "X")) == 0
    )


def test_evaluate_gradient_2d():
    grid = StructuredGrid2D()
    grid.update_property("Y", grid.nodes[:, 1])
    vector = np.mean(grid.evaluate_gradient(grid.barycentre, "Y"), axis=0)
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
    lx, ly = grid.position_to_local_coordinates(point)
    assert np.isclose(lx[0], 0.2)
    assert np.isclose(ly[0], 0.5)


def test_get_element_outside2d():
    grid = StructuredGrid2D()
    point = np.array([grid.origin - np.ones(2)])
    idc, inside = grid.position_to_cell_corners(point)
    assert inside[0] == False
