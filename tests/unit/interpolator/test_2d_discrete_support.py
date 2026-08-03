import numpy as np
import pytest

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
    idc, _inside = grid.position_to_cell_corners(point)
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
    _idc, inside = grid.position_to_cell_corners(point)
    assert not inside[0]


def test_structured_grid2d_vtk_assigns_quad_cell_types():
    try:
        import pyvista as pv
    except ImportError:
        pytest.skip("pyvista is required for vtk support")
    grid = StructuredGrid2D(origin=np.zeros(2), nsteps=np.array([4, 4]))

    vtk_grid = grid.vtk(z=0.0)

    assert vtk_grid.n_cells == grid.n_elements
    assert vtk_grid.celltypes.shape == (grid.n_elements,)
    assert np.all(vtk_grid.celltypes == pv.CellType.QUAD)
    assert vtk_grid.get_cell(0).point_ids == [0, 1, 5, 4]

    triangulated = vtk_grid.triangulate()

    assert triangulated.n_cells == grid.n_elements * 2


def test_evaluate_gradient_2d_world_units():
    """
    Gradient must be returned in world units, independent of the cell size.

    Regression test for #289: get_element_gradient_for_location returned the
    shape-function derivatives with respect to local (0-1) cell coordinates,
    so gradients were inflated by step_vector on each axis.
    """
    grid = StructuredGrid2D(
        origin=np.zeros(2), nsteps=np.array([10, 10]), step_vector=np.array([2.5, 0.5])
    )
    # f(x, y) = x and f(x, y) = y have unit gradients regardless of cell size
    gradient_x = np.mean(grid.evaluate_gradient(grid.barycentre, grid.nodes[:, 0]), axis=0)
    gradient_y = np.mean(grid.evaluate_gradient(grid.barycentre, grid.nodes[:, 1]), axis=0)
    assert np.allclose(gradient_x, np.array([1.0, 0.0]))
    assert np.allclose(gradient_y, np.array([0.0, 1.0]))


def test_fdi_2d_gradient_resolution_independent():
    """
    The interpolated gradient magnitude must not depend on nelements (#289).
    """
    from LoopStructural.geometry import BoundingBox
    from LoopStructural.interpolators import InterpolatorFactory

    sqrt2 = np.sqrt(2.0)
    bbox = BoundingBox(dimensions=2, origin=np.array([0.0, 0.0]), maximum=np.array([100.0, 100.0]))
    points = np.random.default_rng(0).uniform(10, 90, size=(60, 2))
    for nelements in (1e3, 4e3):
        interpolator = InterpolatorFactory.create_interpolator(
            interpolatortype="FDI", boundingbox=bbox, nelements=nelements
        )
        interpolator.set_value_constraints(
            np.column_stack([points, (points[:, 0] + points[:, 1]) / sqrt2])
        )
        interpolator.setup_interpolator()
        interpolator.solve_system(solver="cg")
        gradient_norm = np.linalg.norm(interpolator.evaluate_gradient(np.array([[50.0, 50.0]]))[0])
        assert np.isclose(gradient_norm, 1.0, atol=0.05)
