from LoopStructural.interpolators import TetMesh
from LoopStructural.interpolators import StructuredGrid, StructuredGrid2D
import numpy as np

## structured grid tests
def test_create_support(support):
    """
    support is a fixture that returns a support object created
    with the default constructor. Ensure that it is not none and
    make sure the origin and maximum are correct.
    """
    assert support is not None
    assert np.sum(support.origin - np.zeros(3)) == 0
    assert np.sum(support.maximum - np.ones(3) * 10) == 0


def test_create_support_origin_nsteps(support_class):
    grid = support_class(
        origin=np.zeros(3),
        nsteps=np.array([10, 10, 10]),
        step_vector=np.array([0.1, 0.1, 0.1]),
    )
    assert np.sum(grid.step_vector - np.array([0.1, 0.1, 0.1])) == 0
    assert np.sum(grid.maximum - np.ones(3)) == 0


def test_evaluate_value(support):
    assert (
        np.sum(
            support.barycentre[:, 0]
            - support.evaluate_value(support.barycentre, support.nodes[:, 0])
        )
        == 0
    )


def test_evaluate_gradient(support_class):
    support = support_class()
    # test by setting the scalar field to the y coordinate
    vector = support.evaluate_gradient(support.barycentre, support.nodes[:, 1])
    assert np.sum(vector - np.array([0, 1, 0])) == 0

    # same test but for a bigger grid, making sure scaling for cell is ok
    support = support_class(step_vector=np.array([100, 100, 100]))
    vector = support.evaluate_gradient(support.barycentre, support.nodes[:, 1])
    assert np.sum(vector - np.array([0, 1, 0])) == 0


def test_outside_box(support):
    # test by setting the scalar field to the y coordinate
    inside = support.inside(support.barycentre + 5)
    assert np.all(~inside == np.any((support.barycentre + 5) > support.maximum, axis=1))
    inside = support.inside(support.barycentre - 5)
    assert np.all(~inside == np.any((support.barycentre - 5) < support.origin, axis=1))

    cix, ciy, ciz = support.position_to_cell_index(support.barycentre - 5)
    assert np.all(cix[inside] < support.nsteps_cells[0])
    assert np.all(ciy[inside] < support.nsteps_cells[1])
    assert np.all(ciz[inside] < support.nsteps_cells[2])
    cornersx, cornersy, cornersz = support.cell_corner_indexes(cix, ciy, ciz)
    assert np.all(cornersx[inside] < support.nsteps[0])
    assert np.all(cornersy[inside] < support.nsteps[1])
    assert np.all(cornersz[inside] < support.nsteps[2])
    globalidx = support.global_node_indicies(
        np.dstack([cornersx, cornersy, cornersz]).T
    )
    # print(globalidx[inside],grid.n_nodes,inside)
    assert np.all(globalidx[inside] < support.n_nodes)
    inside = support.inside(support.barycentre - 5)
    # inside, support.position_to_cell_corne    rs(support.barycentre - 5)
    vector = support.evaluate_gradient(support.barycentre - 5, support.nodes[:, 1])
    assert np.sum(np.mean(vector[inside, :], axis=0) - np.array([0, 1, 0])) == 0
    vector = support.evaluate_gradient(support.nodes, support.nodes[:, 1])


def test_evaluate_gradient2(support_class):
    # this test is the same as above but we will use a random vector
    np.random.seed(0)
    for i in range(10):
        step = np.random.uniform(0, 100)
        grid = support_class(step_vector=np.array([step, step, step]))

        # define random vector
        n = np.random.random(3)
        n /= np.linalg.norm(n)
        distance = (
            n[0] * grid.nodes[:, 0] + n[1] * grid.nodes[:, 1] + n[2] * grid.nodes[:, 2]
        )
        vector = grid.evaluate_gradient(
            np.random.uniform(1, 8, size=(100, 3)), distance
        )
        assert np.all(np.isclose(np.sum(vector - n[None, :], axis=1), 0)) == True
    assert i == 9


def test_get_element(support):
    point = support.barycentre[[0], :]
    vertices, dof, idc, inside = support.get_element_for_location(point)
    # vertices = vertices.reshape(-1, 3)
    bary = np.mean(vertices, axis=1)
    assert np.isclose(np.sum(point - bary), 0)


def test_global_to_local_coordinates():
    grid = StructuredGrid()
    point = np.array([[1.2, 1.5, 1.7]])
    lx, ly, lz = grid.position_to_local_coordinates(point)
    assert np.isclose(lx[0], 0.2)
    assert np.isclose(ly[0], 0.5)
    assert np.isclose(lz[0], 0.7)


def test_get_element_outside(support):
    point = np.array([support.origin - np.ones(3)])
    idc, inside = support.position_to_cell_corners(point)
    assert inside[0] == False


def test_node_index_to_position(support):
    assert np.sum(support.node_indexes_to_position(0, 0, 0) - np.array([0, 0, 0])) == 0
    for i in range(10):
        for j in range(10):
            for k in range(10):
                assert (
                    np.sum(
                        support.node_indexes_to_position(i, j, k)
                        - np.array([i, j, k]) * support.step_vector
                    )
                    == 0
                )
    # assert np.sum(support.node_indexes_to_position(0, 0, 0) - np.array([0, 0, 0])) == 0


def test_global_index_to_cell_index(support):
    assert np.sum(support.global_index_to_cell_index(0) - np.array([0, 0, 0])) == 0


def test_get_elements(support):
    pass
