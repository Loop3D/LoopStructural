from LoopStructural.interpolators import StructuredGrid
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


def test_inside(support):
    assert np.all(support.inside(support.barycentre))


def test_evaluate_value(support):
    # print()
    # print(support.evaluate_value(support.barycentre, support.nodes[:, 0]))
    # print(support.barycentre[:, 0])
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

    cell_indexes, inside = support.position_to_cell_index(support.barycentre - 5)
    assert np.all(cell_indexes[inside, 0] < support.nsteps_cells[0])
    assert np.all(cell_indexes[inside, 1] < support.nsteps_cells[1])
    assert np.all(cell_indexes[inside, 2] < support.nsteps_cells[2])
    corners = support.cell_corner_indexes(cell_indexes)
    assert np.all(corners[inside, 0] < support.nsteps[0])
    assert np.all(corners[inside, 1] < support.nsteps[1])
    assert np.all(corners[inside, 2] < support.nsteps[2])
    globalidx = support.global_node_indices(corners)
    # print(globalidx[inside],grid.n_nodes,inside)
    assert np.all(globalidx[inside] < support.n_nodes)
    inside = support.inside(support.barycentre - 5)
    # inside, support.position_to_cell_corne    rs(support.barycentre - 5)
    vector = support.evaluate_gradient(support.barycentre - 5, support.nodes[:, 1])
    assert np.sum(np.mean(vector[inside, :], axis=0) - np.array([0, 1, 0])) == 0
    vector = support.evaluate_gradient(support.nodes, support.nodes[:, 1])


def test_evaluate_gradient2(support_class):
    # this test is the same as above but we will use a random vector
    rng = np.random.default_rng(10)
    _i = 0
    for _i in range(10):
        step = rng.uniform(0, 100)
        grid = support_class(step_vector=np.array([step, step, step]))

        # define random vector
        n = rng.random(3)
        n /= np.linalg.norm(n)
        distance = n[0] * grid.nodes[:, 0] + n[1] * grid.nodes[:, 1] + n[2] * grid.nodes[:, 2]
        vector = grid.evaluate_gradient(rng.uniform(1, 8, size=(100, 3)), distance)
        assert np.all(np.isclose(np.sum(vector - n[None, :], axis=1), 0, atol=1e-3, rtol=1e-3))
    assert _i == 9


def test_get_element(support):
    point = support.barycentre[[0], :]
    # point[0, 0] += 0.1
    vertices, dof, idc, inside = support.get_element_for_location(point)
    # vertices = vertices.reshape(-1, 3)
    bary = np.mean(vertices, axis=1)
    assert np.isclose(np.sum(point - bary), 0)


def test_global_to_local_coordinates():
    grid = StructuredGrid()
    point = np.array([[1.2, 1.5, 1.7]])
    local_coords = grid.position_to_local_coordinates(point)
    assert np.isclose(local_coords[0, 0], 0.2)
    assert np.isclose(local_coords[0, 1], 0.5)
    assert np.isclose(local_coords[0, 2], 0.7)


def test_get_element_outside(support):
    point = np.array([support.origin - np.ones(3)])
    idc, inside = support.position_to_cell_corners(point)
    assert not inside[0]


def test_node_index_to_position(support):
    assert (
        np.sum(support.node_indexes_to_position(np.array([[0, 0, 0]])) - np.array([0, 0, 0])) == 0
    )
    for i in range(10):
        for j in range(10):
            for k in range(10):
                assert (
                    np.sum(
                        support.node_indexes_to_position(np.array([[i, j, k]]))
                        - np.array([i, j, k]) * support.step_vector
                    )
                    == 0
                )
    assert np.sum(support.node_indexes_to_position(np.array([0, 0, 0])) - np.array([0, 0, 0])) == 0


def test_global_index_to_cell_index(support):
    assert np.sum(support.global_index_to_cell_index(np.array([0])) - np.array([0, 0, 0])) == 0


def test_global_index(support):
    indexes = np.array(
        np.meshgrid(
            np.arange(0, support.nsteps[0]),
            np.arange(0, support.nsteps[1]),
            np.arange(0, support.nsteps[2]),
        )
    ).reshape(-1, 3)
    global_node_index = support.global_node_indices(indexes)
    assert np.all(global_node_index >= 0)
    assert np.all(global_node_index < support.n_nodes)

    indexes = np.array(
        np.meshgrid(
            np.arange(0, 3),
            np.arange(0, 1),
            np.arange(0, 1),
        )
    ).reshape(-1, 3)
    global_node_index = support.global_node_indices(indexes)
    assert np.all(global_node_index >= 0)
    assert np.all(global_node_index < support.n_nodes)
