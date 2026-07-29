import numpy as np

from loop_common.supports import StructuredGrid


def test_evaluate_value_on_domain_vertices_no_nan():
    grid = StructuredGrid(origin=np.zeros(3), nsteps=np.array([10, 10, 10]), step_vector=np.ones(3))
    values = grid.nodes[:, 0] + 2.0 * grid.nodes[:, 1] + 3.0 * grid.nodes[:, 2]

    corner_points = np.array(
        [
            [grid.origin[0], grid.origin[1], grid.origin[2]],
            [grid.maximum[0], grid.origin[1], grid.origin[2]],
            [grid.origin[0], grid.maximum[1], grid.origin[2]],
            [grid.origin[0], grid.origin[1], grid.maximum[2]],
            [grid.maximum[0], grid.maximum[1], grid.maximum[2]],
        ]
    )
    expected = corner_points[:, 0] + 2.0 * corner_points[:, 1] + 3.0 * corner_points[:, 2]
    evaluated = grid.evaluate_value(corner_points, values)

    assert np.all(np.isfinite(evaluated))
    assert np.allclose(evaluated, expected)


def test_evaluate_value_on_domain_face_no_nan():
    grid = StructuredGrid(origin=np.zeros(3), nsteps=np.array([10, 10, 10]), step_vector=np.ones(3))
    values = grid.nodes[:, 0] - 0.5 * grid.nodes[:, 1] + 0.25 * grid.nodes[:, 2]

    # Points on the x=max face, including one collocated with a boundary vertex.
    points = np.array(
        [
            [grid.maximum[0], 2.5, 4.5],
            [grid.maximum[0], 0.0, 0.0],
            [grid.maximum[0], grid.maximum[1], 7.25],
        ]
    )
    expected = points[:, 0] - 0.5 * points[:, 1] + 0.25 * points[:, 2]
    evaluated = grid.evaluate_value(points, values)

    assert np.all(np.isfinite(evaluated))
    assert np.allclose(evaluated, expected)
