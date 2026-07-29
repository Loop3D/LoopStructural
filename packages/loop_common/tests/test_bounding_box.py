import numpy as np
import pytest

from loop_common.geometry import BoundingBox


def test_origin_and_maximum_are_world_coordinates():
    bbox = BoundingBox(origin=[10.0, 10.0, 10.0], maximum=[20.0, 20.0, 20.0])

    assert np.allclose(bbox.origin, [10.0, 10.0, 10.0])
    assert np.allclose(bbox.maximum, [20.0, 20.0, 20.0])


def test_default_projection_is_identity():
    bbox = BoundingBox(origin=[1.0, 2.0, 3.0], maximum=[4.0, 5.0, 6.0])
    pts = np.array([[1.5, 2.5, 3.5], [3.5, 4.5, 5.5]])

    local = bbox.project(pts)
    assert np.allclose(local, pts)
    assert np.allclose(bbox.reproject(local), pts)


def test_fit_local_coordinate_keeps_world_bounds_and_sets_local_origin():
    locations = np.array([[2.0, 3.0, 4.0], [5.0, 7.0, 11.0]])
    bbox = BoundingBox()
    bbox.fit(locations, local_coordinate=True)

    assert np.allclose(bbox.origin, [2.0, 3.0, 4.0])
    assert np.allclose(bbox.maximum, [5.0, 7.0, 11.0])
    assert np.allclose(bbox.local_origin, [2.0, 3.0, 4.0])

    local = bbox.project(locations)
    assert np.allclose(local.min(axis=0), [0.0, 0.0, 0.0])
    assert np.allclose(local.max(axis=0), [3.0, 4.0, 7.0])
    assert np.allclose(bbox.reproject(local), locations)


def test_rotation_transform_applies_to_points_and_vectors():
    bbox = BoundingBox(origin=[0.0, 0.0, 0.0], maximum=[10.0, 10.0, 10.0])
    rotation = np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    bbox.set_local_transform(local_origin=[0.0, 0.0, 0.0], rotation_matrix=rotation)

    point = np.array([1.0, 0.0, 0.0])
    vector = np.array([1.0, 0.0, 0.0])

    assert np.allclose(bbox.project(point), [0.0, 1.0, 0.0])
    assert np.allclose(bbox.project_vectors(vector), [0.0, 1.0, 0.0])
    assert np.allclose(bbox.reproject(bbox.project(point)), point)
    assert np.allclose(bbox.reproject_vectors(bbox.project_vectors(vector)), vector)


def test_matrix_matches_world_to_local_transform():
    bbox = BoundingBox(origin=[10.0, 0.0, 0.0], maximum=[20.0, 5.0, 2.0])
    bbox.set_local_transform(local_origin=[10.0, 0.0, 0.0])

    assert np.allclose(bbox.matrix(normalise=False), bbox.world_to_local_matrix)
    assert np.allclose(bbox.project([12.0, 1.0, 1.0]), [2.0, 1.0, 1.0])


def test_legacy_global_arguments_are_rejected():
    with pytest.raises(TypeError):
        BoundingBox(global_origin=[10.0, 10.0, 10.0], global_maximum=[20.0, 20.0, 20.0])
