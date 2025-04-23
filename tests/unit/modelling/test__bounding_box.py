import numpy as np
import pytest

from LoopStructural.datatypes._bounding_box import BoundingBox


def test_bounding_box_creation():
    origin = np.array([0, 0, 0])
    maximum = np.array([1, 1, 1])
    nsteps = np.array([10, 10, 10])
    step_vector = (maximum - origin) / nsteps

    bbox = BoundingBox(origin=origin, maximum=maximum, nsteps=nsteps, step_vector=step_vector)
    assert np.all(np.isclose(bbox.origin, origin))
    assert np.all(np.isclose(bbox.maximum, maximum))
    assert np.all(np.isclose(bbox.nsteps, nsteps))
    assert np.all(np.isclose(bbox.step_vector, step_vector))


def test_bounding_box_fit():
    locations = np.array([[0.5, 0.5, 0.5], [0.2, 0.3, 0.4], [0.8, 0.9, 1.0]])
    expected_origin = np.array([0.2, 0.3, 0.4])
    expected_maximum = np.array([0.8, 0.9, 1.0])

    bbox = BoundingBox()
    bbox.fit(locations)
    assert np.all(np.isclose(bbox.origin, expected_origin))
    assert np.all(np.isclose(bbox.maximum, expected_maximum))
    assert np.all(np.isclose(bbox.maximum, expected_maximum))
    assert np.all(np.isclose(bbox.global_origin, np.zeros(3)))
    bbox.fit(locations, local_coordinate=True)
    assert np.all(np.isclose(bbox.origin, np.zeros(3)))
    assert np.all(np.isclose(bbox.maximum, expected_maximum - expected_origin))
    assert np.all(np.isclose(bbox.global_origin, expected_origin))


def test_bounding_box_volume():
    origin = np.array([0, 0, 0])
    maximum = np.array([1, 1, 1])
    nsteps = np.array([10, 10, 10])
    step_vector = (maximum - origin) / nsteps

    bbox = BoundingBox(origin=origin, maximum=maximum, nsteps=nsteps, step_vector=step_vector)

    expected_volume = 1.0
    assert bbox.volume == expected_volume


def test_bounding_box_is_inside():
    origin = np.array([0, 0, 0])
    maximum = np.array([1, 1, 1])
    nsteps = np.array([10, 10, 10])
    step_vector = (maximum - origin) / nsteps

    bbox = BoundingBox(origin=origin, maximum=maximum, nsteps=nsteps, step_vector=step_vector)

    inside_points = np.array([[0.5, 0.5, 0.5], [0.2, 0.3, 0.4]])
    outside_points = np.array([[1.5, 1.5, 1.5], [-0.1, -0.2, -0.3]])

    assert np.all(bbox.is_inside(inside_points))
    assert not np.any(bbox.is_inside(outside_points))


def test_local_and_global_origin():
    origin = np.array([0, 0, 0])
    maximum = np.array([1, 1, 1])
    nsteps = np.array([10, 10, 10])
    step_vector = (maximum - origin) / nsteps

    bbox = BoundingBox(origin=origin, maximum=maximum, nsteps=nsteps, step_vector=step_vector)
    assert np.all(np.isclose(bbox.global_origin, origin))
    assert np.all(np.isclose(bbox.global_maximum, maximum))
    assert np.all(np.isclose(bbox.origin, np.zeros(3)))
    assert np.all(np.isclose(bbox.maximum, maximum - origin))


def test_buffer():
    origin = np.array([0, 0, 0])
    maximum = np.array([1, 1, 1])
    nsteps = np.array([10, 10, 10])
    step_vector = (maximum - origin) / nsteps

    bbox = BoundingBox(origin=origin, maximum=maximum, nsteps=nsteps, step_vector=step_vector)
    bbox2 = bbox.with_buffer(0.1)
    assert np.all(np.isclose(bbox2.origin, np.array([-0.1, -0.1, -0.1])))
    assert np.all(np.isclose(bbox2.maximum, np.array([1.1, 1.1, 1.1])))


if __name__ == "__main__":
    pytest.main([__file__])
