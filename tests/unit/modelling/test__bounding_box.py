import numpy as np
import pytest
from LoopStructural.geometry import BoundingBox


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
    # origin/maximum are always world coordinates; without local_coordinate=True
    # the local interpolation frame is anchored at zero (no shift).
    assert np.all(np.isclose(bbox.local_origin, np.zeros(3)))

    bbox.fit(locations, local_coordinate=True)
    # origin/maximum stay in world coordinates; only the local interpolation
    # frame's anchor moves to the fitted origin.
    assert np.all(np.isclose(bbox.origin, expected_origin))
    assert np.all(np.isclose(bbox.maximum, expected_maximum))
    assert np.all(np.isclose(bbox.local_origin, expected_origin))
    assert np.all(np.isclose(bbox.project(expected_origin), np.zeros(3)))


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


def test_origin_and_maximum_are_world_coordinates():
    origin = np.array([10, 20, 30])
    maximum = np.array([11, 21, 31])
    nsteps = np.array([10, 10, 10])
    step_vector = (maximum - origin) / nsteps

    bbox = BoundingBox(origin=origin, maximum=maximum, nsteps=nsteps, step_vector=step_vector)
    # No automatic local shift -- origin/maximum are always world coordinates.
    assert np.all(np.isclose(bbox.origin, origin))
    assert np.all(np.isclose(bbox.maximum, maximum))
    assert np.all(np.isclose(bbox.local_origin, np.zeros(3)))

    # Setting a local transform anchors project()/reproject() at that origin,
    # without changing origin/maximum themselves.
    bbox.set_local_transform(local_origin=origin)
    assert np.all(np.isclose(bbox.origin, origin))
    assert np.all(np.isclose(bbox.project(origin), np.zeros(3)))
    assert np.all(np.isclose(bbox.reproject(np.zeros(3)), origin))


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
