from LoopStructural.datatypes import BoundingBox
import numpy as np


def test_create_bounding_box():
    bbox = BoundingBox(origin=[0, 0, 0], maximum=[1, 1, 1])
    assert bbox.origin[0] == 0
    assert bbox.origin[1] == 0
    assert bbox.origin[2] == 0
    assert bbox.maximum[0] == 1
    assert bbox.maximum[1] == 1
    assert bbox.maximum[2] == 1
    assert np.all(bbox.bb == np.array([[0, 0, 0], [1, 1, 1]]))
    assert bbox.valid is True


def test_create_bounding_box_from_points():
    bbox = BoundingBox()
    bbox.fit(np.array([[0, 0, 0], [1, 1, 1]]))
    assert bbox.origin[0] == 0
    assert bbox.origin[1] == 0
    assert bbox.origin[2] == 0
    assert bbox.maximum[0] == 1
    assert bbox.maximum[1] == 1
    assert bbox.maximum[2] == 1
    assert np.all(bbox.bb == np.array([[0, 0, 0], [1, 1, 1]]))
    assert bbox.valid is True


def test_create_with_buffer():
    bbox = BoundingBox(origin=[0, 0, 0], maximum=[1, 1, 1])
    bbox = bbox.with_buffer(0.2)
    assert bbox.origin[0] == -0.2
    assert bbox.origin[1] == -0.2
    assert bbox.origin[2] == -0.2
    assert bbox.maximum[0] == 1.2
    assert bbox.maximum[1] == 1.2
    assert bbox.maximum[2] == 1.2
    assert np.all(bbox.bb == np.array([[-0.2, -0.2, -0.2], [1.2, 1.2, 1.2]]))
    assert bbox.valid is True


def test_is_inside():
    bbox = BoundingBox(origin=[0, 0, 0], maximum=[1, 1, 1])
    assert np.all(bbox.is_inside(np.array([0.5, 0.5, 0.5])))
    assert not np.all(bbox.is_inside(np.array([0.5, 0.5, 1.5])))
    assert not np.all(bbox.is_inside(np.array([0.5, 0.5, -0.5])))
    assert not np.all(bbox.is_inside(np.array([0.5, 1.5, 0.5])))
    assert not np.all(bbox.is_inside(np.array([0.5, -0.5, 0.5])))
    assert not np.all(bbox.is_inside(np.array([1.5, 0.5, 0.5])))
    assert not np.all(bbox.is_inside(np.array([-0.5, 0.5, 0.5])))


if __name__ == "__main__":
    test_create_bounding_box()
    test_create_bounding_box_from_points()
    test_create_with_buffer()
    test_is_inside()
