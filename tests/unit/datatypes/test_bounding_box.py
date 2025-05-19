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
    assert bbox.dimensions == 3
    assert np.all(bbox.bb == np.array([[0, 0, 0], [1, 1, 1]]))
    assert bbox.valid is True


def test_create_bounding_box_2d():
    bbox = BoundingBox(origin=[0, 0], maximum=[1, 1], dimensions=2)
    assert bbox.origin[0] == 0
    assert bbox.origin[1] == 0
    assert bbox.maximum[0] == 1
    assert bbox.maximum[1] == 1
    assert bbox.dimensions == 2
    assert np.all(bbox.bb == np.array([[0, 0], [1, 1]]))
    assert bbox.valid is True


def test_create_bounding_box_from_points():
    bbox = BoundingBox(dimensions=3)
    bbox.fit(np.array([[0, 0, 0], [1, 1, 1]]))
    assert bbox.origin[0] == 0
    assert bbox.origin[1] == 0
    assert bbox.origin[2] == 0
    assert bbox.maximum[0] == 1
    assert bbox.maximum[1] == 1
    assert bbox.maximum[2] == 1
    assert np.all(bbox.bb == np.array([[0, 0, 0], [1, 1, 1]]))
    assert bbox.valid is True


def test_create_bounding_box_from_points_2d():
    bbox = BoundingBox(dimensions=2)
    bbox.fit(np.array([[0, 0], [1, 1]]))
    assert bbox.origin[0] == 0
    assert bbox.origin[1] == 0
    assert bbox.maximum[0] == 1
    assert bbox.maximum[1] == 1
    assert np.all(bbox.bb == np.array([[0, 0], [1, 1]]))
    assert bbox.valid is True


def test_create_3d_bounding_box_from_2d_points():
    bbox = BoundingBox(dimensions=3)
    try:
        bbox.fit(np.array([[0, 0], [1, 1]]))
    except Exception as e:
        assert str(e) == "locations array is 2D but bounding box is 3"
    else:
        assert False


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


def test_regular_grid_3d():
    bbox = BoundingBox(origin=[0, 0, 0], maximum=[1, 1, 1])
    print(bbox.dimensions)
    grid = bbox.regular_grid((10, 10, 10))
    assert grid.shape == (10 * 10 * 10, 3)


def test_regular_grid_2d():
    bbox = BoundingBox(origin=[0, 0], maximum=[1, 1], dimensions=2)
    print(bbox.dimensions)
    grid = bbox.regular_grid((10, 10))
    assert grid.shape == (10 * 10, 2)

def test_project_to_local():
    bbox = BoundingBox(global_origin=[10,10,10], global_maximum=[20,20,20])
    point = np.array([15, 15, 15])
    local_point = bbox.project(point)
    assert np.all(local_point == np.array([5, 5, 5]))
if __name__ == "__main__":
    test_create_bounding_box()
    test_create_bounding_box_from_points()
    test_create_with_buffer()
    test_is_inside()
    test_regular_grid_3d()
    test_regular_grid_2d()
