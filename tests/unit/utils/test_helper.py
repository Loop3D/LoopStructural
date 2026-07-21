import numpy as np
import pandas as pd

from LoopStructural.geometry import BoundingBox
from LoopStructural.utils.helper import (
    get_data_bounding_box,
    get_data_bounding_box_map,
    create_surface,
    create_box,
    xyz_names,
    normal_vec_names,
    tangent_vec_names,
    gradient_vec_names,
    weight_name,
    val_name,
    coord_name,
    interface_name,
    inequality_name,
    feature_name,
    polarity_name,
    pairs_name,
    all_heading,
    empty_dataframe,
)


def _cube_points():
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [0.5, 0.5, 0.5],
        ]
    )


def test_get_data_bounding_box_buffer_scaled_by_extent():
    xyz = _cube_points()
    bb, region = get_data_bounding_box(xyz, 0.1)
    # length of the cube is 1 in each direction, buffer is 10% of that
    expected = np.array([[-0.1, -0.1, -0.1], [1.1, 1.1, 1.1]])
    assert np.allclose(bb, expected)
    # all of the original points should be inside the buffered region
    assert np.all(region(xyz))
    # a point outside the buffered box should be excluded
    outside = np.array([[-1.0, -1.0, -1.0]])
    assert not np.any(region(outside))


def test_get_data_bounding_box_region_checks_all_axes():
    xyz = _cube_points()
    bb, region = get_data_bounding_box(xyz, 0.0)
    # z just above the box should be excluded because get_data_bounding_box
    # applies the mask on all three axes
    outside_z = np.array([[0.5, 0.5, 2.0]])
    assert not np.any(region(outside_z))


def test_get_data_bounding_box_map_absolute_buffer():
    xyz = _cube_points()
    bb, region = get_data_bounding_box_map(xyz, 0.5)
    # get_data_bounding_box_map uses an absolute buffer (not scaled by extent)
    expected = np.array([[-0.5, -0.5, -0.5], [1.5, 1.5, 1.5]])
    assert np.allclose(bb, expected)
    assert np.all(region(xyz))


def test_get_data_bounding_box_map_region_ignores_z():
    xyz = _cube_points()
    # buffer of 0 means the region mask boundary sits exactly on the data extent
    bb, region = get_data_bounding_box_map(xyz, 0.0)
    # region() from get_data_bounding_box_map only thresholds x and y, not z
    # so a point far outside in z but within x/y bounds is still "inside"
    far_z_but_within_xy = np.array([[0.5, 0.5, 100.0]])
    assert np.all(region(far_z_but_within_xy))


def test_create_surface_grid_shapes():
    bounding_box = np.array([[0.0, 0.0], [1.0, 1.0]])
    tri, xx, yy = create_surface(bounding_box, [3, 3])
    # 3x3 grid of points
    assert xx.shape == (9,)
    assert yy.shape == (9,)
    assert np.isclose(xx.min(), 0.0)
    assert np.isclose(xx.max(), 1.0)
    assert np.isclose(yy.min(), 0.0)
    assert np.isclose(yy.max(), 1.0)
    # 2 * (nstep0 - 1) * (nstep1 - 1) triangles
    assert tri.shape == (8, 3)
    # triangle indices must be valid indices into the point arrays
    assert tri.max() < xx.shape[0]
    assert tri.min() >= 0


def test_create_box_returns_closed_hexahedral_mesh():
    bbox = BoundingBox(origin=[0, 0, 0], maximum=[1, 1, 1])
    points, tri = create_box(bbox, np.array([3, 3, 3]))
    assert points.shape[1] == 3
    # 6 faces each built from a 3x3 grid => 6 * 9 points
    assert points.shape[0] == 6 * 9
    # triangle indices should reference valid points
    assert tri.max() < points.shape[0]
    # points should be bound within the (unbuffered) bounding box
    assert np.all(points[:, 0] >= bbox.origin[0] - 1e-9)
    assert np.all(points[:, 0] <= bbox.maximum[0] + 1e-9)
    assert np.all(points[:, 2] >= bbox.origin[2] - 1e-9)
    assert np.all(points[:, 2] <= bbox.maximum[2] + 1e-9)


def test_name_helper_functions():
    assert xyz_names() == ["X", "Y", "Z"]
    assert normal_vec_names() == ["nx", "ny", "nz"]
    assert tangent_vec_names() == ["tx", "ty", "tz"]
    assert gradient_vec_names() == ["gx", "gy", "gz"]
    assert weight_name() == ["w"]
    assert val_name() == ["val"]
    assert coord_name() == ["coord"]
    assert interface_name() == ["interface"]
    assert inequality_name() == ["l", "u"]
    assert feature_name() == ["feature_name"]
    assert polarity_name() == ["polarity"]
    assert pairs_name() == ["pair_id"]


def test_all_heading_concatenates_all_names():
    heading = all_heading()
    expected = (
        xyz_names()
        + normal_vec_names()
        + tangent_vec_names()
        + gradient_vec_names()
        + weight_name()
        + val_name()
        + coord_name()
        + feature_name()
        + interface_name()
        + polarity_name()
        + inequality_name()
        + pairs_name()
    )
    assert heading == expected
    # every expected column name should be present exactly once
    assert len(heading) == len(set(heading))


def test_empty_dataframe_has_expected_number_of_columns():
    df = empty_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 0
    # NOTE: empty_dataframe() constructs the DataFrame with
    # `columns=[all_heading()]`, i.e. a *list containing one list*, rather
    # than `columns=all_heading()`. Pandas therefore builds a MultiIndex of
    # 1-tuples instead of a flat Index of plain column-name strings. This
    # looks like a bug: df["X"] does not return a Series as one would expect
    # for a normal dataframe with an "X" column, it returns a DataFrame
    # (partial MultiIndex selection).
    assert df.shape[1] == len(all_heading())
    assert isinstance(df.columns, pd.MultiIndex)
    flat_names = [c[0] for c in df.columns]
    assert flat_names == all_heading()
