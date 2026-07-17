import numpy as np
import pytest

from LoopStructural.utils.regions import (
    RegionEverywhere,
    RegionFunction,
    PositiveRegion,
    NegativeRegion,
)


class PlaneFeature:
    """A simple mock GeologicalFeature: a scalar field equal to the x
    coordinate, everywhere defined (no NaNs)."""

    def evaluate_value(self, xyz):
        xyz = np.asarray(xyz)
        return xyz[:, 0].astype(float)

    def evaluate_gradient(self, xyz):
        xyz = np.asarray(xyz)
        g = np.zeros((xyz.shape[0], 3))
        g[:, 0] = 1
        return g


class PartiallyUndefinedPlaneFeature:
    """A mock feature whose scalar field is NaN outside of |x| <= 5, to
    exercise the distance-based fallback branch in BaseSignRegion."""

    def evaluate_value(self, xyz):
        xyz = np.asarray(xyz)
        v = xyz[:, 0].astype(float).copy()
        v[np.abs(xyz[:, 0]) > 5] = np.nan
        return v

    def evaluate_gradient(self, xyz):
        xyz = np.asarray(xyz)
        g = np.zeros((xyz.shape[0], 3))
        g[:, 0] = 1
        return g


class AllPositiveFeature:
    """A mock feature whose scalar field never goes negative, used to
    trigger the "cannot find point on surface" error path."""

    def evaluate_value(self, xyz):
        xyz = np.asarray(xyz)
        return np.ones(xyz.shape[0])

    def evaluate_gradient(self, xyz):
        xyz = np.asarray(xyz)
        return np.tile([1.0, 0.0, 0.0], (xyz.shape[0], 1))


def test_region_everywhere_cannot_be_constructed_bug():
    """Documents a real bug: RegionEverywhere.__init__ calls
    `super().__init__()` with no arguments, but BaseRegion.__init__ requires
    a `feature` positional argument. As written, RegionEverywhere() always
    raises TypeError and the class cannot actually be used, despite being
    part of the public `LoopStructural.utils` API
    (`from .regions import RegionEverywhere`).
    """
    with pytest.raises(TypeError):
        RegionEverywhere()


def test_region_function_cannot_be_constructed_bug():
    """Documents the same bug as test_region_everywhere_cannot_be_constructed_bug
    for RegionFunction: it also calls `super().__init__()` with no arguments
    so it always raises TypeError, regardless of the function passed in.
    """
    with pytest.raises(TypeError):
        RegionFunction(lambda xyz: xyz[:, 0] > 0)


def test_positive_region_matches_sign_of_scalar_field():
    feature = PlaneFeature()
    region = PositiveRegion(feature)
    xyz = np.array([[-1.0, 0, 0], [1.0, 0, 0], [0.5, 0, 0], [-0.5, 0, 0]])

    result = region(xyz)

    assert result.dtype == bool
    assert np.array_equal(result, xyz[:, 0] > 0)


def test_negative_region_matches_sign_of_scalar_field():
    feature = PlaneFeature()
    region = NegativeRegion(feature)
    xyz = np.array([[-1.0, 0, 0], [1.0, 0, 0], [0.5, 0, 0], [-0.5, 0, 0]])

    result = region(xyz)

    assert np.array_equal(result, xyz[:, 0] < 0)


def test_positive_region_caches_point_and_vector():
    feature = PlaneFeature()
    region = PositiveRegion(feature)
    assert region.point is None
    assert region.vector is None

    xyz = np.array([[-1.0, 0, 0], [1.0, 0, 0]])
    region(xyz)

    # point/vector should now be cached on the region so subsequent calls
    # don't need to re-derive them
    assert region.point is not None
    assert region.vector is not None
    assert np.allclose(region.vector, [1.0, 0.0, 0.0])


def test_positive_region_uses_precomputed_val():
    feature = PlaneFeature()
    region = PositiveRegion(feature, vector=np.array([1.0, 0, 0]), point=np.array([0.0, 0, 0]))
    xyz = np.array([[1.0, 0, 0], [-1.0, 0, 0]])
    precomputed = np.array([5.0, -5.0])

    result = region(xyz, precomputed_val=precomputed)

    assert np.array_equal(result, precomputed > 0)


def test_region_raises_when_no_point_below_zero_found():
    feature = AllPositiveFeature()
    region = PositiveRegion(feature)
    xyz = np.array([[1.0, 0, 0], [2.0, 0, 0]])

    with pytest.raises(ValueError, match="Cannot find point on surface"):
        region(xyz)


def test_region_falls_back_to_distance_for_nan_values():
    feature = PartiallyUndefinedPlaneFeature()
    region = PositiveRegion(feature, vector=np.array([1.0, 0.0, 0.0]), point=np.array([0.0, 0.0, 0.0]))
    # x = 10 and x = -10 are outside of the feature's support (NaN), so the
    # region must fall back to using signed distance from the cached point
    # along the cached vector to decide in/out. NOTE: the distance is
    # computed as `(centre - xyz) . vector`, i.e. the *opposite* sign
    # convention to the in-support `val > 0` test (which would classify
    # x=10 as "positive" since val=x there). This looks like a possible
    # sign inconsistency between the two branches, but this test documents
    # the actual current behaviour rather than the possibly-intended one.
    xyz = np.array([[10.0, 0, 0], [-10.0, 0, 0], [1.0, 0, 0], [-1.0, 0, 0]])

    result = region(xyz)

    assert np.array_equal(result, np.array([False, True, True, False]))


def test_negative_region_falls_back_to_distance_for_nan_values():
    feature = PartiallyUndefinedPlaneFeature()
    region = NegativeRegion(feature, vector=np.array([1.0, 0.0, 0.0]), point=np.array([0.0, 0.0, 0.0]))
    xyz = np.array([[10.0, 0, 0], [-10.0, 0, 0], [1.0, 0, 0], [-1.0, 0, 0]])

    result = region(xyz)

    assert np.array_equal(result, np.array([True, False, False, True]))


def test_positive_and_negative_regions_are_complementary_away_from_zero():
    feature = PlaneFeature()
    xyz = np.array([[1.0, 0, 0], [-1.0, 0, 0], [2.5, 0, 0], [-2.5, 0, 0]])

    positive = PositiveRegion(feature)(xyz)
    negative = NegativeRegion(feature)(xyz)

    assert np.array_equal(positive, ~negative)
