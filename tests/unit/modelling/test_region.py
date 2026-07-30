import numpy as np

from LoopStructural.modelling.features._region import Region
from LoopStructural.modelling.features._analytical_feature import (
    AnalyticalGeologicalFeature,
)


class PlaneFeature:
    """Minimal mock GeologicalFeature: scalar field equal to x."""

    name = "plane"

    def evaluate_value(self, xyz):
        xyz = np.asarray(xyz)
        return xyz[:, 0].astype(float)


def test_region_positive_sign_selects_positive_values():
    feature = PlaneFeature()
    region = Region(feature, value=0.0, sign=True)

    xyz = np.array([[-1.0, 0, 0], [1.0, 0, 0], [0.0, 0, 0]])
    result = region(xyz)

    assert np.array_equal(result, np.array([False, True, False]))


def test_region_negative_sign_selects_negative_values():
    feature = PlaneFeature()
    region = Region(feature, value=0.0, sign=False)

    xyz = np.array([[-1.0, 0, 0], [1.0, 0, 0], [0.0, 0, 0]])
    result = region(xyz)

    assert np.array_equal(result, np.array([True, False, False]))


def test_region_positive_and_negative_are_complementary_away_from_zero():
    feature = PlaneFeature()
    xyz = np.array([[-2.0, 0, 0], [-0.5, 0, 0], [0.5, 0, 0], [3.0, 0, 0]])

    positive = Region(feature, value=0.0, sign=True)(xyz)
    negative = Region(feature, value=0.0, sign=False)(xyz)

    assert np.array_equal(positive, ~negative)


def test_region_to_json():
    feature = PlaneFeature()
    region = Region(feature, value=0.5, sign=True)

    json = region.to_json()

    assert json == {"feature": "plane", "value": 0.5, "sign": True}


def test_region_call_ignores_value_and_always_thresholds_at_zero():
    """Region.__call__ only ever compares the feature's scalar value against
    zero (`> 0` or `< 0`); `self.value` is stored (and used by `to_json`)
    but is never actually used as the threshold in `__call__`. This test
    documents that behaviour explicitly: even with `value=0.5`, a point
    with scalar value 0.25 (which is < 0.5 but > 0) is still classified as
    "positive".
    """
    feature = PlaneFeature()
    region = Region(feature, value=0.5, sign=True)

    xyz = np.array([[0.25, 0, 0]])
    result = region(xyz)

    assert result[0]


def test_region_stores_constructor_arguments():
    feature = PlaneFeature()
    region = Region(feature, value=1.5, sign=False)

    assert region.feature is feature
    assert region.value == 1.5
    assert region.sign is False


def test_base_feature_regions_not_shared_between_instances():
    """Regression test: BaseFeature used to default `regions`/`faults` to a
    single mutable list shared across all instances (a classic mutable
    default argument bug). This was fixed by defaulting to None and copying
    into a fresh list per-instance in BaseFeature.__init__. Confirm here that
    mutating one instance's `.regions` does not leak into another instance
    that was also constructed with the default (no explicit regions passed).
    """
    feature_a = AnalyticalGeologicalFeature(
        name="feature_a", vector=np.array([1.0, 0.0, 0.0]), origin=np.array([0.0, 0.0, 0.0])
    )
    feature_b = AnalyticalGeologicalFeature(
        name="feature_b", vector=np.array([0.0, 1.0, 0.0]), origin=np.array([0.0, 0.0, 0.0])
    )

    assert feature_a.regions == []
    assert feature_b.regions == []
    assert feature_a.regions is not feature_b.regions

    region = Region(PlaneFeature(), value=0.0, sign=True)
    feature_a.regions.append(region)

    assert feature_a.regions == [region]
    assert feature_b.regions == []
