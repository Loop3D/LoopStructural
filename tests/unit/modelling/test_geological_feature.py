import sys

import numpy as np

from LoopStructural.modelling.features import (
    AnalyticalGeologicalFeature,
    FeatureType,
    GeologicalFeature,
)


def test_constructors():
    # test constructors work and that the types are set correctly
    # base_feature = GeologicalFeature("test", None, [], [], None)
    # assert base_feature.type == FeatureType.BASE
    # assert base_feature.name == "test"
    feature = GeologicalFeature("test", None, [], [], None)
    assert feature.type == FeatureType.INTERPOLATED
    assert feature.name == "test"
    feature = AnalyticalGeologicalFeature("test", [0, 0, 1], [0, 0, 0], [], [], None, None)
    # for analytical feature check that the evaluated value is correct.
    # this should be the distance from origin to the point in the direction
    # of the direction vector
    assert feature.type == FeatureType.ANALYTICAL
    assert feature.name == "test"
    assert feature.evaluate_value([0, 0, 0]) == 0
    assert np.all(feature.evaluate_gradient([0, 0, 0]) - np.array([0, 0, 1]) == 0)
    assert feature.evaluate_value([0, 0, 1]) == 1
    assert feature.evaluate_value([0, 0, -1]) == -1
    assert feature.evaluate_value([0, 1, 0]) == 0


def test_toggle_faults():
    base_feature = GeologicalFeature("test", None, [], [], None)
    assert base_feature.faults_enabled is True
    base_feature.toggle_faults()
    assert base_feature.faults_enabled is False
    base_feature.toggle_faults()
    assert base_feature.faults_enabled is True


def test_tojson():
    base_feature = GeologicalFeature("test", None, [], [], None)
    import json

    from LoopStructural.utils import LoopJSONEncoder

    json.dumps(base_feature, cls=LoopJSONEncoder)


if __name__ == "__main__":
    test_constructors()
    test_toggle_faults()
    test_tojson()
    print("All tests passed")
    sys.exit(0)


def test_fillnan_nearest_keeps_nan_for_points_with_a_nan_coordinate():
    """fillnan='nearest' fills a finite point and leaves a NaN-coordinate point NaN.

    Since a point with a NaN coordinate gives NaN, the KD-tree query got NaN
    points and raised "'x' must be finite" (seen in fault gradients near faults).
    """

    class LinearInterpolator:
        def evaluate_value(self, points):
            return points[:, 0] + 2.0 * points[:, 1]

    feature = GeologicalFeature(
        "field",
        builder=None,
        # The region removes x > 0.5, so that point needs a nearest value.
        regions=[lambda points: ~(points[:, 0] > 0.5)],
        interpolator=LinearInterpolator(),
    )
    points = np.array([[0.0, 1.0, 0.0], [0.6, 1.0, 0.0], [np.nan, 1.0, 0.0]])
    values = feature.evaluate_value(points, fillnan="nearest")
    np.testing.assert_allclose(values, [2.0, 2.0, np.nan], equal_nan=True)


def test_fillnan_nearest_without_known_values_gives_nan():
    class NanInterpolator:
        def evaluate_value(self, points):
            return np.full(points.shape[0], np.nan)

    feature = GeologicalFeature("field", builder=None, regions=[], interpolator=NanInterpolator())
    values = feature.evaluate_value(np.zeros((3, 3)), fillnan="nearest")
    assert np.all(np.isnan(values))
