from LoopStructural.modelling.features import (
    GeologicalFeature,
    AnalyticalGeologicalFeature,
    FeatureType,
)
import numpy as np


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
    exit(0)
