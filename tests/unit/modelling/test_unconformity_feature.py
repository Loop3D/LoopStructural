import numpy as np
import pytest

from LoopStructural import GeologicalModel
from LoopStructural.modelling.features import FeatureType
from LoopStructural.modelling.features._unconformity_feature import UnconformityFeature


@pytest.fixture()
def strati_feature(horizontal_data):
    model = GeologicalModel([0, 0, 0], [1, 1, 1])
    model.data = horizontal_data
    return model.create_and_add_foliation("strati")


def test_unconformity_feature_name_and_type(strati_feature):
    uc = UnconformityFeature(strati_feature, 0.15, sign=True)

    assert uc.name == "__strati_unconformity"
    assert uc.type == FeatureType.UNCONFORMITY
    assert uc.sign is True
    assert uc.value == 0.15
    assert uc.parent is strati_feature


def test_unconformity_feature_onlap_sets_onlap_type(strati_feature):
    uc = UnconformityFeature(strati_feature, 0.15, sign=True, onlap=True)

    assert uc.type == FeatureType.ONLAPUNCONFORMITY


def test_unconformity_feature_faults_delegates_to_parent(strati_feature):
    uc = UnconformityFeature(strati_feature, 0.15, sign=True)

    assert uc.faults is strati_feature.faults


def test_unconformity_feature_evaluate_sign_true_is_less_equal(strati_feature):
    uc = UnconformityFeature(strati_feature, 0.15, sign=True)

    # from horizontal_data: z=0.25 -> val ~0, z=0.55 -> val ~0.3
    points = np.array([[0.5, 0.5, 0.25], [0.5, 0.5, 0.55]])
    result = uc.evaluate(points)

    assert result.dtype == bool
    # value at z=0.25 (~0) is <= 0.15 -> True (above unconformity)
    # value at z=0.55 (~0.3) is > 0.15 -> False (below unconformity)
    assert np.array_equal(result, np.array([True, False]))


def test_unconformity_feature_evaluate_sign_false_is_greater_equal(strati_feature):
    uc = UnconformityFeature(strati_feature, 0.15, sign=False)

    points = np.array([[0.5, 0.5, 0.25], [0.5, 0.5, 0.55]])
    result = uc.evaluate(points)

    assert np.array_equal(result, np.array([False, True]))


def test_unconformity_feature_call_matches_evaluate(strati_feature):
    uc = UnconformityFeature(strati_feature, 0.15, sign=True)

    points = np.array([[0.5, 0.5, 0.25], [0.5, 0.5, 0.55]])
    assert np.array_equal(uc(points), uc.evaluate(points))


def test_unconformity_feature_inverse_flips_sign_and_keeps_parent(strati_feature):
    uc = UnconformityFeature(strati_feature, 0.15, sign=True)
    inv = uc.inverse()

    assert inv.sign is False
    assert inv.parent is strati_feature
    assert inv.value == uc.value
    assert inv.name == uc.name + "_inverse"
    assert inv.type == FeatureType.UNCONFORMITY


def test_unconformity_feature_inverse_preserves_onlap_type(strati_feature):
    uc = UnconformityFeature(strati_feature, 0.15, sign=True, onlap=True)
    inv = uc.inverse()

    assert inv.type == FeatureType.ONLAPUNCONFORMITY


def test_unconformity_feature_to_json(strati_feature):
    uc = UnconformityFeature(strati_feature, 0.15, sign=True)
    json = uc.to_json()

    assert json["value"] == 0.15
    assert json["sign"] is True
    assert json["parent"] == strati_feature.name
