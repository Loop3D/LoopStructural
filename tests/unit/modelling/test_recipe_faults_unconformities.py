"""Recipes rebuild faults and unconformities, not just foliations."""
import json

import numpy as np
import pandas as pd
import pytest

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius
from LoopStructural.modelling.core._model_serializer import ModelSerializer
from LoopStructural.utils import LoopValueError


def fault_data():
    return pd.DataFrame(
        [
            [5, 5, 5, 0, 0.70710678, 0.0, 0.70710678, 0, "fault"],
            [5, 5, 5, 0, -0.70710678, 0.0, 0.70710678, 1, "fault"],
            [8, 5, 5, 0, 0, 0, 1, np.nan, "strati"],
        ],
        columns=["X", "Y", "Z", "val", "nx", "ny", "nz", "coord", "feature_name"],
    )


def fault_model(**foliation_parameters):
    model = GeologicalModel([0, 0, 0], [10, 10, 10])
    model.data = fault_data()
    model.create_and_add_fault(
        "fault", 1, nelements=1000, buffer=0.5, major_axis=10, minor_axis=3, intermediate_axis=10
    )
    model.create_and_add_foliation("strati", nelements=1000, **foliation_parameters)
    model.update()
    return model


def two_group_model(onlap):
    data, bb = load_claudius()
    data = data.reset_index()
    data.loc[:, "val"] *= -1
    data.loc[:, ["nx", "ny", "nz"]] *= -1
    data.loc[792, "feature_name"] = "strati2"
    data.loc[792, ["nx", "ny", "nz"]] = [0, 0, 1]
    data.loc[792, "val"] = 0
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.data = data
    upper = model.create_and_add_foliation("strati2", nelements=1000)
    if onlap:
        model.add_onlap_unconformity(upper, 1)
    else:
        model.add_unconformity(upper, 1)
    model.create_and_add_foliation("strati", nelements=1000)
    model.update()
    return model, bb


def rebuild(model):
    rebuilt = GeologicalModel.from_recipe_json(model.to_recipe_json())
    rebuilt.update()
    return rebuilt


def same_field(a, b, name, points):
    first, second = a.evaluate_feature_value(name, points), b.evaluate_feature_value(name, points)
    assert np.array_equal(np.isnan(first), np.isnan(second))
    np.testing.assert_allclose(second, first, atol=1e-6, equal_nan=True)


# ------------------------------------------------------------------------------- faults
def test_a_fault_is_rebuilt_as_a_fault():
    model = fault_model()
    rebuilt = rebuild(model)
    assert [(f.name, type(f)) for f in rebuilt.features] == [(f.name, type(f)) for f in model.features]
    assert type(rebuilt["fault"]).__name__ == "FaultSegment"


def test_the_faulted_field_and_the_displacement_are_reproduced():
    model = fault_model()
    rebuilt = rebuild(model)
    points = np.random.default_rng(0).uniform(0, 10, (300, 3))
    same_field(model, rebuilt, "strati", points)
    np.testing.assert_allclose(
        rebuilt.evaluate_fault_displacements(points),
        model.evaluate_fault_displacements(points),
        atol=1e-8,
    )
    # the fault does something here, so the test is not vacuous
    assert np.abs(model.evaluate_fault_displacements(points)).max() > 0.1


def test_the_fault_parameters_are_stored():
    model = fault_model()
    entry = json.loads(model.to_recipe_json())["model"]["features"][0]
    assert entry["kind"] == "fault"
    parameters = entry["build"]["parameters"]
    assert entry["build"]["feature_type"] == "fault"
    assert parameters["displacement"] == 1
    assert parameters["major_axis"] == 10
    assert parameters["faults"] == []  # a fault with no faults is not "all faults"


def test_a_foliation_given_no_faults_is_not_given_the_faults_when_rebuilt():
    """faults=[] means none; faults=None (the default) means every fault in the model."""
    without = fault_model(faults=[])
    with_faults = fault_model()
    points = np.random.default_rng(1).uniform(0, 10, (300, 3))
    # the two models differ, so the test can tell the cases apart
    assert not np.allclose(
        without.evaluate_feature_value("strati", points),
        with_faults.evaluate_feature_value("strati", points),
    )
    same_field(without, rebuild(without), "strati", points)
    same_field(with_faults, rebuild(with_faults), "strati", points)


def test_a_fault_that_refers_to_a_later_feature_is_refused():
    model = fault_model()
    recipe = model.to_recipe_dict()
    recipe["model"]["features"][1]["build"]["parameters"]["faults"] = ["not_there"]
    with pytest.raises(LoopValueError, match="not earlier entries"):
        GeologicalModel.from_recipe_dict(recipe)


# ------------------------------------------------------------------------ unconformities
@pytest.mark.parametrize("onlap", [False, True], ids=["erode", "onlap"])
def test_an_unconformity_is_rebuilt_and_gives_the_same_model(onlap):
    model, bb = two_group_model(onlap)
    rebuilt = rebuild(model)
    assert [(f.name, int(f.type)) for f in rebuilt.features] == [
        (f.name, int(f.type)) for f in model.features
    ]
    points = np.random.default_rng(0).uniform(bb[0], bb[1], (400, 3))
    same_field(model, rebuilt, "strati2", points)
    same_field(model, rebuilt, "strati", points)
    np.testing.assert_array_equal(rebuilt.evaluate_model(points), model.evaluate_model(points))


@pytest.mark.parametrize("onlap", [False, True], ids=["erode", "onlap"])
def test_the_unconformity_entry(onlap):
    model, _ = two_group_model(onlap)
    entry = json.loads(model.to_recipe_json())["model"]["features"][1]
    assert entry["kind"] == "unconformity"
    assert entry["build"] == {
        "feature_type": "unconformity",
        "parameters": {"feature": "strati2", "value": 1, "onlap": onlap},
    }
    assert ModelSerializer.feature_recipe_kind(model.features[1]) == "unconformity"


def test_an_unconformity_that_refers_to_a_missing_feature_is_refused():
    model, _ = two_group_model(False)
    recipe = model.to_recipe_dict()
    recipe["model"]["features"][1]["build"]["parameters"]["feature"] = "not_there"
    with pytest.raises(LoopValueError, match="not an earlier entry"):
        GeologicalModel.from_recipe_dict(recipe)
