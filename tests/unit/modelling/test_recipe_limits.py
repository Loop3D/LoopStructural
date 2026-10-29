"""What a recipe can and cannot store, and that it says so instead of rebuilding something wrong."""
import json

import numpy as np
import pandas as pd
import pytest

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius
from LoopStructural.modelling.core._model_serializer import ModelSerializer
from LoopStructural.utils import LoopValueError


def claudius_model(**build_parameters):
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.data = data
    model.stratigraphic_column.clear(basement=False)
    model.create_and_add_foliation("strati", **{"nelements": 1000, **build_parameters})
    return model, data, bb


def fault_model():
    data = pd.DataFrame(
        [
            [5, 5, 5, 0, 0.70710678, 0.0, 0.70710678, 0, "fault"],
            [5, 5, 5, 0, -0.70710678, 0.0, 0.70710678, 1, "fault"],
            [8, 5, 5, 0, 0, 0, 1, np.nan, "strati"],
        ],
        columns=["X", "Y", "Z", "val", "nx", "ny", "nz", "coord", "feature_name"],
    )
    model = GeologicalModel([0, 0, 0], [10, 10, 10])
    model.data = data
    model.create_and_add_fault(
        "fault", 1, nelements=1000, buffer=0.5, major_axis=10, minor_axis=3, intermediate_axis=10
    )
    model.create_and_add_foliation("strati", nelements=1000, faults=[model["fault"]])
    return model


# ------------------------------------------------------------------ what is stored
def test_the_requested_parameters_are_stored_not_the_realised_ones():
    model, _, _ = claudius_model(nelements=2000, regularisation=10.0)
    entry = json.loads(model.to_recipe_json())["model"]["features"][0]
    parameters = entry["build"]["parameters"]
    assert entry["build"]["feature_type"] == "foliation"
    assert parameters["nelements"] == 2000
    assert parameters["regularisation"] == 10.0
    assert parameters["interpolatortype"] == "FDI"
    # the realised element count is different, and storing it would not reproduce the grid
    model.update()
    assert model["strati"].builder.build_arguments["nelements"] != 2000


def test_data_and_other_features_are_not_stored_as_parameters():
    model, _, _ = claudius_model()
    parameters = json.loads(model.to_recipe_json())["model"]["features"][0]["build"]["parameters"]
    assert not {"data", "faults", "index"} & set(parameters)


def test_group_mapping_survives_a_round_trip():
    model, _, _ = claudius_model()
    model.stratigraphic_column.add_unit("unit_0", thickness=60, id=0)
    model.stratigraphic_column.group_mapping["Group_0"] = "strati"
    rebuilt = GeologicalModel.from_recipe_json(model.to_recipe_json())
    assert rebuilt.stratigraphic_column.group_mapping == {"Group_0": "strati"}


def test_a_column_without_a_group_mapping_still_loads():
    """Recipes written before group_mapping was stored have no such key."""
    model, _, _ = claudius_model()
    model.stratigraphic_column.add_unit("unit_0", thickness=60, id=0)
    recipe = model.to_recipe_dict()
    del recipe["model"]["stratigraphic_column"]["group_mapping"]
    rebuilt = GeologicalModel.from_recipe_dict(recipe)
    assert rebuilt.stratigraphic_column.group_mapping == {}


# ------------------------------------------------------------------ kind detection
def test_a_fault_is_recognised_as_a_fault():
    model = fault_model()
    assert ModelSerializer.feature_recipe_kind(model["fault"]) == "fault"
    assert ModelSerializer.feature_recipe_kind(model["strati"]) == "foliation"


def test_an_unconformity_is_recognised_as_an_unconformity():
    model, _, _ = claudius_model()
    model.add_unconformity(model["strati"], 0)
    kinds = {f.name: ModelSerializer.feature_recipe_kind(f) for f in model.features}
    assert kinds["__strati_unconformity"] == "unconformity"


# ------------------------------------------------------------------ what is refused
def test_a_model_with_a_fault_is_refused_not_rebuilt_without_the_fault():
    with pytest.raises(LoopValueError, match=r"fault.*foliations only"):
        fault_model().to_recipe_json()


def test_a_model_with_an_unconformity_is_refused():
    model, _, _ = claudius_model()
    model.add_unconformity(model["strati"], 0)
    with pytest.raises(LoopValueError, match=r"__strati_unconformity.*unconformity"):
        model.to_recipe_json()


def test_a_feature_built_with_its_own_data_is_refused():
    model, data, _ = claudius_model()
    model.create_and_add_foliation(
        "strati_2", data=data[data["feature_name"] == "strati"].iloc[::2], nelements=1000
    )
    with pytest.raises(LoopValueError, match=r"strati_2.*its own data"):
        model.to_recipe_json()


def test_a_parameter_json_cannot_hold_is_refused():
    model, _, _ = claudius_model(solver=lambda a, b: None)
    with pytest.raises(LoopValueError, match=r"parameter 'solver' of feature 'strati'"):
        model.to_recipe_json()


def test_a_feature_with_no_creation_record_is_refused():
    model, _, _ = claudius_model()
    del model["strati"]._recipe_creation
    with pytest.raises(LoopValueError, match="no record of how it was created"):
        model.to_recipe_json()


# ------------------------------------------------------------------ older recipes
def test_an_older_recipe_without_build_parameters_loads_with_a_warning(caplog):
    import logging

    model, _, _ = claudius_model()
    recipe = model.to_recipe_dict()
    del recipe["model"]["features"][0]["build"]
    with caplog.at_level(logging.WARNING):
        rebuilt = GeologicalModel.from_recipe_dict(recipe)
    assert [f.name for f in rebuilt.features] == ["strati"]
    assert any("no build parameters" in record.getMessage() for record in caplog.records)


def test_an_older_recipe_entry_that_is_not_a_foliation_is_refused():
    model, _, _ = claudius_model()
    recipe = model.to_recipe_dict()
    entry = recipe["model"]["features"][0]
    del entry["build"]
    entry["kind"] = "structural_frame"
    with pytest.raises(LoopValueError, match="foliations only"):
        GeologicalModel.from_recipe_dict(recipe)


# ------------------------------------------------------- parameters are copied, not aliased
def test_changing_a_parameter_after_creation_does_not_change_the_recipe():
    model, _, _ = claudius_model()
    values = [1, 2, 3]
    weights = {"a": [1.0, 2.0]}
    ModelSerializer.record_creation(
        model["strati"], "foliation", {"nelements": 1000, "v": values, "w": weights}
    )
    before = model.to_recipe_json()
    values.append(4)
    weights["a"].append(3.0)
    weights["b"] = 1
    assert model.to_recipe_json() == before
    parameters = json.loads(before)["model"]["features"][0]["build"]["parameters"]
    assert parameters["v"] == [1, 2, 3]
    assert parameters["w"] == {"a": [1.0, 2.0]}


def test_a_numpy_array_parameter_is_stored_as_a_list_and_copied():
    model, _, _ = claudius_model()
    vector = np.array([0.0, 1.0, 0.5])
    ModelSerializer.record_creation(model["strati"], "foliation", {"nelements": 1000, "vec": vector})
    vector[0] = 99.0
    parameters = json.loads(model.to_recipe_json())["model"]["features"][0]["build"]["parameters"]
    assert parameters["vec"] == [0.0, 1.0, 0.5]
