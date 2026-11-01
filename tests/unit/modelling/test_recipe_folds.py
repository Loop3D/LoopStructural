"""Recipes rebuild fold frames and folded foliations, including refolded folds."""
import json

import numpy as np
import pandas as pd
import pytest

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_laurent2016, load_noddy_single_fold
from LoopStructural.modelling.core._model_serializer import ModelSerializer
from LoopStructural.utils import LoopValueError

FOLD_AXIS = [-6.51626577e-06, -5.00013645e-01, -8.66017526e-01]


def single_fold_model(explicit_frame=True):
    data, boundary = load_noddy_single_fold()
    data = data.sort_values("random")
    model = GeologicalModel(boundary[0, :], boundary[1, :])
    model.data = pd.concat([data[:20], data[data["feature_name"] == "s1"]])
    frame = model.create_and_add_fold_frame(
        "s1", interpolatortype="PLI", nelements=3000, buffer=0.5
    )
    model.create_and_add_folded_foliation(
        "s0",
        **({"fold_frame": frame} if explicit_frame else {}),
        nelements=3000,
        fold_axis=FOLD_AXIS,
        limb_wl=12000,
        buffer=0.5,
    )
    model.update()
    return model, boundary


def refolded_model():
    data, bb = load_laurent2016()
    extra = pd.DataFrame(
        [[5923.504395, 4748.135254, 3588.621094, "s2", 1.0]],
        columns=["X", "Y", "Z", "feature_name", "val"],
    )
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.data = pd.concat([data, extra], sort=False)
    s2 = model.create_and_add_fold_frame("s2", nelements=3000, buffer=0.5, solver="lu", damp=True)
    s1 = model.create_and_add_folded_fold_frame(
        "s1", fold_frame=s2, av_fold_axis=True, nelements=5000, buffer=0.3, limb_wl=4
    )
    model.create_and_add_folded_foliation(
        "s0", fold_frame=s1, av_fold_axis=True, nelements=5000, buffer=0.2
    )
    model.update()
    return model, bb


def rebuild(model):
    rebuilt = GeologicalModel.from_recipe_json(model.to_recipe_json())
    rebuilt.update()
    return rebuilt


def assert_same(model, rebuilt, boundary, names_and_frames):
    points = np.random.default_rng(0).uniform(boundary[0], boundary[1], (300, 3))
    for name, is_frame in names_and_frames:
        if is_frame:
            for i in range(3):
                np.testing.assert_allclose(
                    rebuilt[name][i].evaluate_value(points),
                    model[name][i].evaluate_value(points),
                    atol=1e-6,
                    equal_nan=True,
                )
        else:
            np.testing.assert_allclose(
                rebuilt.evaluate_feature_value(name, points),
                model.evaluate_feature_value(name, points),
                atol=1e-6,
                equal_nan=True,
            )


# ----------------------------------------------------------------------------- single fold
def test_a_fold_frame_and_folded_foliation_are_rebuilt_exactly():
    model, boundary = single_fold_model()
    rebuilt = rebuild(model)
    assert [(f.name, type(f)) for f in rebuilt.features] == [(f.name, type(f)) for f in model.features]
    assert_same(model, rebuilt, boundary, [("s1", True), ("s0", False)])


def test_the_fold_parameters_are_stored():
    model, _ = single_fold_model()
    frame, foliation = json.loads(model.to_recipe_json())["model"]["features"]
    assert (frame["kind"], frame["build"]["feature_type"]) == ("fold_frame", "fold_frame")
    assert frame["build"]["parameters"]["interpolatortype"] == "PLI"
    assert (foliation["kind"], foliation["build"]["feature_type"]) == ("foliation", "folded_foliation")
    parameters = foliation["build"]["parameters"]
    assert parameters["fold_frame"] == "s1"
    assert parameters["fold_axis"] == FOLD_AXIS
    assert parameters["limb_wl"] == 12000
    assert parameters["nelements"] == 3000


def test_a_folded_foliation_that_used_the_last_feature_as_its_frame_rebuilds_exactly():
    """fold_frame=None means 'the last feature'; the recipe records the frame that was used."""
    model, boundary = single_fold_model(explicit_frame=False)
    foliation = json.loads(model.to_recipe_json())["model"]["features"][1]
    assert foliation["build"]["parameters"]["fold_frame"] == "s1"
    assert_same(model, rebuild(model), boundary, [("s1", True), ("s0", False)])


def test_a_model_built_with_index_records_the_frame_it_used_and_rebuilds_identically():
    """A frame inserted between a folded foliation and its frame must not become its frame."""
    data, boundary = load_noddy_single_fold()
    data = data.sort_values("random")
    s1_rows = data[data["feature_name"] == "s1"]
    other = s1_rows.copy()
    other["feature_name"] = "sB"
    model = GeologicalModel(boundary[0, :], boundary[1, :])
    model.data = pd.concat([data[:20], s1_rows, other])
    model.create_and_add_fold_frame("s1", interpolatortype="PLI", nelements=3000, buffer=0.5)
    model.create_and_add_folded_foliation(
        "s0", nelements=3000, fold_axis=FOLD_AXIS, limb_wl=12000, buffer=0.5
    )
    model.create_and_add_fold_frame(
        "sB", interpolatortype="PLI", nelements=3000, buffer=0.5, index=1
    )
    model.update()
    assert [f.name for f in model.features] == ["s1", "sB", "s0"]
    entries = json.loads(model.to_recipe_json())["model"]["features"]
    assert entries[2]["build"]["parameters"]["fold_frame"] == "s1"
    assert_same(
        model, rebuild(model), boundary, [("s1", True), ("sB", True), ("s0", False)]
    )


def test_an_array_fold_axis_is_stored_as_a_list_and_rebuilds_exactly():
    data, boundary = load_noddy_single_fold()
    data = data.sort_values("random")
    model = GeologicalModel(boundary[0, :], boundary[1, :])
    model.data = pd.concat([data[:20], data[data["feature_name"] == "s1"]])
    frame = model.create_and_add_fold_frame(
        "s1", interpolatortype="PLI", nelements=3000, buffer=0.5
    )
    axis = np.array(FOLD_AXIS)
    model.create_and_add_folded_foliation(
        "s0", fold_frame=frame, nelements=3000, fold_axis=axis, limb_wl=12000, buffer=0.5
    )
    model.update()
    axis[:] = 0.0  # changing the array afterwards must not change the recipe
    entries = json.loads(model.to_recipe_json())["model"]["features"]
    assert entries[1]["build"]["parameters"]["fold_axis"] == FOLD_AXIS
    assert_same(model, rebuild(model), boundary, [("s1", True), ("s0", False)])


def test_a_fold_frame_is_recognised_as_a_fold_frame():
    model, _ = single_fold_model()
    assert ModelSerializer.feature_recipe_kind(model["s1"]) == "fold_frame"
    assert ModelSerializer.feature_recipe_kind(model["s0"]) == "foliation"


# ------------------------------------------------------------------------------- refolded
def test_a_refolded_fold_is_rebuilt_exactly():
    model, bb = refolded_model()
    recipe = json.loads(model.to_recipe_json())["model"]["features"]
    assert [(f["name"], f["build"]["feature_type"]) for f in recipe] == [
        ("s2", "fold_frame"),
        ("s1", "folded_fold_frame"),
        ("s0", "folded_foliation"),
    ]
    assert recipe[1]["build"]["parameters"]["fold_frame"] == "s2"
    assert recipe[2]["build"]["parameters"]["fold_frame"] == "s1"
    rebuilt = rebuild(model)
    assert [(f.name, type(f)) for f in rebuilt.features] == [(f.name, type(f)) for f in model.features]
    assert_same(model, rebuilt, bb, [("s2", True), ("s1", True), ("s0", False)])


# -------------------------------------------------------------------------------- refusals
def test_a_fold_frame_that_is_not_an_earlier_entry_is_refused():
    model, _ = single_fold_model()
    recipe = model.to_recipe_dict()
    recipe["model"]["features"][1]["build"]["parameters"]["fold_frame"] = "not_there"
    with pytest.raises(LoopValueError, match="fold frame 'not_there'.*not an earlier entry"):
        GeologicalModel.from_recipe_dict(recipe)


def test_a_fold_with_a_function_parameter_is_refused():
    model, _ = single_fold_model()
    model["s0"]._recipe_creation["parameters"]["fold_limb_rotation"] = lambda x: x
    with pytest.raises(LoopValueError, match=r"parameter 'fold_limb_rotation' of feature 's0'"):
        model.to_recipe_json()
