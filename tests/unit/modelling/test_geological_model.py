import json

import numpy as np
import pandas as pd
import pytest

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius


@pytest.mark.parametrize("origin, maximum", [([0, 0, 0], [5, 5, 5]), ([10, 10, 10], [15, 15, 15])])
def test_create_geological_model(origin, maximum):
    model = GeologicalModel(origin, maximum)
    # origin/maximum are world coordinates; the interpolation frame is
    # anchored at `origin` internally via set_local_transform, without
    # shifting the public origin/maximum themselves.
    assert (model.bounding_box.origin - np.array(origin)).sum() == 0
    assert (model.bounding_box.maximum - np.array(maximum)).sum() == 0
    assert (model.bounding_box.local_origin - np.array(origin)).sum() == 0


def test_prepare_data_keeps_world_coordinates():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    # Data is kept in world coordinates end-to-end; the interpolator projects
    # into its local frame when constraints are set, not at data ingestion.
    expected = data[['X', 'Y', 'Z']].values
    actual = model.prepare_data(model.data)[['X', 'Y', 'Z']].values
    assert np.allclose(actual, expected, atol=1e-6)


def test_access_feature_model():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation("strati")
    assert s0 == model["strati"]


def test_model_recipe_roundtrip_inline_data():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data.iloc[:3].copy())
    model.stratigraphic_column.add_unit("strati", thickness=1.0)

    recipe = model.to_recipe_dict()
    restored = GeologicalModel.from_recipe_dict(recipe)

    assert restored.bounding_box.to_dict() == model.bounding_box.to_dict()
    pd.testing.assert_frame_equal(restored.data, model.data)
    assert restored.stratigraphic_column.to_dict() == model.stratigraphic_column.to_dict()


def test_model_recipe_roundtrip_data_reference(tmp_path):
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data.iloc[:3].copy())

    data_path = tmp_path / "model-data.csv"
    model.data.to_csv(data_path, index=False)

    recipe = model.to_recipe_dict(data_reference=data_path)
    restored = GeologicalModel.from_recipe_dict(recipe)

    assert recipe["model"]["data_source"]["kind"] == "reference"
    assert recipe["model"]["data_source"]["path"] == str(data_path)
    pd.testing.assert_frame_equal(restored.data, model.data)


def test_model_recipe_roundtrip_state_includes_features_and_faults():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    feature_data = data.iloc[:10].copy()
    feature_data.loc[:, "feature_name"] = "strati"
    feature_data_2 = feature_data.copy()
    feature_data_2.loc[:, "feature_name"] = "strati_2"
    model.set_model_data(pd.concat([feature_data, feature_data_2], ignore_index=True))

    feature_a = model.create_and_add_foliation("strati")
    feature_b = model.create_and_add_foliation("strati_2")
    feature_a.faults = [feature_b]

    recipe = model.to_recipe_dict()
    restored = GeologicalModel.from_recipe_dict(recipe)

    assert [feature["name"] for feature in recipe["model"]["features"]] == [
        feature_a.name,
        feature_b.name,
    ]
    assert [feature.name for feature in restored.features] == [feature_a.name, feature_b.name]
    assert [fault.name for fault in restored.features[0].faults] == [feature_b.name]


def test_recipe_to_json_string():
    """Test that to_recipe_json returns a valid JSON string."""
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data.iloc[:3].copy())
    model.stratigraphic_column.add_unit("strati", thickness=1.0)

    json_str = model.to_recipe_json()

    # Verify it's a valid JSON string
    assert isinstance(json_str, str)
    recipe = json.loads(json_str)
    assert recipe["schema"] == "LoopStructural.GeologicalModelRecipe"
    assert recipe["version"] == 1
    assert "model" in recipe


def test_recipe_from_json_string():
    """Test that from_recipe_json can parse JSON and reconstruct the model."""
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data.iloc[:3].copy())
    model.stratigraphic_column.add_unit("strati", thickness=1.0)

    json_str = model.to_recipe_json()
    restored = GeologicalModel.from_recipe_json(json_str)

    assert restored.bounding_box.to_dict() == model.bounding_box.to_dict()
    pd.testing.assert_frame_equal(restored.data, model.data)
    assert restored.stratigraphic_column.to_dict() == model.stratigraphic_column.to_dict()


def test_recipe_json_roundtrip_with_features():
    """Test JSON roundtrip preserves feature relationships."""
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    feature_data = data.iloc[:10].copy()
    feature_data.loc[:, "feature_name"] = "strati"
    feature_data_2 = feature_data.copy()
    feature_data_2.loc[:, "feature_name"] = "strati_2"
    model.set_model_data(pd.concat([feature_data, feature_data_2], ignore_index=True))

    feature_a = model.create_and_add_foliation("strati")
    feature_b = model.create_and_add_foliation("strati_2")
    feature_a.faults = [feature_b]

    json_str = model.to_recipe_json()
    restored = GeologicalModel.from_recipe_json(json_str)

    assert len(restored.features) == 2
    assert [f.name for f in restored.features] == ["strati", "strati_2"]
    assert [fault.name for fault in restored.features[0].faults] == ["strati_2"]


def test_save_recipe_inline_data(tmp_path):
    """Test saving recipe with inline data to JSON file."""
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data.iloc[:3].copy())

    recipe_file = tmp_path / "model_recipe.json"
    model.save_recipe(recipe_file)

    # Verify file exists and is valid JSON
    assert recipe_file.exists()
    with open(recipe_file) as f:
        recipe = json.load(f)
    assert recipe["schema"] == "LoopStructural.GeologicalModelRecipe"
    assert recipe["model"]["data_source"]["kind"] == "inline"


def test_load_recipe_inline_data(tmp_path):
    """Test loading recipe with inline data from JSON file."""
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data.iloc[:3].copy())

    recipe_file = tmp_path / "model_recipe.json"
    model.save_recipe(recipe_file)

    restored = GeologicalModel.load_recipe(recipe_file)

    assert restored.bounding_box.to_dict() == model.bounding_box.to_dict()
    pd.testing.assert_frame_equal(restored.data, model.data)


def test_save_load_recipe_with_data_reference(tmp_path):
    """Test save/load roundtrip with external data reference."""
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data.iloc[:3].copy())

    data_file = tmp_path / "model_data.csv"
    model.data.to_csv(data_file, index=False)

    recipe_file = tmp_path / "model_recipe.json"
    model.save_recipe(recipe_file, data_reference=data_file)

    # Verify recipe references the data file
    with open(recipe_file) as f:
        recipe = json.load(f)
    assert recipe["model"]["data_source"]["kind"] == "reference"

    # Load and verify
    restored = GeologicalModel.load_recipe(recipe_file)
    pd.testing.assert_frame_equal(restored.data, model.data)


def test_recipe_json_error_handling():
    """Test error handling for invalid JSON input."""
    with pytest.raises(TypeError, match="json_str must be a string"):
        GeologicalModel.from_recipe_json(123)

    with pytest.raises(TypeError, match="json_str is not valid JSON"):
        GeologicalModel.from_recipe_json("not valid json")


def test_load_recipe_file_not_found(tmp_path):
    """Test error handling for missing recipe file."""
    recipe_file = tmp_path / "nonexistent_recipe.json"
    with pytest.raises(FileNotFoundError, match="Recipe file not found"):
        GeologicalModel.load_recipe(recipe_file)


def test_recipe_json_formatting():
    """Test that JSON is properly formatted with indentation."""
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data.iloc[:3].copy())

    json_str = model.to_recipe_json(indent=2)
    # Check that it's properly indented (has newlines and spaces)
    assert "\n" in json_str
    assert "  " in json_str

    # Verify custom indent works
    json_str_no_indent = model.to_recipe_json(indent=None)
    assert len(json_str_no_indent) < len(json_str)


def test_from_file_default_loads_pickled_model(tmp_path):
    """Default behaviour (allow_pickle=True) should still load a valid
    pickled model file successfully, preserving backward compatibility."""
    pytest.importorskip("dill")
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data.iloc[:3].copy())

    model_file = tmp_path / "model.pkl"
    model.to_file(model_file)

    restored = GeologicalModel.from_file(model_file)

    assert restored is not None
    assert isinstance(restored, GeologicalModel)
    assert restored.bounding_box.to_dict() == model.bounding_box.to_dict()


def test_from_file_allow_pickle_false_raises_without_loading(tmp_path):
    """allow_pickle=False must refuse to unpickle the file and raise a clear,
    actionable error rather than attempting to load it."""
    from LoopStructural.utils import LoopValueError

    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data.iloc[:3].copy())

    model_file = tmp_path / "model.pkl"
    model.to_file(model_file)

    with pytest.raises(LoopValueError, match="allow_pickle"):
        GeologicalModel.from_file(model_file, allow_pickle=False)


def test_intrusions_lists_intrusion_features():
    from LoopStructural.modelling.intrusions import IntrusionFeature

    model = GeologicalModel(np.zeros(3), np.ones(3))
    intrusion = IntrusionFeature(None, None, name="intrusion", model=model)
    model.features.append(intrusion)
    assert model.intrusions == [intrusion]


def test_save_with_default_arguments(tmp_path):
    from tests.fixtures.fault_models import OLD, add_fault, add_layers

    model = GeologicalModel(np.zeros(3), np.ones(3))
    model.nsteps = [5, 5, 5]
    fault = add_fault(
        model,
        "fault",
        OLD["center"],
        OLD["dip"],
        OLD["dip_azimuth"],
        OLD["displacement"],
        nelements=500,
    )
    add_layers(model, [fault], nelements=500)
    model.stratigraphic_column.add_unit("a", thickness=0.3)
    model.stratigraphic_column.add_unit("b", thickness=0.3)
    model.stratigraphic_column.group_mapping["Group_0"] = "strat"

    model.save(str(tmp_path / "m.pkl"))

    names = {p.name for p in tmp_path.iterdir()}
    assert "m_block_model.pkl" in names
    # one data file per data set: the fault frame data do not overwrite each other
    fault_data = [d.name for d in model["fault"].get_data()]
    assert len(fault_data) > 1
    for d in fault_data:
        assert f"m_{d}.pkl" in names
    for d in model["strat"].get_data():
        assert f"m_{d.name}.pkl" in names


def test_get_stratigraphic_surfaces_returns_only_the_requested_units():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.data = data
    model.nsteps = [5, 5, 5]
    model.create_and_add_foliation("strati", nelements=500)
    for name, value in (("a", 0.0), ("b", 60.0), ("c", 250.0)):
        model.stratigraphic_column.add_unit(name, thickness=value or 60.0)
    model.stratigraphic_column.group_mapping["Group_0"] = "strati"
    surfaces = model.get_stratigraphic_surfaces(units=["b"])
    assert [s.name for s in surfaces] == ["b"]
def test_set_stratigraphic_column_keeps_group_names():
    """Regression test: with two groups, each unit got the name of the other group."""
    model = GeologicalModel([0, 0, 0], [1, 1, 1])
    with pytest.warns(DeprecationWarning):
        model.set_stratigraphic_column(
            {
                "A": {"a1": {"min": 0, "max": 1}, "a2": {"min": 1, "max": 2}},
                "B": {"b1": {"min": 0, "max": 1}},
            }
        )
    column = model.stratigraphic_column
    assert column.get_group_for_unit_name("a1").name == "A"
    assert column.get_group_for_unit_name("b1").name == "B"


def test_evaluate_model_returns_unit_ids():
    """Regression test: evaluate_model gave a running counter, not unit.id."""
    from LoopStructural.modelling.features._analytical_feature import (
        AnalyticalGeologicalFeature,
    )

    model = GeologicalModel([0, 0, 0], [10, 10, 10])
    column = model.stratigraphic_column
    column.clear(basement=False)
    for i in (1, 2, 3):
        column.add_unit(f"u{i}", thickness=1.0, id=i)
    column.group_mapping["Group_0"] = "strati"
    model._add_feature(
        AnalyticalGeologicalFeature(
            name="strati", vector=np.array([0.0, 0.0, 1.0]), origin=np.zeros(3)
        )
    )
    points = np.array([[1.0, 1.0, z] for z in (0.5, 1.5, 2.5)])
    assert model.evaluate_model(points).tolist() == [1, 2, 3]


if __name__ == "__main__":
    test_prepare_data_keeps_world_coordinates()
