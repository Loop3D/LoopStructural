from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius
import numpy as np
import pandas as pd
import pytest

@pytest.mark.parametrize("origin, maximum", [([0,0,0],[5,5,5]), ([10,10,10],[15,15,15])])
def test_create_geological_model(origin, maximum):
    model = GeologicalModel(origin, maximum)
    assert (model.bounding_box.global_origin - np.array(origin)).sum() == 0
    assert (model.bounding_box.global_maximum - np.array(maximum)).sum() == 0
    assert (model.bounding_box.origin - np.zeros(3)).sum() == 0
    assert (model.bounding_box.maximum - np.ones(3)*5).sum() == 0

def test_rescale_model_data():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    # Check that the model data is rescaled to local coordinates
    expected = data[['X', 'Y', 'Z']].values - bb[None, 0, :]
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

if __name__ == "__main__":
    test_rescale_model_data()