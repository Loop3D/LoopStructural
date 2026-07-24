import numpy as np
import pytest

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius
from LoopStructural.modelling.core._feature_registry import FeatureBuilderRegistry


def test_builtin_feature_types_registered():
    assert FeatureBuilderRegistry.registered_types() == sorted(
        [
            "domain_fault",
            "fault",
            "fold_frame",
            "folded_fold_frame",
            "folded_foliation",
            "foliation",
            "intrusion",
        ]
    )


def test_create_and_add_feature_matches_convenience_method():
    data, bb = load_claudius()

    model_a = GeologicalModel(bb[0, :], bb[1, :])
    model_a.set_model_data(data)
    via_wrapper = model_a.create_and_add_foliation("strati")

    model_b = GeologicalModel(bb[0, :], bb[1, :])
    model_b.set_model_data(data)
    via_generic = model_b.create_and_add_feature("foliation", "strati")

    assert via_wrapper.name == via_generic.name == "strati"
    xyz = model_a.regular_grid(shuffle=False)
    assert np.allclose(
        via_wrapper.evaluate_value(xyz),
        via_generic.evaluate_value(xyz),
        equal_nan=True,
    )


def test_convert_feature_to_structural_frame_returns_frame():
    from LoopStructural.modelling.features import StructuralFrame

    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    model.create_and_add_foliation("strati")

    frame = model.convert_feature_to_structural_frame("strati")

    assert isinstance(frame, StructuralFrame)
    assert model["strati"] is frame


def test_add_fold_to_feature_rejects_non_fold_frame():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    model.create_and_add_foliation("strati")

    with pytest.raises(ValueError):
        model.add_fold_to_feature("strati", fold_frame="not a fold frame")
