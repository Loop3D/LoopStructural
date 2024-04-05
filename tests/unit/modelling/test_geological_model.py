from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius
import numpy as np


def test_create_geological_model():
    model = GeologicalModel([0, 0, 0], [5, 5, 5])
    assert (model.origin - np.array([0, 0, 0])).sum() == 0
    assert (model.maximum - np.array([5, 5, 5])).sum() == 0


def test_access_feature_model():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation("strati")
    assert s0 == model["strati"]
