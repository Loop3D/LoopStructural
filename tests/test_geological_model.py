from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius

def test_create_scale_factor_model():
    model = GeologicalModel([0,0,0],[5,5,5])
    assert model.scale_factor == 5

def test_access_feature_model():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation('s0',
                                        interpolatortype='FDI',
                                        nelements=1000,
                                        solver='cg',
                                        damp=False)
    assert s0 == model['s0']

