from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius

def test_create_model():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0,:],bb[1,:])

def test_add_data():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0,:],bb[1,:])
    model.set_model_data(data)

def test_create_stratigraphy_FDI_cg():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation('s0',
                                        interpolatortype='FDI',
                                        nelements=1000,
                                        solver='cg',
                                        damp=False)


def test_remove_constraints_PLI():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation('s0',
                                        interpolatortype='FDI',
                                        nelements=1000,
                                        solver='cg',
                                        damp=False)

def test_create_stratigraphy_FDI_lu():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation('s0',
                                        interpolatortype='FDI',
                                        nelements=1000,
                                        solver='lu',
                                        damp=True)


def test_create_stratigraphy_FDI_pyamg():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation('s0',
                                        interpolatortype='FDI',
                                        nelements=1000,
                                        solver='pyamg',
                                        damp=True)

def test_create_stratigraphy_PLI_cg():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation('s0',
                                        interpolatortype='PLI',
                                        nelements=1000,
                                        solver='cg',
                                        damp=False)

def test_create_stratigraphy_PLI_lu():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation('s0',
                                        interpolatortype='PLI',
                                        nelements=1000,
                                        solver='lu',
                                        damp=True)


def test_create_stratigraphy_PLI_pyamg():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation('s0',
                                        interpolatortype='PLI',
                                        nelements=1000,
                                        solver='pyamg',
                                        damp=True)

def test_model_with_data_outside_of_bounding_box():
    pass

