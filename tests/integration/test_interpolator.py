from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius, load_horizontal
import numpy as np


def model_fit(model, data):
    diff = (
        model.evaluate_feature_value("strati", data.loc[data["val"].notna(), ["X", "Y", "Z"]])
        - data.loc[data["val"].notna(), "val"]
    )
    assert np.std(diff) < 30


def test_create_model():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])


def test_add_data():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)


def test_create_stratigraphy_FDI_cg():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation(
        "strati", interpolatortype="FDI", nelements=1000, solver="cg", damp=False
    )
    model.update()
    model_fit(model, data)


def test_remove_constraints_PLI():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation(
        "strati", interpolatortype="FDI", nelements=1000, solver="cg", damp=False
    )
    model.update()
    model_fit(model, data)


def test_create_stratigraphy_FDI_lu():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation(
        "strati", interpolatortype="FDI", nelements=1000, solver="lu", damp=True
    )
    model.update()
    model_fit(model, data)


def test_create_stratigraphy_FDI_pyamg():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation(
        "strati", interpolatortype="FDI", nelements=1000, solver="pyamg", damp=True
    )
    model.update()
    model_fit(model, data)


def test_create_stratigraphy_PLI_cg():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation(
        "strati", interpolatortype="PLI", nelements=1000, solver="cg", damp=False
    )
    model.update()
    model_fit(model, data)


def test_create_stratigraphy_PLI_lu():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation(
        "strati", interpolatortype="PLI", nelements=1000, solver="lu", damp=True
    )
    model.update()
    model_fit(model, data)


def test_create_stratigraphy_PLI_pyamg():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation(
        "strati", interpolatortype="PLI", nelements=1000, solver="pyamg", damp=True
    )
    model.update()
    model_fit(model, data)


def test_model_with_data_outside_of_bounding_box():
    pass


def test_horizontal_layers(interpolatortype, nelements):
    data, bb = load_horizontal()
    model = GeologicalModel(bb[0, :], bb[1, :])

    model.data = data
    model.create_and_add_foliation("strati", interpolatortype=interpolatortype, nelements=1e4)

    assert np.all(np.isclose(model["strati"].evaluate_value(data[["X", "Y", "Z"]]), data["val"]))


def test_horizontal_layers(interpolatortype, nelements):
    data, bb = load_horizontal()
    model = GeologicalModel(bb[0, :], bb[1, :])

    model.data = data
    model.create_and_add_foliation("strati", interpolatortype=interpolatortype, nelements=1e4)

    assert np.all(np.isclose(model["strati"].evaluate_value(data[["X", "Y", "Z"]]), data["val"]))


if __name__ == "__main__":
    test_create_model()
    test_add_data()
    test_create_stratigraphy_FDI_cg()
    test_remove_constraints_PLI()
    test_create_stratigraphy_FDI_lu()
    test_create_stratigraphy_FDI_pyamg()
    test_create_stratigraphy_PLI_cg()
    test_create_stratigraphy_PLI_lu()
    test_create_stratigraphy_PLI_pyamg()
    test_model_with_data_outside_of_bounding_box()
    test_horizontal_layers("FDI", 1000)
    test_horizontal_layers("PLI", 1000)
    test_horizontal_layers("FDI", 1000)
    test_horizontal_layers("PLI", 1000)
    print("ok")
