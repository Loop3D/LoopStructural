import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel
from LoopStructural.datatypes import BoundingBox
from LoopStructural.modelling.features import (
    GeologicalFeature,
    StructuralFrame,
)


def test_structural_frame():
    coordinate_0 = GeologicalFeature("coord0", None)
    coordinate_1 = GeologicalFeature("coord1", None)
    coordinate_2 = GeologicalFeature("coord2", None)
    frame = StructuralFrame("structural_frame", [coordinate_0, coordinate_1, coordinate_2])
    assert frame is not None
    assert frame.name == "structural_frame"


def set_model():
    model = GeologicalModel(np.zeros(3), np.ones(3))
    coordinate_0 = GeologicalFeature("coord0", None)
    coordinate_1 = GeologicalFeature("coord1", None)
    coordinate_2 = GeologicalFeature("coord2", None)
    frame = StructuralFrame("structural_frame", [coordinate_0, coordinate_1, coordinate_2])
    frame.set_model(model)
    assert frame.model == model
    assert frame[0].model == model
    assert frame[1].model == model
    assert frame[2].model == model


def get_item():
    coordinate_0 = GeologicalFeature("coord0", None)
    coordinate_1 = GeologicalFeature("coord1", None)
    coordinate_2 = GeologicalFeature("coord2", None)
    frame = StructuralFrame("structural_frame", [coordinate_0, coordinate_1, coordinate_2])
    assert frame[0] == coordinate_0
    assert frame[1] == coordinate_0
    assert frame[2] == coordinate_0


def test_structural_frame_cross_product():
    pass


def test_create_structural_frame_pli():
    data = pd.DataFrame(
        [
            [5.1, 5.1, 5, 0, 0, 1, 0, 0],
            [5, 5.1, 5, 0, 1, 0, 1, 0],
            [5.1, 5, 5, 1, 0, 0, 2, 0],
        ],
        columns=["X", "Y", "Z", "nx", "ny", "nz", "coord", "val"],
    )
    data["feature_name"] = "fault"

    bb = BoundingBox(origin=np.zeros(3), maximum=np.ones(3) * 10)
    model = GeologicalModel(bb.origin, bb.maximum)

    model.data = data

    fault = model.create_and_add_fault(
        "fault", 10, nelements=2000, steps=4, interpolatortype="PLI", buffer=2
    )
    model.update()

    assert np.all(
        np.isclose(fault[0].evaluate_gradient(np.array([[5, 5, 5]])), [0, 0, 1], atol=1e-1)
    )
    assert np.all(
        np.isclose(fault[1].evaluate_gradient(np.array([[5, 5, 5]])), [0, 1, 0], atol=1e-1)
    )
    assert np.all(
        np.isclose(fault[2].evaluate_gradient(np.array([[5, 5, 5]])), [1, 0, 0], atol=1e-1)
    )


def test_create_structural_frame_fdi():
    data = pd.DataFrame(
        [
            [5.1, 5.1, 5, 0, 0, 1, 0, 0],
            [5, 5.1, 5, 0, 1, 0, 1, 0],
            [5.1, 5, 5, 1, 0, 0, 2, 0],
        ],
        columns=["X", "Y", "Z", "nx", "ny", "nz", "coord", "val"],
    )
    data["feature_name"] = "fault"

    bb = BoundingBox(origin=np.zeros(3), maximum=np.ones(3) * 10)
    model = GeologicalModel(bb.origin, bb.maximum)

    model.data = data

    fault = model.create_and_add_fault(
        "fault", 10, nelements=2000, steps=4, interpolatortype="FDI", buffer=2
    )
    model.update()

    assert np.all(
        np.isclose(fault[0].evaluate_gradient(np.array([[5, 5, 5]])), [0, 0, 1], atol=1e-1)
    )
    assert np.all(
        np.isclose(fault[1].evaluate_gradient(np.array([[5, 5, 5]])), [0, 1, 0], atol=1e-1)
    )
    assert np.all(
        np.isclose(fault[2].evaluate_gradient(np.array([[5, 5, 5]])), [1, 0, 0], atol=1e-1)
    )


if __name__ == "__main__":
    test_create_structural_frame_pli()
