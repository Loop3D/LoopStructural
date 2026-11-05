"""Tests that consumers of model data show and export world locations."""

import itertools

import numpy as np
import pandas as pd
import pytest

from tests.fixtures.coordinate_frames import (
    ANCHORS,
    ROTATIONS_DEG,
    tp,
    world_foliation_model,
)

LENGTH = 12.0
WIDTH = 6500000.0 + 12.0


def _input_rows(data: pd.DataFrame):
    """Return the world locations of the value, gradient and norm rows."""
    values = data[data["val"].notna()]
    gradients = data[data["gx"].notna()]
    norms = data[data["nx"].notna()]
    return values, gradients, norms


def _split_objects(objects):
    """Return the ValuePoints, the gradient VectorPoints and the norm VectorPoints."""
    by_name = {obj.name: obj for obj in objects}
    value = next(o for o in objects if hasattr(o, "values"))
    gradient = next(o for n, o in by_name.items() if n.endswith("_gradient"))
    norm = next(o for n, o in by_name.items() if n.endswith("_norm"))
    return value, gradient, norm


@pytest.mark.parametrize(
    "anchor_index, rotation",
    list(itertools.product(range(len(ANCHORS)), ROTATIONS_DEG)),
)
def test_feature_data_locations_are_world(anchor_index, rotation):
    """The data of a feature has the same world locations as the input rows."""
    anchor = ANCHORS[anchor_index]
    _, feature, data = world_foliation_model(anchor, rotation)
    values, gradients, norms = _input_rows(data)
    value_obj, gradient_obj, norm_obj = _split_objects(feature.get_data())
    tol = tp(LENGTH, max(WIDTH, np.abs(anchor).max()))
    np.testing.assert_allclose(
        value_obj.locations, values[["X", "Y", "Z"]].to_numpy(float), rtol=0, atol=tol
    )
    np.testing.assert_allclose(
        gradient_obj.locations, gradients[["X", "Y", "Z"]].to_numpy(float), rtol=0, atol=tol
    )
    np.testing.assert_allclose(
        norm_obj.locations, norms[["X", "Y", "Z"]].to_numpy(float), rtol=0, atol=tol
    )
    np.testing.assert_allclose(
        gradient_obj.vectors, gradients[["gx", "gy", "gz"]].to_numpy(float), rtol=0, atol=1e-12
    )
    np.testing.assert_allclose(
        norm_obj.vectors, norms[["nx", "ny", "nz"]].to_numpy(float), rtol=0, atol=1e-12
    )


@pytest.mark.parametrize("rotation", ROTATIONS_DEG)
def test_spec_example_observation_is_not_shifted(rotation):
    """The observation at anchor + (2, 3, 4) is at that point, not shifted by the origin."""
    anchor = ANCHORS[1]
    _, feature, _ = world_foliation_model(anchor, rotation)
    _, gradient_obj, _ = _split_objects(feature.get_data())
    expected = anchor + np.array([2.0, 3.0, 4.0])
    shifted = expected + anchor
    tol = tp(LENGTH, WIDTH)
    distance = np.linalg.norm(gradient_obj.locations - expected, axis=1)
    shifted_distance = np.linalg.norm(gradient_obj.locations - shifted, axis=1)
    assert distance.min() <= tol
    assert shifted_distance.min() > 1.0


@pytest.mark.parametrize("anchor_index", range(len(ANCHORS)))
def test_exported_data_keeps_world_locations(anchor_index, tmp_path):
    """The CSV files of the feature data hold the world locations."""
    anchor = ANCHORS[anchor_index]
    _, feature, data = world_foliation_model(anchor, 30.0)
    values, gradients, norms = _input_rows(data)
    expected = {
        "value": values[["X", "Y", "Z"]].to_numpy(float),
        "gradient": gradients[["X", "Y", "Z"]].to_numpy(float),
        "norm": norms[["X", "Y", "Z"]].to_numpy(float),
    }
    tol = tp(LENGTH, max(WIDTH, np.abs(anchor).max()))
    objects = feature.get_data()
    assert len(objects) == 3
    for index, obj in enumerate(objects):
        path = tmp_path / f"object_{index}.csv"
        obj.save(str(path))
        frame = pd.read_csv(path)
        kind = "value" if hasattr(obj, "values") else obj.name.rsplit("_", 1)[-1]
        np.testing.assert_allclose(
            frame[["x", "y", "z"]].to_numpy(float), expected[kind], rtol=0, atol=tol
        )


@pytest.mark.parametrize(
    "anchor_index, rotation",
    list(itertools.product(range(len(ANCHORS)), ROTATIONS_DEG)),
)
def test_intrusion_classification_uses_world_points(anchor_index, rotation):
    """The intrusion builder evaluates the model at the world points of the data."""
    from LoopStructural.modelling.intrusions.intrusion_frame_builder import (
        IntrusionFrameBuilder,
    )
    from tests.fixtures.coordinate_frames import model_with_frame

    anchor = ANCHORS[anchor_index]
    model = model_with_frame(anchor, anchor + np.array([20.0, 20.0, 20.0]), rotation)
    received = []

    def spy(xyz):
        received.append(np.array(xyz, copy=True))
        return np.zeros(len(xyz))

    model.evaluate_model = spy
    builder = object.__new__(IntrusionFrameBuilder)
    builder.model = model
    offsets = np.array(
        [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [2.0, 9.0, 4.0], [8.0, 3.0, 5.0]]
    )
    points = anchor + offsets
    data = pd.DataFrame({"X": points[:, 0], "Y": points[:, 1], "Z": points[:, 2]})
    builder._network_model_values(data)
    assert len(received) == 1
    assert np.array_equal(received[0], data[["X", "Y", "Z"]].to_numpy(float))
