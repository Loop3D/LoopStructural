"""Tests that consumers of model data show and export world locations."""

import itertools

import numpy as np
import pandas as pd
import pytest

from tests.fixtures.coordinate_frames import (
    ANCHORS,
    ROTATIONS_DEG,
    model_with_frame,
    tp,
    world_foliation_model,
)
from tests.unit.modelling.test_fault_frame_analytic import _fault

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


@pytest.mark.parametrize(
    "anchor_index, rotation",
    list(itertools.product((0, 3), (0, 30))),
)
def test_2d_viewer_fault_ellipse_centre_is_world(anchor_index, rotation):
    """The 2D viewer draws the fault ellipse at the world centre of the fault."""
    anchor = ANCHORS[anchor_index]
    # Skip the test if the viewer package is not installed.
    viewer = pytest.importorskip("loopstructuralvisualisation._2d_viewer")
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fault = _fault(anchor, rotation, "PLI")
    fig, ax = plt.subplots()
    try:
        view = type("View", (), {})()
        view.ax = ax
        view.model = None
        viewer.Loop2DView.add_fault_ellipse(view, faults=[fault])
        assert len(ax.patches) == 1
        tol = tp(LENGTH, max(WIDTH, np.abs(anchor).max()))
        np.testing.assert_allclose(ax.patches[0].center, fault.fault_centre[:2], rtol=0, atol=tol)
    finally:
        plt.close(fig)


_OFFSETS = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [2.0, 9.0, 4.0]])


def _deprecated_call_cell(anchor_index, rotation, method_name, target_name, match):
    """Check that a deprecated model method warns and keeps its result."""
    anchor = ANCHORS[anchor_index]
    model = model_with_frame(anchor, anchor + np.array([20.0, 20.0, 20.0]), rotation)
    tol = tp(LENGTH, max(WIDTH, np.abs(anchor).max()))
    points = anchor + _OFFSETS
    expected = getattr(model.bounding_box, target_name)(points.copy())
    with pytest.warns(DeprecationWarning, match=match) as record:
        result = getattr(model, method_name)(points.copy())
    assert len(record) == 1
    assert record[0].filename == __file__
    np.testing.assert_allclose(result, expected, rtol=0, atol=tol)
    # With inplace=True the method changes the array that the caller passes.
    expected_inplace = points.copy()
    getattr(model.bounding_box, target_name)(expected_inplace, inplace=True)
    passed = points.copy()
    with pytest.warns(DeprecationWarning, match=match):
        getattr(model, method_name)(passed, inplace=True)
    np.testing.assert_allclose(passed, expected_inplace, rtol=0, atol=tol)


@pytest.mark.parametrize(
    "anchor_index, rotation",
    list(itertools.product(range(len(ANCHORS)), ROTATIONS_DEG)),
)
def test_rescale_warns_and_keeps_behaviour(anchor_index, rotation):
    """model.rescale gives a DeprecationWarning and the result of reproject."""
    _deprecated_call_cell(anchor_index, rotation, "rescale", "reproject", "bounding_box.reproject")


@pytest.mark.parametrize(
    "anchor_index, rotation",
    list(itertools.product(range(len(ANCHORS)), ROTATIONS_DEG)),
)
def test_scale_warns_and_keeps_behaviour(anchor_index, rotation):
    """model.scale gives a DeprecationWarning and the result of project."""
    _deprecated_call_cell(anchor_index, rotation, "scale", "project", "bounding_box.project")


def test_internal_callers_do_not_warn(tmp_path):
    """No library caller uses the deprecated rescale or scale methods."""
    import warnings

    from LoopStructural.modelling.core._model_exporter import ModelExporter
    from LoopStructural.modelling.intrusions.intrusion_frame_builder import (
        IntrusionFrameBuilder,
    )

    anchor = ANCHORS[1]
    model, feature, _ = world_foliation_model(anchor, 30.0)
    # The model is not solved. Use a stub for the model evaluation.
    model.evaluate_model = lambda xyz: np.zeros(len(xyz))
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=".*(rescale|scale) is deprecated.*")
        feature.get_data()
        ModelExporter.get_block_model(model)
        builder = object.__new__(IntrusionFrameBuilder)
        builder.model = model
        points = anchor + _OFFSETS
        data = pd.DataFrame({"X": points[:, 0], "Y": points[:, 1], "Z": points[:, 2]})
        builder._network_model_values(data)
        # Call the exporter functions in the way that tests/unit/io/test_exporters.py does.
        pytest.importorskip("pyevtk")
        from LoopStructural.export import exporters
        from LoopStructural.export.file_formats import FileFormat
        from tests.unit.io import test_exporters as fakes

        class _Surfaces(fakes._FakeModelForSurfaces):
            def rescale(self, points, *, inplace=False):
                return model.rescale(points, inplace=inplace)

        class _Volume(fakes._FakeModelForVolume):
            def rescale(self, points, *, inplace=False):
                return model.rescale(points, inplace=inplace)

        exporters.write_feat_surfs(_Surfaces(), "strati", file_format=FileFormat.NUMPY)
        exporters.write_cubeface(
            _Volume(), str(tmp_path / "cube"), "label", np.array([5, 5, 5]), FileFormat.VTK
        )
        exporters.write_vol(
            _Volume(), str(tmp_path / "vol"), "label", np.array([5, 5, 5]), FileFormat.VTK
        )
