"""Regression tests for fold and fault audit fixes."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from LoopStructural import GeologicalModel
from LoopStructural.modelling.features.fault import FaultSegment
from LoopStructural.modelling.features.fold import FoldFrame
from LoopStructural.modelling.features.fold.fold_function._fourier_series_fold_rotation_angle import (
    FourierSeriesFoldRotationAngleProfile,
)


class _PlaneFeature:
    def __init__(self, gradient):
        self.gradient = np.array(gradient, dtype=float)

    def evaluate_value(self, points):
        return points @ self.gradient

    def evaluate_gradient(self, points):
        return np.tile(self.gradient, (points.shape[0], 1))


@pytest.mark.parametrize("limb", [0.0, 45.0, 70.0])
def test_fold_axis_rotation_not_biased_by_limb_dip(limb):
    theta, phi = np.deg2rad(60.0), np.deg2rad(limb)
    lineation = [np.sin(theta) * np.cos(phi), np.cos(theta) * np.cos(phi), np.sin(phi)]
    builder = SimpleNamespace(
        get_gradient_constraints=lambda: np.zeros((0, 6)),
        get_norm_constraints=lambda: np.array([[0.0, 0.0, 0.0, *lineation]]),
    )
    frame = SimpleNamespace(features=[_PlaneFeature([0, 0, 1]), _PlaneFeature([0, 1, 0])])
    far, _ = FoldFrame.calculate_fold_axis_rotation(frame, builder)
    np.testing.assert_allclose(far, [60.0 - 90.0], atol=1e-8)


def test_fold_profile_fit_ignores_one_nan_row():
    s = np.linspace(0, 1000, 50)
    angle = 20 + 30 * np.cos(2 * np.pi / 500 * s)
    angle[10] = np.nan
    profile = FourierSeriesFoldRotationAngleProfile(angle, s)
    assert profile.fit(params={"guess": [0, 1, 1, 500]})


def _fault_model(**kwargs):
    model = GeologicalModel([0, 0, 0], [1, 1, 1])
    model.data = pd.DataFrame(
        [
            [0.5, 0.5, 0.5, 1, 0, 0, 0, "fault", 0],
            [0.5, 0.5, 0.5, 0, 0, 1, 1, "fault", 0],
            [0.5, 0.5, 0.5, 0, 1, 0, 2, "fault", 0],
            [0.2, 0.5, 0.5, 0, 0, 1, 0, "strati", 0.5],
            [0.8, 0.5, 0.5, 0, 0, 1, 0, "strati", 0.5],
        ],
        columns=["X", "Y", "Z", "nx", "ny", "nz", "coord", "feature_name", "val"],
    )
    fault = model.create_and_add_fault("fault", 0.1, nelements=1e3, **kwargs)
    return model, fault


def test_fault_pitch_rotates_strike_in_fault_plane():
    # The plane dips 45 degrees to +x: strike is along y, down-dip is (1, 0, -1).
    normal = np.array([1.0, 0.0, 1.0]) / np.sqrt(2)
    _, fault = _fault_model(fault_normal_vector=normal, fault_pitch=90)
    down_dip = np.array([1.0, 0.0, -1.0]) / np.sqrt(2)
    np.testing.assert_allclose(fault.builder.fault_slip_vector, down_dip, atol=1e-8)
    _, fault = _fault_model(fault_normal_vector=normal, fault_pitch=0)
    np.testing.assert_allclose(np.abs(fault.builder.fault_slip_vector), [0, 1, 0], atol=1e-8)


def test_faulted_gradient_on_fault_is_nan_and_finishes():
    model, _ = _fault_model()
    strati = model.create_and_add_foliation("strati", nelements=1e3)
    points = np.array([[0.5, 0.5, 0.5], [0.2, 0.5, 0.5]])
    g = strati.evaluate_gradient(points)
    assert np.all(np.isnan(g[0]))
    assert np.all(np.isfinite(g[1]))


def test_faultfunction_setter_accepts_callable_rejects_object():
    segment = FaultSegment.__new__(FaultSegment)
    segment.faultfunction = lambda x: x
    assert callable(segment.faultfunction)
    segment.faultfunction = "BaseFault"
    with pytest.raises(ValueError):
        segment.faultfunction = object()
