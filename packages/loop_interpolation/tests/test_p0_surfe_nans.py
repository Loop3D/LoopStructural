"""Regression tests for SurfeRBFInterpolator NaN-handling bugs (P0 fixes)."""

import pytest
import numpy as np

pytest.importorskip("surfe", minversion=None)

from loop_interpolation import SurfeRBFInterpolator


def test_surfe_evaluate_value_masks_nan_coordinates():
    """Test that evaluate_value properly masks out NaN coordinates instead of passing them to surfe."""
    # Create a simple interpolator with dummy data
    # (This is a minimal test to verify the NaN-masking logic; full integration test would need surfe setup)
    interpolator = SurfeRBFInterpolator()

    # Mock out surfe methods to verify masking behavior
    call_count = [0]

    def mock_evaluate(points):
        call_count[0] += 1
        return np.ones(points.shape[0])

    interpolator.surfe.EvaluateInterpolantAtPoints = mock_evaluate

    # Evaluate with some NaN coordinates
    eval_points = np.array([
        [0.0, 0.0, 0.0],  # valid
        [1.0, 1.0, 1.0],  # valid
        [np.nan, np.nan, np.nan],  # invalid
    ])

    result = interpolator.evaluate_value(eval_points)

    # Should have called surfe only with 2 valid points (not 3)
    assert call_count[0] == 1
    # Result should have NaN for the invalid point
    assert np.isnan(result[2])
    # Valid points should have surfe's return value
    assert result[0] == 1.0
    assert result[1] == 1.0


def test_surfe_evaluate_gradient_returns_array_not_none():
    """Test that evaluate_gradient returns the evaluated array, not None."""
    interpolator = SurfeRBFInterpolator()

    # Mock surfe's vector evaluation
    def mock_evaluate_vector(points):
        return np.ones((points.shape[0], 3)) * 0.5

    interpolator.surfe.EvaluateVectorInterpolantAtPoints = mock_evaluate_vector

    eval_points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0],
    ])

    result = interpolator.evaluate_gradient(eval_points)

    # Result should be an array, not None
    assert result is not None
    assert isinstance(result, np.ndarray)
    assert result.shape == (2, 3)
    # Valid points should have surfe's return value
    assert np.allclose(result[0], 0.5)
    assert np.allclose(result[1], 0.5)


def test_surfe_evaluate_gradient_nan_coordinates_become_nan_output():
    """Test that NaN input coordinates produce NaN output in evaluate_gradient."""
    interpolator = SurfeRBFInterpolator()

    # Mock surfe's vector evaluation
    def mock_evaluate_vector(points):
        return np.ones((points.shape[0], 3)) * 0.5

    interpolator.surfe.EvaluateVectorInterpolantAtPoints = mock_evaluate_vector

    eval_points = np.array([
        [0.0, 0.0, 0.0],  # valid
        [np.nan, 1.0, 1.0],  # invalid (one NaN)
        [np.nan, np.nan, np.nan],  # invalid (all NaN)
    ])

    result = interpolator.evaluate_gradient(eval_points)

    # Valid point should have surfe's return value
    assert np.allclose(result[0], 0.5)
    # Invalid points should be NaN
    assert np.all(np.isnan(result[1]))
    assert np.all(np.isnan(result[2]))
