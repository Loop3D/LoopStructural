"""Comprehensive tests for SurfeRBFInterpolator."""


import numpy as np
import pytest

try:
    import surfepy  # noqa: F401
    from loop_interpolation import SurfeRBFInterpolator
    HAS_SURFE = True
except ImportError:
    HAS_SURFE = False


@pytest.mark.skipif(not HAS_SURFE, reason="surfe not installed")
class TestSurfeRBFInterpolatorBasics:
    """Test basic SurfeRBFInterpolator functionality."""

    def test_surfe_rbf_creation(self):
        """Test creating a SurfeRBFInterpolator."""
        interpolator = SurfeRBFInterpolator()
        assert interpolator is not None

    def test_surfe_rbf_has_surfe_library(self):
        """Test that interpolator has access to surfe library."""
        interpolator = SurfeRBFInterpolator()
        assert hasattr(interpolator, "surfe")
        assert interpolator.surfe is not None


@pytest.mark.skipif(not HAS_SURFE, reason="surfe not installed")
class TestSurfeRBFEvaluateValue:
    """Test SurfeRBFInterpolator.evaluate_value method."""

    def test_evaluate_value_with_valid_points(self):
        """Test evaluate_value with valid coordinate points."""
        interpolator = SurfeRBFInterpolator()

        # Mock the surfe evaluation
        def mock_evaluate(points):
            return np.ones(points.shape[0]) * 42.0

        interpolator.surfe.EvaluateInterpolantAtPoints = mock_evaluate

        points = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]])
        result = interpolator.evaluate_value(points)

        assert result is not None
        assert result.shape == (3,)
        assert np.all(result == 42.0)

    def test_evaluate_value_with_nan_coordinates(self):
        """Test evaluate_value properly handles NaN coordinates."""
        interpolator = SurfeRBFInterpolator()

        def mock_evaluate(points):
            return np.ones(points.shape[0]) * 10.0

        interpolator.surfe.EvaluateInterpolantAtPoints = mock_evaluate

        points = np.array([
            [0.0, 0.0, 0.0],
            [np.nan, 1.0, 1.0],
            [2.0, 2.0, 2.0],
        ])

        result = interpolator.evaluate_value(points)

        assert result is not None
        # Valid points should have values from surfe
        assert result[0] == 10.0
        assert result[2] == 10.0
        # NaN input should produce NaN output
        assert np.isnan(result[1])

    def test_evaluate_value_all_nan_input(self):
        """Test evaluate_value with all NaN coordinates."""
        interpolator = SurfeRBFInterpolator()

        def mock_evaluate(points):
            if points.shape[0] == 0:
                return np.array([])
            return np.ones(points.shape[0]) * 5.0

        interpolator.surfe.EvaluateInterpolantAtPoints = mock_evaluate

        points = np.array([[np.nan, np.nan, np.nan]])
        result = interpolator.evaluate_value(points)

        assert result is not None
        assert np.isnan(result[0])

    def test_evaluate_value_empty_input(self):
        """Test evaluate_value with empty array."""
        interpolator = SurfeRBFInterpolator()

        def mock_evaluate(points):
            if points.shape[0] == 0:
                return np.array([])
            return np.ones(points.shape[0]) * 5.0

        interpolator.surfe.EvaluateInterpolantAtPoints = mock_evaluate

        points = np.array([]).reshape(0, 3)
        result = interpolator.evaluate_value(points)

        assert result is not None
        assert result.shape == (0,)

    def test_evaluate_value_single_point(self):
        """Test evaluate_value with single point."""
        interpolator = SurfeRBFInterpolator()

        def mock_evaluate(points):
            return np.array([99.0])

        interpolator.surfe.EvaluateInterpolantAtPoints = mock_evaluate

        points = np.array([[1.0, 2.0, 3.0]])
        result = interpolator.evaluate_value(points)

        assert result.shape == (1,)
        assert result[0] == 99.0


@pytest.mark.skipif(not HAS_SURFE, reason="surfe not installed")
class TestSurfeRBFEvaluateGradient:
    """Test SurfeRBFInterpolator.evaluate_gradient method."""

    def test_evaluate_gradient_returns_array(self):
        """Test that evaluate_gradient returns an array, not None."""
        interpolator = SurfeRBFInterpolator()

        def mock_evaluate_vector(points):
            return np.ones((points.shape[0], 3)) * 2.5

        interpolator.surfe.EvaluateVectorInterpolantAtPoints = mock_evaluate_vector

        points = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
        result = interpolator.evaluate_gradient(points)

        # Bug fix: should return array, not None
        assert result is not None
        assert isinstance(result, np.ndarray)

    def test_evaluate_gradient_shape(self):
        """Test that evaluate_gradient returns correct shape."""
        interpolator = SurfeRBFInterpolator()

        def mock_evaluate_vector(points):
            return np.ones((points.shape[0], 3)) * 3.0

        interpolator.surfe.EvaluateVectorInterpolantAtPoints = mock_evaluate_vector

        points = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]])
        result = interpolator.evaluate_gradient(points)

        assert result.shape == (3, 3), "Should return (n_points, 3)"

    def test_evaluate_gradient_with_nan_coordinates(self):
        """Test evaluate_gradient handles NaN coordinates."""
        interpolator = SurfeRBFInterpolator()

        def mock_evaluate_vector(points):
            return np.ones((points.shape[0], 3)) * 0.5

        interpolator.surfe.EvaluateVectorInterpolantAtPoints = mock_evaluate_vector

        points = np.array([
            [0.0, 0.0, 0.0],
            [np.nan, 1.0, 1.0],
            [np.nan, np.nan, np.nan],
        ])

        result = interpolator.evaluate_gradient(points)

        assert result is not None
        # Valid point should have values
        assert np.allclose(result[0], 0.5)
        # NaN inputs should produce NaN outputs
        assert np.all(np.isnan(result[1]))
        assert np.all(np.isnan(result[2]))

    def test_evaluate_gradient_single_point(self):
        """Test evaluate_gradient with single point."""
        interpolator = SurfeRBFInterpolator()

        gradient_value = np.array([[1.0, 2.0, 3.0]])

        def mock_evaluate_vector(points):
            return gradient_value

        interpolator.surfe.EvaluateVectorInterpolantAtPoints = mock_evaluate_vector

        points = np.array([[5.0, 6.0, 7.0]])
        result = interpolator.evaluate_gradient(points)

        assert result.shape == (1, 3)
        assert np.allclose(result[0], [1.0, 2.0, 3.0])

    def test_evaluate_gradient_empty_input(self):
        """Test evaluate_gradient with empty array."""
        interpolator = SurfeRBFInterpolator()

        def mock_evaluate_vector(points):
            if points.shape[0] == 0:
                return np.empty((0, 3))
            return np.ones((points.shape[0], 3))

        interpolator.surfe.EvaluateVectorInterpolantAtPoints = mock_evaluate_vector

        points = np.array([]).reshape(0, 3)
        result = interpolator.evaluate_gradient(points)

        assert result is not None
        assert result.shape == (0, 3)


@pytest.mark.skipif(not HAS_SURFE, reason="surfe not installed")
class TestSurfeRBFSetup:
    """Test SurfeRBFInterpolator setup and configuration."""

    def test_interpolator_type(self):
        """Test that interpolator has correct type."""
        interpolator = SurfeRBFInterpolator()
        assert interpolator.type is not None

    def test_interpolator_is_rbf(self):
        """Test that SurfeRBFInterpolator is recognized as RBF type."""
        interpolator = SurfeRBFInterpolator()
        # The type should indicate it's a surfe/RBF interpolator
        assert "Surfe" in str(type(interpolator).__name__) or "RBF" in str(
            interpolator.type
        )


@pytest.mark.skipif(not HAS_SURFE, reason="surfe not installed")
class TestSurfeRBFNaNMaskingCorrectness:
    """Test NaN masking implementation correctness."""

    def test_nan_masking_uses_isnan_not_comparison(self):
        """Test that NaN detection uses np.isnan, not direct comparison."""
        interpolator = SurfeRBFInterpolator()

        call_log = []

        def mock_evaluate(points):
            call_log.append(points.copy())
            return np.ones(points.shape[0])

        interpolator.surfe.EvaluateInterpolantAtPoints = mock_evaluate

        # Points with various NaN patterns
        points = np.array([
            [0.0, 0.0, 0.0],
            [np.nan, 0.0, 0.0],
            [0.0, np.nan, 0.0],
            [0.0, 0.0, np.nan],
            [np.nan, np.nan, 0.0],
            [np.nan, np.nan, np.nan],
        ])

        result = interpolator.evaluate_value(points)

        # Surfe should only be called with valid (non-NaN) rows
        surfe_input = call_log[0]
        # Row with all NaN should be filtered
        assert surfe_input.shape[0] == 1
        assert np.allclose(surfe_input[0], [0.0, 0.0, 0.0])

    def test_nan_gradient_masking_uses_isnan(self):
        """Test that gradient NaN masking uses np.isnan."""
        interpolator = SurfeRBFInterpolator()

        call_log = []

        def mock_evaluate_vector(points):
            call_log.append(points.copy())
            return np.ones((points.shape[0], 3))

        interpolator.surfe.EvaluateVectorInterpolantAtPoints = mock_evaluate_vector

        points = np.array([
            [0.0, 0.0, 0.0],
            [np.nan, 0.0, 0.0],
        ])

        result = interpolator.evaluate_gradient(points)

        # Surfe should only receive valid point
        surfe_input = call_log[0]
        assert surfe_input.shape[0] == 1
        # The result should have NaN for the second point
        assert np.all(np.isnan(result[1]))


@pytest.mark.skipif(not HAS_SURFE, reason="surfe not installed")
class TestSurfeRBFEdgeCases:
    """Test edge cases and error handling."""

    def test_very_large_coordinate_values(self):
        """Test with very large coordinate values."""
        interpolator = SurfeRBFInterpolator()

        def mock_evaluate(points):
            return np.ones(points.shape[0]) * 1e-10

        interpolator.surfe.EvaluateInterpolantAtPoints = mock_evaluate

        points = np.array([[1e10, 1e10, 1e10], [1e-10, 1e-10, 1e-10]])
        result = interpolator.evaluate_value(points)

        assert result is not None
        assert result.shape == (2,)

    def test_mixed_valid_and_invalid_coordinates(self):
        """Test with mixed valid and invalid coordinates."""
        interpolator = SurfeRBFInterpolator()

        def mock_evaluate(points):
            return np.arange(points.shape[0]) * 1.0

        interpolator.surfe.EvaluateInterpolantAtPoints = mock_evaluate

        points = np.array([
            [0.0, 0.0, 0.0],
            [np.inf, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ])

        result = interpolator.evaluate_value(points)

        # Should handle infinity as a real value (surfe will handle it)
        assert result is not None

    def test_partial_nan_rows(self):
        """Test rows with partial NaN values."""
        interpolator = SurfeRBFInterpolator()

        def mock_evaluate(points):
            return np.ones(points.shape[0])

        interpolator.surfe.EvaluateInterpolantAtPoints = mock_evaluate

        points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, np.nan, 2.0],  # One NaN in the row
            [3.0, 4.0, 5.0],
        ])

        result = interpolator.evaluate_value(points)

        # Point with partial NaN should still be detected as invalid
        assert np.isnan(result[1])
        assert not np.isnan(result[0])
        assert not np.isnan(result[2])


@pytest.mark.skipif(not HAS_SURFE, reason="surfe not installed")
class TestSurfeRBFGradientConsistency:
    """Test consistency between value and gradient methods."""

    def test_gradient_and_value_handle_nans_same_way(self):
        """Test that NaN handling is consistent between value and gradient."""
        interpolator = SurfeRBFInterpolator()

        def mock_evaluate(points):
            return np.ones(points.shape[0]) * 42.0

        def mock_evaluate_vector(points):
            return np.ones((points.shape[0], 3)) * 10.0

        interpolator.surfe.EvaluateInterpolantAtPoints = mock_evaluate
        interpolator.surfe.EvaluateVectorInterpolantAtPoints = mock_evaluate_vector

        points = np.array([
            [0.0, 0.0, 0.0],
            [np.nan, 1.0, 1.0],
            [2.0, 2.0, 2.0],
        ])

        value_result = interpolator.evaluate_value(points)
        gradient_result = interpolator.evaluate_gradient(points)

        # Both should have NaN at same indices
        assert np.isnan(value_result[1])
        assert np.all(np.isnan(gradient_result[1]))

        # Both should be valid for same valid points
        assert not np.isnan(value_result[0])
        assert not np.all(np.isnan(gradient_result[0]))
