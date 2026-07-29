"""Regression test for NaN constraints being skipped (P0 fix)."""

import pytest
import numpy as np
from unittest.mock import Mock, patch

from loop_interpolation._discrete_interpolator import DiscreteInterpolator


def test_add_constraints_to_least_squares_skips_nan_constraints():
    """Test that constraints with NaN values are actually skipped, not added to the system."""
    # Create a mock DiscreteInterpolator (can't instantiate abstract class directly)
    interpolator = Mock(spec=DiscreteInterpolator)
    interpolator.constraints = {}
    interpolator.dof = 100
    interpolator.n_nodes = 10

    # Get the real method (not mocked)
    from loop_interpolation._discrete_interpolator import DiscreteInterpolator as RealDI
    add_method = RealDI.add_constraints_to_least_squares.__get__(interpolator, type(interpolator))

    # Test case 1: constraint with NaN in the data points.
    # idc must match A's shape - it gives the global dof index of every
    # entry of A (e.g. one column per local node of an element), not a
    # single index per row.
    idc = np.array(
        [[0, 1, 2], [3, 4, 5], [6, 7, np.nan], [9, 10, 11]], dtype=float
    )
    A = np.ones((4, 3))
    B = np.array([1.0, 2.0, 3.0, 4.0])

    # add_constraints_to_least_squares logs via the module-level `logger`,
    # not a per-instance attribute, so patch that directly.
    with patch("loop_interpolation._discrete_interpolator.logger") as mock_logger:
        add_method(A, B, idc, w=1.0, name="test_nan_constraint")

        # Constraint should NOT be added to self.constraints due to NaN
        assert "test_nan_constraint" not in interpolator.constraints
        assert mock_logger.warning.called  # Should log the warning

        # Test case 2: valid constraint (no NaN)
        interpolator.constraints.clear()
        mock_logger.reset_mock()

        idc_valid = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]], dtype=float)
        A_valid = np.ones((4, 3))
        B_valid = np.array([1.0, 2.0, 3.0, 4.0])

        add_method(A_valid, B_valid, idc_valid, w=1.0, name="test_valid_constraint")

        # Constraint SHOULD be added for valid data
        assert "test_valid_constraint" in interpolator.constraints
        assert not mock_logger.warning.called  # No warning for valid data


def test_add_constraints_to_least_squares_skips_nan_in_matrix():
    """Test that constraints with NaN in matrix A are skipped."""
    interpolator = Mock(spec=DiscreteInterpolator)
    interpolator.constraints = {}
    interpolator.dof = 100
    interpolator.n_nodes = 10

    from loop_interpolation._discrete_interpolator import DiscreteInterpolator as RealDI
    add_method = RealDI.add_constraints_to_least_squares.__get__(interpolator, type(interpolator))

    # idc must match A's shape (see comment in the previous test)
    idc = np.array([[0, 1], [2, 3], [4, 5], [6, 7]], dtype=float)
    A_with_nan = np.array([[1.0, 2.0], [3.0, np.nan], [5.0, 6.0], [7.0, 8.0]])  # NaN in matrix
    B = np.array([1.0, 2.0, 3.0, 4.0])

    with patch("loop_interpolation._discrete_interpolator.logger") as mock_logger:
        add_method(A_with_nan, B, idc, w=1.0, name="test_nan_matrix")

        # Constraint should NOT be added
        assert "test_nan_matrix" not in interpolator.constraints
        assert mock_logger.warning.called


def test_add_constraints_to_least_squares_skips_nan_in_b():
    """Test that constraints with NaN in vector B are skipped."""
    interpolator = Mock(spec=DiscreteInterpolator)
    interpolator.constraints = {}
    interpolator.dof = 100
    interpolator.n_nodes = 10

    from loop_interpolation._discrete_interpolator import DiscreteInterpolator as RealDI
    add_method = RealDI.add_constraints_to_least_squares.__get__(interpolator, type(interpolator))

    # idc must match A's shape (see comment in the first test)
    idc = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]], dtype=float)
    A = np.ones((4, 3))
    B_with_nan = np.array([1.0, np.nan, 3.0, 4.0])  # NaN in B

    with patch("loop_interpolation._discrete_interpolator.logger") as mock_logger:
        add_method(A, B_with_nan, idc, w=1.0, name="test_nan_b_vector")

        # Constraint should NOT be added
        assert "test_nan_b_vector" not in interpolator.constraints
        assert mock_logger.warning.called
