"""Regression tests: NaN constraint rows are dropped, the rest of the family is kept."""

from unittest.mock import Mock, patch

import numpy as np
from loop_interpolation._discrete_interpolator import DiscreteInterpolator


def test_add_constraints_to_least_squares_drops_nan_rows():
    """Test that rows with NaN values are dropped and the other rows are added."""
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

        # the nan row is dropped, the other three are added
        assert interpolator.constraints["test_nan_constraint"]["matrix"].shape[0] == 3
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


def test_add_constraints_to_least_squares_drops_nan_row_in_matrix():
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

        # the nan row is dropped
        assert interpolator.constraints["test_nan_matrix"]["matrix"].shape[0] == 3
        assert mock_logger.warning.called


def test_add_constraints_to_least_squares_drops_nan_row_in_b():
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

        # the nan row is dropped
        assert interpolator.constraints["test_nan_b_vector"]["matrix"].shape[0] == 3
        assert mock_logger.warning.called


def test_nan_row_filter_keeps_weights_aligned():
    interpolator = Mock(spec=DiscreteInterpolator)
    interpolator.constraints = {}
    interpolator.dof = 100
    add_method = DiscreteInterpolator.add_constraints_to_least_squares.__get__(
        interpolator, type(interpolator)
    )
    idc = np.arange(8).reshape(4, 2)
    A = np.array([[1.0, 0.0], [np.nan, 1.0], [0.0, 2.0], [3.0, 0.0]])
    add_method(A, np.array([1.0, 2.0, 3.0, 4.0]), idc, w=np.array([1.0, 2.0, 3.0, 4.0]), name="fam")
    fam = interpolator.constraints["fam"]
    assert fam["matrix"].shape == (3, 100)
    np.testing.assert_array_equal(fam["w"], [1.0, 3.0, 4.0])
    np.testing.assert_allclose(fam["b"], [1.0, 3.0 / 2.0, 4.0 / 3.0])


def test_fold_constraints_drop_only_nan_elements():
    """Fold families go through the same filter: one nan element drops one row."""
    from tests.golden import engines, geomodels

    case = geomodels.build("fold_sparse")
    interp = engines.make_interpolator("DFI", case)
    engines.add_data(interp, case)
    original = interp.fold.get_deformed_orientation

    def with_nan(points):
        ori, axis, dgz = (np.array(a, dtype=float) for a in original(points))
        ori[0, :] = np.nan
        return ori, axis, dgz

    interp.fold.get_deformed_orientation = with_nan
    interp.setup_interpolator()
    n = interp.support.n_elements
    assert interp.constraints["fold orientation"]["matrix"].shape[0] == n - 1
    assert interp.constraints["fold axis"]["matrix"].shape[0] == n
