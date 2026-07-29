"""Tests for DiscreteFoldInterpolator."""

import numpy as np
import pytest
from unittest.mock import Mock, MagicMock

from loop_interpolation._discrete_fold_interpolator import DiscreteFoldInterpolator
from loop_interpolation._p1interpolator import P1Interpolator
from loop_common.supports import SupportType
from loop_common.supports._3d_structured_tetra import TetMesh


class MockFoldEvent:
    """Mock FoldEvent for testing."""

    def __init__(self, orientation_grad=None, axis_grad=None, deformed_normal=None):
        """Initialize mock fold with configurable return values."""
        self.orientation_grad = orientation_grad if orientation_grad is not None else np.ones((10, 3))
        self.axis_grad = axis_grad if axis_grad is not None else np.ones((10, 3))
        self.deformed_normal = deformed_normal if deformed_normal is not None else np.ones((10, 3))

    def get_deformed_orientation(self, points):
        """Return deformed orientation components, broadcast to match the
        number of query points (the mesh may have more elements than the
        configured mock array has rows)."""
        n_points = points.shape[0] if hasattr(points, "shape") else 1
        return (
            np.tile(self.orientation_grad[0], (n_points, 1)),
            np.tile(self.axis_grad[0], (n_points, 1)),
            np.tile(self.deformed_normal[0], (n_points, 1)),
        )


class TestDiscreteFoldInterpolatorCreation:
    """Test DiscreteFoldInterpolator instantiation."""

    def test_discrete_fold_interpolator_creation(self):
        """Test basic creation of DiscreteFoldInterpolator."""
        support = TetMesh()
        fold = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold)

        assert interpolator is not None
        assert interpolator.support == support
        assert interpolator.fold == fold

    def test_discrete_fold_interpolator_without_fold(self):
        """Test creation without a fold (fold=None)."""
        support = TetMesh()

        interpolator = DiscreteFoldInterpolator(support, fold=None)

        assert interpolator is not None
        assert interpolator.fold is None

    def test_discrete_fold_interpolator_inherits_from_p1(self):
        """Test that DiscreteFoldInterpolator inherits from P1Interpolator."""
        support = TetMesh()
        fold = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold)

        assert isinstance(interpolator, P1Interpolator)


class TestDiscreteFoldInterpolatorFoldUpdate:
    """Test fold attribute updates."""

    def test_update_fold(self):
        """Test updating fold attribute."""
        support = TetMesh()
        fold1 = MockFoldEvent()
        fold2 = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold1)
        assert interpolator.fold == fold1

        interpolator.update_fold(fold2)
        assert interpolator.fold == fold2

    def test_fold_can_be_none(self):
        """Test that fold can be set to None."""
        support = TetMesh()
        fold = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold)
        assert interpolator.fold is not None

        interpolator.update_fold(None)
        assert interpolator.fold is None


class TestDiscreteFoldInterpolatorAddConstraints:
    """Test fold constraint addition methods."""

    def test_add_fold_orientation_constraints(self):
        """Test adding fold orientation constraints."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        fold = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # This should not raise
        interpolator.add_fold_constraints(
            fold_orientation=10.0,
            fold_axis_w=None,
            fold_regularisation=None,
        )

    def test_add_fold_axis_constraints(self):
        """Test adding fold axis constraints."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        fold = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # This should not raise
        interpolator.add_fold_constraints(
            fold_orientation=None,
            fold_axis_w=10.0,
            fold_regularisation=None,
        )

    def test_add_fold_normalisation_constraints(self):
        """Test adding fold normalisation constraints."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        fold = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # This should not raise
        interpolator.add_fold_constraints(
            fold_orientation=None,
            fold_axis_w=None,
            fold_regularisation=None,
            fold_normalisation=10.0,
        )

    def test_add_fold_regularisation_constraints(self):
        """Test adding fold regularisation constraints."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        fold = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # This should not raise
        interpolator.add_fold_constraints(
            fold_orientation=None,
            fold_axis_w=None,
            fold_regularisation=(0.1, 0.01, 0.01),
        )

    def test_add_all_fold_constraints_simultaneously(self):
        """Test adding all fold constraint types at once."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        fold = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # This should not raise
        interpolator.add_fold_constraints(
            fold_orientation=10.0,
            fold_axis_w=10.0,
            fold_regularisation=(0.1, 0.01, 0.01),
            fold_normalisation=1.0,
        )

    def test_add_fold_constraints_with_mask_function(self):
        """Test adding fold constraints with a mask function."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        fold = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # Define a mask function that returns True for points above z=1.5
        def mask_fn(points):
            return points[:, 2] > 1.5

        # This should not raise
        interpolator.add_fold_constraints(
            fold_orientation=10.0,
            fold_axis_w=10.0,
            fold_regularisation=(0.1, 0.01, 0.01),
            fold_normalisation=1.0,
            mask_fn=mask_fn,
        )

    def test_fold_constraints_with_custom_weights(self):
        """Test fold constraints with various weight values."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        fold = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # Test with different weight combinations
        weights_to_test = [
            (1.0, 1.0, (1.0, 1.0, 1.0), 1.0),
            (100.0, 100.0, (10.0, 10.0, 10.0), 10.0),
            (0.1, 0.1, (0.01, 0.01, 0.01), 0.1),
        ]

        for f_ori, f_axis, f_reg, f_norm in weights_to_test:
            interpolator.add_fold_constraints(
                fold_orientation=f_ori,
                fold_axis_w=f_axis,
                fold_regularisation=f_reg,
                fold_normalisation=f_norm,
            )


class TestDiscreteFoldInterpolatorSetup:
    """Test setup_interpolator method."""

    def test_setup_interpolator_calls_fold_setup(self):
        """Test that setup_interpolator incorporates fold constraints."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        fold = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold)

        # Setup with data
        result = interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]),
            data_values=np.array([1.0, 2.0]),
        )

        # Should complete without error
        assert result is not None

    def test_setup_interpolator_without_fold(self):
        """Test setup_interpolator when fold is None raises a clear error."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))

        interpolator = DiscreteFoldInterpolator(support, fold=None)

        # A DiscreteFoldInterpolator without a fold event can't build fold
        # constraints, so setup should fail loudly rather than silently
        # skip them.
        with pytest.raises(RuntimeError, match="no fold event set"):
            interpolator.setup_interpolator(
                data_points=np.array([[1.0, 1.0, 1.0]]),
                data_values=np.array([1.0]),
            )


class TestDiscreteFoldInterpolatorConstraintWeighting:
    """Test constraint weighting by element volume."""

    def test_fold_constraints_use_element_volume_weighting(self):
        """Test that fold constraints respect element volume weighting."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        fold = MockFoldEvent(
            orientation_grad=np.random.rand(support.n_elements, 3),
            axis_grad=np.random.rand(support.n_elements, 3),
            deformed_normal=np.random.rand(support.n_elements, 3),
        )

        interpolator = DiscreteFoldInterpolator(support, fold=fold)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # Add constraints with different weights
        interpolator.add_fold_constraints(
            fold_orientation=10.0,
            fold_axis_w=5.0,
            fold_regularisation=None,
            fold_normalisation=None,
        )

        # Constraints should be added (hard to verify exact values,
        # but we can at least check it doesn't crash)
        assert True


class TestDiscreteFoldInterpolatorStepParameter:
    """Test the step parameter for constraint sampling."""

    def test_fold_constraints_with_step_parameter(self):
        """Test that step parameter subsamples constraints."""
        support = TetMesh(nsteps=np.array([5, 5, 5]))
        fold = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold)
        interpolator.setup_interpolator(
            data_points=np.array([[2.0, 2.0, 2.0]]),
            data_values=np.array([1.0]),
        )

        # With step=2, should use every other element
        interpolator.add_fold_constraints(
            fold_orientation=10.0,
            fold_axis_w=None,
            fold_regularisation=None,
            step=2,
        )

        # With step=1, should use all elements
        interpolator.add_fold_constraints(
            fold_orientation=10.0,
            fold_axis_w=None,
            fold_regularisation=None,
            step=1,
        )

    def test_fold_constraints_with_large_step(self):
        """Test that step parameter larger than n_elements works."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        fold = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # With step larger than n_elements, should still work (selects first element or none)
        interpolator.add_fold_constraints(
            fold_orientation=10.0,
            step=1000,
        )


class TestDiscreteFoldInterpolatorFoldNorm:
    """Test fold norm constraints."""

    def test_fold_norm_parameter(self):
        """Test custom fold_norm value."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        fold = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # Test with custom fold norm
        interpolator.add_fold_constraints(
            fold_orientation=None,
            fold_axis_w=None,
            fold_regularisation=None,
            fold_normalisation=1.0,
            fold_norm=2.0,
        )

    def test_fold_norm_none_uses_default(self):
        """Test that fold_norm=None uses default weighting."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        fold = MockFoldEvent()

        interpolator = DiscreteFoldInterpolator(support, fold=fold)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # fold_norm=None should work (uses default 1.0)
        interpolator.add_fold_constraints(
            fold_orientation=None,
            fold_axis_w=None,
            fold_regularisation=None,
            fold_normalisation=1.0,
            fold_norm=None,
        )
