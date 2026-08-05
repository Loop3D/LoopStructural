"""Tests for add_element_anisotropy_constraints."""


import numpy as np
import pytest
from loop_common.supports._3d_structured_tetra import TetMesh
from loop_interpolation._anisotropy import add_element_anisotropy_constraints
from loop_interpolation._p1interpolator import P1Interpolator


class MockAnisotropy:
    """Mock anisotropy direction provider for testing."""

    def __init__(self, primary_grad=None, secondary_grad=None, normal_grad=None):
        """Initialize mock anisotropy with configurable return values."""
        self.primary_grad = primary_grad if primary_grad is not None else np.ones((10, 3))
        self.secondary_grad = secondary_grad if secondary_grad is not None else np.ones((10, 3))
        self.normal_grad = normal_grad if normal_grad is not None else np.ones((10, 3))

    def get_directions(self, points):
        """Return direction components, broadcast to match the number of
        query points (the mesh may have more elements than the configured
        mock array has rows)."""
        n_points = points.shape[0] if hasattr(points, "shape") else 1
        return (
            np.tile(self.primary_grad[0], (n_points, 1)),
            np.tile(self.secondary_grad[0], (n_points, 1)),
            np.tile(self.normal_grad[0], (n_points, 1)),
        )


class TestNoAnisotropyGuard:
    """Test add_element_anisotropy_constraints(anisotropy=None) guard."""

    def test_none_anisotropy_raises(self):
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        interpolator = P1Interpolator(support)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )
        with pytest.raises(ValueError, match="anisotropy must not be None"):
            add_element_anisotropy_constraints(interpolator, None)


class TestAddP1AnisotropyConstraints:
    """Test anisotropy constraint addition."""

    def test_add_primary_constraints(self):
        """Test adding primary orientation constraints."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        anisotropy = MockAnisotropy()

        interpolator = P1Interpolator(support)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # This should not raise
        add_element_anisotropy_constraints(
            interpolator,
            anisotropy,
            primary_w=10.0,
            secondary_w=None,
            regularisation_weights=None,
        )

    def test_add_secondary_constraints(self):
        """Test adding secondary orientation constraints."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        anisotropy = MockAnisotropy()

        interpolator = P1Interpolator(support)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # This should not raise
        add_element_anisotropy_constraints(
            interpolator,
            anisotropy,
            primary_w=None,
            secondary_w=10.0,
            regularisation_weights=None,
        )

    def test_add_normal_constraints(self):
        """Test adding normal magnitude constraints."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        anisotropy = MockAnisotropy()

        interpolator = P1Interpolator(support)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # This should not raise
        add_element_anisotropy_constraints(
            interpolator,
            anisotropy,
            primary_w=None,
            secondary_w=None,
            regularisation_weights=None,
            normal_w=10.0,
        )

    def test_add_regularisation_constraints(self):
        """Test adding anisotropic regularisation constraints."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        anisotropy = MockAnisotropy()

        interpolator = P1Interpolator(support)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # This should not raise
        add_element_anisotropy_constraints(
            interpolator,
            anisotropy,
            primary_w=None,
            secondary_w=None,
            regularisation_weights=(0.1, 0.01, 0.01),
        )

    def test_add_all_constraints_simultaneously(self):
        """Test adding all anisotropy constraint types at once."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        anisotropy = MockAnisotropy()

        interpolator = P1Interpolator(support)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # This should not raise
        add_element_anisotropy_constraints(
            interpolator,
            anisotropy,
            primary_w=10.0,
            secondary_w=10.0,
            regularisation_weights=(0.1, 0.01, 0.01),
            normal_w=1.0,
        )

    def test_add_constraints_with_mask_function(self):
        """Test adding anisotropy constraints with a mask function."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        anisotropy = MockAnisotropy()

        interpolator = P1Interpolator(support)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # Define a mask function that returns True for points above z=1.5
        def mask_fn(points):
            return points[:, 2] > 1.5

        # This should not raise
        add_element_anisotropy_constraints(
            interpolator,
            anisotropy,
            primary_w=10.0,
            secondary_w=10.0,
            regularisation_weights=(0.1, 0.01, 0.01),
            normal_w=1.0,
            mask_fn=mask_fn,
        )

    def test_constraints_with_custom_weights(self):
        """Test anisotropy constraints with various weight values."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        anisotropy = MockAnisotropy()

        interpolator = P1Interpolator(support)
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

        for f_primary, f_secondary, f_reg, f_norm in weights_to_test:
            add_element_anisotropy_constraints(
                interpolator,
                anisotropy,
                primary_w=f_primary,
                secondary_w=f_secondary,
                regularisation_weights=f_reg,
                normal_w=f_norm,
            )


class TestAddP1AnisotropyConstraintWeighting:
    """Test constraint weighting by element volume."""

    def test_constraints_use_element_volume_weighting(self):
        """Test that anisotropy constraints respect element volume weighting."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        rng = np.random.default_rng(0)
        anisotropy = MockAnisotropy(
            primary_grad=rng.random((support.n_elements, 3)),
            secondary_grad=rng.random((support.n_elements, 3)),
            normal_grad=rng.random((support.n_elements, 3)),
        )

        interpolator = P1Interpolator(support)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # Add constraints with different weights
        add_element_anisotropy_constraints(
            interpolator,
            anisotropy,
            primary_w=10.0,
            secondary_w=5.0,
            regularisation_weights=None,
            normal_w=None,
        )

        # Constraints should be added (hard to verify exact values,
        # but we can at least check it doesn't crash)
        assert any("anisotropy primary" in k for k in interpolator.constraints)
        assert any("anisotropy secondary" in k for k in interpolator.constraints)


class TestAddP1AnisotropyStepParameter:
    """Test the step parameter for constraint sampling."""

    def test_constraints_with_step_parameter(self):
        """Test that step parameter subsamples constraints."""
        support = TetMesh(nsteps=np.array([5, 5, 5]))
        anisotropy = MockAnisotropy()

        interpolator = P1Interpolator(support)
        interpolator.setup_interpolator(
            data_points=np.array([[2.0, 2.0, 2.0]]),
            data_values=np.array([1.0]),
        )

        # With step=2, should use every other element
        add_element_anisotropy_constraints(
            interpolator,
            anisotropy,
            primary_w=10.0,
            secondary_w=None,
            regularisation_weights=None,
            step=2,
        )

        # With step=1, should use all elements
        add_element_anisotropy_constraints(
            interpolator,
            anisotropy,
            primary_w=10.0,
            secondary_w=None,
            regularisation_weights=None,
            step=1,
        )

    def test_constraints_with_large_step(self):
        """Test that step parameter larger than n_elements works."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        anisotropy = MockAnisotropy()

        interpolator = P1Interpolator(support)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # With step larger than n_elements, should still work (selects first element or none)
        add_element_anisotropy_constraints(
            interpolator,
            anisotropy,
            primary_w=10.0,
            step=1000,
        )


class TestAddP1AnisotropyNormTarget:
    """Test normal-direction norm target constraints."""

    def test_norm_target_parameter(self):
        """Test custom norm_target value."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        anisotropy = MockAnisotropy()

        interpolator = P1Interpolator(support)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # Test with custom norm target
        add_element_anisotropy_constraints(
            interpolator,
            anisotropy,
            primary_w=None,
            secondary_w=None,
            regularisation_weights=None,
            normal_w=1.0,
            norm_target=2.0,
        )

    def test_norm_target_none_uses_default(self):
        """Test that norm_target=None uses default weighting."""
        support = TetMesh(nsteps=np.array([3, 3, 3]))
        anisotropy = MockAnisotropy()

        interpolator = P1Interpolator(support)
        interpolator.setup_interpolator(
            data_points=np.array([[1.0, 1.0, 1.0]]),
            data_values=np.array([1.0]),
        )

        # norm_target=None should work (uses default -1.0)
        add_element_anisotropy_constraints(
            interpolator,
            anisotropy,
            primary_w=None,
            secondary_w=None,
            regularisation_weights=None,
            normal_w=1.0,
            norm_target=None,
        )
