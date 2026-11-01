"""Regression tests for gradients evaluated through faults."""

import numpy as np
import pytest

from LoopStructural.modelling.features import GeologicalFeature
from LoopStructural.modelling.features._base_geological_feature import BaseFeature


class IdentityFault(BaseFeature):
    """A fault double that leaves coordinates unchanged."""

    def __init__(self, name="identity", coordinate=None):
        super().__init__(name)
        self.coordinate = coordinate or self

    def evaluate_value(self, points, **kwargs):
        return self.coordinate_values(points)

    def coordinate_values(self, points):
        return points[:, 0]

    def evaluate_gradient(self, points, **kwargs):
        return np.broadcast_to([1.0, 0.0, 0.0], points.shape).copy()

    def __getitem__(self, key):
        assert key == 0
        return self.coordinate

    def apply_to_points(self, points, reverse=False):
        return points

    def copy(self, name=None):
        return IdentityFault(name or self.name, self.coordinate)

    def get_data(self, value_map=None):
        return []


class ControlledFault(IdentityFault):
    """Return mixed corner signs once, then positive signs."""

    def __init__(self):
        super().__init__("controlled")
        self.value_calls = 0

    def evaluate_value(self, points, **kwargs):
        self.value_calls += 1
        if self.value_calls == 1:
            return np.where(np.arange(len(points)) % 4 == 0, 1.0, -1.0)
        return np.ones(len(points))


class LinearField:
    """A linear interpolator with a constant gradient."""

    gradient = np.array([2.0, -3.0, 0.5])

    def evaluate_value(self, points):
        return points @ self.gradient + 7.0

    def evaluate_gradient(self, points):
        return np.broadcast_to(self.gradient, points.shape).copy()


def make_feature(faults=(), regions=None):
    return GeologicalFeature(
        "linear",
        builder=None,
        regions=regions,
        faults=list(faults),
        interpolator=LinearField(),
    )


@pytest.mark.parametrize("fault_count", [1, 2])
def test_faulted_gradient_has_three_components(fault_count):
    points = np.array([[1.0, 2.0, 3.0], [4.0, -2.0, 1.0]])
    feature = make_feature(IdentityFault(str(i)) for i in range(fault_count))

    gradient = feature.evaluate_gradient(points, element_scale_parameter=1)

    expected = np.tile([2.0, -3.0, 0.5], (2, 1))
    assert gradient.shape == (2, 3)
    np.testing.assert_allclose(gradient, expected, rtol=0, atol=1e-12)


def test_faulted_gradient_keeps_region_masked_point_nan():
    points = np.array([[1.0, 2.0, 3.0], [-1.0, 2.0, 3.0]])
    feature = make_feature([IdentityFault()], regions=[lambda p: p[:, 0] > 0])

    gradient = feature.evaluate_gradient(points, element_scale_parameter=1)

    np.testing.assert_allclose(gradient[0], [2.0, -3.0, 0.5], rtol=0, atol=1e-12)
    assert np.isnan(gradient[1]).all()


def test_faulted_gradient_refines_once_when_corner_signs_are_mixed():
    points = np.array([[1.0, 2.0, 3.0], [4.0, -2.0, 1.0]])
    fault = ControlledFault()
    feature = make_feature([fault])

    gradient = feature.evaluate_gradient(points, element_scale_parameter=1)

    assert fault.value_calls == 2
    expected = np.tile([2.0, -3.0, 0.5], (2, 1))
    assert gradient.shape == (2, 3)
    np.testing.assert_allclose(gradient, expected, rtol=0, atol=1e-12)
