import pytest
import numpy as np
from LoopStructural.datatypes import BoundingBox
from LoopStructural.interpolators._interpolator_builder import InterpolatorBuilder
from LoopStructural.interpolators import InterpolatorType


@pytest.fixture
def setup_builder():
    bounding_box = BoundingBox(np.array([0, 0, 0]), np.array([1, 1, 1]))
    interpolatortype = InterpolatorType.FINITE_DIFFERENCE
    nelements = 1000
    buffer = 0.2
    builder = InterpolatorBuilder(
        interpolatortype=interpolatortype,
        bounding_box=bounding_box,
        nelements=nelements,
        buffer=buffer,
    )
    return builder


def test_create_interpolator(setup_builder):
    builder = setup_builder
    builder.build()
    assert builder.interpolator is not None, "Interpolator should be created"


def test_set_value_constraints(setup_builder):
    builder = setup_builder
    builder.build()
    value_constraints = np.array([[0.5, 0.5, 0.5, 1.0, 1.0]])
    builder.add_value_constraints(value_constraints)
    assert np.array_equal(
        builder.interpolator.data["value"], value_constraints
    ), "Value constraints should be set correctly"


def test_set_gradient_constraints(setup_builder):
    builder = setup_builder
    gradient_constraints = np.array([[0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0]])
    builder.add_gradient_constraints(gradient_constraints)
    assert np.array_equal(
        builder.interpolator.data["gradient"], gradient_constraints
    ), "Gradient constraints should be set correctly"


def test_set_normal_constraints(setup_builder):
    builder = setup_builder
    normal_constraints = np.array([[0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0]])
    builder.add_normal_constraints(normal_constraints)
    assert np.array_equal(
        builder.interpolator.data["normal"], normal_constraints
    ), "Normal constraints should be set correctly"


def test_setup_interpolator(setup_builder):
    builder = setup_builder
    builder.build()
    value_constraints = np.array([[0.5, 0.5, 0.5, 1.0, 1.0]])
    interpolator = builder.add_value_constraints(value_constraints).setup_interpolator().build()
    assert interpolator is not None, "Interpolator should be set up"
    assert np.array_equal(
        interpolator.data["value"], value_constraints
    ), "Value constraints should be set correctly after setup"


def test_evaluate_scalar_value(setup_builder):
    builder = setup_builder
    builder.build()
    value_constraints = np.array([[0.5, 0.5, 0.5, 1.0]])
    interpolator = builder.add_value_constraints(value_constraints).setup_interpolator().build()
    locations = np.array([[0.5, 0.5, 0.5]])
    values = interpolator.evaluate_value(locations)
    assert values is not None, "Evaluation should return values"
    assert values.shape == (1,), "Evaluation should return correct shape"
