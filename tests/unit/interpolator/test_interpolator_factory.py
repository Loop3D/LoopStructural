import numpy as np
import pytest

from LoopStructural.geometry import BoundingBox
from LoopStructural.interpolators import (
    FiniteDifferenceInterpolator,
    InterpolatorFactory,
    InterpolatorType,
    P1Interpolator,
    StructuredGridSupport,
    TetMesh,
)


@pytest.fixture
def bounding_box():
    return BoundingBox(np.array([0, 0, 0]), np.array([1, 1, 1]))


def test_create_interpolator_with_string_fdi(bounding_box):
    interpolator = InterpolatorFactory.create_interpolator("FDI", bounding_box, 1000)
    assert isinstance(interpolator, FiniteDifferenceInterpolator)
    assert isinstance(interpolator.support, StructuredGridSupport)


def test_create_interpolator_with_string_pli(bounding_box):
    interpolator = InterpolatorFactory.create_interpolator("PLI", bounding_box, 1000)
    assert isinstance(interpolator, P1Interpolator)
    assert isinstance(interpolator.support, TetMesh)


def test_create_interpolator_with_enum(bounding_box):
    interpolator = InterpolatorFactory.create_interpolator(
        InterpolatorType.FINITE_DIFFERENCE, bounding_box, 1000
    )
    assert isinstance(interpolator, FiniteDifferenceInterpolator)


def test_create_interpolator_with_explicit_support(bounding_box):
    support = StructuredGridSupport(
        origin=bounding_box.origin, nsteps=np.array([5, 5, 5]), step_vector=np.array([0.2, 0.2, 0.2])
    )
    interpolator = InterpolatorFactory.create_interpolator(
        "FDI", bounding_box, nelements=None, support=support
    )
    assert interpolator.support is support


def test_create_interpolator_missing_type_raises(bounding_box):
    with pytest.raises(ValueError, match="No interpolator type specified"):
        InterpolatorFactory.create_interpolator(None, bounding_box, 1000)


def test_create_interpolator_missing_bounding_box_raises():
    with pytest.raises(ValueError, match="No bounding box specified"):
        InterpolatorFactory.create_interpolator("FDI", None, 1000)


def test_from_dict_missing_type_raises(bounding_box):
    with pytest.raises(ValueError, match="No interpolator type specified"):
        InterpolatorFactory.from_dict({"boundingbox": bounding_box, "nelements": 1000})


def test_from_dict_builds_interpolator(bounding_box):
    d = {"type": "FDI", "boundingbox": bounding_box, "nelements": 1000}
    interpolator = InterpolatorFactory.from_dict(d)
    assert isinstance(interpolator, FiniteDifferenceInterpolator)
    # from_dict should not mutate the caller's dictionary
    assert "type" in d


def test_get_supported_interpolators_contains_fdi_and_pli():
    supported = InterpolatorFactory.get_supported_interpolators()
    assert InterpolatorType.FINITE_DIFFERENCE in supported
    assert InterpolatorType.PIECEWISE_LINEAR in supported


def test_create_interpolator_with_data_sets_value_constraints(bounding_box):
    value_constraints = np.array([[0.5, 0.5, 0.5, 1.0]])
    interpolator = InterpolatorFactory.create_interpolator_with_data(
        "FDI",
        bounding_box,
        1000,
        value_constraints=value_constraints,
    )
    assert np.array_equal(interpolator.get_value_constraints()[:, :4], value_constraints)


def test_create_interpolator_with_data_sets_gradient_constraints(bounding_box):
    gradient_constraints = np.array([[0.5, 0.5, 0.5, 0.0, 0.0, 1.0]])
    interpolator = InterpolatorFactory.create_interpolator_with_data(
        "FDI",
        bounding_box,
        1000,
        gradient_constraints=gradient_constraints,
    )
    assert np.array_equal(interpolator.get_gradient_constraints()[:, :6], gradient_constraints)


def test_create_interpolator_with_data_sets_normal_constraints(bounding_box):
    gradient_norm_constraints = np.array([[0.5, 0.5, 0.5, 0.0, 0.0, 1.0]])
    interpolator = InterpolatorFactory.create_interpolator_with_data(
        "FDI",
        bounding_box,
        1000,
        gradient_norm_constraints=gradient_norm_constraints,
    )
    assert np.array_equal(interpolator.get_norm_constraints()[:, :6], gradient_norm_constraints)
