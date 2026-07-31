import numpy as np
import pytest
from loop_interpolation import _validation

ValidationError = _validation.ValidationError
ShapeError = _validation.ShapeError
DtypeError = _validation.DtypeError
VectorError = _validation.VectorError
WeightError = _validation.WeightError


def test_value_constraint_valid_and_additional_weight_column_supported():
    pts = np.array([[0, 0, 0, 1.0], [1, 1, 1, 2.0]])
    out = _validation.validate_value_constraint(pts)
    assert out.shape == (2, 4)
    assert out.dtype == np.float64


def test_value_constraint_non_finite_rejected():
    pts = np.array([[0.0, 0.0, 0.0, np.nan]])
    out = _validation.validate_value_constraint(pts)
    assert out.shape == (0, 4)


def test_value_constraint_non_numeric_rejected():
    pts = np.array([["x", "y", "z", "v"]])
    with pytest.raises(DtypeError):
        _validation.validate_value_constraint(pts)


def test_gradient_constraint_zero_vector_rejected():
    pts = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
    with pytest.raises(VectorError):
        _validation.validate_gradient_constraint(pts)


def test_gradient_constraint_non_finite_rejected():
    pts = np.array([[0.0, 0.0, 0.0, 1.0, np.inf, 0.0]])
    out = _validation.validate_gradient_constraint(pts)
    assert out.shape == (0, 6)


def test_normal_constraint_zero_vector_rejected():
    pts = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
    with pytest.raises(VectorError):
        _validation.validate_normal_constraint(pts)


def test_tangent_constraint_zero_vector_rejected():
    pts = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
    with pytest.raises(VectorError):
        _validation.validate_tangent_constraint(pts)


def test_interface_constraint_bad_shape_rejected():
    pts = np.array([[0.0, 0.0, 0.0]])
    with pytest.raises(ShapeError):
        _validation.validate_interface_constraint(pts)


def test_inequality_value_constraint_invalid_bounds_rejected():
    pts = np.array([[0.0, 0.0, 0.0, 2.0, 1.0]])
    with pytest.raises(ValidationError):
        _validation.validate_inequality_value_constraint(pts)


def test_inequality_pairs_constraint_bad_shape_rejected():
    pts = np.array([[0.0, 0.0, 0.0]])
    with pytest.raises(ShapeError):
        _validation.validate_inequality_pairs_constraint(pts)


def test_validate_weights_scalar_and_array():
    assert _validation.validate_weights(1.5, n_constraints=3) == 1.5
    arr = np.array([1.0, 2.0, 3.0])
    out = _validation.validate_weights(arr, n_constraints=3)
    assert np.array_equal(out, arr)


def test_validate_weights_non_positive_rejected():
    with pytest.raises(WeightError):
        _validation.validate_weights(0.0, n_constraints=1)
    with pytest.raises(WeightError):
        _validation.validate_weights(np.array([1.0, -1.0]), n_constraints=2)
