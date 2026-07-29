import numpy as np
import pytest
from loop_interpolation.constraints import (
    ValueConstraint,
    GradientConstraint,
    InequalityConstraint,
    InequalityPair,
    InterfaceConstraint,
)


def test_value_constraint():
    points = np.array([[0, 0, 0], [1, 1, 1]])
    values = np.array([10, 20])
    weights = np.array([1.0, 0.5])
    constraint = ValueConstraint(points=points, values=values, weights=weights)

    assert np.array_equal(constraint.points, points)
    assert np.array_equal(constraint.values, values)
    assert np.array_equal(constraint.weights, weights)


def test_gradient_constraint():
    points = np.array([[0, 0, 0], [1, 1, 1]])
    vectors = np.array([[1, 0, 0], [0, 1, 0]])
    weights = np.array([1.0, 0.5])
    constraint = GradientConstraint(points=points, vectors=vectors, weights=weights, is_normal=True)

    assert np.array_equal(constraint.points, points)
    assert np.array_equal(constraint.vectors, vectors)
    assert np.array_equal(constraint.weights, weights)
    assert constraint.is_normal


def test_inequality_constraint():
    points = np.array([[0, 0, 0], [1, 1, 1]])
    bounds = np.array([[0, 10], [5, 15]])
    weights = np.array([1.0, 0.5])
    constraint = InequalityConstraint(points=points, bounds=bounds, weights=weights)

    assert np.array_equal(constraint.points, points)
    assert np.array_equal(constraint.bounds, bounds)
    assert np.array_equal(constraint.weights, weights)


def test_inequality_pair():
    points = np.array([[0, 0, 0], [1, 1, 1]])
    pair_ids = np.array([0, 1])
    weights = np.array([1.0, 0.5])
    constraint = InequalityPair(points=points, pair_ids=pair_ids, weights=weights)

    assert np.array_equal(constraint.points, points)
    assert np.array_equal(constraint.pair_ids, pair_ids)
    assert np.array_equal(constraint.weights, weights)


def test_interface_constraint():
    points = np.array([[0, 0, 0], [1, 1, 1]])
    interface_ids = np.array([10.0, 20.0])
    weights = np.array([1.0, 0.5])
    constraint = InterfaceConstraint(points=points, interface_ids=interface_ids, weights=weights)

    assert np.array_equal(constraint.points, points)
    assert np.array_equal(constraint.interface_ids, interface_ids)
    assert np.array_equal(constraint.weights, weights)


def test_constraint_json_round_trip():
    points = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    values = np.array([10.0, 20.0])
    constraint = ValueConstraint(points=points, values=values, weights=np.array([1.0, 0.5]))

    payload = constraint.model_dump_json()
    restored = ValueConstraint.model_validate_json(payload)

    assert np.array_equal(restored.points, points)
    assert np.array_equal(restored.values, values)
    assert np.array_equal(restored.weights, np.array([1.0, 0.5]))


def test_value_constraint_drops_non_finite_rows_and_repairs_nan_weights():
    constraint = ValueConstraint.from_array(
        np.array(
            [
                [0.0, 0.0, 0.0, 1.0, np.nan],
                [1.0, 1.0, 1.0, np.nan, 2.0],
            ]
        )
    )

    assert constraint.points.shape == (1, 3)
    assert np.array_equal(constraint.values, np.array([1.0]))
    assert np.array_equal(constraint.weights, np.array([1.0]))


def test_gradient_constraint_rejects_zero_vector_via_object_validation():
    with pytest.raises(Exception) as excinfo:
        GradientConstraint(points=np.array([[0.0, 0.0, 0.0]]), vectors=np.array([[0.0, 0.0, 0.0]]))

    assert "zero or near-zero magnitude" in str(excinfo.value)
