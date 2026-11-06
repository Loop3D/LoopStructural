"""Tests for `EuclideanTransformation` in `loop_common.math`."""

import importlib
import sys

import numpy as np
import pytest

ANCHOR = np.array([500000.0, 6500000.0, 100.0])


def _tp(length, width):
    """Return a scale-aware floating point tolerance."""
    return 1e-10 * length + 32 * np.finfo(float).eps * max(width, length)


def _load_class():
    """Return the class from its module, so that each test fails alone."""
    module = importlib.import_module("loop_common.math._transformation")
    return module.EuclideanTransformation


def test_module_imports():
    module = importlib.import_module("loop_common.math._transformation")
    assert hasattr(module, "EuclideanTransformation")


def test_export_from_loop_common_math():
    import loop_common.math

    assert loop_common.math.EuclideanTransformation is _load_class()


@pytest.mark.parametrize("n", [1, 4, 50])
@pytest.mark.parametrize("dimensions", [2, 3])
def test_round_trip_without_fit(dimensions, n):
    """The inverse of the forward transform returns the input points."""
    rng = np.random.default_rng(0)
    points = ANCHOR + rng.uniform(-50, 50, size=(n, 3))
    original = points.copy()
    length = max(float(np.ptp(points, axis=0).max()), 1.0)
    tol = _tp(length, float(np.abs(points).max()))
    transformation_class = _load_class()
    t = transformation_class(dimensions=dimensions, angle=0.7, translation=ANCHOR.copy())

    forward = t.transform(points)

    assert np.abs(forward - points).max() > 0.1
    if dimensions == 3:
        assert np.abs(forward[:, 2] - (points[:, 2] - ANCHOR[2])).max() <= 1e-12 * length
    back = t.inverse_transform(forward)
    assert np.abs(back - original).max() <= tol
    assert np.array_equal(points, original)


def test_two_dimensional_rotation_of_known_angle():
    transformation_class = _load_class()
    t = transformation_class(dimensions=2, angle=np.pi / 6, translation=np.zeros(2))

    result = t.transform(np.array([[2.0, 0.0, 7.0]]))

    assert np.allclose(result[0, :2], [np.sqrt(3.0), 1.0], atol=1e-12, rtol=0)
    assert result[0, 2] == 7.0


def test_three_dimensional_rotation_keeps_the_z_direction():
    transformation_class = _load_class()
    t = transformation_class(dimensions=3, angle=np.pi / 2, translation=np.array([1.0, 2.0, 3.0]))
    points = np.array([[2.0, 2.0, 5.0], [1.0, 3.0, 3.0]])
    expected = np.array([[0.0, 1.0, 2.0], [-1.0, 0.0, 0.0]])

    result = t.transform(points)

    assert np.abs(result - expected).max() <= 1e-12
    assert np.isclose(np.linalg.det(t.rotation[:3, :3]), 1.0, atol=1e-12, rtol=0)
    assert np.abs(t.inverse_transform(result) - points).max() <= 1e-12


def test_fit_without_scikit_learn_raises_import_error(monkeypatch):
    monkeypatch.setitem(sys.modules, "sklearn", None)
    transformation_class = _load_class()
    t = transformation_class(dimensions=2, angle=0.3, translation=np.array([1.0, 2.0]))
    points = np.array([[0.0, 0.0, 5.0], [2.0, 4.0, 6.0], [3.0, 1.0, 7.0]])

    with pytest.raises(ImportError, match="scikit-learn"):
        t.fit(points)

    assert np.array_equal(t.translation, [1.0, 2.0])
    assert t.angle == 0.3


def test_fit_without_rotation_needs_no_scikit_learn(monkeypatch):
    monkeypatch.setitem(sys.modules, "sklearn", None)
    transformation_class = _load_class()
    t = transformation_class(dimensions=2, fit_rotation=False)

    t.fit(np.array([[0.0, 0.0, 5.0], [2.0, 4.0, 6.0]]))

    assert np.allclose(t.translation, [1.0, 2.0])
    assert t.angle == 0
