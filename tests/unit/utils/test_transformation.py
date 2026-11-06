import pickle

import numpy as np
import pytest
from tests.fixtures.coordinate_frames import ANCHORS, tp

from LoopStructural.utils import EuclideanTransformation


def _line_points(n=50, angle_deg=30.0, centre=(5.0, 5.0)):
    """Points scattered along a line through `centre` at `angle_deg` to x,
    with an unrelated z coordinate so we can check that dimensions=2
    leaves z untouched."""
    t_param = np.linspace(-5, 5, n)
    angle = np.deg2rad(angle_deg)
    pts = np.zeros((n, 3))
    pts[:, 0] = centre[0] + t_param * np.cos(angle)
    pts[:, 1] = centre[1] + t_param * np.sin(angle)
    pts[:, 2] = np.arange(n) * 0.1
    return pts


def _tolerance(points):
    """Return the scale-aware tolerance `Tp` for a point cloud.

    L is the largest extent of the cloud (at least 1.0).
    W is the largest absolute coordinate."""
    length = max(float(np.ptp(points[:, :3], axis=0).max()), 1.0)
    width = float(np.abs(points[:, :3]).max())
    return tp(length, width), length


def _angle_error(angle, expected):
    """Return the difference of two angles, modulo pi."""
    return abs((angle - expected + np.pi / 2) % np.pi - np.pi / 2)


def test_default_construction():
    t = EuclideanTransformation()
    assert t.dimensions == 2
    assert t.angle == 0
    assert t.fit_rotation is True
    assert np.allclose(t.translation, [0, 0])


def test_fit_finds_rotation_that_aligns_main_axis_with_x():
    pts = _line_points(angle_deg=30.0)
    t = EuclideanTransformation(dimensions=2)
    t.fit(pts)

    # translation should recover the centre of the point cloud
    assert np.allclose(t.translation, [5.0, 5.0])
    # the fitted angle should align the 30 degree line with the x axis:
    # -30 degrees in radians
    assert np.isclose(t.angle, np.deg2rad(-30.0))


def test_transform_aligns_variance_with_x_axis():
    pts = _line_points(angle_deg=30.0)
    t = EuclideanTransformation(dimensions=2)
    transformed = t.fit_transform(pts)

    # after alignment nearly all variance should be along x, none along y
    assert np.var(transformed[:, 0]) > 1.0
    assert np.var(transformed[:, 1]) < 1e-20
    # z (untouched dimension) must be preserved exactly
    assert np.allclose(transformed[:, 2], pts[:, 2])


def test_fit_rotation_false_keeps_angle_zero():
    pts = _line_points(angle_deg=30.0)
    t = EuclideanTransformation(dimensions=2, fit_rotation=False)
    t.fit(pts)

    assert t.angle == 0
    # translation is still fitted even when rotation fitting is disabled
    assert np.allclose(t.translation, [5.0, 5.0])


def test_transform_raises_if_points_have_too_few_columns():
    t = EuclideanTransformation(dimensions=3)
    pts_2d = np.zeros((5, 2))
    with pytest.raises(ValueError):
        t.transform(pts_2d)


def test_fit_raises_if_points_have_too_few_columns():
    t = EuclideanTransformation(dimensions=3)
    pts_2d = np.zeros((5, 2))
    with pytest.raises(ValueError):
        t.fit(pts_2d)


def test_rotation_and_inverse_rotation_are_transposes_in_plane():
    t = EuclideanTransformation(dimensions=2, angle=np.pi / 4)
    rot = t.rotation
    inv_rot = t.inverse_rotation
    # for the 2D (x, y) block, rotating by -angle should be the transpose
    # (inverse) of rotating by +angle
    assert np.allclose(rot[:2, :2].T, inv_rot[:2, :2])


def test_call_is_equivalent_to_transform():
    pts = _line_points(angle_deg=10.0)
    t = EuclideanTransformation(dimensions=2)
    t.fit(pts)

    assert np.allclose(t(pts), t.transform(pts))


def test_repr_html_contains_translation_and_angle():
    t = EuclideanTransformation(dimensions=2)
    pts = _line_points(angle_deg=15.0)
    t.fit(pts)

    html = t._repr_html_()
    assert isinstance(html, str)
    assert "Translation" in html
    assert "Rotation Angle" in html


@pytest.mark.parametrize("angle_deg", [30.0, -75.0])
@pytest.mark.parametrize("columns", [3, 5])
@pytest.mark.parametrize("anchor", ANCHORS, ids=lambda a: f"anchor{a[0]:g}")
@pytest.mark.parametrize("n", [1, 4, 50])
@pytest.mark.parametrize("dimensions", [2, 3])
def test_round_trip_explicit_parameters(dimensions, n, anchor, columns, angle_deg):
    """The inverse of the forward transform returns the input points.

    The forward result must also match the closed form of the rotation."""
    rng = np.random.default_rng(0)
    points = rng.uniform(-50, 50, size=(n, columns))
    points[:, :3] += anchor
    original = points.copy()
    tol, length = _tolerance(points)
    angle = np.deg2rad(angle_deg)
    translation = anchor + np.array([1.0, 2.0, 3.0])
    t = EuclideanTransformation(dimensions=dimensions, angle=angle, translation=translation)

    forward = t.transform(points)

    # the transform must really move the points
    assert np.abs(forward[:, :dimensions] - points[:, :dimensions]).max() > 0.1
    cos, sin = np.cos(angle), np.sin(angle)
    dx = points[:, 0] - translation[0]
    dy = points[:, 1] - translation[1]
    assert np.abs(forward[:, 0] - (cos * dx - sin * dy)).max() <= 1e-12 * length
    assert np.abs(forward[:, 1] - (sin * dx + cos * dy)).max() <= 1e-12 * length
    if dimensions == 3:
        assert np.abs(forward[:, 2] - (points[:, 2] - translation[2])).max() <= 1e-12 * length
    # the columns after the transformed ones are not changed
    assert np.array_equal(forward[:, dimensions:], points[:, dimensions:])

    back = t.inverse_transform(forward)

    assert back.shape == points.shape
    assert np.abs(back - original).max() <= tol
    assert np.array_equal(back[:, 3:], original[:, 3:])
    # the input is not changed and the results do not share memory with it
    assert np.array_equal(points, original)
    assert not np.shares_memory(forward, points)
    assert not np.shares_memory(back, forward)
    assert not np.shares_memory(back, points)


@pytest.mark.parametrize("fit_rotation", [True, False])
@pytest.mark.parametrize("n", [1, 4, 50])
@pytest.mark.parametrize("anchor", ANCHORS, ids=lambda a: f"anchor{a[0]:g}")
def test_round_trip_after_fit(anchor, n, fit_rotation):
    """A fitted transformation is inverted by `inverse_transform`."""
    points = _line_points(50, 30.0, centre=(anchor[0], anchor[1]))[:n]
    points[:, 2] += anchor[2]
    original = points.copy()
    tol, _ = _tolerance(points)
    t = EuclideanTransformation(dimensions=2, fit_rotation=fit_rotation)

    forward = t.fit_transform(points)

    assert np.abs(t.translation - np.mean(original[:, :2], axis=0)).max() <= tol
    back = t.inverse_transform(forward)
    assert np.abs(back - original).max() <= tol
    assert np.array_equal(points, original)


def test_integer_input_is_not_truncated():
    """An integer array must give the float result, not a truncated one."""
    points = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]], dtype=np.int64)
    tol, _ = _tolerance(points.astype(float))
    angle = np.deg2rad(30.0)
    translation = np.array([0.5, 0.25])
    t = EuclideanTransformation(dimensions=2, angle=angle, translation=translation)

    result = t.transform(points)

    dx = points[:, 0] - translation[0]
    dy = points[:, 1] - translation[1]
    expected = np.column_stack(
        [
            np.cos(angle) * dx - np.sin(angle) * dy,
            np.sin(angle) * dx + np.cos(angle) * dy,
            points[:, 2].astype(float),
        ]
    )
    assert np.abs(result - expected).max() <= tol


@pytest.mark.parametrize("method", ["transform", "inverse_transform", "fit"])
def test_one_dimensional_input_raises_value_error(method):
    t = EuclideanTransformation(dimensions=2)
    with pytest.raises(ValueError, match="2D array"):
        getattr(t, method)(np.zeros(3))


@pytest.mark.parametrize("method", ["transform", "fit"])
def test_too_few_columns_raises_value_error(method):
    t = EuclideanTransformation(dimensions=3)
    with pytest.raises(ValueError, match="at least 3 dimensions"):
        getattr(t, method)(np.zeros((5, 2)))


def test_entry_points_are_one_class():
    """All four import paths give the same class object."""
    import loop_common.math as common_math
    import loop_common.math._transformation as new_module

    import LoopStructural.utils as utils
    import LoopStructural.utils._transformation as old_module

    assert utils.EuclideanTransformation is new_module.EuclideanTransformation
    assert old_module.EuclideanTransformation is new_module.EuclideanTransformation
    assert common_math.EuclideanTransformation is new_module.EuclideanTransformation


def _fitted_common_object():
    """Return a fitted object of the class in `loop_common.math`."""
    from loop_common.math import EuclideanTransformation as CommonClass

    fitted = CommonClass(dimensions=2)
    fitted.fit(_line_points(angle_deg=30.0))
    return fitted


def test_pickle_with_the_old_module_path_loads(monkeypatch):
    """A pickle that names the old module path still loads."""
    fitted = _fitted_common_object()
    with monkeypatch.context() as patch:
        patch.setattr(type(fitted), "__module__", "LoopStructural.utils._transformation")
        data = pickle.dumps(fitted)
    assert b"LoopStructural.utils._transformation" in data

    loaded = pickle.loads(data)

    assert np.array_equal(loaded.translation, fitted.translation)
    assert loaded.angle == fitted.angle
    assert loaded.dimensions == fitted.dimensions
    assert loaded.fit_rotation == fitted.fit_rotation


def test_new_pickle_names_the_loop_common_path():
    data = pickle.dumps(_fitted_common_object())
    assert b"loop_common.math._transformation" in data


@pytest.mark.parametrize("slope_deg", [-30.0, -10.0])
def test_fit_aligns_a_line_with_negative_slope(slope_deg):
    """A line with a negative slope must end up along the x axis."""
    pts = _line_points(angle_deg=slope_deg)
    t = EuclideanTransformation(dimensions=2)
    transformed = t.fit_transform(pts)

    assert _angle_error(t.angle, np.deg2rad(-slope_deg)) < 1e-9
    assert np.var(transformed[:, 1]) < 1e-20


@pytest.mark.parametrize("slope_deg", [30.0, -30.0])
def test_fit_three_dimensions(slope_deg):
    """The fit works with `dimensions=3` and keeps z as z - mean(z)."""
    pts = _line_points(angle_deg=slope_deg)
    tol, _ = _tolerance(pts)
    t = EuclideanTransformation(dimensions=3)
    transformed = t.fit_transform(pts)

    assert np.allclose(t.translation, [5.0, 5.0, 2.45], atol=1e-9, rtol=0)
    assert _angle_error(t.angle, np.deg2rad(-slope_deg)) < 1e-9
    assert np.abs(transformed[:, 2] - (pts[:, 2] - 2.45)).max() <= 1e-12
    assert np.var(transformed[:, 1]) < 1e-20
    assert np.abs(t.inverse_transform(transformed) - pts).max() <= tol
