import numpy as np
import pytest

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


def test_inverse_transform_is_broken_for_normal_point_clouds_bug():
    """Documents a real bug in EuclideanTransformation.inverse_transform().

    The implementation slices `points[: self.dimensions]` which slices the
    first `self.dimensions` *rows* of the array (not columns, unlike every
    other method on this class which uses `points[:, : self.dimensions]`).
    Combined with an einsum contracting over the last axis against the
    (dimensions x dimensions) rotation matrix, this raises a ValueError for
    any input whose number of columns does not equal `self.dimensions.`
    In practice this means `inverse_transform` cannot be used to invert the
    output of `transform`/`fit_transform` for ordinary xyz point arrays.
    """
    pts = _line_points(angle_deg=30.0)
    t = EuclideanTransformation(dimensions=2)
    transformed = t.fit_transform(pts)

    with pytest.raises(ValueError):
        t.inverse_transform(transformed)


def test_repr_html_contains_translation_and_angle():
    t = EuclideanTransformation(dimensions=2)
    pts = _line_points(angle_deg=15.0)
    t.fit(pts)

    html = t._repr_html_()
    assert isinstance(html, str)
    assert "Translation" in html
    assert "Rotation Angle" in html
