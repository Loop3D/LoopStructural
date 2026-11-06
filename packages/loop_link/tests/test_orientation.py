"""Tests for the vector conversion functions in loop_link.orientation.

Convention: x is east, y is north, z is up. All angles are in degrees.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from loop_link.orientation import (
    directions_to_lineation,
    lineation_to_directions,
    normals_to_orientation,
    orientation_to_normals,
)

from loop_link.errors import InvalidInputError

AZIMUTHS = [0, 45, 90, 135, 180, 270, 359]
DIPS = [10, 45, 80]


def _upward_normal(azimuth: float, dip: float) -> np.ndarray:
    """Return the expected upward normal of a plane.

    A plane that dips to azimuth a with dip d has the upward normal
    (sin d sin a, sin d cos a, cos d).
    """
    a = math.radians(azimuth)
    d = math.radians(dip)
    return np.array([math.sin(d) * math.sin(a), math.sin(d) * math.cos(a), math.cos(d)])


def _direction(trend: float, plunge: float) -> np.ndarray:
    """Return the expected direction of a lineation.

    The direction is (cos p sin t, cos p cos t, -sin p).
    """
    t = math.radians(trend)
    p = math.radians(plunge)
    return np.array([math.cos(p) * math.sin(t), math.cos(p) * math.cos(t), -math.sin(p)])


@pytest.mark.parametrize("dip", DIPS)
@pytest.mark.parametrize("azimuth", AZIMUTHS)
def test_normal_to_orientation(azimuth, dip):
    normal = _upward_normal(azimuth, dip)
    dip_direction, got_dip = normals_to_orientation(normal[None, :])
    assert dip_direction[0] == pytest.approx(azimuth, abs=1e-9)
    assert got_dip[0] == pytest.approx(dip, abs=1e-9)


@pytest.mark.parametrize("dip", DIPS)
@pytest.mark.parametrize("azimuth", AZIMUTHS)
def test_downward_normal_gives_same_result(azimuth, dip):
    # A downward normal and its upward twin describe the same plane.
    normal = _upward_normal(azimuth, dip)
    up = normals_to_orientation(normal[None, :])
    down = normals_to_orientation(-normal[None, :])
    assert down[0][0] == pytest.approx(up[0][0], abs=1e-9)
    assert down[1][0] == pytest.approx(up[1][0], abs=1e-9)


def test_vertical_plane_has_dip_90():
    # A vertical plane has n_z = 0.
    normal = np.array([[1.0, 0.0, 0.0]])
    dip_direction, dip = normals_to_orientation(normal)
    assert dip[0] == pytest.approx(90.0, abs=1e-9)
    assert dip_direction[0] == pytest.approx(90.0, abs=1e-9)


def test_horizontal_plane_has_dip_0_and_dip_direction_0():
    dip_direction, dip = normals_to_orientation(np.array([[0.0, 0.0, 1.0]]))
    assert dip[0] == pytest.approx(0.0, abs=1e-9)
    assert dip_direction[0] == 0.0


def test_dip_direction_is_in_range():
    normals = np.array([_upward_normal(a, d) for a in AZIMUTHS for d in DIPS])
    dip_direction, _ = normals_to_orientation(normals)
    assert np.all(dip_direction >= 0.0)
    assert np.all(dip_direction < 360.0)


@pytest.mark.parametrize("plunge", [0, 10, 45, 80])
@pytest.mark.parametrize("trend", AZIMUTHS)
def test_lineation_to_directions(trend, plunge):
    got = lineation_to_directions(np.array([trend], dtype=float), np.array([plunge], dtype=float))
    assert got.shape == (1, 3)
    np.testing.assert_allclose(got[0], _direction(trend, plunge), atol=1e-9)


@pytest.mark.parametrize("plunge", [0, 10, 45, 80])
@pytest.mark.parametrize("trend", AZIMUTHS)
def test_directions_to_lineation(trend, plunge):
    direction = _direction(trend, plunge)
    got_trend, got_plunge = directions_to_lineation(direction[None, :])
    assert got_trend[0] == pytest.approx(trend, abs=1e-9)
    assert got_plunge[0] == pytest.approx(plunge, abs=1e-9)


def test_upward_direction_flips():
    # A lineation is an axis. A direction with z > 0 flips to point down.
    direction = -_direction(45.0, 30.0)
    assert direction[2] > 0.0
    trend, plunge = directions_to_lineation(direction[None, :])
    assert trend[0] == pytest.approx(45.0, abs=1e-9)
    assert plunge[0] == pytest.approx(30.0, abs=1e-9)


def test_vertical_line_has_trend_0_and_plunge_90():
    trend, plunge = directions_to_lineation(np.array([[0.0, 0.0, -1.0]]))
    assert trend[0] == 0.0
    assert plunge[0] == pytest.approx(90.0, abs=1e-9)


def test_trend_is_in_range():
    directions = np.array([_direction(t, p) for t in AZIMUTHS for p in [10, 45, 80]])
    trend, _ = directions_to_lineation(directions)
    assert np.all(trend >= 0.0)
    assert np.all(trend < 360.0)


@pytest.mark.parametrize("dip", DIPS)
@pytest.mark.parametrize("azimuth", AZIMUTHS)
def test_orientation_round_trip(azimuth, dip):
    normals = orientation_to_normals(np.array([azimuth], dtype=float), np.array([dip], dtype=float))
    assert normals.shape == (1, 3)
    np.testing.assert_allclose(normals[0], _upward_normal(azimuth, dip), atol=1e-9)
    dip_direction, got_dip = normals_to_orientation(normals)
    assert dip_direction[0] == pytest.approx(azimuth, abs=1e-9)
    assert got_dip[0] == pytest.approx(dip, abs=1e-9)


@pytest.mark.parametrize("plunge", [10, 45, 80])
@pytest.mark.parametrize("trend", AZIMUTHS)
def test_lineation_round_trip(trend, plunge):
    directions = lineation_to_directions(
        np.array([trend], dtype=float), np.array([plunge], dtype=float)
    )
    got_trend, got_plunge = directions_to_lineation(directions)
    assert got_trend[0] == pytest.approx(trend, abs=1e-9)
    assert got_plunge[0] == pytest.approx(plunge, abs=1e-9)


@pytest.mark.parametrize("bad", [[0.0, 0.0, 0.0], [float("nan"), 0.0, 1.0]], ids=["zero", "nan"])
def test_bad_normal_raises(bad):
    with pytest.raises(InvalidInputError):
        normals_to_orientation(np.array([bad]))


@pytest.mark.parametrize("bad", [[0.0, 0.0, 0.0], [float("nan"), 0.0, 1.0]], ids=["zero", "nan"])
def test_bad_direction_raises(bad):
    with pytest.raises(InvalidInputError):
        directions_to_lineation(np.array([bad]))
