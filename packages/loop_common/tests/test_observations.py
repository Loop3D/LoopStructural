import numpy as np
import pytest

from loop_common.observations.pointset import PointSet
from loop_common.observations.orientation import (
    OrientationObservation,
    OrientationType,
)
from loop_common.observations.lineset import LineSet


def test_pointset_coords_and_json_roundtrip():
    pts = PointSet(coords=np.array([[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]]))
    assert pts.coords.shape == (2, 3)

    j = pts.to_json()
    pts2 = PointSet.from_json(j)
    assert np.allclose(pts2.coords, pts.coords)


def test_orientation_from_strike_dip_and_dimensions():
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    strike = np.array([0.0, 90.0])
    dip = np.array([30.0, 45.0])
    polarity = np.array([1, -1])

    o = OrientationObservation.from_strike_dip(coords, strike, dip, polarity)
    assert o.type == OrientationType.PLANE
    assert o.coords.shape == (2, 3)
    assert o.vector.shape == (2, 3)
    assert o.magnitude.shape[0] == 2


def test_orientation_from_dip_direction_and_plunge_variants():
    coords = np.array([[0.0, 0.0, 0.0]])
    dip_direction = np.array([45.0])
    dip = np.array([10.0])
    polarity = np.array([1])

    o2 = OrientationObservation.from_dip_direction_and_dip(coords, dip_direction, dip, polarity)
    assert o2.coords.shape == (1, 3)

    plunge_dir = np.array([120.0])
    plunge = np.array([5.0])
    o3 = OrientationObservation.from_plunge_and_plunge_direction(
        coords, plunge_dir, plunge, polarity
    )
    assert o3.coords.shape == (1, 3)


def test_lineset_to_pointset_and_tangents():
    # create two lines: first 3 points, second 4 points -> offsets [0,3,7]
    vertices = np.vstack(
        [
            np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            np.array([[0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0], [3.0, 1.0, 0.0]]),
        ]
    )
    offsets = np.array([0, 3, 7])

    ls = LineSet(vertices=vertices, offsets=offsets)

    ps = ls.to_point_set()
    assert isinstance(ps, PointSet)
    assert ps.coords.shape[1] == 3

    tangents = ls.to_tangent_vectors()
    # Expect two Orientation objects (one per segment list)
    assert isinstance(tangents, list)
    assert tangents[0].type == OrientationType.TANGENT
    # check vectors lengths
    total_vectors = sum([t.vector.shape[0] for t in tangents])
    # first segment (3 pts) -> 2 vectors, second (4 pts) -> 3 vectors
    assert total_vectors == 5
