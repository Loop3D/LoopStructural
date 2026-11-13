"""Kinematics of fault restoration: forward movement and displaced faults."""

import numpy as np
import pandas as pd
import pytest

from LoopStructural import GeologicalModel
from tests.fixtures.fault_models import (
    COLUMNS,
    OLD,
    YOUNG,
    add_fault,
    dip_vectors,
    interior_points,
    plane_points,
    young_forward,
    young_side,
)


@pytest.fixture(scope="module")
def planar_fault():
    """Planar normal fault, dip 60 toward +x, constant slip 0.15 down dip."""
    model = GeologicalModel(np.zeros(3), np.ones(3))
    return add_fault(model, "f1", (0.5, 0.5, 0.5), 60, 90, 0.15)


def test_reverse_moves_the_hanging_wall_forward(planar_fault):
    _, down_dip, _ = dip_vectors(60, 90)
    pts = interior_points(n=200, lo=0.3, hi=0.7)
    hanging = planar_fault.evaluate_frame_value(pts) > 0.05
    moved = planar_fault.apply_to_points(pts, reverse=True)
    # Forward: the hanging wall moves 0.15 down dip; the footwall stays.
    expected = np.tile(0.15 * down_dip, (int(hanging.sum()), 1))
    np.testing.assert_allclose(moved[hanging] - pts[hanging], expected, atol=2e-3)
    footwall = planar_fault.evaluate_frame_value(pts) < 0
    np.testing.assert_array_equal(moved[footwall], pts[footwall])


def _restore_exactly(points):
    """Restore present-day points through young, then old (both translations)."""
    points = np.array(points, dtype=float)
    _, young_dip, _ = dip_vectors(YOUNG["dip"], YOUNG["dip_azimuth"])
    points[young_side(points) > 0] -= YOUNG["displacement"] * young_dip
    old_normal, old_dip, _ = dip_vectors(OLD["dip"], OLD["dip_azimuth"])
    old_side = (points - np.asarray(OLD["center"])) @ old_normal
    points[old_side > 0] -= OLD["displacement"] * old_dip
    return points, old_side


@pytest.fixture(scope="module")
def displaced_series():
    """A series cut by an old fault that a young fault displaces.

    The values of the data are the analytic field z of each point restored
    once through each fault, the kinematics that the fault frames describe.
    """
    model = GeologicalModel(np.zeros(3), np.ones(3))
    young = add_fault(model, "young", YOUNG["center"], YOUNG["dip"], YOUNG["dip_azimuth"],
                      YOUNG["displacement"])
    plane, _, _ = plane_points(OLD["center"], OLD["dip"], OLD["dip_azimuth"], n=15)
    keep = np.abs(young_side(plane)) > 0.04
    center = young_forward(np.asarray([OLD["center"]]))[0]
    old = add_fault(model, "old", center, OLD["dip"], OLD["dip_azimuth"], OLD["displacement"],
                    points=young_forward(plane)[keep], faults=[young])
    pts = np.random.default_rng(4).uniform(0.0, 1.0, size=(1500, 3))
    restored, _ = _restore_exactly(pts)
    df = pd.DataFrame(
        {"X": pts[:, 0], "Y": pts[:, 1], "Z": pts[:, 2], "val": restored[:, 2],
         "nx": np.nan, "ny": np.nan, "nz": np.nan, "coord": np.nan, "feature_name": "strat"}
    )[COLUMNS]
    model.create_and_add_foliation("strat", data=df, faults=[old, young], nelements=8000,
                                   buffer=0.3)
    return model


def test_displaced_fault_restores_a_series_once(displaced_series):
    pts = interior_points(n=3000, lo=0.1, hi=0.9, seed=8)
    exact, old_side = _restore_exactly(pts)
    value = displaced_series.evaluate_feature_value("strat", pts)
    # Away from both faults, the field is z of the point restored once through each fault.
    far = (np.abs(young_side(pts)) > 0.03) & (np.abs(old_side) > 0.03)
    np.testing.assert_allclose(value[far], exact[far, 2], atol=5e-3)


def test_reverse_then_restore_returns_the_points(planar_fault):
    pts = interior_points(n=200, lo=0.3, hi=0.7)
    hanging = planar_fault.evaluate_frame_value(pts) > 0.05
    back = planar_fault.apply_to_points(planar_fault.apply_to_points(pts, reverse=True))
    np.testing.assert_allclose(back[hanging], pts[hanging], atol=2e-3)
