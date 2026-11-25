"""An unconformity truncates the faults that are older than it (LS-D2)."""

import numpy as np
import pandas as pd
import pytest

from LoopStructural import GeologicalModel
from LoopStructural.modelling.features import UnconformityFeature
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

TILT = np.radians(20.0)
UNCONFORMITY_Z = 0.6


def cover_field(p):
    return p[:, 2] - UNCONFORMITY_Z


def basement_field(p):
    return -np.sin(TILT) * p[:, 0] + np.cos(TILT) * p[:, 2] - 0.1


def series(model, name, field, faults=None, n=500, seed=5):
    """A series whose data are restored through its faults (given oldest first).

    ``faults=None`` leaves the fault list to the model: every fault that the
    model has at this moment.
    """
    pts = np.random.default_rng(seed).uniform(0.0, 1.0, size=(n, 3))
    use = faults
    if use is None:
        use = [f for f in model.features if f.type.name == "FAULT"]
    restored = pts.copy()
    for f in reversed(use):
        restored = f.apply_to_points(restored)
    df = pd.DataFrame(
        {
            "X": pts[:, 0], "Y": pts[:, 1], "Z": pts[:, 2], "val": field(restored),
            "nx": np.nan, "ny": np.nan, "nz": np.nan, "coord": np.nan, "feature_name": name,
        }
    )[COLUMNS]
    kwargs = {} if faults is None else {"faults": faults}
    return model.create_and_add_foliation(name, data=df, nelements=4000, buffer=0.3, **kwargs)


def fault(model, name, spec, **kwargs):
    return add_fault(model, name, spec["center"], spec["dip"], spec["dip_azimuth"],
                     spec["displacement"], **kwargs)


def new_model():
    return GeologicalModel(np.zeros(3), np.ones(3))


def truncations(fault_feature):
    return [r for r in fault_feature[0].regions if isinstance(r, UnconformityFeature)]


def probes():
    """Points in the hanging wall of the OLD fault, above and below the unconformity."""
    normal, _, _ = dip_vectors(OLD["dip"], OLD["dip_azimuth"])
    rng = np.random.default_rng(1)
    pts = rng.uniform(0.15, 0.85, size=(4000, 3))
    side = (pts - np.asarray(OLD["center"])) @ normal
    pts = pts[side > 0.05]
    above = pts[pts[:, 2] > UNCONFORMITY_Z + 0.1]
    below = pts[pts[:, 2] < UNCONFORMITY_Z - 0.2]
    return above, below


def assert_truncated_by_cover(model, name):
    f = model[name]
    regions = truncations(f)
    assert [r.parent.name for r in regions] == ["cover"]
    assert regions[0].value == 0.0 and regions[0].sign is True
    # the three frame features share the regions of the fault
    assert all(f[i].regions is f[0].regions for i in range(3))
    above, below = probes()
    # the fault surface and frame end on the unconformity
    assert np.all(np.isnan(f[0].evaluate_value(above)))
    assert np.all(np.isfinite(f[0].evaluate_value(below)))
    # the fault moves no point above the unconformity
    np.testing.assert_array_equal(f.apply_to_points(above), above)
    # below it the hanging wall is restored by the displacement
    moved = np.linalg.norm(f.apply_to_points(below) - below, axis=1)
    np.testing.assert_allclose(moved, f.displacement, rtol=0.05)
    relation = model.get_feature_relations()["faults"][name]["truncated_by"]
    assert relation == [
        {"unconformity": "__cover_unconformity", "feature": "cover", "value": 0.0,
         "keep": "below"}
    ]


def assert_not_truncated(model, name):
    f = model[name]
    assert truncations(f) == []
    above, below = probes()
    for pts in (above, below):
        assert np.all(np.isfinite(f[0].evaluate_value(pts)))
        moved = np.linalg.norm(f.apply_to_points(pts) - pts, axis=1)
        np.testing.assert_allclose(moved, f.displacement, rtol=0.05)


def test_fault_left_out_of_the_cover_is_older():
    """All faults first: the fault list of the cover says which fault is older."""
    model = new_model()
    post = fault(model, "f_post", YOUNG)
    pre = fault(model, "f_pre", OLD, faults=[post])
    cover = series(model, "cover", cover_field, faults=[post])
    model.add_unconformity(cover, 0.0)
    basement = series(model, "basement", basement_field, faults=[pre, post])
    model.update()
    assert_truncated_by_cover(model, "f_pre")
    # the younger fault displaces the cover and is not cut
    assert truncations(model["f_post"]) == []
    assert [f.name for f in cover.faults] == ["f_post"]
    # the basement is still displaced by the older fault
    assert sorted(f.name for f in basement.faults) == ["f_post", "f_pre"]


def test_fault_added_after_the_cover_and_before_the_call_is_older():
    model = new_model()
    cover = series(model, "cover", cover_field)
    fault(model, "f_pre", OLD)
    model.add_unconformity(cover, 0.0)
    series(model, "basement", basement_field)
    model.update()
    assert_truncated_by_cover(model, "f_pre")
    assert cover.faults == []
    # the cover is flat on both sides of the fault trace above the unconformity
    above, _ = probes()
    np.testing.assert_allclose(cover.evaluate_value(above), cover_field(above), atol=2e-2)


def test_fault_added_after_the_call_is_older():
    """The behaviour of the released versions: the new fault takes the region."""
    model = new_model()
    cover = series(model, "cover", cover_field)
    model.add_unconformity(cover, 0.0)
    fault(model, "f_pre", OLD)
    series(model, "basement", basement_field)
    model.update()
    assert_truncated_by_cover(model, "f_pre")


def test_fault_that_displaces_the_cover_is_younger():
    model = new_model()
    fault(model, "f_young", OLD)
    cover = series(model, "cover", cover_field)
    model.add_unconformity(cover, 0.0)
    series(model, "basement", basement_field)
    model.update()
    assert [f.name for f in cover.faults] == ["f_young"]
    assert_not_truncated(model, "f_young")


def test_model_without_an_unconformity():
    model = new_model()
    fault(model, "f", OLD)
    series(model, "basement", basement_field)
    model.update()
    assert_not_truncated(model, "f")


def test_fault_that_displaces_a_younger_series_is_not_truncated():
    """Left out of the cover, but it displaces a series above an earlier unconformity."""
    model = new_model()
    f = fault(model, "f", OLD)
    young = series(model, "young", lambda p: p[:, 2] - 0.9)  # displaced by f (model default)
    model.add_unconformity(young, 0.0)
    cover = series(model, "cover", cover_field, faults=[])
    model.add_unconformity(cover, 0.0)
    assert truncations(f) == []


def test_fault_that_displaces_only_truncated_series_is_older():
    """Basement first: the unconformity truncates the basement and the fault that cuts it."""
    model = new_model()
    f = fault(model, "f_pre", OLD)
    series(model, "basement", basement_field)  # displaced by f_pre (model default)
    cover = series(model, "cover", cover_field, faults=[], seed=6)
    model.add_unconformity(cover, 0.0)
    model.update()
    assert_truncated_by_cover(model, "f_pre")


def test_older_fault_is_truncated_by_every_younger_unconformity():
    model = new_model()
    old = fault(model, "f_old", OLD)
    top = series(model, "cover", cover_field, faults=[])
    model.add_unconformity(top, 0.0)
    middle = series(model, "middle", lambda p: p[:, 2] - 0.4, faults=[])
    model.add_unconformity(middle, 0.0)
    late = fault(model, "f_late", YOUNG)  # added after both: older than both
    series(model, "basement", basement_field, faults=[old])
    for f in (old, late):
        assert sorted(r.parent.name for r in truncations(f)) == ["cover", "middle"]


def test_onlap_unconformity_does_not_truncate_faults():
    """Not changed by this draft: no fault takes the region of an onlap unconformity."""
    model = new_model()
    pre = fault(model, "f_pre", OLD)
    cover = series(model, "cover", cover_field, faults=[])
    model.add_onlap_unconformity(cover, 0.0)
    late = fault(model, "f_late", YOUNG)
    assert truncations(pre) == []
    assert truncations(late) == []


def test_the_region_is_added_once():
    model = new_model()
    pre = fault(model, "f_pre", OLD)
    cover = series(model, "cover", cover_field, faults=[])
    model.add_unconformity(cover, 0.0)
    model.add_unconformity(cover, 0.0)
    assert len(truncations(pre)) == 1


def test_features_built_before_the_call_are_built_again():
    """A series restored through the fault before the fault was truncated is rebuilt."""
    model = new_model()
    pre = fault(model, "f_pre", OLD)
    basement = series(model, "basement", basement_field, faults=[pre])
    model.update()
    assert basement.builder._up_to_date
    cover = series(model, "cover", cover_field, faults=[], seed=6)
    model.add_unconformity(cover, 0.0)
    assert not basement.builder._up_to_date


@pytest.mark.parametrize("value", [0.0, 0.1])
def test_the_region_uses_the_value_of_the_unconformity(value):
    model = new_model()
    pre = fault(model, "f_pre", OLD)
    cover = series(model, "cover", cover_field, faults=[])
    model.add_unconformity(cover, value)
    assert [r.value for r in truncations(pre)] == [value]


def _restore_once(points):
    """Present-day points restored through the young fault, then through the old
    fault where the restored point is below the unconformity (both translations)."""
    points = np.array(points, dtype=float)
    _, young_dip, _ = dip_vectors(YOUNG["dip"], YOUNG["dip_azimuth"])
    points[young_side(points) > 0] -= YOUNG["displacement"] * young_dip
    once = points.copy()
    old_normal, old_dip, _ = dip_vectors(OLD["dip"], OLD["dip_azimuth"])
    old_side = (points - np.asarray(OLD["center"])) @ old_normal
    points[(old_side > 0) & (cover_field(once) <= 0)] -= OLD["displacement"] * old_dip
    return points, once, old_side


def test_truncated_fault_that_a_younger_fault_displaces():
    """The old fault ends on the unconformity in the space restored by the young fault.

    The frame of the old fault is built in that space, so its region is read
    there too: at the point restored once through the young fault.
    """
    model = new_model()
    young = fault(model, "young", YOUNG)
    plane, _, _ = plane_points(OLD["center"], OLD["dip"], OLD["dip_azimuth"], n=15)
    keep = np.abs(young_side(plane)) > 0.04
    center = young_forward(np.asarray([OLD["center"]]))[0]
    old = add_fault(model, "old", center, OLD["dip"], OLD["dip_azimuth"], OLD["displacement"],
                    points=young_forward(plane)[keep], faults=[young])

    def add_series(name, values, pts, faults):
        df = pd.DataFrame(
            {
                "X": pts[:, 0], "Y": pts[:, 1], "Z": pts[:, 2], "val": values,
                "nx": np.nan, "ny": np.nan, "nz": np.nan, "coord": np.nan, "feature_name": name,
            }
        )[COLUMNS]
        return model.create_and_add_foliation(name, data=df, faults=faults, nelements=8000,
                                              buffer=0.3)

    pts = np.random.default_rng(4).uniform(0.0, 1.0, size=(1500, 3))
    restored, once, _ = _restore_once(pts)
    cover = add_series("cover", cover_field(once), pts, [young])
    model.add_unconformity(cover, 0.0)
    add_series("basement", basement_field(restored), pts, [old, young])
    model.update()
    assert [r.parent.name for r in truncations(old)] == ["cover"]

    probe = interior_points(n=4000, lo=0.1, hi=0.9, seed=8)
    exact, once, old_side = _restore_once(probe)
    value = model.evaluate_feature_value("basement", probe)
    # below the unconformity, away from it and from both faults
    far = (
        (cover_field(once) < -0.03)
        & (np.abs(young_side(probe)) > 0.03)
        & (np.abs(old_side) > 0.03)
    )
    assert far.sum() > 500
    np.testing.assert_allclose(value[far], basement_field(exact[far]), atol=5e-3)


def test_onlap_unconformity_keeps_the_regions_the_feature_already_has():
    """An onlap unconformity on a series must not remove the erosional region of an older one."""
    model = new_model()
    old = series(model, "old", basement_field, faults=[])
    model.add_unconformity(old, 0.0)
    young = series(model, "young", cover_field, faults=[], seed=6)
    before = [r.name for r in young.regions]
    assert before
    model.add_onlap_unconformity(young, 0.0)
    assert set(before) <= {r.name for r in young.regions}
