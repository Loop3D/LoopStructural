"""Tests for the shared grid lookup of the fault frame (plan D6, task T026).

The private method ``FaultSegment._evaluate_frame(points, need_gradient=False)``
gives the three frame values and, if asked, the gradient of coordinate 1.
It is a fast path for frames with equal StructuredGrid supports. The fast path
makes one grid lookup for the three frame members. Every other frame uses the
old calls of the features (the legacy fallback).

Assumption: the method returns a tuple. The first three items are the values
of coordinates 0, 1 and 2. When ``need_gradient`` is true, the last item is the
(N, 3) gradient of coordinate 1. D6 does not define the exact return type.

The reference code in this file only uses public feature calls. It does not
use ``_evaluate_frame``. The lookups are counted with a wrapper on
``StructuredGrid._lookup_points``.
"""

from __future__ import annotations

import copy

import numpy as np
import pytest
from loop_common.supports import RectilinearGrid, StructuredGrid, TetMesh

from tests.fault_evaluation.models import build_model

POINT_COUNT = 250
# The golden tolerance of D1 is 1e-8 times the box size.
BOX_SIZE = 1000.0
ATOL = 1e-8 * BOX_SIZE


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
class LookupCounter:
    """Count the calls of ``StructuredGrid._lookup_points``."""

    def __init__(self):
        self.count = 0

    def reset(self):
        self.count = 0


@pytest.fixture
def lookups(monkeypatch):
    """Wrap the lookup of the grid. The wrapper calls the original method."""
    counter = LookupCounter()
    original = StructuredGrid._lookup_points

    def counting(self, pos):
        counter.count += 1
        return original(self, pos)

    monkeypatch.setattr(StructuredGrid, "_lookup_points", counting)
    return counter


def _fault_of(model):
    return next(feature for feature in model.features if hasattr(feature, "steps"))


def _fault_by_name(model, name):
    return next(f for f in model.features if hasattr(f, "steps") and f.name == name)


def _build(name="f1_u0", fault_name=None):
    """Return a model and one of its faults. The persistent cache is off."""
    model = build_model(name)
    cache = getattr(model, "evaluation_cache", None)
    if cache is not None:
        cache.enabled = False
    fault = _fault_by_name(model, fault_name) if fault_name else _fault_of(model)
    return model, fault


def _points(model, seed=5):
    """Return points inside and outside the box, near the fault, and one NaN row."""
    rng = np.random.default_rng(seed)
    box = model.bounding_box
    size = box.maximum - box.origin
    inside = box.origin + rng.random((POINT_COUNT, 3)) * size
    outside = box.origin + (1.2 + rng.random((20, 3))) * size
    centre = box.origin + 0.5 * size
    near = centre + (rng.random((40, 3)) - 0.5) * np.array([0.02, 1.0, 1.0]) * size
    nan_row = np.array([[np.nan, centre[1], centre[2]]])
    return np.vstack((inside, outside, near, nan_row))


def _unpack(result, need_gradient):
    """Split the result of ``_evaluate_frame`` (see the assumption above)."""
    items = tuple(result)
    gx, gy, gz = items[:3]
    gradient = items[-1] if need_gradient else None
    if need_gradient:
        assert len(items) == 4, "need_gradient=True must add the gradient as last item"
    return gx, gy, gz, gradient


def _reference_frame(fault, points, need_gradient):
    """Evaluate the frame with the public calls of the features."""
    gx = fault[0].evaluate_value(points)
    gy = fault[1].evaluate_value(points)
    gz = fault[2].evaluate_value(points)
    gradient = fault[1].evaluate_gradient(points, ignore_regions=True) if need_gradient else None
    return gx, gy, gz, gradient


def _displacement(fault, gx, gy, gz):
    """Copy of the displacement code of the legacy step loop."""
    d = np.zeros(gx.shape)
    valid = ~np.isnan(gx) & ~np.isnan(gy) & ~np.isnan(gz)
    positive = np.zeros_like(valid, dtype=bool)
    positive[valid] = gx[valid] > 0
    d[positive] = 1.0
    d[valid] = fault.faultfunction(gx[valid], gy[valid], gz[valid])
    return d * fault.displacement, valid


def _legacy_restore(fault, points, reverse=False, precomputed_gx=None):
    """Copy of the legacy ``apply_to_points`` that uses only public feature calls."""
    steps = fault.steps
    newp = np.copy(points).astype(float)
    if precomputed_gx is not None:
        gx = precomputed_gx
    else:
        gx = fault[0].evaluate_value(newp)
    gy = fault[1].evaluate_value(newp)
    gz = fault[2].evaluate_value(newp)
    d, _ = _displacement(fault, gx, gy, gz)
    # ``reverse`` only changes the sign of the start displacement. The mask is the same.
    mask = np.abs(d) > 0.0
    for i in range(steps):
        if i == 0:
            gx_m, gy_m, gz_m = gx[mask], gy[mask], gz[mask]
        else:
            gx_m = fault[0].evaluate_value(newp[mask, :])
            gy_m = fault[1].evaluate_value(newp[mask, :])
            gz_m = fault[2].evaluate_value(newp[mask, :])
        g = fault[1].evaluate_gradient(newp[mask, :], ignore_regions=True)
        # The legacy loop does not apply ``reverse`` to the step displacement.
        d, valid = _displacement(fault, gx_m, gy_m, gz_m)
        magnitude = np.zeros(g.shape[0])
        magnitude[valid] = np.linalg.norm(g[valid], axis=1)
        g[magnitude > 0.0] /= magnitude[magnitude > 0, None]
        g *= (1.0 / steps) * d[:, None]
        newp[mask, :] += g
    return newp


def _assert_close(actual, expected, atol=ATOL):
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=atol, equal_nan=True)


def _count_legacy(fault, points, lookups, **kwargs):
    """Return the lookup count of the reference restoration."""
    lookups.reset()
    _legacy_restore(fault, points, **kwargs)
    return lookups.count


# ---------------------------------------------------------------------------
# Fast path: values, gradient and shared lookup
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("need_gradient", [False, True])
def test_frame_values_and_coordinate_one_gradient_match_feature_calls(need_gradient):
    """The frame has the same values and gradient as the three features."""
    model, fault = _build()
    points = _points(model)
    expected = _reference_frame(fault, points, need_gradient)

    gx, gy, gz, gradient = _unpack(fault._evaluate_frame(points, need_gradient), need_gradient)

    for actual, wanted in zip((gx, gy, gz), expected[:3]):
        assert actual.shape == (len(points),)
        _assert_close(actual, wanted, atol=1e-12)
    # Points outside the grid and the NaN row must stay NaN.
    assert np.isnan(gx).any()
    assert np.isfinite(gx).sum() > POINT_COUNT // 2
    if need_gradient:
        assert gradient.shape == (len(points), 3)
        _assert_close(gradient, expected[3], atol=1e-12)


def test_frame_default_does_not_return_a_gradient():
    """Without ``need_gradient`` the result has no gradient item."""
    model, fault = _build()
    result = tuple(fault._evaluate_frame(_points(model)))
    assert len(result) == 3


def test_frame_accepts_empty_points():
    """No points give empty arrays."""
    _, fault = _build()
    gx, gy, gz, gradient = _unpack(fault._evaluate_frame(np.empty((0, 3)), True), True)
    assert gx.shape == gy.shape == gz.shape == (0,)
    assert gradient.shape == (0, 3)


@pytest.mark.parametrize(
    ("need_gradient", "legacy_lookups"),
    [(False, 3), (True, 4)],
)
def test_frame_makes_one_lookup_for_all_members(lookups, need_gradient, legacy_lookups):
    """The three values and the gradient share one grid lookup."""
    model, fault = _build()
    points = _points(model)
    lookups.reset()
    _reference_frame(fault, points, need_gradient)
    # This checks the counter: the legacy path does one lookup for each call.
    assert lookups.count == legacy_lookups

    lookups.reset()
    fault._evaluate_frame(points, need_gradient)
    assert lookups.count == 1


# ---------------------------------------------------------------------------
# apply_to_points: equivalence with the legacy path, all steps, masks
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("steps", [10, 3])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("use_precomputed_gx", [False, True])
def test_restoration_matches_legacy_path(steps, reverse, use_precomputed_gx):
    """Every step, the mask and the optional gx give the legacy points."""
    model, fault = _build()
    if steps != fault.steps:
        fault.steps = steps
    points = _points(model)
    gx = fault[0].evaluate_value(points) if use_precomputed_gx else None

    expected = _legacy_restore(fault, points, reverse=reverse, precomputed_gx=gx)
    actual = fault.apply_to_points(points, reverse=reverse, precomputed_gx=gx)

    _assert_close(actual, expected)
    moved = np.any(np.abs(expected - points) > 1e-9, axis=1)
    # The test must contain moved points and points that the mask keeps fixed.
    assert moved.any()
    assert (~moved).any()
    # Points that are not finite or outside the grid do not move.
    np.testing.assert_array_equal(actual[~moved], points[~moved])


def test_restoration_without_displacement_keeps_points():
    """A fault with zero displacement gives an empty mask and no change."""
    model, fault = _build()
    points = _points(model)
    fault.displacement = 0.0
    np.testing.assert_array_equal(fault.apply_to_points(points), points)


def test_precomputed_gx_gives_the_same_points_as_no_gx():
    """The initial gx can come from the caller. The result must not change."""
    model, fault = _build()
    points = _points(model)
    gx = fault[0].evaluate_value(points)
    _assert_close(
        fault.apply_to_points(points, precomputed_gx=gx), fault.apply_to_points(points), atol=ATOL
    )


# ---------------------------------------------------------------------------
# apply_to_points: lookup counts and step count
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("steps", [10, 3])
@pytest.mark.parametrize("use_precomputed_gx", [False, True])
def test_restoration_shares_one_lookup_for_each_step(lookups, steps, use_precomputed_gx):
    """Each step makes one lookup for the frame members and the gradient."""
    model, fault = _build()
    if steps != fault.steps:
        fault.steps = steps
    points = _points(model)
    gx = fault[0].evaluate_value(points) if use_precomputed_gx else None

    legacy = _count_legacy(fault, points, lookups, precomputed_gx=gx)
    # This checks the counter: the legacy path does about four lookups each step.
    assert legacy >= 3 * steps

    lookups.reset()
    fault.apply_to_points(points, precomputed_gx=gx)
    # One lookup for the initial frame, one for the gradient of step 0 and one
    # for each of the other steps. No lookup is kept after the points move.
    assert lookups.count <= steps + 2
    assert lookups.count < legacy


def test_integration_step_count_is_unchanged(monkeypatch):
    """The fast path does one displacement call for the start and for every step."""
    model, fault = _build()
    points = _points(model)
    calls = []
    original = fault.faultfunction

    def counting(gx, gy, gz):
        calls.append(len(gx))
        return original(gx, gy, gz)

    monkeypatch.setattr(fault, "_faultfunction", counting)

    _legacy_restore(fault, points)
    legacy_calls = len(calls)
    assert legacy_calls == fault.steps + 1

    calls.clear()
    fault.apply_to_points(points)
    assert len(calls) == legacy_calls == 11
    assert fault.steps == 10


# ---------------------------------------------------------------------------
# Legacy fallback
# ---------------------------------------------------------------------------
class _SubstituteInterpolator:
    """Interpolator that reports another support but evaluates like the inner one."""

    def __init__(self, inner, support):
        self._inner = inner
        self.support = support

    def __getattr__(self, name):
        return getattr(self._inner, name)


class _NonUniformGrid(StructuredGrid):
    """Grid with its own local coordinates. It cannot use the uniform lookup."""

    def position_to_local_coordinates(self, pos):
        return super().position_to_local_coordinates(pos)


def _case_geometry():
    model, fault = _build()
    support = copy.deepcopy(fault[1].interpolator.support)
    support._origin = support.origin + 0.3 * support.step_vector
    support.onGeometryChange()
    fault[1].interpolator.support = support
    return model, fault


def _case_transform():
    model, fault = _build()
    box = copy.deepcopy(fault[2].interpolator.bounding_box)
    box.set_local_transform(local_origin=np.asarray(box.local_origin) + [5.0, 5.0, 0.0])
    fault[2].interpolator.bounding_box = box
    return model, fault


def _case_rectilinear():
    model, fault = _build()
    inner = fault[1].interpolator
    nodes = [np.linspace(0.0, 1.0, 4)] * 3
    fault[1].interpolator = _SubstituteInterpolator(inner, RectilinearGrid(*nodes))
    return model, fault


def _case_tetra():
    model, fault = _build()
    inner = fault[2].interpolator
    mesh = TetMesh(origin=np.zeros(3), nsteps=np.array([3, 3, 3]), step_vector=np.ones(3))
    fault[2].interpolator = _SubstituteInterpolator(inner, mesh)
    return model, fault


def _case_non_uniform_subclass():
    model, fault = _build()
    fault[0].interpolator.support.__class__ = _NonUniformGrid
    return model, fault


def _case_child_faults():
    return _build("cutting_faults", "fault_1")


def _case_opaque_region():
    model, fault = _build()
    centre_y = float(np.mean([model.bounding_box.origin[1], model.bounding_box.maximum[1]]))
    fault[0].regions.append(lambda pts: np.asarray(pts)[:, 1] > 0.2 * centre_y)
    return model, fault


def _case_domain_region():
    return _build("onlap_domain", "fault_0")


FALLBACK_CASES = {
    "different_geometry": _case_geometry,
    "different_transform": _case_transform,
    "rectilinear_support": _case_rectilinear,
    "tetrahedral_support": _case_tetra,
    "non_uniform_grid_subclass": _case_non_uniform_subclass,
    "child_faults": _case_child_faults,
    "opaque_region_on_one_member": _case_opaque_region,
    "domain_region": _case_domain_region,
}


@pytest.fixture(params=sorted(FALLBACK_CASES))
def fallback_case(request):
    return FALLBACK_CASES[request.param]()


def test_fallback_restoration_equals_legacy_with_legacy_lookup_count(lookups, fallback_case):
    """A frame that cannot use the fast path gives the legacy points and lookups."""
    model, fault = fallback_case
    points = _points(model)
    legacy_count = _count_legacy(fault, points, lookups)
    expected = _legacy_restore(fault, points)

    lookups.reset()
    actual = fault.apply_to_points(points)

    _assert_close(actual, expected)
    assert lookups.count == legacy_count
    # The lookups must be real: the frame is evaluated for a child fault too.
    assert legacy_count > 0


@pytest.mark.parametrize("need_gradient", [False, True])
def test_fallback_frame_equals_legacy_with_legacy_lookup_count(
    lookups, fallback_case, need_gradient
):
    """``_evaluate_frame`` uses the feature calls when the frame is not eligible."""
    model, fault = fallback_case
    # The gradient of a feature with child faults needs finite points inside the box.
    points = _points(model)[:POINT_COUNT]
    lookups.reset()
    expected = _reference_frame(fault, points, need_gradient)
    legacy_count = lookups.count

    lookups.reset()
    gx, gy, gz, gradient = _unpack(fault._evaluate_frame(points, need_gradient), need_gradient)

    for actual, wanted in zip((gx, gy, gz), expected[:3]):
        _assert_close(actual, wanted, atol=1e-12)
    if need_gradient:
        _assert_close(gradient, expected[3], atol=1e-12)
    assert lookups.count == legacy_count
