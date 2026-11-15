"""The displacement of a block on the other side of its fault (LS-A2, decision D1).

A block field is evaluated on both sides of the fault. On the other side the
displacement is the smooth continuation of the profile of the block's side
across gx = 0, so the block field is C2 across the fault surface.
"""

import numpy as np
import pytest

from LoopStructural import GeologicalModel
from LoopStructural.modelling.features.fault._fault_function import (
    BaseFault,
    BaseFault3D,
    Composite,
    CubicFunction,
    FaultDisplacement,
    Ones,
    block_displacement,
)
from LoopStructural.utils import _fault_block
from LoopStructural.utils._fault_block import (
    BLOCK_EXTENSION_EXACT,
    BLOCK_EXTENSION_LIMIT,
    force_side,
)
from tests.fixtures.fault_models import add_fault, add_layers, constant_profile

GY, GZ = 0.3, -0.2


def cubic(value, slope, end):
    """A cubic with a value and a slope at the fault that is flat and zero at ``end``."""
    f = CubicFunction()
    f.add_cstr(0, value)
    f.add_grad(0, slope)
    f.add_cstr(end, 0)
    f.add_grad(end, 0)
    f.set_lim(min(0, end), max(0, end))
    return f


def sloped_cubic_profile():
    """A user profile with a slope at the fault on both sides."""
    return FaultDisplacement(hw=cubic(1.0, -1.5, 1), fw=cubic(-0.6, -0.8, -1), scale=1.0)


def sloped_callable_profile():
    """A user gx profile that is not a polynomial: the second-order expansion is used."""

    def gx(v):
        v = np.asarray(v, dtype=float)
        return np.where(v > 0, np.exp(-2.0 * v), -0.5 * np.exp(3.0 * v))

    return FaultDisplacement(gx=gx, gy=Ones(), gz=Ones(), scale=1.0)


def plain_function(gx, gy, gz):
    """A fault function that is not a FaultDisplacement."""
    return np.where(gx > 0, 1.0 / (1.0 + 2.0 * gx), 0.2 * np.cos(gx) - 0.3 * gx) * (
        1.0 - 0.1 * gy**2
    )


PROFILES = {
    "constant": constant_profile,
    "BaseFault": lambda: BaseFault().fault_displacement,
    "BaseFault3D": lambda: BaseFault3D().fault_displacement,
    "user cubic with a slope": sloped_cubic_profile,
    "user callable with a slope": sloped_callable_profile,
    "plain function": lambda: plain_function,
}
SLOPED = ("user cubic with a slope", "user callable with a slope", "plain function")


def displacement(profile, x, side):
    x = np.asarray(x, dtype=float)
    return block_displacement(
        profile, x, force_side(x, side), np.full(x.shape, GY), np.full(x.shape, GZ), side
    )


def one_sided(profile, side, direction, h=5e-4):
    """Value, slope and curvature of the block displacement at the fault, from one side."""
    x = direction * (1e-9 + h * np.arange(4))
    f = displacement(profile, x, side)
    slope = direction * (-3 * f[0] + 4 * f[1] - f[2]) / (2 * h)
    curvature = (2 * f[0] - 5 * f[1] + 4 * f[2] - f[3]) / h**2
    return f[0], slope, curvature


def test_library_cubics_at_the_fault():
    """Value, slope, curvature and third derivative of the library gx profiles."""
    for fault in (BaseFault, BaseFault3D):
        for piece, expected in ((fault.hw, [2, -3, 0, 1]), (fault.fw, [2, 3, 0, -1])):
            piece.solve()
            np.testing.assert_allclose(piece.w, expected, atol=1e-12)
            value, slope, curvature, third = piece.w[3], piece.w[2], 2 * piece.w[1], 6 * piece.w[0]
            assert slope == pytest.approx(0.0, abs=1e-12)
            assert abs(curvature) == pytest.approx(6.0)
            assert third == pytest.approx(12.0)
            assert abs(value) == pytest.approx(1.0)


@pytest.mark.parametrize("name", PROFILES)
@pytest.mark.parametrize("side", [1, -1])
def test_block_displacement_is_c2_at_the_fault(name, side):
    """Value (1e-8), slope (1e-5) and curvature (1e-3) agree from both sides.

    The two values are taken 1e-9 from the fault on each side.
    """
    profile = PROFILES[name]()
    own = one_sided(profile, side, side)
    other = one_sided(profile, side, -side)
    assert other[0] == pytest.approx(own[0], abs=1e-8)
    assert other[1] == pytest.approx(own[1], abs=1e-5)
    assert other[2] == pytest.approx(own[2], abs=1e-3)


@pytest.mark.parametrize("name", SLOPED)
def test_the_mirror_image_has_a_kink_where_the_profile_has_a_slope(name, monkeypatch):
    profile = PROFILES[name]()
    own = one_sided(profile, 1, 1)
    assert abs(own[1]) > 0.5
    monkeypatch.setattr(_fault_block, "BLOCK_EXTENSION", "mirror")
    mirrored = one_sided(profile, 1, -1)
    # an even extension: the slope changes sign at the fault
    assert mirrored[1] == pytest.approx(-own[1], abs=1e-5)


@pytest.mark.parametrize("name", PROFILES)
@pytest.mark.parametrize("side", [1, -1])
def test_the_block_side_is_not_changed(name, side, monkeypatch):
    """On the block's own side the displacement is the profile, bit for bit."""
    profile = PROFILES[name]()
    x = side * np.concatenate([[0.0, 1e-12], np.linspace(1e-6, 1.5, 400)])
    y, z = np.full(x.shape, GY), np.full(x.shape, GZ)
    expected = profile(force_side(x, side), y, z)
    np.testing.assert_array_equal(displacement(profile, x, side), expected)
    monkeypatch.setattr(_fault_block, "BLOCK_EXTENSION", "mirror")
    np.testing.assert_array_equal(displacement(profile, x, side), expected)


@pytest.mark.parametrize("side", [1, -1])
def test_cubic_continuation_is_the_cubic_itself(side, monkeypatch):
    """The library cubic at the signed coordinate; the mirror image differs by 4|x|^3."""
    profile = BaseFault().fault_displacement
    x = -side * np.linspace(1e-6, BLOCK_EXTENSION_EXACT, 50)
    y, z = np.full(x.shape, GY), np.full(x.shape, GZ)
    w = (profile.gx.positive if side > 0 else profile.gx.negative).w
    cubic_value = w[0] * x**3 + w[1] * x**2 + w[2] * x + w[3]
    rest = profile.scale * profile.gy(y) * profile.gz(z)
    smooth = displacement(profile, x, side)
    np.testing.assert_allclose(smooth, rest * cubic_value, atol=1e-14)
    monkeypatch.setattr(_fault_block, "BLOCK_EXTENSION", "mirror")
    mirror = displacement(profile, x, side)
    np.testing.assert_allclose(smooth - mirror, rest * 4 * x**3, atol=1e-13)


@pytest.mark.parametrize("name", PROFILES)
@pytest.mark.parametrize("side", [1, -1])
def test_the_continuation_is_constant_beyond_the_limit(name, side):
    profile = PROFILES[name]()
    far = -side * np.array([BLOCK_EXTENSION_LIMIT, 0.7, 1.0, 5.0, 100.0])
    values = displacement(profile, far, side)
    np.testing.assert_allclose(values, values[0], rtol=0, atol=1e-12)
    # it stays between the values that the exact part can take: no growth without limit
    near = displacement(profile, -side * np.linspace(1e-6, BLOCK_EXTENSION_LIMIT, 200), side)
    assert np.all(np.isfinite(near))
    assert abs(values[0]) <= np.max(np.abs(near)) + 1e-12


@pytest.mark.parametrize("name", PROFILES)
def test_the_blend_is_c2(name):
    """No jump of the value, the slope or the curvature at the two ends of the blend."""
    profile = PROFILES[name]()
    h = 1e-4
    for joint in (BLOCK_EXTENSION_EXACT, BLOCK_EXTENSION_LIMIT):
        x = -(joint + h * np.arange(-3, 4))
        f = displacement(profile, x, 1)
        left = (2 * f[3] - 5 * f[2] + 4 * f[1] - f[0]) / h**2
        right = (2 * f[3] - 5 * f[4] + 4 * f[5] - f[6]) / h**2
        assert right == pytest.approx(left, abs=5e-2)
        assert (f[4] - f[3]) / h == pytest.approx((f[3] - f[2]) / h, abs=1e-2)


# ---------------------------------------------------------------------------
# a model: the block field across the fault surface


@pytest.fixture(scope="module")
def sloped_model():
    """Planar normal fault (dip 60) whose profile has a slope at the fault."""
    model = GeologicalModel(np.zeros(3), np.ones(3))
    fault = add_fault(model, "f1", (0.5, 0.5, 0.5), 60, 90, 0.15,
                      profile=sloped_cubic_profile(), axes=(4.0, 4.0, 1.0))
    add_layers(model, [fault])
    model.update()
    return model


def crossing(model, y=0.5, z=0.5):
    """The x at which the line (x, y, z) crosses the fault."""
    lo, hi = 0.2, 0.8
    fault = model["f1"]

    def gx(x):
        return fault.evaluate_frame_value(np.array([[x, y, z]]))[0]

    assert gx(lo) * gx(hi) < 0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if gx(lo) * gx(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def line(x, y=0.5, z=0.5):
    x = np.atleast_1d(x)
    return np.column_stack([x, np.full(x.shape, y), np.full(x.shape, z)])


def slope_jump(model, block, x0, d=1e-3):
    """Jump of the slope of the block field along the line, across the fault."""
    v = model.evaluate_feature_block_value(
        "strat", line(x0 + d * np.array([-2.0, -1.0, 1.0, 2.0])), block
    )
    return abs((v[3] - v[2]) / d - (v[1] - v[0]) / d)


@pytest.mark.parametrize("side", [1, -1])
def test_block_field_has_a_continuous_gradient_across_the_fault(sloped_model, side, monkeypatch):
    x0 = crossing(sloped_model)
    block = {"f1": side}
    assert slope_jump(sloped_model, block, x0) < 5e-3
    # the test sees a kink: the mirror image gives one
    monkeypatch.setattr(_fault_block, "BLOCK_EXTENSION", "mirror")
    assert slope_jump(sloped_model, block, x0) > 5e-2


@pytest.mark.parametrize("side", [1, -1])
def test_block_gradient_uses_the_same_continuation(sloped_model, side, monkeypatch):
    """evaluate_feature_block_gradient is the gradient of the block value on both sides."""
    x0 = crossing(sloped_model)
    block = {"f1": side}
    d = 0.01
    points = line(x0 + d * np.array([-3.0, -2.0, -1.0, 1.0, 2.0, 3.0]))
    gradient = sloped_model.evaluate_feature_block_gradient("strat", points, block)
    h = 1e-4
    finite_difference = np.empty(points.shape)
    for axis in range(3):
        step = np.zeros(3)
        step[axis] = h
        finite_difference[:, axis] = (
            sloped_model.evaluate_feature_block_value("strat", points + step, block)
            - sloped_model.evaluate_feature_block_value("strat", points - step, block)
        ) / (2 * h)
    np.testing.assert_allclose(gradient, finite_difference, atol=2e-2)
    # continuous across the fault surface, and not with the mirror image
    near = line(x0 + 1e-3 * np.array([-1.0, 1.0]))
    smooth = sloped_model.evaluate_feature_block_gradient("strat", near, block)
    assert np.max(np.abs(smooth[1] - smooth[0])) < 1e-2
    monkeypatch.setattr(_fault_block, "BLOCK_EXTENSION", "mirror")
    mirror = sloped_model.evaluate_feature_block_gradient("strat", near, block)
    assert np.max(np.abs(mirror[1] - mirror[0])) > 5e-2


@pytest.mark.parametrize("side", [1, -1])
def test_block_field_inside_its_block_is_not_changed(sloped_model, side, monkeypatch):
    """Bit-equal to the mirror image (the first draft) and to the plain field."""
    points = np.random.default_rng(2).uniform(0.1, 0.9, size=(3000, 3))
    gx = sloped_model["f1"].evaluate_frame_value(points)
    inside = side * gx > 0
    block = {"f1": side}
    smooth = sloped_model.evaluate_feature_block_value("strat", points, block)
    plain = sloped_model.evaluate_feature_value("strat", points)
    np.testing.assert_array_equal(smooth[inside], plain[inside])
    monkeypatch.setattr(_fault_block, "BLOCK_EXTENSION", "mirror")
    mirror = sloped_model.evaluate_feature_block_value("strat", points, block)
    np.testing.assert_array_equal(smooth[inside], mirror[inside])
    assert np.max(np.abs(smooth[~inside] - mirror[~inside])) > 1e-3


def test_constant_profile_is_the_same_with_both_extensions(monkeypatch):
    """A constant-slip block: the continuation and the mirror image are equal."""
    model = GeologicalModel(np.zeros(3), np.ones(3))
    fault = add_fault(model, "f1", (0.5, 0.5, 0.5), 60, 90, 0.15)
    add_layers(model, [fault])
    points = np.random.default_rng(2).uniform(0.1, 0.9, size=(2000, 3))
    smooth = [model.evaluate_feature_block_value("strat", points, {"f1": s}) for s in (1, -1)]
    monkeypatch.setattr(_fault_block, "BLOCK_EXTENSION", "mirror")
    mirror = [model.evaluate_feature_block_value("strat", points, {"f1": s}) for s in (1, -1)]
    for a, b in zip(smooth, mirror):
        np.testing.assert_array_equal(a, b)
