"""FaultDisplacement.tip_ends: which frame coordinates carry a tip."""

import numpy as np

from LoopStructural.modelling.features.fault._fault_function import (
    BaseFault,
    BaseFault3D,
    Composite,
    CubicFunction,
    FaultDisplacement,
    Ones,
    Zeros,
    smooth_peak,
)


def test_library_profiles_declare_their_tips():
    assert BaseFault.fault_displacement.tip_ends == {1: (), 2: (1, -1)}
    assert BaseFault3D.fault_displacement.tip_ends == {1: (1, -1), 2: (1, -1)}
    constant = FaultDisplacement(gx=Composite(Ones(), Zeros()), gy=Ones(), gz=Ones())
    assert constant.tip_ends == {1: (), 2: ()}
    assert smooth_peak.tip_ends == (1, -1)


def test_cubic_profile_with_a_zero_limit_has_a_tip_at_that_end():
    f = CubicFunction()
    f.add_cstr(-1, 1)
    f.add_cstr(0, 1)
    f.add_cstr(1, 0)
    f.add_grad(1, 0)
    f.set_lim(-1, 1)
    assert f.tip_ends == (1,)
    assert np.all(f(np.array([1.0, 1.5, 3.0])) == 0.0)


def test_base_fault_displacement_is_exactly_zero_beyond_its_decay():
    # The hanging-wall and footwall cubics end at 0 at |gx| = 1: exactly 0 beyond.
    gx = np.array([1.0, 1.2, 5.0, -1.0, -1.2, -5.0])
    d = BaseFault.fault_displacement(gx, np.zeros(6), np.zeros(6))
    np.testing.assert_array_equal(d, 0.0)


def test_undeclared_profile_is_tested():
    step = FaultDisplacement(gx=Ones(), gy=lambda v: (np.abs(v) < 1) * 1.0, gz=Ones())
    assert step.tip_ends == {1: (1, -1), 2: ()}
