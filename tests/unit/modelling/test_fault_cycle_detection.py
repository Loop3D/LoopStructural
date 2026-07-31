import pytest

from LoopStructural.modelling.features import LambdaGeologicalFeature
from LoopStructural.utils import LoopValueError


def _feature(name):
    return LambdaGeologicalFeature(name=name)


def test_fault_chain_without_cycle_is_allowed():
    a, b, c = _feature("a"), _feature("b"), _feature("c")
    a.faults = [b]
    b.faults = [c]
    assert a.faults == [b]
    assert b.faults == [c]


def test_feature_cannot_be_its_own_fault():
    a = _feature("a")
    with pytest.raises(LoopValueError):
        a.faults = [a]


def test_mutual_fault_cycle_is_rejected():
    a, b = _feature("a"), _feature("b")
    a.faults = [b]
    with pytest.raises(LoopValueError):
        b.faults = [a]
    # the valid assignment made before the cycle was attempted should be unaffected
    assert a.faults == [b]
    assert b.faults == []


def test_transitive_fault_cycle_is_rejected():
    x, y, z = _feature("x"), _feature("y"), _feature("z")
    x.faults = [y]
    y.faults = [z]
    with pytest.raises(LoopValueError):
        z.faults = [x]
