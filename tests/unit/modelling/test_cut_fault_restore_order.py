"""A fault faulted by a younger fault (``faults=[A]``) is restored after ``A``.

Features keep their faults newest first. A fault created with ``faults=[A]``
needs ``A`` to exist, so ``A`` comes later in that list, and was restored after
the fault it cuts: the frame then restored every point through ``A`` inside each
step (a cost that multiplies along chains of cutting faults) and gave a
different result near the intersection.
"""

import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel
from LoopStructural.modelling.features.fault._fault_segment import FaultSegment
from LoopStructural.utils._fault_block import cutters_first, restore_through
from tests.fixtures.fault_models import COLUMNS, OLD, YOUNG, add_fault, interior_points


class _Fault:
    """A stand-in with the frame builders that hold the cutting faults."""

    def __init__(self, name, cut_by=()):
        self.name = name
        self.builder = type("B", (), {"builders": [type("S", (), {"faults": list(cut_by)})()]})()

    def __getitem__(self, key):
        raise IndexError(key)


def test_cutting_fault_comes_first_and_other_faults_keep_their_order():
    a = _Fault("a")
    i = _Fault("i")
    b = _Fault("b", cut_by=[a])
    c = _Fault("c", cut_by=[b])
    assert [f.name for f in cutters_first([c, i, b, a])] == ["i", "a", "b", "c"]
    assert [f.name for f in cutters_first([i, a, b])] == ["i", "a", "b"]
    unrelated = [_Fault("x"), _Fault("y"), _Fault("z")]
    assert cutters_first(unrelated) == unrelated


def _model_with_cut_pair():
    model = GeologicalModel(np.zeros(3), np.ones(3))
    young = add_fault(model, "young", **YOUNG, steps=3)
    old = add_fault(model, "old", **OLD, steps=3, faults=[young])
    rng = np.random.default_rng(1)
    pts = rng.uniform(0, 1, (200, 3))
    data = pd.DataFrame({"X": pts[:, 0], "Y": pts[:, 1], "Z": pts[:, 2], "val": pts[:, 2],
                         "nx": np.nan, "ny": np.nan, "nz": np.nan, "coord": np.nan,
                         "feature_name": "strat"})[COLUMNS]
    strat = model.create_and_add_foliation("strat", data=data, nelements=2000, buffer=0.3)
    model.update()
    return strat, young, old


def test_feature_restores_the_cutting_fault_first(monkeypatch):
    strat, young, old = _model_with_cut_pair()
    # The feature keeps the faults newest first: the cut fault comes before its cutter.
    assert [f.name for f in strat.faults] == ["old", "young"]
    points = interior_points(n=300)
    expected = old.apply_to_points(young.apply_to_points(points.copy()), restored_through=(young,))

    calls = []
    original = FaultSegment.apply_to_points

    def counted(self, pts, *args, **kwargs):
        calls.append(self.name)
        return original(self, pts, *args, **kwargs)

    monkeypatch.setattr(FaultSegment, "apply_to_points", counted)
    restored = restore_through(points.copy(), strat.faults)
    np.testing.assert_allclose(restored, expected, rtol=0.0, atol=1e-12)
    # One application of each fault: the frame of old does not restore through young again.
    assert calls == ["young", "old"]


def test_forward_movement_keeps_the_list_order_without_relations():
    a, b = _Fault("a"), _Fault("b")
    seen = []
    for f in (a, b):
        f.apply_to_points = lambda pts, f=f, **kw: seen.append(f.name) or pts
    restore_through(np.zeros((1, 3)), [a, b], reverse=True)
    assert seen == ["a", "b"]


def test_a_chain_of_cuts_is_restored_in_one_pass(monkeypatch):
    """oldest is cut by old, old by young: no fault restores through its cutters alone.

    Restoring the frame of oldest through old alone made the frame of old restore
    every point through young again at each step (one more level of nesting for
    each link of the chain, so the cost multiplied). The chain is now restored in
    one pass, youngest first, so nesting never goes deeper than 2.
    """
    model = GeologicalModel(np.zeros(3), np.ones(3))
    young = add_fault(model, "young", **YOUNG, steps=3)
    old = add_fault(model, "old", **OLD, steps=3, faults=[young])
    oldest = add_fault(model, "oldest", center=(0.5, 0.5, 0.7), dip=70, dip_azimuth=0,
                       displacement=0.08, steps=3, faults=[old])
    model.update()

    depth, deepest = [0], [0]
    original = FaultSegment.apply_to_points

    def tracked(self, pts, *args, **kwargs):
        depth[0] += 1
        deepest[0] = max(deepest[0], depth[0])
        try:
            return original(self, pts, *args, **kwargs)
        finally:
            depth[0] -= 1

    monkeypatch.setattr(FaultSegment, "apply_to_points", tracked)
    restored = oldest.apply_to_points(interior_points(n=200))
    assert not np.isnan(restored).any()
    assert deepest[0] == 2
