"""The builder rotates data vectors through the faults in the order points are restored.

Each vector is rotated by each fault at its position restored so far, youngest fault
first, and each fault is told which faults are already restored. Rotating at the
present-day position made a fault frame displaced by younger faults restore the
vector through its whole chain of cutting faults at each step.
"""

import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel
from LoopStructural.modelling.features.fault._fault_segment import FaultSegment
from tests.fixtures.fault_models import COLUMNS, OLD, YOUNG, add_fault


def _track_depth(monkeypatch):
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
    return deepest


def test_building_a_frame_in_a_chain_of_cuts_does_not_nest(monkeypatch):
    model = GeologicalModel(np.zeros(3), np.ones(3))
    young = add_fault(model, "young", **YOUNG, steps=3)
    old = add_fault(model, "old", **OLD, steps=3, faults=[young])
    oldest = add_fault(model, "oldest", center=(0.5, 0.5, 0.7), dip=70, dip_azimuth=0,
                       displacement=0.08, steps=3, faults=[old])
    young.builder.up_to_date()
    old.builder.up_to_date()
    deepest = _track_depth(monkeypatch)
    oldest.builder.up_to_date()
    # The data of oldest (with its normal vectors) go through young then old, and the
    # frame of old is read at points young already restored: no fault inside another.
    assert deepest[0] == 1


def test_one_fault_rotates_vectors_at_their_present_position(monkeypatch):
    """With one fault the vectors are rotated as before: at the data positions."""
    model = GeologicalModel(np.zeros(3), np.ones(3))
    fault = add_fault(model, "young", **YOUNG, steps=3)
    fault.builder.up_to_date()
    rng = np.random.default_rng(4)
    xyz = rng.uniform(0.1, 0.9, (40, 3))
    data = pd.DataFrame({"X": xyz[:, 0], "Y": xyz[:, 1], "Z": xyz[:, 2], "val": np.nan,
                         "nx": 0.0, "ny": 0.0, "nz": 1.0, "coord": np.nan,
                         "feature_name": "strat"})[COLUMNS]
    strat = model.create_and_add_foliation("strat", data=data, nelements=2000, buffer=0.3,
                                           faults=[fault])
    seen = []
    original = FaultSegment.apply_to_vectors

    def recorded(self, vector, *args, **kwargs):
        seen.append((self.name, np.array(vector[:, :3]), kwargs.get("restored_through", ())))
        return original(self, vector, *args, **kwargs)

    monkeypatch.setattr(FaultSegment, "apply_to_vectors", recorded)
    strat.builder.up_to_date()
    assert [name for name, _, _ in seen] == ["young"]
    _, positions, restored = seen[0]
    np.testing.assert_array_equal(positions, xyz)
    assert tuple(restored) == ()
