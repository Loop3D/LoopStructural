"""One up_to_date check looks at each builder once.

A feature depends on the builders of its faults, and a fault frame on the builders
of the faults that cut it. The check walked every path through that graph, so in a
network of cutting faults it grew with the number of paths and ran on every
evaluation (36 million calls to build one fault of a 94-fault model).
"""

from collections import Counter

import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel
from LoopStructural.interpolators import DiscreteInterpolator
from tests.fixtures.fault_models import COLUMNS, OLD, YOUNG, add_fault


def _model():
    model = GeologicalModel(np.zeros(3), np.ones(3))
    young = add_fault(model, "young", **YOUNG, steps=3)
    old = add_fault(model, "old", **OLD, steps=3, faults=[young])
    add_fault(model, "oldest", center=(0.5, 0.5, 0.7), dip=70, dip_azimuth=0,
              displacement=0.08, steps=3, faults=[old, young])
    xyz = np.random.default_rng(1).uniform(0, 1, (100, 3))
    data = pd.DataFrame({"X": xyz[:, 0], "Y": xyz[:, 1], "Z": xyz[:, 2], "val": xyz[:, 2],
                         "nx": np.nan, "ny": np.nan, "nz": np.nan, "coord": np.nan,
                         "feature_name": "strat"})[COLUMNS]
    strat = model.create_and_add_foliation("strat", data=data, nelements=2000, buffer=0.3)
    model.update()
    return strat


def test_each_interpolator_is_checked_once(monkeypatch):
    strat = _model()
    reads = Counter()
    flag = DiscreteInterpolator.up_to_date

    def counted(self):
        reads[id(self)] += 1
        return flag.fget(self)

    monkeypatch.setattr(DiscreteInterpolator, "up_to_date", property(counted, flag.fset))
    strat.builder.up_to_date()
    # strat and the 3 coordinates of each of the 3 fault frames
    assert len(reads) == 10
    assert max(reads.values()) == 1


def test_a_stale_fault_is_still_rebuilt():
    strat = _model()
    young = strat.faults[-1] if strat.faults[-1].name == "young" else strat.faults[0]
    young.builder.builders[0]._up_to_date = False
    strat.builder.up_to_date()
    assert young.builder.builders[0]._up_to_date
