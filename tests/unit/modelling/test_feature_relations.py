"""GeologicalModel.get_feature_relations: the relations of a model as plain data."""

import numpy as np
import pandas as pd
import pytest

from LoopStructural import GeologicalModel
from tests.fixtures.fault_models import (
    COLUMNS,
    OLD,
    YOUNG,
    add_fault,
    plane_points,
    young_forward,
    young_side,
)


def _series(model, name, faults, value_of_z, seed):
    pts = np.random.default_rng(seed).uniform(0.0, 1.0, size=(400, 3))
    df = pd.DataFrame(
        {"X": pts[:, 0], "Y": pts[:, 1], "Z": pts[:, 2], "val": value_of_z(pts[:, 2]),
         "nx": np.nan, "ny": np.nan, "nz": np.nan, "coord": np.nan, "feature_name": name}
    )[COLUMNS]
    return model.create_and_add_foliation(name, data=df, faults=faults, nelements=3000)


@pytest.fixture(scope="module")
def model():
    """young displaces old; splay abuts young on its footwall; cover above an unconformity."""
    model = GeologicalModel(np.zeros(3), np.ones(3))
    young = add_fault(model, "young", YOUNG["center"], YOUNG["dip"], YOUNG["dip_azimuth"],
                      YOUNG["displacement"], nelements=2000)
    plane, _, _ = plane_points(OLD["center"], OLD["dip"], OLD["dip_azimuth"], n=11)
    keep = np.abs(young_side(plane)) > 0.04
    center = young_forward(np.asarray([OLD["center"]]))[0]
    old = add_fault(model, "old", center, OLD["dip"], OLD["dip_azimuth"], OLD["displacement"],
                    points=young_forward(plane)[keep], faults=[young], nelements=2000)
    splay = add_fault(model, "splay", (0.85, 0.5, 0.3), 70, 90, 0.05, nelements=2000)
    splay.add_abutting_fault(young, positive=False)
    cover = _series(model, "cover", [young], lambda z: z - 0.6, seed=1)
    model.add_unconformity(cover, 0.0)
    _series(model, "basement", [old, young, splay], lambda z: z, seed=2)
    return model


def test_relations_before_and_after_update_are_equal(model):
    before = model.get_feature_relations()
    model.update(progressbar=False)
    assert model.get_feature_relations() == before


def test_fault_relations(model):
    rel = model.get_feature_relations()
    assert rel["faults"]["young"]["displaced_by"] == []
    assert rel["faults"]["old"]["displaced_by"] == ["young"]
    assert rel["faults"]["splay"]["abuts"] == [{"fault": "young", "side": -1}]
    assert rel["faults"]["old"]["abuts"] == []
    assert rel["faults"]["young"]["displacement"] == pytest.approx(YOUNG["displacement"])


def test_feature_relations_and_unconformities(model):
    rel = model.get_feature_relations()
    # Restore order: youngest first.
    assert rel["features"]["basement"]["faults"] == ["splay", "young", "old"]
    assert rel["features"]["cover"]["faults"] == ["young"]
    (uc_name,) = rel["unconformities"]
    uc = rel["unconformities"][uc_name]
    assert uc["feature"] == "cover" and uc["value"] == 0.0 and not uc["onlap"]
    assert "basement" in uc["truncates"]
    basement_cut = rel["features"]["basement"]["truncated_by"]
    assert basement_cut == [
        {"unconformity": uc_name, "feature": "cover", "value": 0.0, "keep": "below"}
    ]
    # The series of the unconformity keeps the part above it.
    assert rel["features"]["cover"]["truncated_by"][0]["keep"] == "above"


def test_event_order_is_the_order_of_creation(model):
    rel = model.get_feature_relations()
    names = [e["name"] for e in rel["events"]]
    assert names[:3] == ["young", "old", "splay"]
    assert names.index("cover") < names.index("basement")
    types = {e["name"]: e["type"] for e in rel["events"]}
    assert types["young"] == "FAULT" and types["basement"] == "INTERPOLATED"
