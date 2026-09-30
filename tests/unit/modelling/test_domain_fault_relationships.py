import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel


def _domain_fault_data():
    # a vertical plane at X = 0.5, given as value constraints on both sides
    xyz = np.array(
        np.meshgrid(np.linspace(0, 1, 5), np.linspace(0, 1, 5), np.linspace(0, 1, 5))
    ).T.reshape(-1, 3)
    data = pd.DataFrame(xyz, columns=["X", "Y", "Z"])
    data["val"] = data["X"] - 0.5
    data["feature_name"] = "domain_fault"
    return data


def _region_names(feature):
    return [getattr(region, "name", None) for region in feature.regions]


def test_domain_fault_does_not_crop_the_unconformity_above(horizontal_data):
    """cover / unconformity / middle / domain fault / basin, built youngest
    first. The domain fault separates middle from basin. The unconformity at
    the base of cover is the boundary above middle, so the domain fault must
    stop there and not crop it, else the cover base is removed on one side
    of the fault.
    """
    model = GeologicalModel([0, 0, 0], [1, 1, 1])
    cover = model.create_and_add_foliation("cover", data=horizontal_data)
    cover_base = model.add_unconformity(cover, 0)
    middle = model.create_and_add_foliation("middle", data=horizontal_data)
    model.data = model.prepare_data(_domain_fault_data())
    model.create_and_add_domain_fault("domain_fault", nelements=500)
    basin = model.create_and_add_foliation("basin", data=horizontal_data)

    domain_fault_region = "__domain_fault_unconformity"
    middle_signs = [r.sign for r in middle.regions if r.name == domain_fault_region]
    basin_signs = [r.sign for r in basin.regions if r.name == domain_fault_region]
    assert middle_signs == [False]
    assert basin_signs == [True]
    assert domain_fault_region not in _region_names(cover_base)
    assert domain_fault_region not in _region_names(cover)
