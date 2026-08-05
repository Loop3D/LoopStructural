### This file tests the function parse_structure_map() in map2loop/mapdata.py
### at the moment only tests for DIPDIR values lower than 360 degrees
### TODO: add more tests for this function

import pytest
import geopandas
import shapely
from map2loop.mapdata import MapData
from map2loop.m2l_enums import Datatype


def test_if_m2l_returns_all_sampled_structures_with_DIPDIR_lower_than_360():

    # call the class
    md = MapData()

    # add config definition
    md.config.structure_config = {
        "dipdir_column": "DIPDIR",
        "dip_column": "DIP",
        "description_column": "DESCRIPTION",
        "bedding_text": "Bedding",
        "objectid_column": "ID",
        "overturned_column": "facing",
        "overturned_text": "DOWN",
        "orientation_type": "strike",
    }

    # create mock data
    data = {
        'geometry': [shapely.geometry.Point(1, 1), shapely.geometry.Point(2, 2), shapely.geometry.Point(3, 3)],
        'DIPDIR': [45.0, 370.0, 420.0],
        'DIP': [30.0, 60.0, 50],
        'OVERTURNED': ["False", "True", "True"],
        'DESCRIPTION': ["Bedding", "Bedding", "Bedidng"],
        'ID': [1, 2, 3],
    }

    # build geodataframe to hold the data
    data = geopandas.GeoDataFrame(data)

    # set it as the raw_data
    md.raw_data[Datatype.STRUCTURE] = data

    # make it parse the structure map and raise exception if error in parse_structure_map

    try:
        md.parse_structure_map()
    except Exception as e:
        pytest.fail(f"parse_structure_map raised an exception: {e}")

    # check if all values below 360
    assert (
        md.data[Datatype.STRUCTURE]['DIPDIR'].all() < 360
    ), "MapData.STRUCTURE is producing DIPDIRs > 360 degrees"


#
