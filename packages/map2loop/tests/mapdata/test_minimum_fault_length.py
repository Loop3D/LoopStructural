import pytest
import geopandas as gpd
import shapely
import numpy

from map2loop.mapdata import MapData
from map2loop.m2l_enums import VerboseLevel, Datatype


@pytest.fixture
def setup_map_data():
    # Create a mock MapData object with verbose level set to ALL
    map_data = MapData(verbose_level=VerboseLevel.ALL)

    # Simulate config with no minimum_fault_length set
    map_data.config.fault_config['minimum_fault_length'] = -1.0

    return map_data


# for cases when there is no minimum_fault_length set in the config by user - does it update from map?
def test_update_minimum_fault_length_from_faults(setup_map_data):
    map_data = setup_map_data

    # Define a test bounding box (in meters)
    bounding_box = {
        "minx": 0,
        "miny": 0,
        "maxx": 10000,  # width = 10,000 meters
        "maxy": 10000,  # height = 10,000 meters
    }

    # Set the bounding box in the map data
    map_data.set_bounding_box(bounding_box)

    # update config
    map_data.config.fault_config['name_column'] = 'NAME'
    map_data.config.fault_config['dip_column'] = 'DIP'

    # Define a dummy fault GeoDataFrame with faults of varying lengths
    faults = gpd.GeoDataFrame(
        {
            'geometry': [
                shapely.geometry.LineString(
                    [(0, 0), (50, 50)]
                ),  # Fault 1 (small, length ~70.7 meters)
                shapely.geometry.LineString(
                    [(0, 0), (3000, 3000)]
                ),  # Fault 2 (length ~4242 meters)
                shapely.geometry.LineString(
                    [(0, 0), (7000, 7000)]
                ),  # Fault 3 (length ~9899 meters)
            ],
            'NAME': ['Fault_1', 'Fault_2', 'Fault_3'],
            'DIP': [60, 45, 30],
        }
    )

    faults.crs = "EPSG: 7850"

    # get the cropped fault dataset from parse_fault_map
    map_data.raw_data[Datatype.FAULT] = faults
    map_data.parse_fault_map()
    cropped_faults = map_data.data[Datatype.FAULT]

    # calculate 5% length of the bounding box area
    expected_minimum_fault_length = numpy.sqrt(
        0.05
        * (bounding_box['maxx'] - bounding_box['minx'])
        * (bounding_box['maxy'] - bounding_box['miny'])
    )

    # Verify that the minimum_fault_length was calculated correctly
    assert map_data.minimum_fault_length == pytest.approx(expected_minimum_fault_length, rel=1e-3)

    # There should only be two faults remaining (the second and third ones)
    assert len(cropped_faults) == 2

    # Ensure that the remaining faults are the correct ones
    remaining_lengths = cropped_faults.geometry.length
    assert all(remaining_lengths >= expected_minimum_fault_length)
    assert cropped_faults.geometry.equals(
        faults.iloc[1:]['geometry']
    )  # Faults 2 and 3 geometries should be the same in the faults raw and faults cropped


# are faults with length less than minimum_fault_length removed from the dataset?
def test_cropping_faults_by_minimum_fault_length(setup_map_data):
    map_data = setup_map_data

    # Set minimum_fault_length in the config to 10
    map_data.config.fault_config['minimum_fault_length'] = 10.0

    map_data.config.fault_config['name_column'] = 'NAME'
    map_data.config.fault_config['dip_column'] = 'DIP'

    # Create a mock faults dataset with lengths < 10 and > 10
    faults = gpd.GeoDataFrame(
        {
            'geometry': [
                shapely.geometry.LineString([(0, 0), (2, 2)]),  # Length ~2.83 (should be cropped)
                shapely.geometry.LineString([(0, 0), (5, 5)]),  # Length ~7.07 (should be cropped)
                shapely.geometry.LineString([(0, 0), (10, 10)]),  # Length ~14.14 (should remain)
            ],
            'NAME': ['Fault_1', 'Fault_2', 'Fault_3'],
            'DIP': [60, 45, 30],
            'DIPDIR': [90, 120, 150],
        }
    )

    # Set the raw data in the map_data object to simulate loaded fault data
    map_data.raw_data[Datatype.FAULT] = faults
    map_data.parse_fault_map()
    cropped_faults = map_data.data[Datatype.FAULT]

    # Assert that only 1 fault remains (the one with length ~14.14)
    assert (
        len(cropped_faults) == 1
    ), f"Expected only 1 fault remaining after cropping, but found {len(cropped_faults)}"

    # Optionally, check that the remaining fault has the correct geometry (the long one)
    assert cropped_faults.iloc[0].geometry.length == pytest.approx(
        14.14, 0.01
    ), f"Expected remaining fault length to be 14.14, but got {cropped_faults.iloc[0].geometry.length}"
