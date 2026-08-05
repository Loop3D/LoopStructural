import pytest
import geopandas as gpd
import shapely.geometry
from map2loop.mapdata import MapData
from map2loop.m2l_enums import Datatype
from map2loop.data_checks import check_fault_fields_validity


@pytest.mark.parametrize(
    "fault_data, fault_config, expected_validity, expected_message",
    [
        # Valid data
        (
            {
                "geometry": [
                    shapely.geometry.LineString([(0, 0), (1, 1)]),
                    shapely.geometry.MultiLineString([[(0, 0), (1, 1)], [(1, 1), (2, 2)]]),
                ],
                "FEATURE": ["Fault A", "Fault B"],
                "ID": [1, 2],
            },
            {"structtype_column": "FEATURE", "fault_text": "Fault", "objectid_column": "ID"},
            False,
            "",
        ),
        # Invalid geometry
        (
            {
                "geometry": [
                    shapely.geometry.LineString([(0, 0), (1, 1)]),
                    shapely.geometry.Polygon(
                        [(0, 0), (1, 1), (1, 0), (0, 1), (0, 0)]
                    ),  # Invalid geometry
                ],
                "FEATURE": ["Fault A", "Fault B"],
                "ID": [1, 2],
            },
            {"structtype_column": "FEATURE", "fault_text": "Fault", "objectid_column": "ID"},
            True,
            "Invalid geometry types found in datatype FAULT. All geometries must be LineString, MultiLineString.",
        ),
        # Non-string FEATURE column
        (
            {
                "geometry": [
                    shapely.geometry.LineString([(0, 0), (1, 1)]),
                    shapely.geometry.MultiLineString([[(0, 0), (1, 1)], [(1, 1), (2, 2)]]),
                ],
                "FEATURE": [5, 2],
                "ID": [1, 2],
            },
            {"structtype_column": "FEATURE", "fault_text": "Fault", "objectid_column": "ID"},
            True,
            "Datatype FAULT: Column 'FEATURE' (config key: 'structtype_column') contains non-string values. Please ensure all values in this column are strings.",
        ),
        # Invalid values in DIP estimate column
        (
            {
                "geometry": [
                    shapely.geometry.LineString([(0, 0), (1, 1)]),
                    shapely.geometry.MultiLineString([[(0, 0), (1, 1)], [(1, 1), (2, 2)]]),
                ],
                "FEATURE": ["Fault", "Fault"],
                "NAME": ["Zuleika", "Zuleika"],
                "ID": [1, 2],
                "DIP": [70, 50],
                "STRIKE": [150, None],
                "DEC": ["north_east", "southt"],
            },
            {
                "structtype_column": "FEATURE",
                "fault_text": "Fault",
                "objectid_column": "ID",
                "name_column": "NAME",
                "dip_column": "DIP",
                "dipdir_column": "STRIKE",
                "dip_estimate_column": "DEC",
            },
            True,
            "Datatype FAULT: Column 'DEC' contains invalid values. Allowed values: ['north_east', 'south_east', 'south_west', 'north_west', 'north', 'east', 'south', 'west'].",
        ),
    ],
    ids=[
        "Valid fault data",
        "Invalid geometry",
        "Non-string FEATURE column",
        "Invalid DIP estimate column",
    ],
)
def test_check_fault_fields_validity(fault_data, fault_config, expected_validity, expected_message):
    # Dynamically create the mock config for this test case
    class MockConfig:
        def __init__(self, config):
            self.fault_config = config

    # Create a GeoDataFrame
    fault_gdf = gpd.GeoDataFrame(fault_data, crs="EPSG:4326")

    # Instantiate the MapData class with the dynamic mock config and data
    map_data = MapData()
    map_data.config = MockConfig(fault_config)
    map_data.raw_data = [None] * len(Datatype.__dict__)
    map_data.raw_data[Datatype.FAULT] = fault_gdf

    # Test the check_fault_fields_validity function
    validity_check, message = check_fault_fields_validity(map_data)
    assert validity_check == expected_validity
    assert message == expected_message
