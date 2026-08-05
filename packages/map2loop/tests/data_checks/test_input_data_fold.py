import pytest
import geopandas as gpd
import shapely.geometry
from map2loop.mapdata import MapData
from map2loop.m2l_enums import Datatype
from map2loop.data_checks import check_fold_fields_validity

@pytest.mark.parametrize(
    "fold_data, fold_config, expected_validity, expected_message",
    [
        # Valid data
        (
            {
                "geometry": [
                    shapely.geometry.LineString([(0, 0), (1, 1)]),
                    shapely.geometry.MultiLineString([[(0, 0), (1, 1)], [(1, 1), (2, 2)]])
                ],
                "FEATURE": ["fold A", "fold B"],
                "ID": [1, 2],
                "description": ["desc1", "desc2"]
            },
            {"structtype_column": "FEATURE", "fold_text": "fold", "objectid_column": "ID", "description_column": "description"},
            False,
            ""
        ),
        # Missing geometry
        (
            {
                "geometry": [
                    shapely.geometry.LineString([(0,0), (0,0)]),  # Invalid type
                    shapely.geometry.LineString([(0, 0), (1, 1)])
                ],
                "FEATURE": ["fold A", "fold B"],
                "ID": [1, 2],
                "description": ["desc1", "desc2"]
            },
            {"structtype_column": "FEATURE", "fold_text": "fold", "objectid_column": "ID", "description_column": "description"},
            True,
            "Invalid geometry types found in datatype FOLD. All geometries must be LineString, MultiLineString."
        ),
        # Non-string FEATURE column
        (
            {
                "geometry": [
                    shapely.geometry.LineString([(0, 0), (1, 1)]),
                    shapely.geometry.MultiLineString([[(0, 0), (1, 1)], [(1, 1), (2, 2)]])
                ],
                "FEATURE": [123, 456],  # Invalid type
                "ID": [1, 2],
                "description": ["desc1", "desc2"]
            },
            {"structtype_column": "FEATURE", "fold_text": "fold", "objectid_column": "ID", "description_column": "description"},
            True,
            "Datatype FOLD: Column 'FEATURE' (config key: 'structtype_column') contains non-string values. Please ensure all values in this column are strings."
        ),
        # Missing ID column
        (
            {
                "geometry": [
                    shapely.geometry.LineString([(0, 0), (1, 1)]),
                    shapely.geometry.MultiLineString([[(0, 0), (1, 1)], [(1, 1), (2, 2)]])
                ],
                "FEATURE": ["fold A", "fold B"],
                "description": ["desc1", "desc2"]
            },
            {"structtype_column": "FEATURE", "fold_text": "fold", "objectid_column": "ID", "description_column": "description"},
            False,
            ""
        ),
        # Duplicate ID values
        (
            {
                "geometry": [
                    shapely.geometry.LineString([(0, 0), (1, 1)]),
                    shapely.geometry.MultiLineString([[(0, 0), (1, 1)], [(1, 1), (2, 2)]])
                ],
                "FEATURE": ["fold A", "fold B"],
                "ID": [1, 1],  # Duplicate values
                "description": ["desc1", "desc2"]
            },
            {"structtype_column": "FEATURE", "fold_text": "fold", "objectid_column": "ID", "description_column": "description"},
            True,
            "Datatype FOLD: Column 'ID' (config key: 'objectid_column') contains duplicate values."
        ),
    ],
    ids=[
        "Valid fold data",
        "Invalid geometry",
        "Non-string FEATURE column",
        "Missing ID column",
        "Duplicate ID values"
    ]
)
def test_check_fold_fields_validity(fold_data, fold_config, expected_validity, expected_message):
    # Dynamically create the mock config for this test case
    class MockConfig:
        def __init__(self, config):
            self.fold_config = config

    # Create a GeoDataFrame
    fold_gdf = gpd.GeoDataFrame(fold_data, crs="EPSG:4326")

    # Instantiate the MapData class with the dynamic mock config and data
    map_data = MapData()
    map_data.config = MockConfig(fold_config)
    map_data.raw_data = [None] * len(Datatype.__dict__)
    map_data.raw_data[Datatype.FOLD] = fold_gdf

    # Test the check_fold_fields_validity function
    validity_check, message = check_fold_fields_validity(map_data)
    assert validity_check == expected_validity
    assert message == expected_message
