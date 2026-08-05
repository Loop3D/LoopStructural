import pytest
import geopandas as gpd
import shapely.geometry
from map2loop.mapdata import MapData
from map2loop.data_checks import check_geology_fields_validity

# Datatype Enum
class Datatype:
    GEOLOGY = 0

# Config 
class MockConfig:
    def __init__(self):
        self.geology_config = {
            "unitname_column": "UNITNAME",
            "alt_unitname_column": "CODE",
            "group_column": "GROUP",
            "supergroup_column": "SUPERGROUP",
            "description_column": "DESCRIPTION",
            "rocktype_column": "ROCKTYPE1",
            "alt_rocktype_column": "ROCKTYPE2",
            "minage_column": "MIN_AGE",
            "maxage_column": "MAX_AGE",
            "objectid_column": "ID",
            "ignore_lithology_codes": [],
        }

@pytest.mark.parametrize(
    "geology_data, expected_validity, expected_message",
    [
        # Valid data
        (
            {
                "geometry": [shapely.geometry.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
                "UNITNAME": ["Sandstone"],
                "CODE": ["SST"],
                "GROUP": ["Sedimentary"],
                "SUPERGROUP": ["Mesozoic"],
                "DESCRIPTION": ["A type of sandstone"],
                "ROCKTYPE1": ["Clastic"],
                "ROCKTYPE2": ["Quartz"],
                "MIN_AGE": [150.0],
                "MAX_AGE": [200.0],
                "ID": [1],
            },
            False,
            "",
        ),
        # Invalid geometry
        (
            {
                "geometry": [shapely.geometry.Polygon([(0, 0), (1, 1), (1, 0), (0, 1), (0, 0)])],
                "UNITNAME": ["Sandstone"],
                "CODE": ["SST"],
                "GROUP": ["Sedimentary"],
                "SUPERGROUP": ["Mesozoic"],
                "DESCRIPTION": ["A type of sandstone"],
                "ROCKTYPE1": ["Clastic"],
                "ROCKTYPE2": ["Quartz"],
                "MIN_AGE": [150.0],
                "MAX_AGE": [200.0],
                "ID": [1],
            },
            False,
            "",
        ),
        # Missing required column
        (
            {
                "geometry": [shapely.geometry.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
                "UNITNAME": ["Sandstone"],
                # "CODE": ["SST"],  # Missing required column
                "GROUP": ["Sedimentary"],
                "SUPERGROUP": ["Mesozoic"],
                "DESCRIPTION": ["A type of sandstone"],
                "ROCKTYPE1": ["Clastic"],
                "ROCKTYPE2": ["Quartz"],
                "MIN_AGE": [150.0],
                "MAX_AGE": [200.0],
                "ID": [1],
            },
            True,
            "Datatype GEOLOGY: Required column with config key 'alt_unitname_column' (column: 'CODE')  is missing from the data.",
        ),
        # Non-string value in required column
        (
            {
                "geometry": [shapely.geometry.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
                "UNITNAME": ["Sandstone"],
                "CODE": [2],  # Non-string value
                "GROUP": ["Sedimentary"],
                "SUPERGROUP": ["Mesozoic"],
                "DESCRIPTION": ["A type of sandstone"],
                "ROCKTYPE1": ["Clastic"],
                "ROCKTYPE2": ["Quartz"],
                "MIN_AGE": [150.0],
                "MAX_AGE": [200.0],
                "ID": [1],
            },
            True,
            "Datatype GEOLOGY: Column 'alt_unitname_column' (column: 'CODE') must contain only <class 'str'> values.",
        ),
        # NaN or blank value in required column
        (
            {
                "geometry": [shapely.geometry.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
                "UNITNAME": [""],  # Blank value
                "CODE": ["SST"],
                "GROUP": ["Sedimentary"],
                "SUPERGROUP": ["Mesozoic"],
                "DESCRIPTION": ["A type of sandstone"],
                "ROCKTYPE1": ["Clastic"],
                "ROCKTYPE2": ["Quartz"],
                "MIN_AGE": [150.0],
                "MAX_AGE": [200.0],
                "ID": [1],
            },
            True,
            "Datatype GEOLOGY: Column 'unitname_column' (column: 'UNITNAME') contains blank (empty) values. Please ensure all values are populated.",
        ),
        # Duplicate ID values
        (
            {
                "geometry": [
                    shapely.geometry.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
                    shapely.geometry.Polygon([(0, 0), (10, 0), (1, 1), (0, 10)]),
                ],
                "UNITNAME": ["fr", "df"],
                "CODE": ["SST", "FGH"],
                "GROUP": ["Sedimentary", "Ign"],
                "SUPERGROUP": ["Mesozoic", "Arc"],
                "DESCRIPTION": ["A", "B"],
                "ROCKTYPE1": ["A", "B"],
                "ROCKTYPE2": ["Quartz", "FDS"],
                "MIN_AGE": [150.0, 200],
                "MAX_AGE": [200.0, 250],
                "ID": [1, 1],  # Duplicate ID
            },
            True,
            "Datatype GEOLOGY: Column 'ID' (config key: 'objectid_column') contains duplicate values.",
        ),
        # nan in id
        (
            {
                "geometry": [
                    shapely.geometry.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
                    shapely.geometry.Polygon([(0, 0), (10, 0), (1, 1), (0, 10)]),
                ],
                "UNITNAME": ["fr", "df"],
                "CODE": ["SST", "FGH"],
                "GROUP": ["Sedimentary", "Ign"],
                "SUPERGROUP": ["Mesozoic", "Arc"],
                "DESCRIPTION": ["A", "B"],
                "ROCKTYPE1": ["A", "B"],
                "ROCKTYPE2": ["Quartz", "FDS"],
                "MIN_AGE": [150.0, 200],
                "MAX_AGE": [200.0, 250],
                "ID": [1, None],  
            },
            True,
            "Datatype GEOLOGY: Column 'ID' (config key: 'objectid_column') contains non-numeric or NaN values. Please rectify the values, or remove this key from the config dictionary to let map2loop assign IDs.",
        ),
        # nan in unit name
        (
            {
                "geometry": [
                    shapely.geometry.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
                    shapely.geometry.Polygon([(0, 0), (10, 0), (1, 1), (0, 10)]),
                ],
                "UNITNAME": ["fr", None],
                "CODE": ["SST", "FGH"],
                "GROUP": ["Sedimentary", "Ign"],
                "SUPERGROUP": ["Mesozoic", "Arc"],
                "DESCRIPTION": ["A", "B"],
                "ROCKTYPE1": ["A", "B"],
                "ROCKTYPE2": ["Quartz", "FDS"],
                "MIN_AGE": [150.0, 200],
                "MAX_AGE": [200.0, 250],
                "ID": [1, 1],  # Duplicate ID
            },
            True,
            "Datatype GEOLOGY: Column 'unitname_column' (column: 'UNITNAME') must contain only <class 'str'> values.",
        ),
    ],
)



def test_check_geology_fields_validity(geology_data, expected_validity, expected_message):
    # Create a GeoDataFrame
    geology_gdf = gpd.GeoDataFrame(geology_data, crs="EPSG:4326")

    # Instantiate the MapData class with the mock config and data
    map_data = MapData()
    map_data.config = MockConfig()
    map_data.raw_data = [None] * len(Datatype.__dict__)
    map_data.raw_data[Datatype.GEOLOGY] = geology_gdf

    # Test the check_geology_fields_validity function
    validity_check, message = check_geology_fields_validity(map_data)

    assert validity_check == expected_validity
    assert message == expected_message