import pytest
from map2loop.mapdata import MapData


@pytest.mark.parametrize(
    "projection, expected_projection, bounding_box, expected_warning",
    [
        (4326, "EPSG:4326", None, None),  # happy path with int projection
        ("EPSG:3857", "EPSG:3857", None, None),  # happy path with str projection
        (9999, "EPSG:9999", None, None),  # edge case with high int projection
        ("EPSG:9999", "EPSG:9999", None, None),  # edge case with high str projection
        (None, None, None, True),  # error case with None
        ([], None, None, True),  # error case with list
        ({}, None, None, True),  # error case with dict
    ],
    ids=[
        "int_projection",
        "str_projection",
        "high_int_projection",
        "high_str_projection",
        "none_projection",
        "list_projection",
        "dict_projection",
    ],
)
def test_set_working_projection(projection, expected_projection, bounding_box, expected_warning):
    # Set up MapData object
    map_data = MapData()
    map_data.bounding_box = bounding_box

    # Call the method
    map_data.set_working_projection(projection)

    # Assert the working projection is correctly set
    assert map_data.working_projection == expected_projection, (
        f"Expected working_projection to be {expected_projection}, but got {map_data.working_projection}"
    )

    # Check for the presence of warnings via side effects (if applicable)
    if expected_warning:
        assert map_data.working_projection is None, (
            "Expected working_projection to remain None when an invalid projection is provided"
        )
