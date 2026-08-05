### This file tests the function colour_units() in map2loop/mapdata.py
### Two main test cases are considered: cases in which there is clut file and cases in which there is no clut file

import pandas as pd
from map2loop.mapdata import MapData


# are random colours being assigned to stratigraphic units in cases where no clut file is provided?
def test_colour_units_no_clut_file():
    # Create a sample DataFrame with missing 'colour' values
    data = {"name": ["Unit1", "Unit2", "Unit3"], "colour": [None, None, None]}
    stratigraphic_units = pd.DataFrame(data)

    # Instantiate the class and call the method
    md = MapData()
    md.colour_filename = None  # Ensure no file is used
    result = md.colour_units(stratigraphic_units)

    # check that there are no duplicates in the 'unit' column
    assert result['name'].is_unique, "colour_units() in mapdata.py producing duplicate units"

    # Check that the 'colour' column has been assigned random colors
    assert (
        len(result["colour"].dropna()) == 3
    ), "function MapData.colour_units not assigning the right len of random colours"

    # are the generated colours valid hex colours?
    colours = result["colour"]

    assert colours.apply(
        isinstance, args=(str,)
    ).all(), (
        "function MapData.colour_units without clut file not assigning random hex colours as str"
    )
    assert colours.str.startswith(
        "#"
    ).all(), "function MapData.colour_units without clut file not generating the right hex codes with # at the start"
    assert (
        colours.str.len().isin([7, 4]).all()
    ), "function MapData.colour_units without clut file not generating the right hex codes with 7 or 4 characters"


def test_colour_units_with_colour_file(tmp_path):
    # Create a strati units df with missing 'colour' values
    data = {"name": ["Unit1", "Unit2", "Unit3"], "colour": [None, None, None]}
    stratigraphic_units = pd.DataFrame(data)

    # Create a temp clut file
    colour_data = {"UNITNAME": ["Unit1", "Unit2"], "colour": ["#112233", "#445566"]}
    colour_lookup_df = pd.DataFrame(colour_data)
    colour_filename = tmp_path / "colour_lookup.csv"
    colour_lookup_df.to_csv(colour_filename, index=False)

    # Instantiate the class and call the method
    md = MapData()
    md.colour_filename = str(colour_filename)
    result = md.colour_units(stratigraphic_units)

    # check that there are no duplicates in the 'unit' column
    assert result['name'].is_unique, "colour_units() in mapdata.py producing duplicate units"

    # Check that the 'colour' column has been merged correctly and missing colors are assigned
    expected_colors = ["#112233", "#445566"]
    assert (
        result["colour"].iloc[0] == expected_colors[0]
    ), "function MapData.colour_units with clut file not assigning the right colour from the lookup file"
    assert (
        result["colour"].iloc[1] == expected_colors[1]
    ), "function MapData.colour_units with clut file not assigning the right colour from the lookup file"
    assert isinstance(
        result["colour"].iloc[2], str
    ), "function MapData.colour_units with clut file not assigning random hex colours as str"
    assert (
        result["colour"].iloc[2].startswith("#")
    ), "function MapData.colour_units with clut file not generating the right hex codes with # at the start"
    assert len(result["colour"].iloc[2]) in {
        7,
        4,
    }, "function MapData.colour_units with clut file not generating the right hex codes with 7 or 4 characters"
