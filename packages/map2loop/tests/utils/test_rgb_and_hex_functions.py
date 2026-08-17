### This file tests the function generate_random_hex_colors() and hex_to_rgba() in map2loop/utils.py

import pytest
import re
from map2loop.utils import generate_random_hex_colors, hex_to_rgb

# does it return the right number of colors?
def test_generate_random_hex_colors_length():
    n = 5
    colors = generate_random_hex_colors(n)
    assert (
        len(colors) == n
    ), f"utils function generate_random_hex_colors not returning the right number of hex codes.Expected {n} colors, got {len(colors)}"


# are the returned hex strings the right format?
def test_generate_random_hex_colors_format():
    n = 10
    hex_pattern = re.compile(r'^#[0-9A-Fa-f]{6}$')
    colors = generate_random_hex_colors(n)
    for color in colors:
        assert hex_pattern.match(
            color
        ), f"utils function generate_random_hex_colors not returning hex strings in the right format. Got {color} instead."


# is hex conversion to rgba working as expected?
def test_hex_to_rgba_long_hex():
    hex_color = "#1a2b3c"  # long hex versions
    expected_output = (0.10196078431372549, 0.16862745098039217, 0.23529411764705882, 1.0)
    assert (
        hex_to_rgb(hex_color) == expected_output
    ), f"utils function hex_to_rgba not doing hex to rgba conversion correctly. Expected {expected_output}, got {hex_to_rgb(hex_color)}"


def test_hex_to_rgba_short_hex():
    hex_color = "#abc"  # short hex versions
    expected_output = (0.6666666666666666, 0.7333333333333333, 0.8, 1.0)
    assert hex_to_rgb(hex_color) == expected_output


# does it handle invalid inputs correctly?
def test_hex_to_rgba_invalid_hex():
    with pytest.raises(ValueError):
        hex_to_rgb(
            "12FF456"
        ), "utils function hex_to_rgba is expected to raise a ValueError when an invalid hex string is passed, but it did not."
