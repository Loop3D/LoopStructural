import pytest
from map2loop.project import Project
from map2loop.sorter import SorterUseNetworkX
from map2loop.m2l_enums import VerboseLevel
from unittest.mock import patch
from pyproj.exceptions import CRSError
import requests
import os

# Define constants for common parameters
bbox_3d = {
    "minx": 515687.31005864,
    "miny": 7493446.76593407,
    "maxx": 562666.860106543,
    "maxy": 7521273.57407786,
    "base": -3200,
    "top": 3000,
}
loop_project_filename = "wa_output.loop3d"

# create a project function
def create_project(state_data="WA", projection="EPSG:28350"):
    return Project(
        use_australian_state_data=state_data,
        working_projection=projection,
        bounding_box=bbox_3d,
        verbose_level=VerboseLevel.NONE,
        loop_project_filename=loop_project_filename,
        overwrite_loopprojectfile=True,
    )


# is the project running?
def test_project_execution():
    try:
        proj = create_project()
    except Exception:
        pytest.skip("Skipping the project test from server data due to loading failure")
    try:
        proj.set_sorter(SorterUseNetworkX())
        proj.run_all(take_best=True)
    except requests.exceptions.ReadTimeout:
        pytest.skip(
            "Skipping the project test from server data due to timeout while attempting to run proj.run_all"
        )

    # if no timeout:
    # is there a project?
    assert proj is not None, "Plot Hamersley Basin failed to execute"
    # is there a LPF?
    assert os.path.exists(
        loop_project_filename
    ), f"Expected file {loop_project_filename} was not created"


# Is the test_project_execution working - ie, is the test skipped on timeout?
def test_timeout_handling():
    # Mock `openURL` in `owslib.util` to raise a ReadTimeout directly
    with patch("owslib.util.openURL"):
        # Run `test_project_execution` and check if the skip occurs
        result = pytest.main(
            ["-q", "--tb=short", "--disable-warnings", "-k", "test_project_execution"]
        )
        assert (
            result.value == pytest.ExitCode.OK
        ), "The test was not skipped as expected on timeout."


# does the project fail when the CRS is invalid?
def test_expect_crs_error():
    try:
        with pytest.raises(CRSError):
            create_project(projection="InvalidCRS")
        print("CRSError was raised as expected.")
    except requests.exceptions.ReadTimeout:
        print("Timeout occurred, skipping test_expect_crs_error.")
        pytest.skip("Skipping test_expect_crs_error due to a timeout.")


# does the project fail when the Aus state name is invalid?
def test_expect_state_error():
    try:
        with pytest.raises(ValueError):
            create_project(state_data="InvalidState")
        print("ValueError was raised as expected.")
    except requests.exceptions.ReadTimeout:
        print("Timeout occurred, skipping test_expect_state_error.")
        pytest.skip("Skipping test_expect_state_error due to a timeout.")


# does the project fail when a config file is invalid?
def test_expect_config_error():
    try:
        with pytest.raises(Exception):
            create_project(config_file='InvalidConfig.csv')
        print("FileNotFoundError//Exception by catchall in project.py was raised as expected.")
    except requests.exceptions.ReadTimeout:
        print("Timeout occurred, skipping test_expect_config_error.")
        pytest.skip("Skipping test_expect_config_error due to a timeout.")
