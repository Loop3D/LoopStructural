"""Tests for the pandas adapter."""

from __future__ import annotations

import time

import numpy as np
import pytest

pd = pytest.importorskip("pandas")

from loop_link.adapters.pandas import drillholes_from_pandas
from loop_link.errors import InvalidInputError
from loop_link.objects import Drillhole


@pytest.fixture
def collar() -> pd.DataFrame:
    """Return a collar table with two holes."""
    return pd.DataFrame(
        {
            "hole_id": ["DH1", "DH2"],
            "x": [10.0, 20.0],
            "y": [11.0, 21.0],
            "z": [12.0, 22.0],
        }
    )


@pytest.fixture
def survey() -> pd.DataFrame:
    """Return a survey table with two stations for each hole."""
    return pd.DataFrame(
        {
            "hole_id": ["DH1", "DH1", "DH2", "DH2"],
            "depth": [0.0, 50.0, 0.0, 40.0],
            "azimuth": [90.0, 95.0, 180.0, 185.0],
            "dip": [60.0, 62.0, 70.0, 72.0],
        }
    )


@pytest.fixture
def logs() -> pd.DataFrame:
    """Return a log table with a number column and a string column."""
    return pd.DataFrame(
        {
            "hole_id": ["DH1", "DH1", "DH2"],
            "from": [0.0, 10.0, 0.0],
            "to": [10.0, 30.0, 25.0],
            "grade": [1.5, 2.5, 3.5],
            "lith": ["granite", "basalt", "chert"],
        }
    )


def test_keys_and_order(collar, survey, logs):
    result = drillholes_from_pandas(collar, survey, logs)
    assert list(result) == ["DH1", "DH2"]
    assert all(isinstance(hole, Drillhole) for hole in result.values())


def test_collar_row(collar, survey, logs):
    result = drillholes_from_pandas(collar, survey, logs)
    np.testing.assert_allclose(result["DH1"].collar, [10.0, 11.0, 12.0])
    np.testing.assert_allclose(result["DH2"].collar, [20.0, 21.0, 22.0])


def test_survey_shape_and_order(collar, survey, logs):
    hole = drillholes_from_pandas(collar, survey, logs)["DH1"]
    assert hole.survey.shape == (2, 3)
    np.testing.assert_allclose(hole.survey, [[0.0, 90.0, 60.0], [50.0, 95.0, 62.0]])


def test_survey_is_sorted_by_depth(collar, survey, logs):
    shuffled = survey.iloc[[1, 0, 3, 2]]
    hole = drillholes_from_pandas(collar, shuffled, logs)["DH1"]
    np.testing.assert_allclose(hole.survey[:, 0], [0.0, 50.0])
    np.testing.assert_allclose(hole.survey[:, 2], [60.0, 62.0])


def test_intervals(collar, survey, logs):
    result = drillholes_from_pandas(collar, survey, logs)
    assert result["DH1"].intervals.shape == (2, 2)
    np.testing.assert_allclose(result["DH1"].intervals, [[0.0, 10.0], [10.0, 30.0]])
    assert result["DH2"].intervals.shape == (1, 2)
    np.testing.assert_allclose(result["DH2"].intervals, [[0.0, 25.0]])


def test_interval_properties(collar, survey, logs):
    props = drillholes_from_pandas(collar, survey, logs)["DH1"].interval_properties
    assert props["grade"].dtype.kind == "f"
    np.testing.assert_allclose(props["grade"], [1.5, 2.5])
    assert props["lith"].dtype.kind == "U"
    assert list(props["lith"]) == ["granite", "basalt"]


def test_dip_positive_down_false_changes_sign(collar, survey, logs):
    hole = drillholes_from_pandas(collar, survey, logs, dip_positive_down=False)["DH1"]
    np.testing.assert_allclose(hole.survey[:, 2], [-60.0, -62.0])
    np.testing.assert_allclose(hole.survey[:, 1], [90.0, 95.0])


def test_custom_column_name(collar, survey, logs):
    renamed = collar.rename(columns={"x": "east"})
    result = drillholes_from_pandas(renamed, survey, logs, x="east")
    np.testing.assert_allclose(result["DH2"].collar, [20.0, 21.0, 22.0])


def test_survey_hole_not_in_collar(collar, survey, logs):
    extra = pd.concat(
        [
            survey,
            pd.DataFrame({"hole_id": ["DH9"], "depth": [0.0], "azimuth": [0.0], "dip": [90.0]}),
        ],
        ignore_index=True,
    )
    with pytest.raises(InvalidInputError, match="DH9"):
        drillholes_from_pandas(collar, extra, logs)


def test_object_column_with_non_string_value(collar, survey, logs):
    bad = logs.copy()
    bad["lith"] = pd.Series(["granite", 5, "chert"], dtype=object)
    with pytest.raises(InvalidInputError, match="lith"):
        drillholes_from_pandas(collar, survey, bad)


def test_no_logs_gives_empty_intervals(collar, survey):
    result = drillholes_from_pandas(collar, survey, None)
    assert list(result) == ["DH1", "DH2"]
    for hole in result.values():
        assert hole.intervals.shape == (0, 2)


def test_result_sends_to_viewer(collar, survey, logs, fake_viewer, viewer_factory):
    hole = drillholes_from_pandas(collar, survey, logs)["DH1"]
    viewer = viewer_factory()
    viewer.send(hole, id="DH1")
    viewer.flush()
    deadline = time.monotonic() + 5.0
    while time.monotonic() < deadline and "DH1" not in fake_viewer.viewer.objects:
        time.sleep(0.01)
    assert fake_viewer.viewer.objects["DH1"].header["object_type"] == "Drillhole"


def test_logs_hole_not_in_collar(collar, survey, logs):
    extra = pd.concat(
        [
            logs,
            pd.DataFrame(
                {
                    "hole_id": ["DH7"],
                    "from": [0.0],
                    "to": [1.0],
                    "grade": [1.0],
                    "lith": ["x"],
                }
            ),
        ],
        ignore_index=True,
    )
    with pytest.raises(InvalidInputError, match="DH7"):
        drillholes_from_pandas(collar, survey, extra)


def test_hole_without_log_rows(collar, survey, logs):
    only_dh1 = logs[logs["hole_id"] == "DH1"]
    hole = drillholes_from_pandas(collar, survey, only_dh1)["DH2"]
    assert hole.intervals.shape == (0, 2)
    assert hole.interval_properties["grade"].shape == (0,)
    assert hole.interval_properties["lith"].shape == (0,)


def test_hole_without_survey_rows(collar, survey, logs):
    only_dh1 = survey[survey["hole_id"] == "DH1"]
    hole = drillholes_from_pandas(collar, only_dh1, logs)["DH2"]
    assert hole.survey.shape == (0, 3)


def test_bool_and_int_columns_keep_dtype(collar, survey, logs):
    typed = logs.copy()
    typed["flag"] = [True, False, True]
    typed["count"] = [1, 2, 3]
    props = drillholes_from_pandas(collar, survey, typed)["DH1"].interval_properties
    assert props["flag"].dtype == np.bool_
    assert list(props["flag"]) == [True, False]
    assert np.issubdtype(props["count"].dtype, np.integer)
    assert list(props["count"]) == [1, 2]
