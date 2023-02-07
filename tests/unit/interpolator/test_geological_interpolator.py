import pytest
import numpy as np


def test_get_data_locations(interpolator, data):
    interpolator.set_value_constraints(
        data.loc[~data["val"].isna(), ["X", "Y", "Z", "val", "w"]].to_numpy()
    )
    interpolator.set_normal_constraints(
        data.loc[~data["nx"].isna(), ["X", "Y", "Z", "nx", "ny", "nz", "w"]].to_numpy()
    )
    locations = interpolator.get_data_locations()
    assert np.sum(locations - data[["X", "Y", "Z"]].to_numpy()) == 0


def test_get_value_constraints(interpolator, data):
    interpolator.set_value_constraints(
        data.loc[~data["val"].isna(), ["X", "Y", "Z", "val", "w"]].to_numpy()
    )
    interpolator.set_normal_constraints(
        data.loc[~data["nx"].isna(), ["X", "Y", "Z", "nx", "ny", "nz", "w"]].to_numpy()
    )
    val = interpolator.get_value_constraints()
    assert (
        np.sum(
            val - data.loc[~data["val"].isna(), ["X", "Y", "Z", "val", "w"]].to_numpy()
        )
        == 0
    )


def test_get_norm_constraints(interpolator, data):
    interpolator.set_value_constraints(
        data.loc[~data["val"].isna(), ["X", "Y", "Z", "val", "w"]].to_numpy()
    )
    interpolator.set_normal_constraints(
        data.loc[~data["nx"].isna(), ["X", "Y", "Z", "nx", "ny", "nz", "w"]].to_numpy()
    )
    val = interpolator.get_norm_constraints()
    assert (
        np.sum(
            val
            - data.loc[
                ~data["nx"].isna(), ["X", "Y", "Z", "nx", "ny", "nz", "w"]
            ].to_numpy()
        )
        == 0
    )


def test_reset(interpolator, data):
    interpolator.set_value_constraints(
        data.loc[~data["val"].isna(), ["X", "Y", "Z", "val", "w"]].to_numpy()
    )
    interpolator.set_normal_constraints(
        data.loc[~data["nx"].isna(), ["X", "Y", "Z", "nx", "ny", "nz", "w"]].to_numpy()
    )
    interpolator.clean()
    assert interpolator.get_data_locations().shape[0] == 0
    assert interpolator.up_to_date == False
