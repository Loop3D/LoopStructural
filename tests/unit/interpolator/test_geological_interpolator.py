import pytest


@pytest.mark.skip("not currently working 4-10-2022")
def test_get_data_locations(interpolator, data):
    interpolator.set_value_constraints(
        data.loc[~data["val"].isna(), ["X", "Y", "Z", "val", "w"]].to_numpy()
    )
    interpolator.set_normal_constraints(
        data.loc[~data["nx"].isna(), ["X", "Y", "Z", "nx", "ny", "nz", "w"]].to_numpy()
    )
    interpolator.get_data_locations()
