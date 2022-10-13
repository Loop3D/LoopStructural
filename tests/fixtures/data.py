import numpy as np
import pandas as pd
import pytest


@pytest.fixture(params=[0, 1, 2])
def data(request):
    data_list = []
    if request.param == 0:
        value = True
        gradient = False
    if request.param == 1:
        value = False
        gradient = True

    if request.param == 2:
        value = True
        gradient = True
    if value:
        xy = np.array(
            np.meshgrid(np.linspace(0, 1, 50), np.linspace(0, 1, 50))
        ).T.reshape(-1, 2)
        xyz = np.hstack([xy, np.zeros((xy.shape[0], 1))])
        data = pd.DataFrame(xyz, columns=["X", "Y", "Z"])
        data["val"] = np.sin(data["X"])
        data["w"] = 1
        data["feature_name"] = "strati"
        data_list.append(data)
    if gradient:
        data = pd.DataFrame(
            [[0.5, 0.5, 0.5, 0, 0, 1], [0.75, 0.5, 0.75, 0, 0, 1]],
            columns=["X", "Y", "Z", "nx", "ny", "nz"],
        )
        data["w"] = 1
        data["feature_name"] = "strati"
        data_list.append(data)
    if "nx" not in data:
        data["nx"] = np.nan
        data["ny"] = np.nan
        data["nz"] = np.nan
    if "val" not in data:
        data["val"] = np.nan
    return pd.concat(data_list, ignore_index=True)
