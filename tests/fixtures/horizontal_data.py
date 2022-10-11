import numpy as np
import pandas as pd
import pytest


@pytest.fixture()
def horizontal_data():

    xy = np.array(np.meshgrid(np.linspace(0, 1, 50), np.linspace(0, 1, 50))).T.reshape(
        -1, 2
    )
    df1 = pd.DataFrame(xy, columns=["X", "Y"])
    df2 = pd.DataFrame(xy, columns=["X", "Y"])
    df1["Z"] = 0.25
    df1["val"] = 0
    df2["Z"] = 0.55
    df2["val"] = 0.3
    data = pd.concat([df1, df2], ignore_index=True)
    data["w"] = 1
    data["feature_name"] = "strati"

    return data
