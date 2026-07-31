import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel


def test_outside_box_normal():
    for interpolator in ["PLI", "FDI"]:
        print(f"Running test for {interpolator} with normal constraints")
        model = GeologicalModel(np.zeros(3), np.ones(3))
        data = pd.DataFrame(
            [
                [0.5, 0.5, 0.5, 0, 1.0, 0.0, np.nan, "strati"],
                [1.5, 0.5, 0.5, 0, 1.0, 0.0, np.nan, "strati"],
                [0.5, 1.5, 1.5, 0, 1.0, 0.0, 1.0, "strati"],
            ],
            columns=["X", "Y", "Z", "nx", "ny", "nz", "val", "feature_name"],
        )
        model.data = data
        model.create_and_add_foliation("strati", interpolatortype=interpolator)
        model.update()


def test_outside_box_value():
    for interpolator in ["PLI", "FDI"]:
        print(f"Running test for {interpolator} with normal constraints")

        model = GeologicalModel(np.zeros(3), np.ones(3))
        data = pd.DataFrame(
            [
                [0.5, 0.5, 0.5, 0, "strati"],
                [1.5, 0.5, 0.5, 0, "strati"],
                [0.5, 1.5, 1.5, 0, "strati"],
            ],
            columns=["X", "Y", "Z", "val", "feature_name"],
        )
        model.data = data
        model.create_and_add_foliation("strati", interpolatortype=interpolator)
        model.update()


if __name__ == "__main__":
    test_outside_box_normal()
    test_outside_box_value()
