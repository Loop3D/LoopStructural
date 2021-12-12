from LoopStructural.modelling import ProcessInputData
import pandas as pd
import numpy as np


def test_create_processor():
    df = pd.DataFrame(np.random.rand(10, 3), columns=["X", "Y", "Z"])
    df["name"] = ["unit_{}".format(name % 2) for name in range(10)]
    stratigraphic_order = [("sg", ["unit_0", "unit_1"])]
    thicknesses = {"unit_0": 1.0, "unit_1": 0.5}
    processor = ProcessInputData(
        contacts=df, stratigraphic_order=stratigraphic_order, thicknesses=thicknesses
    )
    assert (processor.data["val"].unique() == np.array([0.5, 0])).all()
