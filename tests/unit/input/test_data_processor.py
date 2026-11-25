import numpy as np
import pandas as pd

from LoopStructural.modelling import ProcessInputData
from LoopStructural.utils import rng


def test_create_processor():
    df = pd.DataFrame(rng.random(size=(10, 3)), columns=["X", "Y", "Z"])
    df["name"] = [f"unit_{name % 2}" for name in range(10)]
    stratigraphic_order = [("sg", ["unit_0", "unit_1", "basement"])]
    thicknesses = {"unit_0": 1.0, "unit_1": 0.5}
    processor = ProcessInputData(
        contacts=df, stratigraphic_order=stratigraphic_order, thicknesses=thicknesses
    )
    assert (processor.data["val"].unique() == np.array([0.5, 0])).all()


def test_from_processor_unit_ranges_match_processor_values():
    """
    Regression test: the basal unit (min=-inf, max=0) gave an infinite
    thickness, so every unit above it got the range (inf, inf).
    """
    import warnings

    from LoopStructural import GeologicalModel
    from LoopStructural.datasets import load_geological_map_data

    contacts, orientations, thickness, order, bbox, *_ = load_geological_map_data()
    processor = ProcessInputData(
        contacts=contacts,
        contact_orientations=orientations.rename({"formation": "name"}, axis=1),
        thicknesses=dict(zip(thickness["name"], thickness["thickness"])),
        stratigraphic_order=[("supergroup_0", list(order["unit name"]))],
        origin=bbox.loc["origin"].to_numpy(),
        maximum=bbox.loc["maximum"].to_numpy(),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        model = GeologicalModel.from_processor(processor)

    column = model.stratigraphic_column
    expected = processor.stratigraphic_column["supergroup_0"]

    def check():
        isovalues = column.get_isovalues()
        for name, values in expected.items():
            unit = column.get_unit_by_name(name)
            assert (unit.min(), unit.max()) == (values["min"], values["max"])
            assert isovalues[name]["value"] == values["min"]

    check()
    # a later recompute of the column must keep the same ranges
    column.add_unit("new_top", thickness=10)
    check()
    unit = column.get_unit_by_name("Bunjinah_Formation")
    unit.thickness = unit.thickness
    check()
