import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel
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


def test_from_processor_populates_stratigraphic_column():
    """
    Regression test: from_processor crashed with a DeprecationWarning raised by
    set_stratigraphic_column when assigning the processor's column to the model.
    """
    df = pd.DataFrame(rng.random(size=(10, 3)), columns=["X", "Y", "Z"])
    df["name"] = [f"unit_{name % 2}" for name in range(10)]
    stratigraphic_order = [("sg", ["unit_0", "unit_1", "basement"])]
    thicknesses = {"unit_0": 1.0, "unit_1": 0.5}
    processor = ProcessInputData(
        contacts=df,
        stratigraphic_order=stratigraphic_order,
        thicknesses=thicknesses,
        origin=np.zeros(3),
        maximum=np.ones(3),
    )

    model = GeologicalModel.from_processor(processor)

    unit_names = [unit.name for unit in model.stratigraphic_column.order if hasattr(unit, "name")]
    assert "unit_0" in unit_names
    assert "unit_1" in unit_names
