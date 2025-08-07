"""
1d. Using Stratigraphic Columns
===============================
We will use the previous example Creating a model with multiple geological features, dealing with unconformities.

"""

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius
from LoopStructural.visualisation import Loop3DView

import numpy as np

data, bb = load_claudius()
data = data.reset_index()

data.loc[:, "val"] *= -1
data.loc[:, ["nx", "ny", "nz"]] *= -1

data.loc[792, "feature_name"] = "strati2"
data.loc[792, ["nx", "ny", "nz"]] = [0, 0, 1]
data.loc[792, "val"] = 0

model = GeologicalModel(bb[0, :], bb[1, :])
model.set_model_data(data)

strati2 = model.create_and_add_foliation(
    "strati2",
    interpolatortype="FDI",
    nelements=1e4,
)
uc = model.add_unconformity(strati2, 1)

strati = model.create_and_add_foliation(
    "strati",
    interpolatortype="FDI",
    nelements=1e4,
)

########################################################################
# Stratigraphic columns
# ~~~~~~~~~~~~~~~~~~~~~~~
# We define the stratigraphic column using a nested dictionary

stratigraphic_column = {}
stratigraphic_column["strati2"] = {}
stratigraphic_column["strati2"]["unit1"] = {"min": 1, "max": 10, "id": 0}
stratigraphic_column["strati"] = {}
stratigraphic_column["strati"]["unit2"] = {"min": -60, "max": 0, "id": 1}
stratigraphic_column["strati"]["unit3"] = {"min": -250, "max": -60, "id": 2}
stratigraphic_column["strati"]["unit4"] = {"min": -330, "max": -250, "id": 3}
stratigraphic_column["strati"]["unit5"] = {"min": -np.inf, "max": -330, "id": 4}

########################################################
# Adding stratigraphic column to the model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The stratigraphic column can be added to the geological model. Allowing
# for the `model.evaluate_model(xyz)` function to be called.

model.set_stratigraphic_column(stratigraphic_column)

viewer = Loop3DView(model)
viewer.plot_block_model(cmap='tab20')
viewer.display()
