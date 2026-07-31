"""
1e. Using Stratigraphic Columns
===============================
The previous example (Multiple groups) built a model from two separately
interpolated scalar fields but stopped short of naming the rock units they
represent. A **stratigraphic column** maps ranges of scalar field value
within each group to named units with an integer id, which is what allows
LoopStructural to evaluate a single "which unit is here" answer at any
point in the model via :code:`model.evaluate_model(xyz)`, and to produce
a labelled block model.

This example reuses the two-group model from the previous tutorial and
defines a stratigraphic column for it.
"""

import numpy as np

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius
from LoopStructural.visualisation import Loop3DView

data, bb = load_claudius()
data = data.reset_index()

data.loc[:, "val"] *= -1
data.loc[:, ["nx", "ny", "nz"]] *= -1

data.loc[792, "feature_name"] = "strati2"
data.loc[792, ["nx", "ny", "nz"]] = [0, 0, 1]
data.loc[792, "val"] = 0

model = GeologicalModel(bb[0, :], bb[1, :])
model.data = data

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
# ~~~~~~~~~~~~~~~~~~~~~
# ``model.stratigraphic_column`` is a :class:`StratigraphicColumn` object.
# Units are added with :code:`add_unit(name, thickness=..., id=...)`, from
# the oldest/deepest unit up to the youngest/shallowest, and
# :code:`add_unconformity(name=...)` marks the boundary between two groups
# (features). Each unit's ``thickness`` sets how much of the group's
# scalar field range it occupies - ranges are assigned automatically,
# resetting to zero at every unconformity - and can be :code:`np.inf` for
# the oldest unit in a group. A unique integer ``id`` is used to label
# each unit in the block model.

# "strati" (oldest group) - four units from shallowest to deepest, the
# last of which extends to infinite thickness
model.stratigraphic_column.add_unit("unit2", thickness=60, id=1)
model.stratigraphic_column.add_unit("unit3", thickness=190, id=2)
model.stratigraphic_column.add_unit("unit4", thickness=80, id=3)
model.stratigraphic_column.add_unit("unit5", thickness=np.inf, id=4)

# mark the boundary between the "strati" and "strati2" groups
model.stratigraphic_column.add_unconformity(name="strati_unconformity")

# "strati2" (youngest group) - a single unit
model.stratigraphic_column.add_unit("unit1", thickness=9, id=0)

model.stratigraphic_column.group_mapping["Group_0"] = "strati2"
model.stratigraphic_column.group_mapping["Group_1"] = "strati"

########################################################
# Adding stratigraphic column to the model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# With the stratigraphic column defined on the model, the
# `model.evaluate_model(xyz)` function can be called.

viewer = Loop3DView(model)
viewer.plot_block_model(cmap='tab20')
viewer.display()
