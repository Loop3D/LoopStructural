"""
1f. Unconformities and faults
==============================
This tutorial builds a model that combines both types of unconformity
supported by LoopStructural with a fault, and shows how to evaluate the
resulting scalar fields directly (without going through a Loop3DView)
for a 2D cross-section plotted with matplotlib.

* :code:`add_unconformity` adds an **erosional** unconformity - the
  surface truncates all older features that were added before it.
* :code:`add_onlap_unconformity` adds an **onlap** unconformity - younger
  features added afterwards only exist on one side of the surface,
  onlapping against it rather than eroding what's below.

Features are added to the model one at a time, and the order in which
they are added matters: unconformities and faults only affect features
that are added *after* them.
"""

import numpy as np
import pandas as pd
from LoopStructural import GeologicalModel
import matplotlib.pyplot as plt

# a single data point (with a normal vector) defines each foliation, plus
# one point on the fault surface with its slip direction (nx, ny, nz) and
# no value constraint (val=nan)
data = pd.DataFrame(
    [
        [100, 100, 150, 0.17, 0, 0.98, 0, "strati"],
        [100, 100, 170, 0, 0, 0.86, 0, "strati3"],
        [100, 100, 100, 0, 0, 1, 0, "strati2"],
        [100, 100, 50, 0, 0, 1, 0, "nconf"],
        [100, 100, 50, 0, 0, 1, 0, "strati4"],
        [700, 100, 190, 1, 0, 0, np.nan, "fault"],
    ],
    columns=["X", "Y", "Z", "nx", "ny", "nz", "val", "feature_name"],
)

model = GeologicalModel(np.zeros(3), np.array([1000, 1000, 200]))
model.data = data

# "strati2" is the oldest package - adding an unconformity on it means any
# feature added afterwards will be eroded/truncated where "strati2" < 0
model.create_and_add_foliation("strati2", buffer=0.0)
model.add_unconformity(model["strati2"], 0)

# the fault is added after the "strati2" unconformity, so it only displaces
# features created from this point onwards ("strati2" itself is unaffected)
model.create_and_add_fault(
    "fault",
    50,
    minor_axis=300,
    major_axis=500,
    intermediate_axis=300,
    fault_center=[700, 500, 0],
)

# "strati" is truncated by its own erosional unconformity in the same way
model.create_and_add_foliation("strati", buffer=0.0)
model.add_unconformity(model["strati"], 0)

# "strati3" is conformable with nothing above/below it - no unconformity added
model.create_and_add_foliation("strati3", buffer=0.0)

# "nconf" introduces an onlap unconformity: "strati4", added next, only
# exists where it onlaps against the "nconf" surface rather than eroding it
model.create_and_add_foliation("nconf", buffer=0.0)
model.add_onlap_unconformity(model["nconf"], 0)
model.create_and_add_foliation("strati4")

######################################################################
# Stratigraphic column
# ~~~~~~~~~~~~~~~~~~~~~
# Units are only defined here for "strati", "strati2" and "strati4" -
# "strati3" is left out deliberately to show that a feature can still be
# evaluated directly even if it isn't part of the final stratigraphic
# column. Groups are added oldest-first via
# :code:`model.stratigraphic_column.add_unit`/:code:`add_unconformity`,
# so that "strati2" (oldest) ends up at the bottom of the column and
# "strati4" (youngest) at the top.

# "strati2" (oldest group)
model.stratigraphic_column.add_unit("series1", thickness=2.0, id=0, colour="red")
model.stratigraphic_column.add_unit("series2", thickness=3.0, id=1, colour="red")
model.stratigraphic_column.add_unit("series3", thickness=5.0, id=2, colour="red")
model.stratigraphic_column.add_unconformity(name="strati2_unconformity")

# "strati"
model.stratigraphic_column.add_unit("series2", thickness=np.inf, id=3, colour="blue")
model.stratigraphic_column.add_unit("series3", thickness=np.inf, id=4, colour="blue")
model.stratigraphic_column.add_unconformity(name="strati_unconformity")

# "strati4" (youngest group)
model.stratigraphic_column.add_unit("series4", thickness=np.inf, id=5)

model.stratigraphic_column.group_mapping["Group_0"] = "strati4"
model.stratigraphic_column.group_mapping["Group_1"] = "strati"
model.stratigraphic_column.group_mapping["Group_2"] = "strati2"

######################################################################
# Evaluating features directly on a cross-section
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rather than using the Loop3DView, each feature is evaluated by hand at a
# regular grid of points on a vertical (X-Z) slice through the model, and
# plotted as a stack of filled contours. This is a useful pattern when you
# want to combine model output with your own custom matplotlib figure.

xx, zz = np.meshgrid(np.linspace(0, 1000, 100), np.linspace(0, 200, 100))
yy = np.zeros_like(xx) + 500
points = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
val = model["strati"].evaluate_value(points)
val2 = model["strati2"].evaluate_value(points)
val3 = model["strati3"].evaluate_value(points)
val4 = model["strati4"].evaluate_value(points)
# .regions[0] is the onlap region mask added by add_onlap_unconformity -
# evaluates to True where "strati4" is present
uf = model["strati4"].regions[0](points)
fval = model['fault'].evaluate_value(points)

fig, ax = plt.subplots(figsize=(10, 3))
ax.contourf(val.reshape((100, 100)), extent=(0, 1000, 0, 200), cmap='viridis')
ax.contourf(val2.reshape((100, 100)), extent=(0, 1000, 0, 200), cmap='Reds')
ax.contourf(val3.reshape((100, 100)), extent=(0, 1000, 0, 200), cmap='Blues')
ax.contourf(val4.reshape((100, 100)), extent=(0, 1000, 0, 200), cmap='Greens')
# overlay the fault surface (0-isovalue of the fault scalar field) as a line
ax.contour(fval.reshape((100, 100)), [0], extent=(0, 1000, 0, 200), colors='k')
ax.set_xlabel("X")
ax.set_ylabel("Z")
plt.show()
