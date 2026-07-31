"""
3b. Modelling a fault network in LoopStructural
===============================================
Real fault networks are rarely made up of isolated faults - they interact
with each other, and the way two faults meet (splaying off one another,
or abutting against each other) affects how displacement is distributed
between them. This tutorial builds a network of two interacting faults
from fault traces digitised from a geological map, using
:code:`ProcessInputData` to turn the traces into a model and
:code:`fault_edge_properties` to control how the faults interact.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_fault_trace
from LoopStructural.modelling import ProcessInputData
from LoopStructural.utils import rng
from LoopStructural.visualisation import Loop3DView

##############################
# Read shapefile
# ~~~~~~~~~~~~~~
# Read the shapefile and create a point for each node of the line
fault_trace = load_fault_trace()
faults = []
for i in range(len(fault_trace)):
    for x, y in zip(fault_trace.loc[i, :].geometry.xy[0], fault_trace.loc[i, :].geometry.xy[1]):
        faults.append(
            [fault_trace.loc[i, "fault_name"], x, y, rng.random() * 0.4]
        )  # better results if points aren't from a single plane
df = pd.DataFrame(faults, columns=["fault_name", "X", "Y", "Z"])

fig, ax = plt.subplots()
ax.scatter(df["X"], df["Y"])
ax.axis("square")
plt.show()

# rescale coordinates so the model is a sensible size for the default
# interpolation settings
scale = np.min([df["X"].max() - df["X"].min(), df["Y"].max() - df["Y"].min()])
df["X"] /= scale
df["Y"] /= scale


##############################
# Orientation data
# ~~~~~~~~~~~~~~~~
# The map only gives the trace (location) of each fault, not its dip - we
# generate a vertical dip vector at the centre of each fault trace, using
# the along-trace tangent (rotated 90 degrees) as the strike direction.

ori = []
for f in df["fault_name"].unique():
    centre = df.loc[df["fault_name"] == f, ["X", "Y", "Z"]].mean().to_numpy().tolist()
    tangent = (
        df.loc[df["fault_name"] == f, ["X", "Y", "Z"]].to_numpy()[0, :]
        - df.loc[df["fault_name"] == f, ["X", "Y", "Z"]].to_numpy()[-1, :]
    )
    norm = tangent / np.linalg.norm(tangent)
    norm = norm.dot(np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]]))
    ori.append([f, *centre, *norm])
ori = pd.DataFrame(ori, columns=["fault_name", "X", "Y", "Z", "gx", "gy", "gz"])

##############################
# Model extent
# ~~~~~~~~~~~~
# Calculate the bounding box for the model using the extent of the fault
# traces, buffered by 20% of the extent in each direction (also used for
# the vertical extent, since the traces carry no depth information).

z = np.max([df["X"].max(), df["Y"].max()]) - np.min([df["X"].min(), df["Y"].min()])
z *= 0.2
origin = [df["X"].min() - z, df["Y"].min() - z, -z]
maximum = [df["X"].max() + z, df["Y"].max() + z, z]


##############################
# Modelling abutting faults
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# ``fault_edges`` declares that "fault_2" interacts with "fault_1", and
# ``fault_edge_properties`` sets the angle between them to :math:`40^\circ`.
# LoopStructural uses this angle to decide the fault relationship: faults
# that meet at a shallow angle are treated as **splay** faults (one
# branches off the other and shares its displacement), while faults that
# meet at a higher angle - as here - are treated as **abutting** (one
# fault truncates against the other, each keeping an independent
# displacement field).

processor = ProcessInputData(
    fault_orientations=ori,
    fault_locations=df,
    origin=origin,
    maximum=maximum,
    fault_edges=[("fault_2", "fault_1")],
    fault_edge_properties=[{"angle": 40}],
)

model = GeologicalModel.from_processor(processor)

view = Loop3DView(model)
for f in model.faults:
    view.plot_surface(f, value=[0])
    view.plot_data(f[0])

view.display()
