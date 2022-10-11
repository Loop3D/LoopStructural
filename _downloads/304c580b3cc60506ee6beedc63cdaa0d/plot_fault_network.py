"""
3b. Modelling a fault network in LoopStructural
===============================================
Uses GeologicalModel, ProcessInputData and LavaVuModelViewer from LoopStructural library. 
Also using geopandas to read a shapefile, pandas, matplotlib and numpy."""

import LoopStructural

LoopStructural.__version__

from LoopStructural import GeologicalModel
from LoopStructural.modelling import ProcessInputData
from LoopStructural.visualisation import LavaVuModelViewer
from LoopStructural.datasets import load_fault_trace
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

##############################
# Read shapefile
# ~~~~~~~~~~~~~~
# Read the shapefile and create a point for each node of the line
fault_trace = load_fault_trace()
faults = []
for i in range(len(fault_trace)):
    for x, y in zip(
        fault_trace.loc[i, :].geometry.xy[0], fault_trace.loc[i, :].geometry.xy[1]
    ):
        faults.append(
            [fault_trace.loc[i, "fault_name"], x, y, np.random.random() * 0.4]
        )  # better results if points aren't from a single plane
df = pd.DataFrame(faults, columns=["fault_name", "X", "Y", "Z"])

fig, ax = plt.subplots()
ax.scatter(df["X"], df["Y"])
ax.axis("square")

scale = np.min([df["X"].max() - df["X"].min(), df["Y"].max() - df["Y"].min()])
df["X"] /= scale
df["Y"] /= scale


##############################
# Orientation data
# ~~~~~~~~~~~~~~~~
# We can generate vertical dip data at the centre of the fault.

ori = []
for f in df["fault_name"].unique():
    centre = df.loc[df["fault_name"] == f, ["X", "Y", "Z"]].mean().to_numpy().tolist()
    tangent = (
        df.loc[df["fault_name"] == f, ["X", "Y", "Z"]].to_numpy()[0, :]
        - df.loc[df["fault_name"] == f, ["X", "Y", "Z"]].to_numpy()[-1, :]
    )
    norm = tangent / np.linalg.norm(tangent)
    norm = norm.dot(np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]]))
    ori.append([f, *centre, *norm])  # .extend(centre.extend(norm.tolist())))
# fault_orientations = pd.DataFrame([[
ori = pd.DataFrame(ori, columns=["fault_name", "X", "Y", "Z", "gx", "gy", "gz"])

##############################
# Model extent
# ~~~~~~~~~~~~
# # Calculate the bounding box for the model using the extent of the shapefiles. We make the Z coordinate 10% of the maximum x/y length.

z = np.max([df["X"].max(), df["Y"].max()]) - np.min([df["X"].min(), df["Y"].min()])
z *= 0.2
origin = [df["X"].min() - z, df["Y"].min() - z, -z]
maximum = [df["X"].max() + z, df["Y"].max() + z, z]

##############################
# Setting up the data
# ~~~~~~~~~~~~~~~~~~~
# The `ProcessInputData` class is used to convert common geological map components to the datastructures required by LoopStructural.#
# To build a fault network we need to provide:# * fault locations - a table of x,y,z, and the fault name
# 1. fault orientations - a table recording the orientation observations of the fault, e.g. strike, dip or normal vector and x,y,z, fault_name
# 2. origin - the origin of the model bounding box
# 3. maximum - the maximum extend of the model bounding box
# 4. fault_edges - list of intersection relationships between faults e.g. [('fault1','fault2')] indicates that there is a intersection between fault1 and fault2
# 5. fault_edge_properties - list of properties for the fault edges - this can be the type of intersection e.g. 'splay' or 'abut' or just the angle between the faults
# 6. fault_properties (*optional*)  - a pandas dataframe with any kwargs for the interpolator where the index is the fault name #
#
#  Below is an example of setting the number of interpolation elements for each fault

##############################
# Modelling splay faults
# ~~~~~~~~~~~~~~~~~~~~~~
# A splay fault relationship is defined for any fault where the angle between the faults is less than :math:`30^\circ`.
# In this example we specify the angle between the faults as :math:`10^\circ`.

processor = ProcessInputData(
    fault_orientations=ori,
    fault_locations=df,
    origin=origin,
    maximum=maximum,
    fault_edges=[("fault_2", "fault_1")],
    fault_edge_properties=[{"angle": 10}],
)

model = GeologicalModel.from_processor(processor)
model.update()

view = LavaVuModelViewer(model)
for f in model.faults:
    view.add_isosurface(f, slices=[0])  #
view.rotation = [-50.92916488647461, -30.319700241088867, -20.521053314208984]
view.display()

##############################
# Modelling abutting faults
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# In this exampe we will use the same faults but specify the angle between the faults as :math:`40^\circ` which will change
# the fault relationship to be abutting rather than splay.

processor = ProcessInputData(
    fault_orientations=ori,
    fault_locations=df,
    origin=origin,
    maximum=maximum,
    fault_edges=[("fault_2", "fault_1")],
    fault_edge_properties=[{"angle": 40}],
)

model = GeologicalModel.from_processor(processor)

view = LavaVuModelViewer(model)
for f in model.faults:
    view.add_isosurface(f, slices=[0])  #
    view.add_data(f[0], vectors=True)

view.rotation = [-50.92916488647461, -30.319700241088867, -20.521053314208984]
view.display()
