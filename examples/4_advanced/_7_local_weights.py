"""
============================
1f. Local data weighting
============================
LoopStructural primarily uses discrete interpolation methods (e.g. finite differences on a regular grid,
or linear/quadratic on tetrahedral meshes). The interpolation is determined by combining a regularisation
term and the data weights. The default behaviour is for every data point to be weighted equally, however
it is also possible to vary these weights per datapoint.

"""

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius
from LoopStructural.visualisation import Loop3DView

##################################################################################################
# Use Cladius case study
# ~~~~~~~~~~~~~~~~~~~~~~~~
#
data, bb = load_claudius()
data.head()
##################################################################################################
# Build model with constant weighting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Build model with weight 1.0 for the control points (cpw) and gradient normal constraints (npw)
model = GeologicalModel(bb[0, :], bb[1, :])
model.data = data
model.create_and_add_foliation("strati", interpolatortype="FDI", cpw=1.0, npw=1.0)
view = Loop3DView(model)
view.add_isosurface(model["strati"], slices=data["val"].dropna().unique())
view.rotation = [84.24493408203125, -77.89686584472656, -171.99795532226562]
view.display()
##################################################################################################
# Change weights to

model = GeologicalModel(bb[0, :], bb[1, :])
model.data = data
model.create_and_add_foliation("strati", interpolatortype="FDI", cpw=10.0, npw=1.0)
view = Loop3DView(model)
view.add_isosurface(model["strati"], slices=data["val"].dropna().unique())
view.rotation = [84.24493408203125, -77.89686584472656, -171.99795532226562]
view.display()

##################################################################################################
# Locally vary weights
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add a weight column to the dataframe and decrease the weighting of the points
# in the North of the model.
data, bb = load_claudius()
data["w"] = 1.0
data.loc[data["Y"] > (bb[1, 1] - bb[0, 1]) * 0.2 + bb[0, 1], "w"] = 0.01
data.sample(10)

# cpw/npw are multipliers for the weight column
model.create_and_add_foliation("strati", cpw=1.0, npw=1)
view = Loop3DView(model)
view.add_isosurface(model["strati"], slices=data["val"].dropna().unique())
view.rotation = [84.24493408203125, -77.89686584472656, -171.99795532226562]
view.display()
