"""
4c. Local data weighting
============================
LoopStructural primarily uses discrete interpolation methods (e.g. finite
differences on a regular grid, or linear/quadratic on tetrahedral
meshes). The interpolation is determined by combining a regularisation
term and the data weights. The default behaviour is for every data point
to be weighted equally, however it is also possible to vary these
weights per-datapoint or uniformly across the whole dataset.
"""

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius
from LoopStructural.visualisation import Loop3DView

##################################################################################################
# Use the Claudius case study
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data, bb = load_claudius()
data.head()

##################################################################################################
# Build model with constant weighting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ``cpw``/``npw`` set the weight for the control (value) points and the
# gradient normal constraints respectively, applied uniformly to all data
# of that type - here both are left at the default of 1.0, weighted
# equally against the regularisation term.
model = GeologicalModel(bb[0, :], bb[1, :])
model.data = data
model.create_and_add_foliation(
    "strati", nelements=10_000, interpolatortype="FDI", cpw=1.0, npw=1.0, regularisation=1.0
)
view = Loop3DView(model)
view.plot_surface(model["strati"], value=data["val"].dropna().unique())
view.display()

##################################################################################################
# Increase the weight of the value constraints
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Raising ``cpw`` relative to ``npw`` and the regularisation term makes
# the interpolator honour the value observations more closely, at the
# cost of a less smooth surface.

model = GeologicalModel(bb[0, :], bb[1, :])
model.data = data
model.create_and_add_foliation("strati", interpolatortype="FDI", cpw=10.0, npw=1.0, regularisation=1.0)
view = Loop3DView(model)
view.plot_surface(model["strati"], value=data["val"].dropna().unique())
view.display()

##################################################################################################
# Locally vary weights
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rather than a single uniform weight, an optional ``w`` column in the
# input data is picked up automatically and used as a per-point weight
# multiplier - here reduced to 1% for points in the northern part of the
# model, so those observations barely constrain the surface at all.

data, bb = load_claudius()
data["w"] = 1.0
data.loc[data["Y"] > (bb[1, 1] - bb[0, 1]) * 0.2 + bb[0, 1], "w"] = 0.01
data.sample(10)

model = GeologicalModel(bb[0, :], bb[1, :])
model.data = data
# cpw/npw are multipliers applied on top of the per-point "w" column
model.create_and_add_foliation("strati", cpw=1.0, npw=1, regularisation=1.0)
view = Loop3DView(model)
view.plot_surface(model["strati"], value=data["val"].dropna().unique())
view.display()
