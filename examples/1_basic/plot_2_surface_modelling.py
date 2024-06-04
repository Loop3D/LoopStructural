"""

1b. Implicit surface modelling
===============================
This tutorial will demonstrate how to create an implicit surface
representation of surfaces from a combination of orientation and
location observations.

Implicit surface representation involves finding an unknown function
where :math:`f(x,y,z)` matches observations of the surface geometry. We
generate a scalar field where the scalar value is the distance away from
a reference horizon. The reference horizon is arbritary and can either
be:

-  a single geological surface where the scalar field would represent
   the signed distance away from this surface. (above the surface
   positive and below negative)
-  Where multiple conformable horizons are observed the same scalar
   field can be used to represent these surfaces and the thickness of
   the layers is used to determine the relative scalar value for each
   surface

This tutorial will demonstrate both of these approaches for modelling a
number of horizons picked from seismic data sets, by following the next
steps: 1. Creation of a geological model, which includes: \*
Presentation and visualization of the data \* Addition of a geological
feature, which in this case is the stratigraphy of the model. 2.
Visualization of the scalar field.

"""

#########################################################################
# Imports
# ~~~~~~~
# Import the required objects from LoopStructural for visualisation and
# model building

from LoopStructural import GeologicalModel
from LoopStructural.visualisation import Loop3DView
from LoopStructural.datasets import load_claudius  # demo data

import numpy as np


######################################################################
# The data for this example can be imported from the example datasets
# module in loopstructural.
#

data, bb = load_claudius()

data["val"].unique()


######################################################################
# GeologicalModel
# ~~~~~~~~~~~~~~~
#
# To create a ``GeologicalModel`` the origin (lower left) and max extent
# (upper right) of the model area need to be specified. In this example
# these are provided in the bb array.
#

model = GeologicalModel(bb[0, :], bb[1, :])


######################################################################
# A pandas dataframe with appropriate columns can be used to link the data
# to the geological model.
#
# * ``X`` is the x coordinate
# * ``Y`` is the y # coordinate
# * ``Z`` is the z coordinate
# * ``feature_name`` is a name to link the data to a model object
# * ``val`` is the value of the scalar field which represents the
# distance from a reference horizon. It is comparable
# to the relative thickness
#
# * ``nx`` is the x component of the normal vector to the surface gradient
# * ``ny`` is the y component of the normal vector to the surface gradient
# * ``nz`` is the z component of the normal vector to the surface gradeint
# * ``strike`` is the strike angle
# * ``dip`` is the dip angle
#
# Having a look at the data for this example by looking at the top of the
# dataframe and then using a 3D plot
#

data["feature_name"].unique()

viewer = Loop3DView(background="white")
viewer.add_points(
    data[~np.isnan(data["val"])][["X", "Y", "Z"]].values,
    scalars=data[~np.isnan(data["val"])]["val"].values,
)
viewer.add_arrows(
    data[~np.isnan(data["nx"])][["X", "Y", "Z"]].values,
    direction=data[~np.isnan(data["nx"])][["nx", "ny", "nz"]].values,
)
viewer.display()


######################################################################
# The pandas dataframe can be linked to the ``GeologicalModel`` using
# ``.set_model_data(dataframe``
#

model.set_model_data(data)


######################################################################
# The ``GeologicalModel`` contains an ordered list of the different
# geological features within a model and how these features interact. This
# controls the topology of the different geological features in the model.
# Different geological features can be added to the geological model such
# as:
# * Foliations
# * Faults
# * Unconformities
# * Folded foliations
# *  Structural Frames
#
# In this example we will only add a foliation using the function
#
# .. code:: python
#
#    model.create_and_add_foliation(name)
#
# where name is the name in the ``feature_name`` field, other parameters we
# specified are the:
# * ``interpolatortype`` - we can either use a
# PiecewiseLinearInterpolator ``PLI``, a FiniteDifferenceInterpolator
# ``FDI`` or a radial basis interpolator ``surfe``
# * ``nelements - int`` is the how many elements are used to discretize the resulting solution
# * ``buffer - float`` buffer percentage around the model area
# * ``solver`` - the algorithm to solve the least squares problem e.g.
# ``lu`` for lower upper decomposition, ``cg`` for conjugate gradient,
# ``pyamg`` for an algorithmic multigrid solver
# * ``damp - bool`` - whether to add a small number to the diagonal of the interpolation
# matrix for discrete interpolators - this can help speed up the solver
# and makes the solution more stable for some interpolators
#

vals = [0, 60, 250, 330, 600]

strat_column = {"strati": {}}
for i in range(len(vals) - 1):
    strat_column["strati"]["unit_{}".format(i)] = {
        "min": vals[i],
        "max": vals[i + 1],
        "id": i,
    }


model.set_stratigraphic_column(strat_column)

strati = model.create_and_add_foliation(
    "strati",
    interpolatortype="FDI",  # try changing this to 'PLI'
    nelements=1e4,  # try changing between 1e3 and 5e4
    buffer=0.3,
    damp=True,
)
######################################################################
# Plot the surfaces
# ------------------------------------

viewer = Loop3DView(model)
viewer.plot_model_surfaces(cmap="tab20")
viewer.display()

######################################################################
# Plot block diagram
# -------------------

viewer = Loop3DView(model)
viewer.plot_block_model(cmap="tab20")
viewer.display()
