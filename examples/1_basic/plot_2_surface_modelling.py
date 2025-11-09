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
# Load Example Data
# ~~~~~~~~~~~~~~~~~
# The data for this example can be imported from the example datasets
# module in LoopStructural. The `load_claudius` function provides both
# the data and the bounding box (bb) for the model.

data, bb = load_claudius()

# Display unique scalar field values in the dataset
print(data["val"].unique())

######################################################################
# Create Geological Model
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Define the model area using the bounding box (bb) and create a
# GeologicalModel instance.

model = GeologicalModel(bb[0, :], bb[1, :])

######################################################################
# Link Data to Geological Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A pandas dataframe with appropriate columns is used to link the data
# to the geological model. Key columns include:
# * `X`, `Y`, `Z`: Coordinates of the observation
# * `feature_name`: Name linking the data to a model object
# * `val`: Scalar field value representing distance from a reference horizon
# * `nx`, `ny`, `nz`: Components of the normal vector to the surface gradient
# * `strike`, `dip`: Strike and dip angles

# Display unique feature names in the dataset
print(data["feature_name"].unique())

# Visualize the data points and orientation vectors in 3D
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

# Link the data to the geological model
model.set_model_data(data)

######################################################################
# Add Geological Features
# ~~~~~~~~~~~~~~~~~~~~~~~~
# The GeologicalModel can include various geological features such as
# foliations, faults, unconformities, and folds. In this example, we
# add a foliation using the `create_and_add_foliation` method.

# Define stratigraphic column with scalar field ranges for each unit
vals = [0, 60, 250, 330, 600]
for i in range(len(vals) - 1):
    model.stratigraphic_column.add_unit(
        f"unit_{i}",
        thickness= vals[i + 1] - vals[i],
        id=i,
    )
model.stratigraphic_column.group_mapping['Group_0'] ='strati'
# Add a foliation to the model
strati = model.create_and_add_foliation(
    "strati",
    interpolatortype="FDI",  # Finite Difference Interpolator
    nelements=int(1e4),  # Number of elements for discretization
    buffer=0.3,  # Buffer percentage around the model area
    damp=True,  # Add damping for stability
)

######################################################################
# Visualize Model Surfaces
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Plot the surfaces of the geological model using a 3D viewer.
viewer = Loop3DView(model)
viewer.plot_model_surfaces(cmap="tab20")
viewer.display()

######################################################################
# Visualize Block Diagram
# ~~~~~~~~~~~~~~~~~~~~~~~
# Plot a block diagram of the geological model to visualize the
# stratigraphic units in 3D.

viewer = Loop3DView(model)
viewer.plot_block_model(cmap="tab20")
viewer.display()
