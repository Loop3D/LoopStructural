"""

1c. Visualising models
===============================
The following tutorial will demonstrate how to use the Loop structural visualisation module. 
This module provides a wrapper for the lavavu model that is written by
Owen Kaluza. 

Lavavu allows for interactive visualisation of 3D models within a jupyter
notebook environment. 

"""

#########################################################################
# Imports
# ~~~~~~~
# Import the required objects from LoopStructural for visualisation and
# model building

from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer

from LoopStructural.datasets import load_claudius  # demo data

import numpy as np

#####################
# Build the model
# ~~~~~~~~~~~~~~~~~
data, bb = load_claudius()
model = GeologicalModel(bb[0, :], bb[1, :])
model.set_model_data(data)
strati = model.create_and_add_foliation("strati")
strat_column = {"strati": {}}
vals = [0, 60, 250, 330, 600]
for i in range(len(vals) - 1):
    strat_column["strati"]["unit_{}".format(i)] = {
        "min": vals[i],
        "max": vals[i + 1],
        "id": i,
    }
model.set_stratigraphic_column(strat_column)

######################################################################
# Visualising results
# ~~~~~~~~~~~~~~~~~~~
#
# The LavaVuModelViewer is an LoopStructural class that provides easy 3D
# plotting options for plotting data points and resulting implicit
# functions.
#
# The implicit function can be visualised by looking at isosurfaces of the
# scalar field.
#
# .. code:: python
#
#    viewer = LavaVuModelViewer()
#    viewer.add_isosurface(feature,**kwargs)
#
# Where optional kwargs can be:
#
# -  ``nslices`` specifying the number of regularly spaced isosurfaces
# -  ``slices`` a numpy array or list of isovalues to slice
# -  ``isovalue`` an isovalue to slice
# -  ``paint_with`` the geological feature to colour the surface with
# -  ``cmap`` colour map for the colouring
# -  ``normals`` to plot the normal vectors to the surface
# -  ``name`` to give the surface
# -  ``colour`` the colour of the surface
# -  other kwargs for passing directly to lavavu
#
# Alternatively the scalar fields can be displayed on a rectangular cuboid.
#
# .. code:: python
#
#    viewer.add_scalar_field(geological_feature)
#
#
# Other possible kwargs are:
#
# -  ``cmap`` colour map for the property
#
# The input data for the model can be visualised by calling either:
#
# .. code:: python
#
#    viewer.add_data(feature,addgrad=True,addvalue=True,**kwargs)
#
# Where both the point and vector data linked to the feature are added to
# the plot or by calling.
#
# .. code:: python
#
#    viewer.add_vector_data(position,vector,name,**kwargs)
#
# Where ``position`` is an array or x,y,z coordinates and vector is a
# similarly sized array of ``vectors``. These can be extracted from a
# geological feature by calling.
# ``feature.support.interpolator.get_gradient_constraint()`` which returns
# a Nx6 matrix of position and vectors.
#
# The value data can be plotted by calling.
#
# .. code:: python
#
#    viewer.add_value_data(position,value,name,**kwargs)
#
# Where ``position`` is an array or x,y,z coordinates and value is a
# similarly sized vector of values. These can be extracted from a
# geological feature by calling.
# ``feature.support.interpolator.get_value_constraint()`` which returns a
# Nx4 matrix of position and values.
#
# Other possible options for plotting are to \* plot point locations.
#
# .. code:: python
#
#    viewer.add_points(position, name, **kwargs)
#
# -  plot a vector field using the gradient of a geological feature
#
# .. code:: python
#
#    viewer.add_vector_field(feature, **kwargs)
#
# Where ``locations`` can be specified to control specific evaluation locations
# It is recommended to visualise
# the vectorfield at a lower resolution than the mesh otherwise it can be
# difficult to see the vectors. You can use numpy stepping along the
# array: ``locations = mesh.barycentre[::20,:]`` which will sample every
# 20th sample in the numpy array.
#

viewer = LavaVuModelViewer(model, background="white")

# determine the number of unique surfaces in the model from
# the input data and then calculate isosurfaces for this
unique = np.unique(strati.interpolator.get_value_constraints()[:, 3])
viewer.add_isosurface(strati, slices=unique, cmap="prism", paint_with=strati)

viewer.add_section(
    strati,
    axis="x",
    value=0.0,
    boundary_points=model.bounding_box,
    nsteps=np.array([30, 30, 30]),
    cmap="prism",
)
viewer.add_scalar_field(strati, cmap="prism")
viewer.add_model(cmap="tab20")

# Add the data addgrad/addvalue arguments are optional
viewer.add_data(strati, addgrad=True, addvalue=True, cmap="prism")
viewer.lv.rotate([-85.18760681152344, 42.93233871459961, 0.8641873002052307])
viewer.display()  # to add an interactive display
