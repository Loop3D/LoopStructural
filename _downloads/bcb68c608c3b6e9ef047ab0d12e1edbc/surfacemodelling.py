"""

test
====
"""


######################################################################
# Implicit surface modelling
# ==========================
# 
# This tutorial will demonstrate how to create an implicit surface
# representation of surfaces from a combination of orientation and
# location observations.
# 
# Implicit surface representation involves finding an unknown function
# where :math:`f(x,y,z)` matches observations of the surface geometry. We
# generate a scalar field where the scalar value is the distance away from
# a reference horizon. The reference horizon is arbritary and can either
# be:
# 
# -  a single geological surface where the scalar field would represent
#    the signed distance away from this surface. (above the surface
#    positive and below negative)
# -  Where multiple conformable horizons are observed the same scalar
#    field can be used to represent these surfaces and the thickness of
#    the layers is used to determine the relative scalar value for each
#    surface
# 
# This tutorial will demonstrate both of these approaches for modelling a
# number of horizons picked from seismic data sets, by following the next
# steps: 1. Creation of a geological model, which includes: \*
# Presentation and visualization of the data \* Addition of a geological
# feature, which in this case is the stratigraphy of the model. 2.
# Visualization of the scalar field.
# 

from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer 

from LoopStructural.datasets import load_claudius #demo data 

import pandas as pd
import glob
import numpy as np
# %matplotlib inline


######################################################################
# The data for this example can be imported from the example datasets
# module in loopstructural.
# 

data, bb = load_claudius()

data['val'].unique()


######################################################################
# GeologicalModel
# ~~~~~~~~~~~~~~~
# 
# To create a ``GeologicalModel`` the origin (lower left) and max extent
# (upper right) of the model area need to be specified. In this example
# these are provided in the bb array.
# 

model = GeologicalModel(bb[0,:],bb[1,:])


######################################################################
# A pandas dataframe with appropriate columns can be used to link the data
# to the geological model. \* ``X`` is the x coordinate \* ``Y`` is the y
# coordinate \* ``Z`` is the z coordinate \* ``type`` is a name to link
# the data to a model object \* ``val`` is the value of the scalar field
# which represents the distance from a reference horizon. It is comparable
# to the relative thickness \* ``nx`` is the x component of the normal
# vector to the surface gradient \* ``ny`` is the y component of the
# normal vector to the surface gradient \* ``nz`` is the z component of
# the normal vector to the surface gradeint \* ``strike`` is the strike
# angle \* ``dip`` is the dip angle
# 
# Having a look at the data for this example by looking at the top of the
# dataframe and then using a 3D plot
# 

data['type'].unique()

viewer = LavaVuModelViewer()
viewer.add_value_data(data[~np.isnan(data['val'])][['X','Y','Z']],data[~np.isnan(data['val'])]['val'],name='value points')
viewer.add_vector_data(data[~np.isnan(data['nx'])][['X','Y','Z']],
                       data[~np.isnan(data['nx'])][['nx','ny','nz']],name='orientation points')

viewer.interactive()


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
# as: \* Foliations \* Faults \* Unconformities \* Folded foliations \*
# Structural Frames
# 
# In this example we will only add a foliation using the function
# 
# .. code:: python
# 
#    model.create_and_add_foliation(name)
# 
# where name is the name in the ``type`` field, other parameters we
# specified are the: \* ``interpolatortype`` - we can either use a
# PiecewiseLinearInterpolator ``PLI``, a FiniteDifferenceInterpolator
# ``FDI`` or a radial basis interpolator ``surfe`` \* ``nelements - int``
# is the how many elements are used to discretize the resulting solution
# \* ``buffer - float`` buffer percentage around the model area \*
# ``solver`` - the algorithm to solve the least squares problem e.g.
# ``lu`` for lower upper decomposition, ``cg`` for conjugate gradient,
# ``pyamg`` for an algorithmic multigrid solver \* ``damp - bool`` -
# whether to add a small number to the diagonal of the interpolation
# matrix for discrete interpolators - this can help speed up the solver
# and makes the solution more stable for some interpolators
# 

vals = [0,60,250,330,600]
strat_column = {'strati':{}}
for i in range(len(vals)-1):
    strat_column['strati']['unit_{}'.format(i)] = {'min':vals[i],'max':vals[i+1],'id':i}
    


print(strat_column)

model.set_stratigraphic_column(strat_column)

strati = model.create_and_add_foliation("strati",
                                           interpolatortype="FDI", # try changing this to 'PLI'
                                           nelements=1e4, # try changing between 1e3 and 5e4
                                           buffer=0.3,
                                           solver='pyamg',
                                           damp=True
                                          )


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
# -  ``voxet`` dict with ``bounding_box=boundary_points`` and
#    ``nsteps = (nx,ny,nz)``
# -  other kwargs for passing directly to lavavu
# 
# Alternatively the scalarfields can be displayed on a rectangular cuboid.
# 
# .. code:: python
# 
#    viewer.add_scalar_field(boundary_points,dimensions,**kwargs)
# 
# Where ``boundary_points`` is a numpy array
# ``[[minx,miny,minz],[maxx,maxy,maxz]]`` and ``dimensions`` corresponds
# to the number of samples along each axis.
# 
# Other possible kwargs are:
# 
# -  ``paint_with`` the geological feature to colour the box with
# -  ``colour`` a single colour to colour the surfaces with
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
#    viewer.add_vector_field(feature, locations, **kwargs)
# 
# Where ``locations`` are an array of points to evaluate the gradient at,
# for example the barycentric coordinates. It is recommended to visualise
# the vectorfield at a lower resolution than the mesh otherwise it can be
# difficult to see the vectors. You can use numpy stepping along the
# array: ``locations = mesh.barycentre[::20,:]`` which will sample every
# 20th sample in the numpy array.
# 

viewer = LavaVuModelViewer(model,background="white")

# determine the number of unique surfaces in the model from 
# the input data and then calculate isosurfaces for this
unique = np.unique(strati['feature'].interpolator.get_value_constraints()[:,3])
viewer.add_isosurface(model.features[0], 
                       slices=unique,  
                       cmap='prism',
                      paint_with=model.features[0],
                     voxet=model.voxet())

viewer.add_section(model.features[0],
                   axis='x',
                   value=0.,
                   boundary_points=model.bounding_box, 
                   nsteps=np.array([30,30,30]),
                  cmap='prism')
viewer.add_scalar_field(model.features[0],
                     cmap='prism')
viewer.add_model()

# Add the data addgrad/addvalue arguments are optional
viewer.add_data(model.features[0],addgrad=True,addvalue=True, cmap='prism')
viewer.lv.rotate([-85.18760681152344, 42.93233871459961, 0.8641873002052307])
viewer.interactive()# to add an interactive display


np.unique(model.evaluate_model(model.regular_grid()))

model.feature_name_index.get('strati')

strat_column

