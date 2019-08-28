# ### Import FME objects for interpolation
from FME.interpolators.finite_difference_interpolator import FiniteDifferenceInterpolator as FDI
from FME.supports.structured_grid import StructuredGrid
from FME.modelling.features.geological_feature import GeologicalFeatureInterpolator
from FME.visualisation.model_visualisation import LavaVuModelViewer
from FME.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from FME.supports.tet_mesh import TetMesh
# other necessary libraries
import numpy as np
import lavavu
import pyvista
import matplotlib.pyplot as plt

# ### Loop3D Forward Modelling Engine
# This is an example using FME to interpolate stratigraphy within a cube. It shows a very basic
# example of using FME.

# ### Defining Model Region
# Define the area to model represented by the domain of the tetrahedral mesh
grid = StructuredGrid(nsteps=(30,30,30),step_vector=np.array([.5,.5,.5]))
# ### Meshing
# Create a TetMesh object and build the mesh using tetgen. The number of tetrahedron can be specified
# by adding the n_tetra flag.


# ### GeologicalFeatureInterpolator
# FME uses an object oriented design so

interpolator = FDI(grid)
feature_builder = GeologicalFeatureInterpolator(interpolator, name='stratigraphy')
for x in range(10):
    feature_builder.add_point(np.array([[x,1., 1.]]),0)
    feature_builder.add_strike_and_dip(np.array([[x, 1., 1.]]), 90,80)

    # feature_builder.add_point(np.array([[x,1., 2.]]),0)
    #
    # feature_builder.add_point(np.array([[x, 5., 5.]]), 0.0)
    # feature_builder.add_point(np.array([[x, 5., 1.]]), 1.0)

    feature_builder.add_strike_and_dip(np.array([[x, 8., 8.]]), 90,10)
    feature_builder.add_strike_and_dip(np.array([[x, 14., 8.]]), 90,30)
    feature_builder.add_point(np.array([[x,14., 8.]]),0)

# feature_builder.add_point([-.9,0,0],.8)
# feature_builder.add_strike_and_dip([0.4,0,0],70,50)
feature_builder.interpolator.add_gradient_orthogonal_constraint(np.arange(grid.n_elements), np.zeros((grid.n_elements, 3)))
# feature_builder.add_strike_and_dip([0,0,0],90,50)
cgw = 6000
# cgw /= mesh.n_elements
feature = feature_builder.build(
    solver='cg',
    dxy=.1,
    dyz=.1,
    dxz= .1,
    dxx= .2,
    dyy= .2,
    dzz= .2,)

