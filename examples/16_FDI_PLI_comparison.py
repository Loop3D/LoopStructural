# ### Import LoopStructural objects for interpolation
from LoopStructural.interpolators.finite_difference_interpolator import FiniteDifferenceInterpolator as FDI
from LoopStructural.supports.structured_grid import StructuredGrid
from LoopStructural.modelling.features.geological_feature import GeologicalFeatureInterpolator
from LoopStructural.visualisation.model_visualisation import LavaVuModelViewer
from LoopStructural.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from LoopStructural.supports.tet_mesh import TetMesh
# other necessary libraries
import numpy as np
import lavavu
import pyvista
import matplotlib.pyplot as plt

# ### Loop3D Forward Modelling Engine
# This is an example using LoopStructural to interpolate stratigraphy within a cube. It shows a very basic
# example of using LoopStructural.

# ### Defining Model Region
# Define the area to model represented by the domain of the tetrahedral mesh
grid = StructuredGrid(nsteps=(30,30,30),step_vector=np.array([.5,.5,.5]))
# ### Meshing
# Create a TetMesh object and build the mesh using tetgen. The number of tetrahedron can be specified
# by adding the n_tetra flag.


# ### GeologicalFeatureInterpolator
# LoopStructural uses an object oriented design so

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
    dzz= .2,
    cgw=cgw)


boundary_points = np.zeros((2,3))
boundary_points[0,0] = 0
boundary_points[0,1] = 0
boundary_points[0,2] = 0
boundary_points[1,0] = 30*0.5
boundary_points[1,1] = 30*0.5
boundary_points[1,2] = 30*0.5

# ### Meshing
# Create a TetMesh object and build the mesh using tetgen. The number of tetrahedron can be specified
# by adding the n_tetra flag.
mesh = TetMesh()
mesh.setup_mesh(boundary_points, n_tetra=30*30*30)

# ### GeologicalFeatureInterpolator
# LoopStructural uses an object oriented design so

plit_interpolator = PLI(mesh)
pli_feature_builder = GeologicalFeatureInterpolator(plit_interpolator, name='stratigraphy2')
for x in range(10):
    pli_feature_builder.add_point(np.array([[x,1., 1.]]),0)
    pli_feature_builder.add_strike_and_dip(np.array([[x, 1., 1.]]), 90,80)

    # feature_builder.add_point(np.array([[x,1., 2.]]),0)
    #
    # feature_builder.add_point(np.array([[x, 5., 5.]]), 0.0)
    # feature_builder.add_point(np.array([[x, 5., 1.]]), 1.0)

    pli_feature_builder.add_strike_and_dip(np.array([[x, 8., 8.]]), 90,10)
    pli_feature_builder.add_strike_and_dip(np.array([[x, 14., 8.]]), 90,30)
    pli_feature_builder.add_point(np.array([[x,14., 8.]]),0)

# feature_builder.add_strike_and_dip([0,0,0],90,50)
cgw = 1000
# cgw /= mesh.n_elements
feature2 = pli_feature_builder.build(
    solver='cg',
    cgw=cgw)
# feature.evaluate_value([[0,0,0]])
# Export the input data and the model domain to a text file so it can be imported into gocad

# ### Visualisation using LavaVu
# LoopStructural uses Lavavu for visualising objects. The LavaVuModelViewer class interfaces between lavavu and LoopStructural
# kwargs can be passed from the wrapper functions to the lavavu objects.

viewer = LavaVuModelViewer(background="white")
viewer.plot_isosurface(
    feature,
    # slices=[0], #specify multiple isosurfaces
    isovalue=0, # a single isosurface
    # nslices=10 #the number of evenly space isosurfaces
)
viewer.plot_isosurface(
    feature2,
    # slices=[0], #specify multiple isosurfaces
    isovalue=0, # a single isosurface
    # nslices=10 #the number of evenly space isosurfaces
    colour='blue'
)
viewer.plot_vector_data(
    feature_builder.interpolator.get_gradient_control()[:,:3],
    feature_builder.interpolator.get_gradient_control()[:,3:],
    "grad" # object name
)
viewer.plot_value_data(
    feature_builder.interpolator.get_control_points()[:,:3],
    feature_builder.interpolator.get_control_points()[:,3:],
    "value",
    pointsize=10,
    colourmap=lavavu.matplotlib_colourmap("Greys"))



viewer.lv.interactive()