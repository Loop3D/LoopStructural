# ### Import LoopStructural objects for interpolation
from LoopStructural.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from LoopStructural.supports.tet_mesh import TetMesh
from LoopStructural.modelling.features.geological_feature import GeologicalFeatureInterpolator
from LoopStructural.visualisation.model_visualisation import LavaVuModelViewer
# other necessary libraries
import numpy as np
import lavavu

# ### Loop3D Forward Modelling Engine
# This is an example using LoopStructural to interpolate stratigraphy within a cube. It shows a very basic
# example of using LoopStructural.

# ### Defining Model Region
# Define the area to model represented by the domain of the tetrahedral mesh
boundary_points = np.zeros((2,3))
boundary_points[0,0] = -1
boundary_points[0,1] = -1
boundary_points[0,2] = -1
boundary_points[1,0] = 1
boundary_points[1,1] = 1
boundary_points[1,2] = 1

# ### Meshing
# Create a TetMesh object and build the mesh using tetgen. The number of tetrahedron can be specified
# by adding the n_tetra flag.
mesh = TetMesh()
mesh.setup_mesh(boundary_points, n_tetra=10000,)

# ### GeologicalFeatureInterpolator
# LoopStructural uses an object oriented design so

interpolator = PLI(mesh)
feature_builder = GeologicalFeatureInterpolator(interpolator, name='stratigraphy')

feature_builder.add_point([0,0,0],0)
feature_builder.add_point([0.5,0,0],1)
# feature_builder.add_point([-.9,0,0],.8)
# feature_builder.add_strike_and_dip([0.4,0,0],70,50)
#
feature_builder.add_strike_and_dip([0,0,0],90,50)
cgw = 6000
# cgw /= mesh.n_elements
feature = feature_builder.build(
    solver='lu',
    cgw=cgw)
feature_builder2 = GeologicalFeatureInterpolator(interpolator, name='stratigraphy2')

# feature_builder2.add_point([0,0,0],0)
feature_builder2.add_point([0.9,0,0],1)
# feature_builder.add_point([-.9,0,0],.8)
# feature_builder.add_strike_and_dip([0.4,0,0],70,50)
#
# feature_builder.add_strike_and_dip([0,0,0],90,50)
mask = mesh.nodes[:,0] < 0
idc = np.arange(0,mesh.n_nodes)[mask]
v = feature.support.get_node_values()[mask] #interpolators[0].

feature_builder2.interpolator.add_equality_constraints(idc,v)
cgw = 6000
# cgw /= mesh.n_elements
feature2 = feature_builder2.build(
    solver='lueq',
    cgw=cgw)

# Export the input data and the model domain to a text file so it can be imported into gocad
np.savetxt("01_gradient.txt", feature_builder.interpolator.get_gradient_control())
np.savetxt("01_value.txt", feature_builder.interpolator.get_control_points())
np.savetxt("01_box_coords.txt", mesh.points)


# ### Visualisation using LavaVu
# LoopStructural uses Lavavu for visualising objects. The LavaVuModelViewer class interfaces between lavavu and LoopStructural
# kwargs can be passed from the wrapper functions to the lavavu objects.

viewer = LavaVuModelViewer(background="white")
viewer.add_isosurface(
    feature,
    colour='green',
    # slices=[0], #specify multiple isosurfaces
    isovalue=0, # a single isosurface
    # nslices=10 #the number of evenly space isosurfaces
)
viewer.add_isosurface(
    feature2,
    colour='blue',
    # slices=[0], #specify multiple isosurfaces
    isovalue=0, # a single isosurface
    # nslices=10 #the number of evenly space isosurfaces
)
viewer.add_vector_data(
    feature_builder.interpolator.get_gradient_control()[:,:3],
    feature_builder.interpolator.get_gradient_control()[:,3:],
    "grad" # object name
)
viewer.add_value_data(
    feature_builder.interpolator.get_control_points()[:,:3],
    feature_builder.interpolator.get_control_points()[:,3:],
    "value",
    pointsize=10,
    colourmap=lavavu.matplotlib_colourmap("Greys"))
viewer.lv.interactive()