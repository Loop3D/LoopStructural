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
scale = 1#0000
boundary_points = np.zeros((2,3))
boundary_points[0,0] = -1
boundary_points[0,1] = -1
boundary_points[0,2] = -1
boundary_points[1,0] = 1
boundary_points[1,1] = 1
boundary_points[1,2] = 1
boundary_points*=scale
# ### Meshing
# Create a TetMesh object and build the mesh using tetgen. The number of tetrahedron can be specified
# by adding the n_tetra flag.
mesh = TetMesh()
mesh.setup_mesh(boundary_points, n_tetra=50000,)

# ### GeologicalFeatureInterpolator
# LoopStructural uses an abstract class representation of geological objects. A **GeologicalFeature** is
# an abstract class that represents a geological object within the model. A **GeologicalFeature**
# can be evaluated for the gradient or scalar value at any location within the model.
# It can be simply a wrapper for a scalar field, or could be a piecewise combination of scalar
# fields.
# The **GeologicalFeature** needs to be built using a builder class, this class manages the data,
# interpolator and any geological structures. The **GeologicalFeatureInterpolator** provides a way
# of linking any LoopStructural interpolator to a **GeologicalFeature**. The GeologicalFeatureInterpolator
# needs to be given an interpolator object which should be a daughter class of a
# GeologicalInterpolator and a name to give the feature that is being built. The name is simply used
# as an identifier for the interpolated properties when exported into file formats or for visualisation.
# There are currently has two interpolator options available , a piecewise linear interpolator using
# a tetrahedral mesh where the roughness of the interpolated surfaces is minimised.
# A finite difference interpolator where the gaussian curvature is minimised.
# The GeologicalFeatureInterpolator can be given value constraints **add_point(position,value)** or
# gradient constraints **add_strike_and_dip**(pos,strike,dip)
# The resulting feature can be interpolated by calling **build**(kwargs) where the kwargs are
# passed to the interpolator you have chosen e.g. which solver to use, what weights etc..
# mesh.regions['test'] = mesh.nodes[:,1]>-.5

interpolator = PLI(mesh)
feature_builder = GeologicalFeatureInterpolator(interpolator, name='stratigraphy')

feature_builder.add_strike_and_dip([0,0,0],0,90)
feature_builder.add_point([0.1,0,0],0)
feature_builder.add_point([0.2,0,0],1)

# for y in np.arange(-0.5,.5,.10):
#     feature_builder.add_point([-0.5*scale,y*scale,0],0)
#     feature_builder.add_point([0.5*scale,y*scale,0],1)
# # feature_builder.add_point([0.5,0,0],1)
# for x in np.arange(-0.5,.5,.10):
#
# # feature_builder.add_point([-.9,0,0],.8)
#     feature_builder.add_strike_and_dip([x*scale,-0.5*scale,0],90,90)
#     feature_builder.add_strike_and_dip([x*scale,0.5*scale,0],90,90)

# feature_builder.add_strike_and_dip([0,0,0],90,50)
# cgw /= mesh.n_elements
feature = feature_builder.build(
    solver='cg',
    tol=1e-10, #need to add tolerance constraint to cg if the scalar values are low
    cgw=0.1)

# Export the input data and the model domain to a text file so it can be imported into gocad
# np.savetxt("01_gradient.txt", feature_builder.interpolator.get_gradient_control())
# np.savetxt("01_value.txt", feature_builder.interpolator.get_control_points())
# np.savetxt("01_box_coords.txt", mesh.points)


# ### Visualisation using LavaVu
# LoopStructural uses Lavavu for visualising objects. The LavaVuModelViewer class interfaces between lavavu and LoopStructural
# kwargs can be passed from the wrapper functions to the lavavu objects.

viewer = LavaVuModelViewer(background="white")
viewer.add_isosurface(
    feature,
    # slices=[0], #specify multiple isosurfaces
    # isovalue=0, # a single isosurface
    # nslices=10 #the number of evenly space isosurfaces
    )
viewer.add_data(feature)
viewer.add_section(feature, axis='y', boundary_points=boundary_points,nsteps=np.array([10,10,10]))
# viewer.add_vector_data(
#     feature_builder.interpolator.get_gradient_constraints()[:, :3],
#     feature_builder.interpolator.get_gradient_constraints()[:, 3:],
#     "grad" # object name
# )
# viewer.add_value_data(
#     feature_builder.interpolator.get_value_constraints()[:, :3],
#     feature_builder.interpolator.get_value_constraints()[:, 3:],
#     "value",
#     pointsize=10,
#     colourmap=lavavu.matplotlib_colourmap("Greys"))
# viewer.add_scalar_field(boundary_points,(38,55,30),
#                       'box',
#                      paint_with=feature,
#                      cmap='prism')
viewer.lv.rotate([-85.18760681152344, 42.93233871459961, 0.8641873002052307])
# viewer.save('01_flat_stratigraphy.png',transparent=True)
viewer.interactive()

#print(feature.support.get_node_values()[np.isnan(feature.support.get_node_values())==True])