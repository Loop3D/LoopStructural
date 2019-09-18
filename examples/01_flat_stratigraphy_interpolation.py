# ### Import FME objects for interpolation
from FME.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from FME.supports.tet_mesh import TetMesh
from FME.modelling.features.geological_feature import GeologicalFeatureInterpolator
from FME.visualisation.model_visualisation import LavaVuModelViewer
# other necessary libraries
import numpy as np
import lavavu

# ### Loop3D Forward Modelling Engine
# This is an example using FME to interpolate stratigraphy within a cube. It shows a very basic
# example of using FME.

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
# FME uses an abstract class representation of geological objects. A **GeologicalFeature** is
# an abstract class that represents a geological object within the model. A **GeologicalFeature**
# can be evaluated for the gradient or scalar value at any location within the model.
# It can be simply a wrapper for a scalar field, or could be a piecewise combination of scalar
# fields.
# The **GeologicalFeature** needs to be built using a builder class, this class manages the data,
# interpolator and any geological structures. The **GeologicalFeatureInterpolator** provides a way
# of linking any FME interpolator to a **GeologicalFeature**. The GeologicalFeatureInterpolator
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


interpolator = PLI(mesh)
feature_builder = GeologicalFeatureInterpolator(interpolator, name='stratigraphy')

feature_builder.add_point([0,0,0],0)
feature_builder.add_point([0.5,0,0],1)
# feature_builder.add_point([-.9,0,0],.8)
# feature_builder.add_strike_and_dip([0.4,0,0],70,50)
#
# feature_builder.add_strike_and_dip([0,0,0],90,50)
cgw = 6000
# cgw /= mesh.n_elements
feature = feature_builder.build(
    solver='lu')

# Export the input data and the model domain to a text file so it can be imported into gocad
np.savetxt("01_gradient.txt", feature_builder.interpolator.get_gradient_control())
np.savetxt("01_value.txt", feature_builder.interpolator.get_control_points())
np.savetxt("01_box_coords.txt", mesh.points)


# ### Visualisation using LavaVu
# FME uses Lavavu for visualising objects. The LavaVuModelViewer class interfaces between lavavu and FME
# kwargs can be passed from the wrapper functions to the lavavu objects.

viewer = LavaVuModelViewer(background="white")
viewer.plot_isosurface(
    feature,
    slices=[0], #specify multiple isosurfaces
    # isovalue=0, # a single isosurface
    # nslices=10 #the number of evenly space isosurfaces
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
