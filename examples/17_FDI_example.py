# ### Import LoopStructural objects for interpolation
# other necessary libraries
import numpy as np

from LoopStructural.interpolators.finite_difference_interpolator import \
    FiniteDifferenceInterpolator as FDI
from LoopStructural.modelling.features.geological_feature import \
    GeologicalFeatureInterpolator
from LoopStructural.supports.structured_grid import StructuredGrid
from LoopStructural.visualisation.model_visualisation import LavaVuModelViewer

# import pyvista

# ### Loop3D Forward Modelling Engine
# This is an example using LoopStructural to interpolate stratigraphy within
# a cube. It shows a very basic
# example of using LoopStructural.

# ### Defining Model Region
# Define the area to model represented by the domain of the tetrahedral mesh
grid = StructuredGrid(nsteps=(30, 30, 30), step_vector=np.array([1, 1, 1]))
# ### Meshing
# Create a TetMesh object and build the mesh using tetgen. The number of
# tetrahedron can be specified
# by adding the n_tetra flag.
boundary_points = np.zeros((2, 3))
boundary_points[1, :] = 30
# ### GeologicalFeatureInterpolator
# LoopStructural uses an object oriented design so

interpolator = FDI(grid)
feature_builder = GeologicalFeatureInterpolator(interpolator,
                                                name='stratigraphy')
for x in range(1, 10):
    feature_builder.add_point(np.array([[x, 1., 1.]]), 0)
    feature_builder.add_strike_and_dip(np.array([[x, 1., 1.]]), 90, 80)

    # feature_builder.add_point(np.array([[x,1., 2.]]),0)
    #
    # feature_builder.add_point(np.array([[x, 5., 5.]]), 0.0)
    # feature_builder.add_point(np.array([[x, 5., 1.]]), 1.0)

    # feature_builder.add_strike_and_dip(np.array([[x, 8., 8.]]), 90,10)
    feature_builder.add_strike_and_dip(np.array([[x, 14., 8.]]), 90, 30)
    feature_builder.add_point(np.array([[x, 14., 8.]]), 0)

# feature_builder.add_point([-.9,0,0],.8)
# feature_builder.add_strike_and_dip([0.4,0,0],70,50)
feature_builder.interpolator.add_gradient_orthogonal_constraint(
    np.arange(grid.n_elements), np.zeros((grid.n_elements, 3)))
# feature_builder.add_strike_and_dip([0,0,0],90,50)
cgw = 6000
# cgw /= mesh.n_elements
feature = feature_builder.build(
    solver='cg',
    # maxiter=500,
    tol=1e-10
)
print(feature.get_node_values())
viewer = LavaVuModelViewer(background="white")
viewer.add_isosurface(
    feature,
    # slices=[0], #specify multiple isosurfaces
    # isovalue=0, # a single isosurface
    # nslices=10 #the number of evenly space isosurfaces
)
viewer.add_data(feature)
# viewer.add_section(feature, axis='y', boundary_points=boundary_points,
# nsteps=np.array([10,10,10]))
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
