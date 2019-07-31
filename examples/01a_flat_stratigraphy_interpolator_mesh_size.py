from FME.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from FME.supports.tet_mesh import TetMesh
from FME.modelling.geological_feature import GeologicalFeature, GeologicalFeatureBuilder
from FME.visualisation.model_visualisation import LavaVuModelViewer
import numpy as np
import lavavu
"""
This is a basic example showing how to use the Piecewise Linear Interpolator for orientation and
value data points. 
"""
boundary_points = np.zeros((2,3))

boundary_points[0,0] = -1
boundary_points[0,1] = -1
boundary_points[0,2] = -1
boundary_points[1,0] = 1
boundary_points[1,1] = 1
boundary_points[1,2] = 1
mesh = TetMesh()
mesh.setup_mesh(boundary_points, nstep=1, n_tetra=100000,)
np.savetxt("01_box_coords.txt",mesh.points)
interpolator = PLI(mesh)
feature_builder = GeologicalFeatureBuilder(interpolator,name='stratigraphy')

feature_builder.add_point([.90,.90,.90],0)
feature_builder.add_point([-.50,0,0],1)
feature_builder.add_point([-.90,0,0],.8)

feature_builder.add_strike_and_dip([0,0,0],90,40)
feature = feature_builder.build(
    solver='chol',
    cgw=0.2)

np.savetxt("01_gradient.txt",feature_builder.interpolator.get_gradient_control())
np.savetxt("01_value.txt",feature_builder.interpolator.get_control_points())
cmap = lavavu.matplotlib_colourmap("Greys")
#
# viewer = LavaVuModelViewer(background="white")
# viewer.plot_isosurface(feature, slices=[0, 1., .8])
#
# viewer.plot_vector_data(feature_builder.interpolator.get_gradient_control()[:,:3],feature_builder.interpolator.get_gradient_control()[:,3:],"grad")
# viewer.plot_value_data(feature_builder.interpolator.get_control_points()[:,:3],feature_builder.interpolator.get_control_points()[:,3:],"value",pointsize=10,colourmap=cmap)
# viewer.lv.interactive()