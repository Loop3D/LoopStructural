from FME.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from FME.supports.tet_mesh import TetMesh
from FME.modelling.geological_feature import GeologicalFeature, GeologicalFeatureBuilder
from FME.visualisation.model_visualisation import LavaVuModelViewer
import numpy as np

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
mesh.setup_mesh(boundary_points, nstep=1, n_tetra=1000000,)

interpolator = PLI(mesh)
feature_builder = GeologicalFeatureBuilder(interpolator,name='stratigraphy')

feature_builder.add_point([1,1,1],0)
feature_builder.add_point([-0.5,0,0],1)
feature_builder.add_point([-1,0,0],.8)

feature_builder.add_strike_and_dip([0,0,0],90,40)
feature = feature_builder.build(solver='cg')



viewer = LavaVuModelViewer()
viewer.plot_isosurface(feature,isovalue=0)

viewer.lv.interactive()