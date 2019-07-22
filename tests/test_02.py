from FME.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from FME.supports.tet_mesh import TetMesh
from FME.modelling.geological_feature import GeologicalFeature
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
mesh.setup_mesh(boundary_points, nstep=1, n_tetra=10000,)

interpolator = PLI(mesh)
interpolator.add_point([0,0,0],0)
a = np.zeros((3,3,3))
interpolator.add_point([0.5,0,0],0.5)
# interpolator.add_strike_and_dip([0,0,0],90,10)
interpolator.setup_interpolator()
interpolator.solve_system(solver='lsmr')

support = interpolator.get_support()

feature = GeologicalFeature(0, 'Name', support)

