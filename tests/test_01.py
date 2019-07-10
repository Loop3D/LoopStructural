from FME.pli_interpolator import PiecewiseLinearInterpolator as PLI
from FME.tet_mesh import TetMesh
import numpy as np
boundary_points = np.zeros((2,3))
boundary_points[0,0] = -1
boundary_points[0,1] = -1
boundary_points[0,2] = -1
boundary_points[1,0] = 1
boundary_points[1,1] = 1
boundary_points[1,2] = 1
mesh = TetMesh()
mesh.setup_mesh(boundary_points,nstep=1,n_tetra=80000,)
interpolator = PLI(mesh)
#interpolator.add_constant_gradient()
interpolator.add_point([0,0,0],0)
interpolator.add_point([0.5,0,0],0.5)
interpolator.setup_interpolator()
interpolator.solve_system(solver='lsmr')
