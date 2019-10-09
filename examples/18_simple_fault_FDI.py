from LoopStructural.interpolators.finite_difference_interpolator import FiniteDifferenceInterpolator as FDI
from LoopStructural.supports.structured_grid import StructuredGrid
from LoopStructural.modelling.features.geological_feature import GeologicalFeatureInterpolator
from LoopStructural.modelling.features.faulted_geological_feature import FaultedGeologicalFeature
from LoopStructural.visualisation.model_visualisation import LavaVuModelViewer
from LoopStructural.modelling.features.structural_frame import StructuralFrameBuilder
from LoopStructural.modelling.fault.fault_segment import FaultSegment
import numpy as np
import timeit

"""
This is a basic example showing how to use the Piecewise Linear Interpolator for orientation and
value data points. 
"""
solver = 'cg'
start = timeit.default_timer()
# boundary_points = np.zeros((2,3))
#
# boundary_points[0,0] = -20
# boundary_points[0,1] = -20
# boundary_points[0,2] = -20
# boundary_points[1,0] = 20
# boundary_points[1,1] = 20
# boundary_points[1,2] = 20
# mesh = TetMesh()
# mesh.setup_mesh(boundary_points, nstep=1, n_tetra=50000,)
grid = StructuredGrid(nsteps=(10,10,10),step_vector=np.array([2,2,2]))
interpolator = FDI(grid)
feature_builder = GeologicalFeatureInterpolator(interpolator, name='stratigraphy')

feature_builder.add_point([1,1,1],0)
feature_builder.add_point([2,2,1],-0.5)
feature_builder.add_strike_and_dip([1,1,1],90,0)
feature = feature_builder.build(
    solver=solver,
    tol=1e-30)


fault_frame_interpolator = FDI(grid)
fault = StructuralFrameBuilder(interpolator=fault_frame_interpolator,support=grid,name='FaultSegment1')

fault.add_point([2.5,.5,1.5],0.,itype='gz')
fault.add_point([2.5,0,1.5],1.,itype='gz')

# fault.add_strike_dip_and_value([18., y, 18], strike=180, dip=35, val=0, itype='gx')

# fault.add_point([10,0,5],1.,itype='gy')

for y in range(0,15):
    fault.add_point([12.,y,12],0,itype='gy')
    fault.add_point([10,y,0],1.,itype='gy')

    fault.add_strike_dip_and_value([12.,y,12],strike=180,dip=35,val=0,itype='gx')

ogw = 300
ogw /= grid.n_elements
cgw = 5000
fault_frame = fault.build(
    solver=solver,
   shape='rectangular',
)
#
fault = FaultSegment(fault_frame, displacement=4)
faulted_feature = FaultedGeologicalFeature(feature, fault)
viewer = LavaVuModelViewer()
# viewer.add_isosurface(faulted_feature, isovalue=0)
viewer.add_isosurface(feature)#, isovalue=0)

# mask = fault_frame.features[0].support.get_node_values() > 0
# mask[grid.elements] = np.any(mask[grid.elements] == True, axis=1)[:, None]
# viewer.plot_points(grid.nodes[mask], "nodes", col="red")
# viewer.add_isosurface(fault_frame.features[0], isovalue=0, colour='blue')
viewer.save('finite_difference_test.png')
