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
boundary_points = np.zeros((2,3))
#
boundary_points[0,0] = 2
boundary_points[0,1] = 2
boundary_points[0,2] = 2
boundary_points[1,0] = 16
boundary_points[1,1] = 16
boundary_points[1,2] = 16
# mesh = TetMesh()
# mesh.setup_mesh(boundary_points, nstep=1, n_tetra=50000,)
grid = StructuredGrid(nsteps=(40,40,20),step_vector=np.array([.5,.5,1]))
interpolator = FDI(grid)
feature_builder = GeologicalFeatureInterpolator(interpolator, name='stratigraphy')
for y in range(0,20,4):

    feature_builder.add_point([10,y+.5,10],0)
    # feature_builder.add_point([2,2,1],-0.5)
    feature_builder.add_strike_and_dip([10,y+.5,10],90,0)
    feature_builder.add_strike_and_dip([15,y+.5,10],0,40)
    feature_builder.add_point([15,y+.5,10],0)

interpolation_weights = {'dxy': 1,
                          'dyz': 1,
                          'dxz': 1,
                          'dxx': 1,
                          'dyy': 1,
                          'dzz': 1,
                          'gpw': 5.,
                          'dy': 1.,
                          'dz': 1.}

feature_builder.interpolator.set_interpolation_weights(interpolation_weights)
feature = feature_builder.build(
    solver='cg',
    tol=1e-10,
    # maxiter=1000
)


# fault_frame_interpolator = FDI(grid)
# fault = StructuralFrameBuilder(interpolator=fault_frame_interpolator,support=grid,name='FaultSegment1')
#
# fault.add_point([2.5,.5,1.5],0.,itype='gz')
# fault.add_point([2.5,0,1.5],1.,itype='gz')
#
# # fault.add_strike_dip_and_value([18., y, 18], strike=180, dip=35, val=0, itype='gx')
#
# # fault.add_point([10,0,5],1.,itype='gy')
#
# for y in range(0,15):
#     fault.add_point([12.,y,12],0,itype='gy')
#     fault.add_point([10,y,0],1.,itype='gy')
#
#     fault.add_strike_dip_and_value([12.,y,12],strike=180,dip=35,val=0,itype='gx')
#
# ogw = 300
# ogw /= grid.n_elements
# cgw = 5000
# fault_frame = fault.build(
#     solver=solver,
#    shape='rectangular',
# )
# #
# fault = FaultSegment(fault_frame, displacement=4)
# faulted_feature = FaultedGeologicalFeature(feature, fault)
viewer = LavaVuModelViewer(background='white')
viewer.add_isosurface(feature,
                      nslices=10,
                      # isovalue=0,
                                            # voxet={'bounding_box':boundary_points,'nsteps':(50,50,25)}
                      )
# viewer.add_isosurface(fault_frame[0],colour='black')
viewer.add_data(feature)
# viewer.add_isosurface(feature)#, isovalue=0)

# mask = fault_frame.features[0].support.get_node_values() > 0
# mask[grid.elements] = np.any(mask[grid.elements] == True, axis=1)[:, None]
# viewer.plot_points(grid.nodes[mask], "nodes", col="red")
# viewer.add_isosurface(fault_frame.features[0], isovalue=0, colour='blue')
# viewer.save('finite_difference_test.png')
viewer.interactive()