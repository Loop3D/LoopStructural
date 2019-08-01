from FME.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from FME.supports.tet_mesh import TetMesh
from FME.modelling.geological_feature import GeologicalFeature, FaultedGeologicalFeature, GeologicalFeatureBuilder
from FME.visualisation.model_visualisation import LavaVuModelViewer
from FME.modelling.structural_frame import StructuralFrameBuilder
from FME.modelling.fault.fault_segment import FaultSegment
import numpy as np
import timeit
import lavavu

"""
This is a basic example showing how to use the Piecewise Linear Interpolator for orientation and
value data points. 
"""
start = timeit.default_timer()
boundary_points = np.zeros((2,3))

boundary_points[0,0] = -20
boundary_points[0,1] = -20
boundary_points[0,2] = -20
boundary_points[1,0] = 20
boundary_points[1,1] = 20
boundary_points[1,2] = 20
mesh = TetMesh()
mesh.setup_mesh(boundary_points, nstep=1, n_tetra=100000,)

interpolator = PLI(mesh)
feature_builder = GeologicalFeatureBuilder(interpolator,name='stratigraphy')

feature_builder.add_point([0,0,0],0)
feature_builder.add_point([0,0,1],-0.5)
feature_builder.add_strike_and_dip([0,0,0],90,0)
feature = feature_builder.build(solver='cg')


fault_frame_interpolator = PLI(mesh)
fault = StructuralFrameBuilder(interpolator=fault_frame_interpolator,mesh=mesh,name='FaultSegment1')

fault.add_point([2.5,.5,1.5],0.,itype='gz')
fault.add_point([2.5,-.5,1.5],1.,itype='gz')


fault.add_point([10,0,-5],1.,itype='gy')

for y in range(-15,15):
    fault.add_point([18.,y,18],0,itype='gy')
    fault.add_point([10,y,-5],1.,itype='gy')

    fault.add_strike_dip_and_value([18.,y,18],strike=180,dip=35,val=0,itype='gx')

ogw = 300
ogw /= mesh.n_elements
cgw = 500
cgw = cgw / mesh.n_elements
fault_frame = fault.build(
    solver='cg',
    guess=None,
   gxxgy=2*ogw,
   gxxgz=2*ogw,
   gyxgz=ogw,
   gxcg=cgw,
   gycg=cgw,
   gzcg=cgw,
   shape='rectangular',
    gx=True,
    gy=True,
    gz=True
)
#
fault = FaultSegment(fault_frame)
faulted_feature = FaultedGeologicalFeature(feature, fault)
viewer = LavaVuModelViewer()
viewer.plot_isosurface(faulted_feature.hw_feature,isovalue=0)
viewer.plot_isosurface(faulted_feature.fw_feature,isovalue=0)
mask = fault_frame.features[0].support.get_node_values() > 0
mask[mesh.elements] = np.any(mask[mesh.elements] == True, axis=1)[:, None]
viewer.plot_points(mesh.nodes[mask], "nodes", col="red")
viewer.plot_isosurface(fault_frame.features[0], isovalue=0, colour='blue')
viewer.lv.interactive()
