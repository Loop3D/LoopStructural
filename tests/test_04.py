from FME.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from FME.supports.tet_mesh import TetMesh
from FME.modelling.geological_feature import GeologicalFeature, FaultedGeologicalFeature, GeologicalFeatureBuilder
from FME.visualisation.model_visualisation import LavaVuModelViewer
from FME.modelling.structural_frame import StructuralFrameBuilder, StructuralFrame
from FME.modelling.fault.fault_segment import FaultSegment
import numpy as np

boundary_points = np.zeros((2,3))

boundary_points[0,0] = -20
boundary_points[0,1] = -20
boundary_points[0,2] = -20
boundary_points[1,0] = 20
boundary_points[1,1] = 20
boundary_points[1,2] = 20
mesh = TetMesh()
mesh.setup_mesh(boundary_points, nstep=1, n_tetra=10000,)
interpolator = PLI(mesh)
stratigraphy_builder = GeologicalFeatureBuilder(interpolator=interpolator,name='stratigraphy')


stratigraphy_builder.add_point([6.1,0.1,1.1],0.)


stratigraphy_builder.add_strike_and_dip([1,1,1],90.,0.)
stratigraphy = stratigraphy_builder.build(solver='lu',cgw=.1)

floor = -6
roof = 6
fault_interpolator = PLI(mesh)
fault = StructuralFrameBuilder(interpolator=fault_interpolator,mesh=mesh,name='FaultSegment1')
for y in range(-20,20,1):
    fault.add_strike_dip_and_value([-15.17,y,floor],strike=0.,dip=0.,val=0,itype='gx')
    fault.add_strike_dip_and_value([-6.17,y,floor],strike=0.,dip=0.,val=0,itype='gx')
    fault.add_strike_dip_and_value([-4.,y,floor],strike=0.,dip=45.,val=0,itype='gx')
    fault.add_strike_dip_and_value([6.17,y,roof],strike=0.,dip=0.,val=0,itype='gx')
    fault.add_strike_dip_and_value([4.,y,roof],strike=0.,dip=0.,val=0,itype='gx')
    fault.add_strike_dip_and_value([8.17,y,roof],strike=0.,dip=0.,val=0,itype='gx')
    fault.add_strike_dip_and_value([15.17,y,roof],strike=0.,dip=0.,val=0,itype='gx')


fault.add_point([2.5,19,1.5],0.,itype='gz')
fault.add_point([2.5,-19,1.5],1.,itype='gz')
for y in range(-20,20,1):
    fault.add_point([-10.5,y,floor],0.,itype='gy')
    fault.add_point([11.56,y,roof],1.,itype='gy')


fault_frame = fault.build(solver='lu', gxxgy=0.1, gxxgz=1, gyxgz=0.05, gycg=5, gzcg=0.1)

fault_operator = FaultSegment(fault_frame)

fault_frame2_interpolator = PLI(mesh)
fault2 = StructuralFrameBuilder(interpolator=fault_frame2_interpolator,mesh=mesh,name='FaultSegment2')
#
fault2.add_point([2.5, .5, 1.5], 0., itype='gz')
fault2.add_point([2.5, -.5, 1.5], 1., itype='gz')
#
#
fault2.add_point([10, 0,-5],1.,itype='gy')
#
for y in range(-15,15):
    fault2.add_point([18.,y,18],0,itype='gy')
    fault2.add_point([10,y,-5],1.,itype='gy')
#
    fault2.add_strike_dip_and_value([18.,y,18],strike=0,dip=-55,val=0,itype='gx')
#
ogw = 300
ogw /= mesh.n_elements
cgw = 500
cgw = cgw / mesh.n_elements
fault_frame2 = fault2.build(
    solver='lu',
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
# #
fault2 = FaultSegment(fault_frame2)
faulted_frame = []
for f in fault_frame.features:
    faulted_frame.append(FaultedGeologicalFeature(f, fault2))
#
structural_frame2 = StructuralFrame("faulted_frame", faulted_frame)

faulted1_op = FaultSegment(fault_frame)
faulted_strati = FaultedGeologicalFeature(stratigraphy, faulted1_op)
viewer = LavaVuModelViewer()
viewer.plot_isosurface(faulted_frame[0].hw_feature, isovalue=0, colour='green')
viewer.plot_isosurface(faulted_frame[0].fw_feature, isovalue=0, colour='grey')

# viewer.plot_isosurface(fault_frame.features[0], isovalue=0, colour='red')
viewer.plot_isosurface(fault_frame2.features[0], isovalue=0, colour='red')

viewer.plot_isosurface(faulted_strati.fw_feature, slices=[-10, -15, -20], colour='blue')
viewer.plot_isosurface(faulted_strati.hw_feature, slices=[-10, -15, -20], colour='red')

# viewer.plot_isosurface(faulted_frame[0].fw_feature, isovalue=0, colour='green')
# viewer.plot_structural_frame_isosurface(fault_frame, 0, isovalue=0, colour='blue')
viewer.lv.interactive()
