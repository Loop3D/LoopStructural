from FME.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from FME.supports.tet_mesh import TetMesh
from FME.modelling.geological_feature import GeologicalFeature, FaultedGeologicalFeature
from FME.visualisation.model_visualisation import LavaVuModelViewer
from FME.modelling.structural_frame import StructuralFrameBuilder
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

stratigraphy = PLI(mesh, region='everywhere', propertyname='Stratigraphy')
# for i in range(-15,15):
#     for j in range(-10,10):
#         stratigraphy.add_point([i,j,3*np.sin(j/5)],0.)
#         stratigraphy.add_point([i,j,3*np.sin(j/5)+1],1.)

stratigraphy.add_point([6.1,0.1,1.1],0.)
# stratigraphy.add_point([-6.1,0.1,3.1],0.)
# stratigraphy.add_point([0.1,0.1,6.1],0.)
#
# stratigraphy.add_point([2.2,1.,2.],0)
# stratigraphy.add_point([2.1,3.1,4.2],1)

stratigraphy.add_strike_and_dip([1,1,1],90.,0.)

stratigraphy.setup_interpolator(cgw=.1)
stratigraphy.solve_system(solver='lu')

support = stratigraphy.get_support()

feature = GeologicalFeature(0, 'Stratigraphy', support)


floor = -6
roof = 6

fault = StructuralFrameBuilder(interpolator=PLI,mesh=mesh,name='FaultSegment1')
for y in range(-20,20,1):
    fault.add_strike_dip_and_value([-15.17,y,floor],strike=0.,dip=0.,val=0,itype='gx')
    fault.add_strike_dip_and_value([-6.17,y,floor],strike=0.,dip=0.,val=0,itype='gx')
    fault.add_strike_dip_and_value([-4.,y,floor],strike=0.,dip=45.,val=0,itype='gx')
    fault.add_strike_dip_and_value([6.17,y,roof],strike=0.,dip=0.,val=0,itype='gx')
    fault.add_strike_dip_and_value([4.,y,roof],strike=0.,dip=0.,val=0,itype='gx')
    fault.add_strike_dip_and_value([8.17,y,roof],strike=0.,dip=0.,val=0,itype='gx')
    fault.add_strike_dip_and_value([15.17,y,roof],strike=0.,dip=0.,val=0,itype='gx')

#fault.add_strike_dip_and_value([0,0.1,floor],strike=0.,dip=90.,val=0,itype='gx')

fault.add_point([2.5,19,1.5],0.,itype='gz')
fault.add_point([2.5,-19,1.5],1.,itype='gz')
for y in range(-20,20,1):
    fault.add_point([-10.5,y,floor],0.,itype='gy')
    fault.add_point([11.56,y,roof],1.,itype='gy')


fault_frame = fault.build(solver='lu',gxxgy=0.1,gxxgz=1,gyxgz=0.05,gycg=5,gzcg=0.1)

fault_operator = FaultSegment(fault_frame)


faulted_feature = FaultedGeologicalFeature(feature, fault_operator)

viewer = LavaVuModelViewer()
viewer.plot_isosurface(faulted_feature.hw_feature, slices=[-20,-15,-13,-10])
viewer.plot_isosurface(faulted_feature.fw_feature, slices=[-20,-15,-13,-10])
viewer.plot_structural_frame_isosurface(fault_frame, 0, isovalue=0, colour='blue')
viewer.lv.interactive()
