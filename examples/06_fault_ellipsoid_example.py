from FME.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from FME.supports.tet_mesh import TetMesh
from FME.modelling.geological_feature import GeologicalFeature, FaultedGeologicalFeature, GeologicalFeatureBuilder
from FME.visualisation.model_visualisation import LavaVuModelViewer
from FME.modelling.structural_frame import StructuralFrameBuilder, StructuralFrame
from FME.modelling.fault.fault_segment import FaultSegment
from FME.modelling.fault.fault_function import CubicFunction, FaultDisplacement, Ones
import numpy as np
import matplotlib.pyplot as plt
"""
This is a simple case study of a duplex fault system. The model is really just 2.5 
"""
boundary_points = np.zeros((2,3))

boundary_points[0,0] = 0
boundary_points[0,1] = -20
boundary_points[0,2] = -15
boundary_points[1,0] = 20
boundary_points[1,1] = 20
boundary_points[1,2] = 1
mesh = TetMesh()
mesh.setup_mesh(boundary_points, nstep=1, n_tetra=100000,)
interpolator = PLI(mesh)
stratigraphy_builder = GeologicalFeatureBuilder(
    interpolator=interpolator,
    name='stratigraphy')

solver = 'cgp'
stratigraphy_builder.add_point([6.1,0.1,1.1],0.)
stratigraphy_builder.add_point([6.1,0.1,2.1],1.)


stratigraphy_builder.add_strike_and_dip([1,1,1],90.,0.)
stratigraphy = stratigraphy_builder.build(solver=solver,cgw=6000)

floor = -4
roof = 4

fault_interpolator = PLI(mesh)
fault = StructuralFrameBuilder(interpolator=fault_interpolator,mesh=mesh,name='FaultSegment1')
# fault.add_strike_dip_and_value([-5,-15, 0], strike=0., dip=90., val=0, itype='gx')
# fault.add_strike_dip_and_value([-5,-10, 0], strike=0., dip=90., val=0, itype='gx')
# fault.add_strike_dip_and_value([-5,-18, 0], strike=0., dip=90., val=0, itype='gx')

fault.add_strike_dip_and_value([5,0, -2], strike=90., dip=50., val=0, itype='gx')
fault.add_strike_dip_and_value([15,0, -2], strike=90., dip=50., val=0, itype='gx')
fault.add_strike_dip_and_value([10,0, -2], strike=90., dip=50., val=0, itype='gx')

# fault.add_strike_dip_and_value([5,15, 0], strike=90., dip=5., val=0, itype='gx')
# fault.add_strike_dip_and_value([15,15, 0], strike=90., dip=5., val=0, itype='gx')
# fault.add_strike_dip_and_value([10,15, 0], strike=90., dip=5., val=0, itype='gx')
#
# fault.add_strike_dip_and_value([5,-10,-10], strike=45,dip=10,val=0,itype='gx')
# fault.add_strike_dip_and_value([15,-10,-10], strike=45,dip=10,val=0,itype='gx')
# fault.add_strike_dip_and_value([10,-10,-10], strike=45,dip=10,val=0,itype='gx')

# for y in range(-5,5,1):
#     fault.add_strike_dip_and_value([-35.17,y,floor],strike=0.,dip=0.,val=0,itype='gx')
#     fault.add_strike_dip_and_value([-25.17,y,floor],strike=0.,dip=0.,val=0,itype='gx')
#     # fault.add_strike_dip_and_value([-15.17,y,floor],strike=0.,dip=0.,val=0,itype='gx')
#     # fault.add_strike_dip_and_value([-4.,y,floor],strike=0.,dip=45.,val=0,itype='gx')
#     # fault.add_strike_dip_and_value([6.17,y,roof],strike=0.,dip=0.,val=0,itype='gx')
#     # fault.add_strike_dip_and_value([4.,y,roof],strike=0.,dip=0.,val=0,itype='gx')
#     fault.add_strike_dip_and_value([-5.17,y,roof],strike=0.,dip=0.,val=0,itype='gx')
#     fault.add_strike_dip_and_value([15.17,y,roof],strike=0.,dip=0.,val=0,itype='gx')
#     fault.add_strike_dip_and_value([18.17,y,roof],strike=0.,dip=0.,val=0,itype='gx')


fault.add_point([-5,-18,0],0.,itype='gz')
fault.add_point([15,10,0],1.,itype='gz')
fault.add_point([-5,-10,0],0.,itype='gy')
fault.add_point([-5,-10,-10],1.,itype='gy')

# fault.add_point([-19,-19,-18],1.,itype='gy')

#
# for y in range(-20,20,1):
#     fault.add_point([-10.5,y,floor],0.,itype='gy')
#     fault.add_point([11.56,y,roof],1.,itype='gy')

ogw = 3000
ogw /= mesh.n_elements
cgw = 6000
# cgw = cgw / mesh.n_elements
solver='cgp'
fault_frame = fault.build(
    solver=solver,
    guess=None,
    gxxgy=2 * ogw,
    gxxgz=2 * ogw,
    gyxgz=ogw,
    gxcg=cgw,
    gycg=cgw,
    gzcg=cgw,
    shape='rectangular')
hw = CubicFunction()
hw.add_cstr(0,1)
hw.add_grad(0,0)
hw.add_cstr(20,0)
hw.add_grad(20,0)
hw.add_max(20)
fw = CubicFunction()
fw.add_cstr(0,-1)
fw.add_grad(0,0)
fw.add_cstr(-20,0)
fw.add_grad(-20,0)
fw.add_min(-20)
gyf = Ones()# CubicFunction()
# gyf.add_cstr(-.5,0)
# gyf.add_cstr(.5,0)
# gyf.add_cstr(-0.2,1)
# gyf.add_cstr(0.2,1)
# gyf.add_grad(0,0)
# gyf.add_min(-.5)
# gyf.add_max(.5)
gzf = CubicFunction()
gzf.add_cstr(-.5,0)
gzf.add_cstr(.5,0)
gzf.add_cstr(-0.2,1)
gzf.add_cstr(0.2,1)
gzf.add_grad(0,0)
gzf.add_min(-.5)
gzf.add_max(.5)

fault_displacement = FaultDisplacement(fw=fw,hw=hw,gy=gyf,gz=gzf)
fault = FaultSegment(fault_frame,
                     displacement=2,
                     # gy_min=0.15,
                     # gy_max=0.6,
                     faultfunction=fault_displacement)
faulted_strat = FaultedGeologicalFeature(stratigraphy,fault)
plt.plot(np.linspace(0,20,100),hw(np.linspace(0,20,100)))
plt.plot(np.linspace(-20,0,100),fw(np.linspace(-20,0,100)))

plt.savefig("fault_functions.png")
viewer = LavaVuModelViewer(background="white")
# viewer.plot_isosurface(faulted_frame[0].hw_feature, isovalue=0, colour='green')
# viewer.plot_isosurface(faulted_frame[0].fw_feature, isovalue=0, colour='grey')

# viewer.plot_isosurface(fault_frame.features[0], isovalue=0, colour='black')
# viewer.plot_isosurface(fault_frame.features[1], isovalue=0, colour='black')
slices = [-4,-2,0]
viewer.plot_vector_data(fault_interpolator.get_gradient_control()[:,:3],fault_interpolator.get_gradient_control()[:,3:],"grd")
viewer.plot_isosurface(fault_frame.features[0], isovalue=0, colour='green')
# viewer.plot_isosurface(fault_frame.features[1], isovalue=0.15, colour='blue')

locations = mesh.barycentre[::20,:]
# viewer.plot_vector_field(fault_frame.features[0], locations=locations, colour='red')
viewer.plot_vector_field(fault_frame.features[1], locations=locations, colour='red')
# viewer.plot_isosurface(
#     faulted_strat.hw_feature,
#     isovalue=-9)
# viewer.plot_isosurface(
#     faulted_strat.fw_feature,
#     isovalue=-9)
viewer.plot_isosurface(
    faulted_strat.feature,
    isosurface=0)
viewer.lv.interactive()