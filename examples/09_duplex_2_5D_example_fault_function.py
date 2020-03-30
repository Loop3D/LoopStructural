# # Example 04 2.5 duplex
# In this example we will create a similar model to example 04 except
# the fault displacement magnitude will be defined by a function of
# the fault frame coordinates.
#
import numpy as np

from LoopStructural.interpolators.piecewiselinear_interpolator import \
    PiecewiseLinearInterpolator as PLI
from LoopStructural.modelling.fault.fault_function import CubicFunction, \
    FaultDisplacement, Ones
from LoopStructural.modelling.fault.fault_segment import FaultSegment
from LoopStructural.core.features.faulted_geological_feature import \
    FaultedGeologicalFeature
from LoopStructural.modelling.features.geological_feature import \
    GeologicalFeatureInterpolator
from LoopStructural.modelling.features.structural_frame import \
    StructuralFrameBuilder, StructuralFrame
from LoopStructural.supports.tet_mesh import TetMesh
from LoopStructural.visualisation.model_visualisation import LavaVuModelViewer

boundary_points = np.zeros((2, 3))

boundary_points[0, 0] = -40
boundary_points[0, 1] = -15
boundary_points[0, 2] = -10
boundary_points[1, 0] = 40
boundary_points[1, 1] = 15
boundary_points[1, 2] = 10
mesh = TetMesh()
mesh.setup_mesh(boundary_points, nstep=1, n_tetra=60000, )
interpolator = PLI(mesh)
stratigraphy_builder = GeologicalFeatureInterpolator(
    interpolator=interpolator,
    name='stratigraphy')

solver = 'lu'
stratigraphy_builder.add_point([6.1, 0.1, 1.1], 0.)
stratigraphy_builder.add_point([6.1, 0.1, 2.1], 1.)

stratigraphy_builder.add_strike_and_dip([1, 1, 1], 90., 0.)
stratigraphy = stratigraphy_builder.build(solver=solver, cgw=.1)

floor = -4
roof = 4

fault_interpolator = PLI(mesh)
fault = StructuralFrameBuilder(interpolator=fault_interpolator, mesh=mesh,
                               name='FaultSegment1')
for y in range(-5, 5, 1):
    fault.add_strike_dip_and_value([-35.17, y, floor], strike=0., dip=0.,
                                   val=0, itype='gx')
    fault.add_strike_dip_and_value([-25.17, y, floor], strike=0., dip=0.,
                                   val=0, itype='gx')
    fault.add_strike_dip_and_value([-5.17, y, roof], strike=0., dip=0., val=0,
                                   itype='gx')
    fault.add_strike_dip_and_value([15.17, y, roof], strike=0., dip=0., val=0,
                                   itype='gx')
    fault.add_strike_dip_and_value([18.17, y, roof], strike=0., dip=0., val=0,
                                   itype='gx')

fault.add_point([2.5, 5, 1.5], 0., itype='gz')
fault.add_point([2.5, -5, 1.5], 1., itype='gz')
for y in range(-20, 20, 1):
    fault.add_point([-10.5, y, floor], 0., itype='gy')
    fault.add_point([11.56, y, roof], 1., itype='gy')

ogw = 3000
ogw /= mesh.n_elements
cgw = 5000
# cgw = cgw / mesh.n_elements

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

fault_frame2_interpolator = PLI(mesh)
fault2 = StructuralFrameBuilder(
    interpolator=fault_frame2_interpolator,
    mesh=mesh,
    name='FaultSegment2'
)
#
for y in range(-5, 5, 1):
    fault2.add_strike_dip_and_value([-18.17, y, floor], strike=0., dip=0.,
                                    val=0, itype='gx')
    fault2.add_strike_dip_and_value([-15.17, y, floor], strike=0., dip=0.,
                                    val=0, itype='gx')
    fault2.add_strike_dip_and_value([-10.17, y, floor], strike=0., dip=0.,
                                    val=0, itype='gx')
    fault2.add_strike_dip_and_value([8.17, y, roof], strike=0., dip=0., val=0,
                                    itype='gx')
    fault2.add_strike_dip_and_value([15.17, y, roof], strike=0., dip=0., val=0,
                                    itype='gx')
    fault2.add_strike_dip_and_value([18.17, y, roof], strike=0., dip=0., val=0,
                                    itype='gx')

    fault_operator = FaultSegment(fault_frame)

fault2.add_point([2.5, 5, 1.5], 0., itype='gz')
fault2.add_point([2.5, -5, 1.5], 1., itype='gz')
for y in range(-5, 5, 1):
    fault2.add_point([-10.5, y, floor], 0., itype='gy')
    fault2.add_point([11.56, y, roof], 1., itype='gy')
#

fault_frame2 = fault2.build(
    solver=solver,
    guess=None,
    gxxgy=2 * ogw,
    gxxgz=2 * ogw,
    gyxgz=ogw,
    gxcg=cgw,
    gycg=cgw,
    gzcg=cgw,
    shape='rectangular',
    gx=True,
    gy=True,
    gz=True
)
# #
# # #
hw = CubicFunction()
hw.add_cstr(0, -1)
hw.add_grad(0, 0)
hw.add_cstr(20, 0)
hw.add_grad(20, 0)
hw.add_max(20)
fw = CubicFunction()
fw.add_cstr(0, 1)
fw.add_grad(0, 0)
fw.add_cstr(-20, 0)
fw.add_grad(-20, 0)
fw.add_min(-20)
gyf = Ones()
# gzf = Ones()
gzf = CubicFunction()
gzf.add_cstr(-1, 0)
gzf.add_cstr(1, 0)
gzf.add_cstr(-0.2, 1)
gzf.add_cstr(0.2, 1)
gzf.add_grad(0, 0)
gzf.add_min(-1)
gzf.add_max(1)
fault_displacement = FaultDisplacement(fw=fw, hw=hw, gy=gyf, gz=gzf)

fault = FaultSegment(
    fault_frame,
    displacement=2.5,
    faultfunction=fault_displacement)

fault2 = FaultSegment(fault_frame2,
                      displacement=2.5,
                      faultfunction=fault_displacement)
faulted_frame = []
for f in fault_frame.features:
    faulted_frame.append(FaultedGeologicalFeature(f, fault2))
#
structural_frame2 = StructuralFrame("faulted_frame", faulted_frame)

faulted1_op = FaultSegment(structural_frame2,
                           displacement=5,
                           faultfunction=fault_displacement)

faulted_strati = FaultedGeologicalFeature(
    FaultedGeologicalFeature(stratigraphy, fault2),
    faulted1_op)  # fault_operator)#

viewer = LavaVuModelViewer(background="white")
# viewer.plot_isosurface(faulted_frame[0].hw_feature, isovalue=0,
# colour='green')
# viewer.plot_isosurface(faulted_frame[0].fw_feature, isovalue=0,
# colour='grey')

# viewer.plot_isosurface(fault_frame.features[0], isovalue=0, colour='black')
# viewer.plot_isosurface(fault_frame.features[1], isovalue=0, colour='black')
slices = [-4, -2, 0]

viewer.add_isosurface(
    fault_frame2.features[0],
    isovalue=0,
    colour='black'
)
# viewer.plot_isosurface(
#     fault_frame.features[0],
#     isovalue=0,
#     colour='pink'
# )
viewer.add_isosurface(
    structural_frame2.features[0],
    isovalue=0,
    colour='black'
)
# viewer.plot_isosurface(fault_frame2.features[0], isovalue=0, colour='red')
# viewer.plot_isosurface(fault_frame2.features[1], isovalue=0, colour='purple')

locations = mesh.barycentre[::10, :]
# viewer.plot_vector_field(fault_frame2.features[2],locations,colour='green')
# viewer.plot_vector_field(fault_frame2.features[1], locations,colour='blue')
# viewer.plot_vector_field(fault_frame2.features[2], locations,colour='green')

# viewer.plot_vector_field(fault_frame2.features[1],locations,colour='red')

# viewer.plot_isosurface(faulted_strati.fw_feature, slices=[-10, -15, -20], colour='blue')
# print(faulted_strati.support.property_name)
# viewer.plot_isosurface(
#     faulted_strati,
#     nslices=10,
#     paint_with=faulted_strati
#     # colour='green'
# )
viewer.add_isosurface(
    faulted_strati,
    slices=slices,
    paint_with=faulted_strati
    # colour='blue'
)

# viewer.plot_isosurface(faulted_frame[0].fw_feature, isovalue=0, colour='green')
# viewer.plot_structural_frame_isosurface(fault_frame, 0, isovalue=0, colour='blue')
viewer.interactive()
#
