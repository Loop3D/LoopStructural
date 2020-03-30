import numpy as np

from LoopStructural.interpolators.piecewiselinear_interpolator import \
    PiecewiseLinearInterpolator as PLI
from LoopStructural.core.fault.fault_segment import FaultSegment
from LoopStructural.modelling.features.faulted_geological_feature import \
    FaultedGeologicalFeature
from LoopStructural.core.features.geological_feature import \
    GeologicalFeatureInterpolator
from LoopStructural.modelling.features.structural_frame import \
    StructuralFrameBuilder, StructuralFrame
from LoopStructural.supports.tet_mesh import TetMesh
from LoopStructural.visualisation.model_visualisation import LavaVuModelViewer

boundary_points = np.zeros((2, 3))
boundary_points[0, 0] = -40
boundary_points[0, 1] = -5
boundary_points[0, 2] = -10
boundary_points[1, 0] = 40
boundary_points[1, 1] = 5
boundary_points[1, 2] = 10
mesh = TetMesh()
mesh.setup_mesh(boundary_points, nstep=1, n_tetra=30000, )
interpolator = PLI(mesh)
stratigraphy_builder = GeologicalFeatureInterpolator(
    interpolator=interpolator,
    name='stratigraphy')
solver = 'lu'
stratigraphy_builder.add_point([6.1, 0.1, 1.1], 0.)
stratigraphy_builder.add_point([6.1, 0.1, 2.1], 1.)

stratigraphy_builder.add_strike_and_dip([1, 1, 1], 90., 0.)
stratigraphy = stratigraphy_builder.build(solver=solver)

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
    # fault.add_strike_dip_and_value([-15.17,y,floor],strike=0.,dip=0.,
    # val=0,itype='gx')
    # fault.add_strike_dip_and_value([-4.,y,floor],strike=0.,dip=45.,val=0,
    # itype='gx')
    # fault.add_strike_dip_and_value([6.17,y,roof],strike=0.,dip=0.,val=0,
    # itype='gx')
    # fault.add_strike_dip_and_value([4.,y,roof],strike=0.,dip=0.,val=0,
    # itype='gx')
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

fault_frame = fault.build(
    solver=solver,
)

fault_frame2_interpolator = PLI(mesh)
fault2 = StructuralFrameBuilder(
    interpolator=fault_frame2_interpolator,
    mesh=mesh,
    name='FaultSegment2'
)

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
    solver='lu',
)
# #
# # #
fault = FaultSegment(fault_frame, displacement=-4.5)

fault2 = FaultSegment(fault_frame2, displacement=-4.5)
faulted_frame = []
for f in fault_frame.features:
    faulted_frame.append(FaultedGeologicalFeature(f, fault2))
#
structural_frame2 = StructuralFrame("faulted_frame", faulted_frame)
faulted1_op = FaultSegment(structural_frame2, displacement=-4.5)
faulted_strati = FaultedGeologicalFeature(
    FaultedGeologicalFeature(stratigraphy, fault2),
    faulted1_op)  # fault_operator)#

viewer = LavaVuModelViewer(background="white")
viewer.add_isosurface(
    fault_frame2.features[0],
    isovalue=0,
    colour='black'
)
viewer.add_isosurface(
    structural_frame2.features[0],
    isovalue=0,
    colour='black'
)
viewer.add_isosurface(
    faulted_strati,
    paint_with=faulted_strati,
    # colour='green',
    nslices=2
)
viewer.save('duplex')
