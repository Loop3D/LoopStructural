# # Example 05 Reverse Fault
# This is an example of a reverse fault with a roll over anticline. There is
# displacement on the hanging wall
# ### Imports
from LoopStructural.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from LoopStructural.supports.tet_mesh import TetMesh
from LoopStructural.modelling.features.geological_feature import GeologicalFeatureInterpolator
from LoopStructural.modelling.features.faulted_geological_feature import FaultedGeologicalFeature
from LoopStructural.visualisation.model_visualisation import LavaVuModelViewer
from LoopStructural.modelling.structural_frame import StructuralFrameBuilder
from LoopStructural.modelling.fault.fault_segment import FaultSegment
import numpy as np

# ### Define model area and build mesh
boundary_points = np.zeros((2,3))
boundary_points[0,0] = 0
boundary_points[0,1] = -20
boundary_points[0,2] = -15
boundary_points[1,0] = 20
boundary_points[1,1] = 20
boundary_points[1,2] = 1
mesh = TetMesh()
mesh.setup_mesh(boundary_points, nstep=1, n_tetra=30000,)

# ### Build a geological feature for the faulted surface
# In this example we are going to cheat and just assume the feature was
# flat lying before faulting. In a real example the data points would
# need to be unfaulted before interpolation occured.
interpolator = PLI(mesh)
stratigraphy_builder = GeologicalFeatureInterpolator(
    interpolator=interpolator,
    name='stratigraphy')
solver = 'lu'
stratigraphy_builder.add_point([6.1,0.1,1.1],0.)
stratigraphy_builder.add_point([6.1,0.1,2.1],1.)
stratigraphy_builder.add_strike_and_dip([1,1,1],90.,0.)
stratigraphy = stratigraphy_builder.build(solver=solver,cgw=6000)

floor = -4
roof = 4

fault_interpolator = PLI(mesh)
fault = StructuralFrameBuilder(interpolator=fault_interpolator,mesh=mesh,name='FaultSegment1')
fault.add_strike_dip_and_value([5,0, -2], strike=90., dip=50., val=0, itype='gx')
fault.add_strike_dip_and_value([15,0, -2], strike=90., dip=50., val=0, itype='gx')
fault.add_strike_dip_and_value([10,0, -2], strike=90., dip=50., val=0, itype='gx')
fault.add_strike_dip_and_value([5,15, 0], strike=90., dip=5., val=0, itype='gx')
fault.add_strike_dip_and_value([15,15, 0], strike=90., dip=5., val=0, itype='gx')
fault.add_strike_dip_and_value([10,15, 0], strike=90., dip=5., val=0, itype='gx')
fault.add_strike_dip_and_value([5,-10,-10], strike=45,dip=10,val=0,itype='gx')
fault.add_strike_dip_and_value([15,-10,-10], strike=45,dip=10,val=0,itype='gx')
fault.add_strike_dip_and_value([10,-10,-10], strike=45,dip=10,val=0,itype='gx')

# To build the other coordintes just specify the extents of the fault
fault.add_point([-5,-18,0],0.,itype='gz')
fault.add_point([15,10,0],1.,itype='gz')
fault.add_point([-5,-10,0],0.,itype='gy')
fault.add_point([-5,-18,0],1.,itype='gy')

ogw = 3000
ogw /= mesh.n_elements
cgw = 6000
solver='lu'
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
fault = FaultSegment(fault_frame, displacement=10.5)
faulted_strat =FaultedGeologicalFeature(stratigraphy,fault)

viewer = LavaVuModelViewer(background="white")
slices = [-4,-2,0]
viewer.add_vector_data(fault_interpolator.get_gradient_control()[:, :3], fault_interpolator.get_gradient_control()[:, 3:], "grd")
viewer.add_isosurface(fault_frame.features[0], isovalue=0, colour='green')
locations = mesh.barycentre[::20,:]
viewer.add_vector_field(fault_frame.features[1], locations=locations, colour='red')
viewer.add_isosurface(faulted_strat, nslices=10)
viewer.add_isosurface(faulted_strat, nslices=10)
viewer.lv.interactive()