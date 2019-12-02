# # Example 06 ellipsoid fault displacement
# This example shows how to create a fault with an ellipsoidal displacement
# field with
# displacements magnitude controlled by three cubic splines
# ### Imports
import numpy as np

from LoopStructural.interpolators.piecewiselinear_interpolator import \
    PiecewiseLinearInterpolator as PLI
from LoopStructural.modelling.fault.fault_function import CubicFunction, \
    FaultDisplacement, Ones
from LoopStructural.modelling.fault.fault_segment import FaultSegment
from LoopStructural.modelling.features.faulted_geological_feature import \
    FaultedGeologicalFeature
from LoopStructural.modelling.features.geological_feature import \
    GeologicalFeatureInterpolator
from LoopStructural.modelling.features.structural_frame import \
    StructuralFrameBuilder
from LoopStructural.supports.tet_mesh import TetMesh
from LoopStructural.visualisation.model_visualisation import LavaVuModelViewer

# ### Define model area and build mesh
boundary_points = np.zeros((2, 3))

boundary_points[0, 0] = 0
boundary_points[0, 1] = -20
boundary_points[0, 2] = -15
boundary_points[1, 0] = 20
boundary_points[1, 1] = 20
boundary_points[1, 2] = 1
mesh = TetMesh()
mesh.setup_mesh(boundary_points, nstep=1, n_tetra=10000, )

# ### Build a geological feature for the faulted surface
# In this example we are going to cheat and just assume the feature was
# flat lying before faulting. In a real example the data points would
# need to be unfaulted before interpolation occured.

interpolator = PLI(mesh)
stratigraphy_builder = GeologicalFeatureInterpolator(
    interpolator=interpolator,
    name='stratigraphy')
solver = 'lu'
# stratigraphy_builder.add_point([6.1,0.1,1.1],0.)
stratigraphy_builder.add_point([6.1, 0.1, 2.1], 1.)
stratigraphy_builder.add_strike_and_dip([1, 1, 1], 90., 0.)
stratigraphy = stratigraphy_builder.build(solver=solver)

# ### Build fault frame
# Create a base interpolator for the fault frame and add the data to the
# builder
fault_interpolator = PLI(mesh)
fault = StructuralFrameBuilder(interpolator=fault_interpolator, mesh=mesh,
                               name='FaultSegment1')
fault.add_strike_dip_and_value([5, 0, -2], strike=90., dip=50., val=0,
                               itype='gx')
fault.add_strike_dip_and_value([15, 0, -2], strike=90., dip=50., val=0,
                               itype='gx')
fault.add_strike_dip_and_value([10, 0, -2], strike=90., dip=50., val=0,
                               itype='gx')
fault.add_point([5, -5, 0], -1, itype='gx')
fault.add_point([5, 5, 0], 1, itype='gx')

fault.add_point([-5, -18, 0], 0., itype='gz')
fault.add_point([15, 10, 0], 1., itype='gz')
fault.add_point([-5, -10, 0], 0., itype='gy')
fault.add_point([-5, -10, -10], 1., itype='gy')
fault.add_strike_and_dip([10, 0, -2], 270., 40., itype='gy')

ogw = 3000
ogw /= mesh.n_elements
cgw = 6000
fault_frame = fault.build(
    solver=solver,

    shape='rectangular')

# ### Define fault displacements
# Hanging wall defines the displacement of the fault for the hanging wall
#

hw = CubicFunction()
hw.add_cstr(0, 1)
hw.add_grad(0, 0)
hw.add_cstr(1, 0)
hw.add_grad(1, 0)
hw.add_max(1)
fw = CubicFunction()
fw.add_cstr(0, -1)
fw.add_grad(0, 0)
fw.add_cstr(-1, 0)
fw.add_grad(-1, 0)
fw.add_min(-1)
gyf = Ones()  # CubicFunction()
# gyf.add_cstr(-.5,0)
# gyf.add_cstr(.5,0)
# gyf.add_cstr(-0.2,1)
# gyf.add_cstr(0.2,1)
# gyf.add_grad(0,0)
# gyf.add_min(-.5)
# gyf.add_max(.5)
gzf = Ones()
# gzf.add_cstr(-.5,0)
# gzf.add_cstr(.5,0)
# gzf.add_cstr(-0.2,1)
# gzf.add_cstr(0.2,1)
# gzf.add_grad(0,0)
# gzf.add_min(-.5)
# gzf.add_max(.5)

# ### Fault displacement field
# Fault displacement can be defined by a volumetric function that will
# return the
# displacemnt of the fault anywhere in the model, using either the hanging wall
# or the footwall functions
fault_displacement = FaultDisplacement(fw=fw, hw=hw, gy=gyf, gz=gzf)
fault = FaultSegment(fault_frame,
                     displacement=2,  # scaling parameter
                     faultfunction=fault_displacement
                     )

# ### Faulted geological feature
# Create a faulted geological feature from the stratigraphy feature and the
# fault
faulted_strat = FaultedGeologicalFeature(stratigraphy, fault)
# plt.plot(np.linspace(0,20,100),hw(np.linspace(0,20,100)))
# plt.plot(np.linspace(-20,0,100),fw(np.linspace(-20,0,100)))

# ### Visualisation
# To visualise faults we can either visualise the whole faulted feature by
# accessing
# **faulted_feature.feature**  or we can look at the hanging wall feature by
# accessing
# ** hw_feature** or the **fw_feature**.
viewer = LavaVuModelViewer(background="white")
slices = [-4, -2, 0]
# viewer.plot_vector_data(fault_interpolator.get_gradient_control()[:,:3],
# fault_interpolator.get_gradient_control()[:,3:],"grd")
# viewer.plot_isosurface(fault_frame.features[0], isovalue=0, colour='green')
# viewer.plot_isosurface(fault_frame.features[1], isovalue=0.15, colour='blue')
locations = mesh.barycentre[::20, :]
# viewer.plot_vector_field(fault_frame.features[1], locations=locations,
# colour='red')
# viewer.plot_isosurface(
#     faulted_strat.hw_feature,
#     isovalue=-9)
viewer.add_isosurface(
    fault,
    slices=[-1, 0, 1],
    paint_with=fault,
    cmap='rainbow')
viewer.add_isosurface(
    fault_frame[1],
    slices=[-1, 0, 1],
    paint_with=fault_frame[1],
    cmap='rainbow')
viewer.add_isosurface(
    fault_frame[2],
    slices=[-1, 0, 1],
    paint_with=fault_frame[2],
    cmap='rainbow')
viewer.add_isosurface(
    faulted_strat,
    nslices=2,
    paint_with=faulted_strat
)
# viewer.add_isosurface(
#     stratigraphy,
#     nslices=10,
#     paint_with=faulted_strat
# )
viewer.lv.rotate([-85.18760681152344, 42.93233871459961, 0.8641873002052307])

# viewer.save('fault_ellipsoid.png')
viewer.interactive()
