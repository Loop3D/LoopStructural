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
mesh.setup_mesh(boundary_points,n_tetra=60000)

stratigraphy = PLI(mesh,region='everywhere',propertyname='stratigraphy')
for y in range(-18,18,1):
    stratigraphy.add_strike_dip_and_value([-10,y,-10],90.,0.,0.)
    # stratigraphy.add_strike_dip_and_value([0,y,-6],90.,40.,0.)
    stratigraphy.add_strike_dip_and_value([5,y,7.1],90.,0.,0.)
#     stratigraphy.add_point([-10,y,-10],0.)
#     stratigraphy.add_point([0,y,-6],0.)
#     stratigraphy.add_point([5,y,7.1],0.)
# stratigraphy.add_point([0.1,0,2.1],1.)
# stratigraphy.add_strike_and_dip([6.1,0.1,-15.1],90.,0.)
# stratigraphy.add_strike_and_dip([1.1,0.1,-15.1],90.,0.)
# stratigraphy.add_strike_and_dip([1.1,1.1,-15.1],90.,0.)
cgw = 100 / mesh.n_elements
stratigraphy.setup_interpolator(cgw=cgw)
# stratigraphy.calculate_constant_gradient_with_element_weighting(w)
stratigraphy.solve_system(solver='chol',clear=True)
print(stratigraphy.c)
support = stratigraphy.get_support()

feature = GeologicalFeature(0, 'Stratigraphy', support)

viewer = LavaVuModelViewer()
viewer.plot_isosurface(feature, isovalue=0)
viewer.plot_points(stratigraphy.get_control_points()[:,0:3],"value")
viewer.plot_vector_data(stratigraphy.get_gradient_control()[:,0:3],stratigraphy.get_gradient_control()[:,3:],"vector")
viewer.lv.control.Panel()
viewer.lv.control.ObjectList()
viewer.lv.window()
viewer.lv.interactive()
