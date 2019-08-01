from FME.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from FME.supports.tet_mesh import TetMesh
from FME.modelling.geological_feature import GeologicalFeature, GeologicalFeatureBuilder
from FME.visualisation.model_visualisation import LavaVuModelViewer
import numpy as np
import lavavu
import matplotlib.pyplot as plt
"""
This is a basic example showing how to use the Piecewise Linear Interpolator for orientation and
value data points. 
"""
boundary_points = np.zeros((2,3))

boundary_points[0,0] = -1
boundary_points[0,1] = -1
boundary_points[0,2] = -1
boundary_points[1,0] = 1
boundary_points[1,1] = 1
boundary_points[1,2] = 1
mesh = TetMesh()
mesh.setup_mesh(boundary_points, nstep=1, n_tetra=100000,)
np.savetxt("01_box_coords.txt",mesh.points)
interpolator = PLI(mesh)
feature_builder = GeologicalFeatureBuilder(interpolator,name='stratigraphy')

feature_builder.add_point([.9,.9,.9],0)
feature_builder.add_point([-0.5,0,0],1)
feature_builder.add_point([-.9,0,0],.8)

# feature_builder.add_strike_and_dip([0,0,0],90,40)
cgw = 300
cgw /= mesh.n_elements
feature = feature_builder.build(
    solver='cg',
    cgw=cgw)

np.savetxt("01_gradient.txt",feature_builder.interpolator.get_gradient_control())
np.savetxt("01_value.txt",feature_builder.interpolator.get_control_points())
cmap = lavavu.matplotlib_colourmap("Greys")

# mesh.save()
# viewer = LavaVuModelViewer(background="white")
# viewer.plot_isosurface(feature, slices=[0, 1., .8])
#
# viewer.plot_vector_data(feature_builder.interpolator.get_gradient_control()[:,:3],feature_builder.interpolator.get_gradient_control()[:,3:],"grad")
# viewer.plot_value_data(feature_builder.interpolator.get_control_points()[:,:3],feature_builder.interpolator.get_control_points()[:,3:],"value",pointsize=10,colourmap=cmap)
# viewer.lv.interactive()
# AA = interpolator.AA.T.dot(interpolator.AA)
# # from scipy.sparse.csgraph import reverse_cuthill_mckee
# # from scipy.sparse import coo_matrix
# # def permute_sparse_matrix(M, order):
# #     print(M.row)
# #     permuted_row = order[M.row]
# #     permuted_col = order[M.col]
# #     new_M = coo_matrix((M.data, (permuted_row, permuted_col)), shape=M.shape)
# #     return new_M
# # order = reverse_cuthill_mckee(AA.tocsr())
# # print(order)
# # permuted = permute_sparse_matrix(AA.tocoo(),order).todense()
# # permuted[permuted==0] = np.nan
# AA = np.linalg.inv(AA.todense())
# AA[AA==0]=np.nan
# plt.imshow(AA)
# # plt.figure()
# # plt.imshow(permuted)
# plt.show()
# print(AA)
