# from scipy.sparse import linalg
# from scipy.sparse import coo_matrix, spdiags
# from .discete_interpolator import DiscreteInterpolator
# class FiniteDifferenceInterpolator(DiscreteInterpolator):
#     """
#     Finite Difference Interpolator
#     """
#     def __init__(self, **kwargs):
#         self.support =
#         #define cartesian grid
#         self.params = {}
#         self.params['shape'] = 'rectangular'
#         self.params['nx'] = 10
#         self.params['ny'] = 10
#         self.params['nz'] = 10
#         self.AA = []
#         self.B = []
#         self.constraints = []
#
#     def add_vaue_constraint(self,pos,v,w=1.):
#         """
#         Add a value constraint to the interpolator
#         :param pos: location of the constraint
#         :param v: vaue to add
#         :param w: weight
#         :return:
#         """
#         a = self.grid.position_to_dof_coefs(pos)
#         #a*=w
#         node_idx = self.grid.position_to_cell_corners(pos)
#         self.add_constraint_to_least_squares(node_idx, a, v)
#     def add_gradient_constraint(self,pos,g,w=1.):
#         """
#         Add a gradient constraint to the interpolator
#         :param pos:
#         :param g:
#         :param w:
#         :return:
#         """
#         node_idx = self.grid.position_to_cell_corners(pos)
#         T = self.grid.calcul_T(pos)
#
#         self.add_constraint_to_least_squares(node_idx,T[:,0,:],g[:,0][0])
#         self.add_constraint_to_least_squares(node_idx,T[:,1,:],g[:,1][0])
#
#     # def add_gradient_orthogonality(self,pos,g,w=1.):
#
#         #do stuff
#
#     def add_regularisation(self,operator,w=0.1):
#         """
#
#         :param operator:
#         :param w:
#         :return:
#         """
#         self.assemble_inner(operator)
#         # self.assemble_borders()
#
#
#
