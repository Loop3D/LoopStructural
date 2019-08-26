from .discete_interpolator import DiscreteInterpolator
from .operator import Operator
import numpy as np


class FiniteDifferenceInterpolator(DiscreteInterpolator):
    """
    Finite Difference Interpolator
    """
    def __init__(self, grid):
        self.grid = grid
        self.nx = self.grid.n
        self.shape = 'rectangular'
        self.region = np.arange(0,self.grid.n).astype(int)#'everywhere'
        self.region_map = np.zeros(self.grid.n).astype(int)
        self.region_map[np.array(range(0, self.grid.n)).astype(int)] = \
            np.array(range(0, self.grid.n)).astype(int)
        DiscreteInterpolator.__init__(self)

    def _setup_interpolator(self, **kwargs):
        """
        adds all of the constraints to the interpolation matrix
        :param kwargs: 'cgw' is the constant gradient weight
        'cpw' control point weight
        'gpw' gradient control point weight
        'tpw' tangent control point weight
        'cg' boolean is cg being used
        :return:
        """
        operator = Operator.Dxy_mask
        self.assemble_inner(operator)
        operator = Operator.Dyz_mask
        self.assemble_inner(operator)
        operator = Operator.Dxz_mask
        self.assemble_inner(operator)
        operator = Operator.Dxx_mask
        self.assemble_inner(operator)
        operator = Operator.Dyy_mask
        self.assemble_inner(operator)
        operator = Operator.Dzz_mask
        self.assemble_inner(operator)

        self.add_gradient_constraint()
        self.add_vaue_constraint()


    def add_vaue_constraint(self, w=1.):
        """
        Add a value constraint to the interpolator
        :param pos: location of the constraint
        :param v: vaue to add
        :param w: weight
        :return:
        """
        points = self.get_control_points()
        # check that we have added some points
        if points.shape[0]>0:
            a = self.grid.position_to_dof_coefs(points[:,:3])
            #a*=w
            node_idx = self.grid.position_to_cell_corners(points[:,:3])
            self.add_constraints_to_least_squares(a.T, points[:,3], node_idx)
    def add_gradient_constraint(self,w=1.):
        """
        Add a gradient constraint to the interpolator
        :param pos:
        :param g:
        :param w:
        :return:
        """

        points = self.get_gradient_control()
        if points.shape[0] > 0:
            node_idx = self.grid.position_to_cell_corners(points[:,:3])
            T = self.grid.calcul_T(points[:,:3])
            self.add_constraints_to_least_squares( T[:, 0, :], points[:,3], node_idx)
            self.add_constraints_to_least_squares( T[:, 1, :], points[:,4], node_idx)
            self.add_constraints_to_least_squares( T[:, 2, :], points[:,5], node_idx)
            # node_idx = self.grid.position_to_cell_corners(pos)
            # T = self.grid.calcul_T(pos)
            #
            # self.add_constraint_to_least_squares(node_idx,T[:,0,:],g[:,0][0])
            # self.add_constraint_to_least_squares(node_idx,T[:,1,:],g[:,1][0])

        # def add_gradient_orthogonality(self,pos,g,w=1.):

            #do stuff

    def add_regularisation(self,operator,w=0.1):
        """

        :param operator:
        :param w:
        :return:
        """
        self.assemble_inner(operator)
        # self.assemble_borders()

    def assemble_inner(self, operator):
        """

        :param operator:
        :return:
        """
        # First get the global indicies of the pairs of neighbours this should be an
        # Nx27 array for 3d and an Nx9 array for 2d

        global_indexes = self.grid.neighbour_global_indexes()  # np.array([ii,jj]))

        a = np.tile(operator.flatten(), (global_indexes.shape[1],1))
        self.add_constraints_to_least_squares(a,np.zeros(global_indexes.shape[1]),
                                              global_indexes)
        return

