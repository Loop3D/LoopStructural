from .discete_interpolator import DiscreteInterpolator
from .operator import Operator
import numpy as np


class FiniteDifferenceInterpolator(DiscreteInterpolator):
    """
    Finite Difference Interpolator
    """
    def __init__(self, grid):
        self.support = grid
        self.nx = self.support.n_nodes
        self.shape = 'rectangular'
        self.region = np.arange(0, self.nx).astype(int)#'everywhere'
        self.region_map = np.zeros(self.nx).astype(int)
        self.region_map[np.array(range(0, self.nx)).astype(int)] = \
            np.array(range(0, self.nx)).astype(int)
        DiscreteInterpolator.__init__(self)
        self.support = grid
        # default weights for the interpolation matrix are 1 in x,y,z and
        # 1/
        self.interpolation_weights = {'dxy': .7,
                                      'dyz': .7,
                                      'dxz': .7,
                                      'dxx': 1.,
                                      'dyy': 1.,
                                      'dzz': 1.,
                                      'dx': 1.,
                                      'dy': 1.,
                                      'dz': 1.,
                                      'cpw':1.,
                                      'gpw':1.}
        self.vol = grid.step_vector[0]*grid.step_vector[1]*grid.step_vector[2]

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
        for key in kwargs:
            self.up_to_date = False
            self.interpolation_weights[key] = kwargs[key]
        # if we want to define the operators manually
        if 'operators' in kwargs:
            for n,o in kwargs['operators'].items():
                self.assemble_inner(o[0], o[1])
        # otherwise just use defaults
        if 'operators' not in kwargs:
            operator = Operator.Dxy_mask
            self.assemble_inner(operator, np.sqrt(2*self.vol)*self.interpolation_weights['dxy'])
            operator = Operator.Dyz_mask
            self.assemble_inner(operator, np.sqrt(2*self.vol)*self.interpolation_weights['dyz'])
            operator = Operator.Dxz_mask
            self.assemble_inner(operator, np.sqrt(2*self.vol)*self.interpolation_weights['dxz'])
            operator = Operator.Dxx_mask
            self.assemble_inner(operator, np.sqrt(self.vol)*self.interpolation_weights['dxx'])
            operator = Operator.Dyy_mask
            self.assemble_inner(operator, np.sqrt(self.vol)*self.interpolation_weights['dyy'])
            operator = Operator.Dzz_mask
            self.assemble_inner(operator, np.sqrt(self.vol)*self.interpolation_weights['dzz'])

        self.add_gradient_constraint(np.sqrt(self.vol)*self.interpolation_weights['gpw'])
        self.add_vaue_constraint(np.sqrt(self.vol)*self.interpolation_weights['cpw'])

    def copy(self):
        return FiniteDifferenceInterpolator(self.support)

    def add_vaue_constraint(self, w=1.):
        """
        Add a value constraint to the interpolator
        :param pos: location of the constraint
        :param v: vaue to add
        :param w: weight
        :return:
        """
        points = self.get_value_constraints()
        # check that we have added some points
        if points.shape[0]>0:
            node_idx, inside = self.support.position_to_cell_corners(points[:, :3])
            #print(points[inside,:].shape)

            a = self.support.position_to_dof_coefs(points[inside, :3])
            #a*=w
            self.add_constraints_to_least_squares(a.T*w, points[inside,3]*w, node_idx[inside,:])

    def add_gradient_constraint(self,w=1.):
        """
        Add a gradient constraint to the interpolator
        :param pos:
        :param g:
        :param w:
        :return:
        """

        points = self.get_gradient_constraints()
        if points.shape[0] > 0:
            # calculate unit vector for orientation data
            # points[:,3:]/=np.linalg.norm(points[:,3:],axis=1)[:,None]

            node_idx, inside = self.support.position_to_cell_corners(points[:, :3])
            # calculate unit vector for node gradients
            # this means we are only constraining direction of grad not the magnitude
            T = self.support.calcul_T(points[inside, :3])
            norm = np.linalg.norm(T,axis=2)
            T /= norm[:,:,None]
            points[inside,3:]/= norm
            # T /= np.linalg.norm(T,axis=1)[:,None,:]
            w /= 3
            self.add_constraints_to_least_squares(T[:, 0, :]*w, points[inside, 3]*w, node_idx[inside, :])
            self.add_constraints_to_least_squares(T[:, 1, :]*w, points[inside, 4]*w, node_idx[inside, :])
            self.add_constraints_to_least_squares(T[:, 2, :]*w, points[inside, 5]*w, node_idx[inside, :])

    def add_gradient_orthogonal_constraint(self, elements, normals, w=1.0, B=0):
        """
        constraints scalar field to be orthogonal to a given vector
        Parameters
        ----------
        elements
        normals
        w
        B

        Returns
        -------

        """
        """
        :param elements: index of elements to apply constraint to
        :param normals: list of normals for elements
        :param w: global weighting per constraint
        :param B: norm value
        :return:
        """

        # get the cell gradient for the global indices
        ix,iy,iz = self.support.global_index_to_cell_index(elements)
        cornerx, cornery, cornerz = self.support.cell_corner_indexes(
            ix,iy,iz
            )
        # posx, posy, posz = self.support.node_indexes_to_position(cornerx, cornery, cornerz)
        pos = np.array(self.support.cell_centres(elements))#np.array([np.mean(posx, axis=1),np.mean(posz, axis=1), np.mean(posz, axis=1)]).T
        T = self.support.calcul_T(pos)
        # find the dot product between the normals and the gradient and add this as a
        # constraint
        A = np.einsum('ij,ijk->ik',normals,T)
        a = np.array([cornerx.T,cornery.T,cornerz.T])
        idc = self.support.global_indicies(a)#elements[elements]
        B = np.zeros(len(elements))
        self.add_constraints_to_least_squares(A*w, B, idc)
    def add_regularisation(self,operator,w=0.1):
        """

        :param operator:
        :param w:
        :return:
        """
        self.assemble_inner(operator)
        # self.assemble_borders()

    def assemble_inner(self, operator, w):
        """

        :param operator:
        :return:
        """
        # First get the global indicies of the pairs of neighbours this should be an
        # Nx27 array for 3d and an Nx9 array for 2d

        global_indexes = self.support.neighbour_global_indexes()  # np.array([ii,jj]))

        a = np.tile(operator.flatten(), (global_indexes.shape[1], 1))

        self.add_constraints_to_least_squares(a*w, np.zeros(global_indexes.shape[1]),
                                              global_indexes.T)
        return

