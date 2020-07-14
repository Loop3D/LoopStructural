"""
FiniteDifference interpolator
"""
import logging

import numpy as np

from LoopStructural.utils.helper import get_vectors
from .discrete_interpolator import DiscreteInterpolator
from .operator import Operator

logger = logging.getLogger(__name__)


class FiniteDifferenceInterpolator(DiscreteInterpolator):
    """

    """

    def __init__(self, grid):
        """
        Finite difference interpolation on a regular cartesian grid

        Parameters
        ----------
        grid : StructuredGrid
        """
        self.shape = 'rectangular'
        DiscreteInterpolator.__init__(self, grid)
        # default weights for the interpolation matrix are 1 in x,y,z and
        # 1/
        self.set_interpolation_weights({'dxy': .7,
                                      'dyz': .7,
                                      'dxz': .7,
                                      'dxx': 1.,
                                      'dyy': 1.,
                                      'dzz': 1.,
                                      'dx': 1.,
                                      'dy': 1.,
                                      'dz': 1.,
                                      'cpw': 1.,
                                      'gpw': 1.,
                                      'npw': 1.,
                                      'tpw': 1.})

        self.vol = grid.step_vector[0] * grid.step_vector[1] * \
                   grid.step_vector[2]

    def _setup_interpolator(self, **kwargs):
        """

        Parameters
        ----------
        kwargs
            possible kwargs are weights for the different masks and masks.

        Notes
        -----
        Default masks are the second derivative in x,y,z direction and the second derivative of x wrt
        y and y wrt z and z wrt x. Custom masks can be used by specifying the operator as a 3d numpy array
        e.g. [ [ [ 0 0 0 ]
                 [ 0 1 0 ]
                 [ 0 0 0 ] ]
                 [ [ 1 1 1 ]
                 [ 1 1 1 ]
                 [ 1 1 1 ] ]
                 [ [ 0 0 0 ]
                 [ 0 1 0 ]
                 [ 0 0 0 ] ]

        Returns
        -------

        """

        for key in kwargs:
            self.up_to_date = False
            if 'regularisation' in kwargs:
                self.interpolation_weights['dxy'] = kwargs[
                                                        'regularisation'] * 0.7
                self.interpolation_weights['dyz'] = kwargs[
                                                        'regularisation'] * 0.7
                self.interpolation_weights['dxz'] = kwargs[
                                                        'regularisation'] * 0.7
                self.interpolation_weights['dxx'] = kwargs[
                                                        'regularisation'] * 1.
                self.interpolation_weights['dyy'] = kwargs[
                                                        'regularisation'] * 1.
                self.interpolation_weights['dzz'] = kwargs[
                                                        'regularisation'] * 1.
            self.interpolation_weights[key] = kwargs[key]
        # if we want to define the operators manually
        if 'operators' in kwargs:
            for n, o in kwargs['operators'].items():
                self.assemble_inner(o[0], o[1])
        # otherwise just use defaults
        if 'operators' not in kwargs:
            operator = Operator.Dxy_mask

            self.assemble_inner(operator, np.sqrt(2 * self.vol) *
                                self.interpolation_weights['dxy'])
            operator = Operator.Dyz_mask
            self.assemble_inner(operator, np.sqrt(2 * self.vol) *
                                self.interpolation_weights['dyz'])
            operator = Operator.Dxz_mask
            self.assemble_inner(operator, np.sqrt(2 * self.vol) *
                                self.interpolation_weights['dxz'])
            operator = Operator.Dxx_mask
            self.assemble_inner(operator,
                                np.sqrt(self.vol) * self.interpolation_weights[
                                    'dxx'])
            operator = Operator.Dyy_mask
            self.assemble_inner(operator,
                                np.sqrt(self.vol) * self.interpolation_weights[
                                    'dyy'])
            operator = Operator.Dzz_mask
            self.assemble_inner(operator,
                                np.sqrt(self.vol) * self.interpolation_weights[
                                    'dzz'])
        self.add_norm_constraint(
            np.sqrt(self.vol) * self.interpolation_weights['npw'])
        self.add_gradient_constraint(
            np.sqrt(self.vol) * self.interpolation_weights['gpw'])
        self.add_vaue_constraint(
            np.sqrt(self.vol) * self.interpolation_weights['cpw'])
        self.add_tangent_ctr_pts(
            np.sqrt(self.vol) * self.interpolation_weights['tpw']
        )

    def copy(self):
        """
        Create a new identical interpolator

        Returns
        -------
        returns a new empy interpolator from the same support
        """
        return FiniteDifferenceInterpolator(self.support)

    def add_vaue_constraint(self, w=1.):
        """

        Parameters
        ----------
        w : double or numpy array

        Returns
        -------

        """

        points = self.get_value_constraints()
        # check that we have added some points
        if points.shape[0] > 0:
            node_idx, inside = self.support.position_to_cell_corners(
                points[:, :3])
            # print(points[inside,:].shape)

            gi = np.zeros(self.support.n_nodes)
            gi[:] = -1
            gi[self.region] = np.arange(0, self.nx)
            idc = np.zeros(node_idx.shape)
            idc[:] = -1

            idc[inside, :] = gi[node_idx[inside, :]]
            inside = np.logical_and(~np.any(idc == -1, axis=1), inside)
            a = self.support.position_to_dof_coefs(points[inside, :3])
            # a*=w

            self.add_constraints_to_least_squares(a.T * w,
                                                  points[inside, 3] * w,
                                                  idc[inside, :])

    def add_gradient_constraint(self, w=1.):
        """

        Parameters
        ----------
        w : double / numpy array

        Returns
        -------

        """

        points = self.get_gradient_constraints()
        if points.shape[0] > 0:
            # calculate unit vector for orientation data
            # points[:,3:]/=np.linalg.norm(points[:,3:],axis=1)[:,None]

            node_idx, inside = self.support.position_to_cell_corners(
                points[:, :3])
            # calculate unit vector for node gradients
            # this means we are only constraining direction of grad not the
            # magnitude
            gi = np.zeros(self.support.n_nodes)
            gi[:] = -1
            gi[self.region] = np.arange(0, self.nx)
            idc = np.zeros(node_idx.shape)
            idc[:] = -1
            idc[inside, :] = gi[node_idx[inside, :]]
            inside = np.logical_and(~np.any(idc == -1, axis=1), inside)

            T = self.support.calcul_T(points[inside, :3])
            strike_vector, dip_vector = get_vectors(points[inside, 3:6])
            A = np.einsum('ij,ijk->ik', strike_vector.T, T)

            B = np.zeros(points[inside, :].shape[0])
            self.add_constraints_to_least_squares(A * w, B, idc[inside, :])
            A = np.einsum('ij,ijk->ik', dip_vector.T, T)
            self.add_constraints_to_least_squares(A * w, B, idc[inside, :])

    def add_norm_constraint(self, w=1.):
        """
        Add constraints to control the norm of the gradient of the scalar field

        Parameters
        ----------
        w : double
            weighting of this constraint (double)

        Returns
        -------

        """
        points = self.get_norm_constraints()
        if points.shape[0] > 0:
            # calculate unit vector for orientation data
            # points[:,3:]/=np.linalg.norm(points[:,3:],axis=1)[:,None]

            node_idx, inside = self.support.position_to_cell_corners(
                points[:, :3])
            gi = np.zeros(self.support.n_nodes)
            gi[:] = -1
            gi[self.region] = np.arange(0, self.nx)
            idc = np.zeros(node_idx.shape)
            idc[:] = -1
            idc[inside, :] = gi[node_idx[inside, :]]
            inside = np.logical_and(~np.any(idc == -1, axis=1), inside)

            # calculate unit vector for node gradients
            # this means we are only constraining direction of grad not the
            # magnitude
            T = self.support.calcul_T(points[inside, :3])

            w /= 3
            self.add_constraints_to_least_squares(T[:, 0, :] * w,
                                                  points[inside, 3] * w,
                                                  idc[inside, :])
            self.add_constraints_to_least_squares(T[:, 1, :] * w,
                                                  points[inside, 4] * w,
                                                  idc[inside, :])
            self.add_constraints_to_least_squares(T[:, 2, :] * w,
                                                  points[inside, 5] * w,
                                                  idc[inside, :])

    def add_gradient_orthogonal_constraint(self, points, vector, w=1.0,
                                           B=0):
        """
        constraints scalar field to be orthogonal to a given vector

        Parameters
        ----------
        elements : np.array
        normals : np.array
        w : double
        B : np.array

        Returns
        -------

        """
        if points.shape[0] > 0:
            # calculate unit vector for orientation data
            # points[:,3:]/=np.linalg.norm(points[:,3:],axis=1)[:,None]

            node_idx, inside = self.support.position_to_cell_corners(
                points[:, :3])
            # calculate unit vector for node gradients
            # this means we are only constraining direction of grad not the
            # magnitude
            gi = np.zeros(self.support.n_nodes)
            gi[:] = -1
            gi[self.region] = np.arange(0, self.nx)
            idc = np.zeros(node_idx.shape)
            idc[:] = -1

            idc[inside, :] = gi[node_idx[inside, :]]
            inside = np.logical_and(~np.any(idc == -1, axis=1), inside)

            T = self.support.calcul_T(points[inside, :3])
            A = np.einsum('ij,ijk->ik', vector[inside, :3], T)

            B = np.zeros(points[inside, :].shape[0])
            self.add_constraints_to_least_squares(A * w, B, idc[inside, :])

    def add_regularisation(self, operator, w=0.1):
        """

        Parameters
        ----------
        operator
        w

        Returns
        -------

        """
        self.assemble_inner(operator)
        # self.assemble_borders()

    def assemble_inner(self, operator, w):
        """

        Parameters
        ----------
        operator : Operator
        w : double

        Returns
        -------

        """
        # First get the global indicies of the pairs of neighbours this should be an
        # Nx27 array for 3d and an Nx9 array for 2d

        global_indexes = self.support.neighbour_global_indexes()  # np.array([ii,jj]))

        a = np.tile(operator.flatten(), (global_indexes.shape[1], 1))
        idc = global_indexes.T

        gi = np.zeros(self.support.n_nodes)
        gi[:] = -1
        gi[self.region] = np.arange(0, self.nx)
        idc = gi[idc]
        inside = ~np.any(idc == -1, axis=1)
        B = np.zeros(global_indexes.shape[1])
        self.add_constraints_to_least_squares(a[inside, :] * w,
                                              B[inside],
                                              idc[inside, :]
                                              )
        return
