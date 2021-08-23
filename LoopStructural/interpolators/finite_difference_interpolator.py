"""
FiniteDifference interpolator
"""
import logging

import numpy as np

from LoopStructural.utils.helper import get_vectors
from .discrete_interpolator import DiscreteInterpolator
from .operator import Operator

from LoopStructural.utils import getLogger
logger = getLogger(__name__)


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
                                      'tpw': 1.,
                                      'ipw': 1.
                                      })

        self.vol = grid.step_vector[0] * grid.step_vector[1] * \
                   grid.step_vector[2]
        self.type = 'FDI'
        
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
                self.interpolation_weights['dxy'] = 0.1*kwargs[
                                                        'regularisation'] 
                self.interpolation_weights['dyz'] = 0.1*kwargs[
                                                        'regularisation'] 
                self.interpolation_weights['dxz'] = 0.1*kwargs[
                                                        'regularisation'] 
                self.interpolation_weights['dxx'] = 0.1*kwargs[
                                                        'regularisation'] 
                self.interpolation_weights['dyy'] = 0.1*kwargs[
                                                        'regularisation'] 
                self.interpolation_weights['dzz'] = 0.1*kwargs[
                                                        'regularisation'] 
            self.interpolation_weights[key] = kwargs[key]
        # if we want to define the operators manually
        if 'operators' in kwargs:
            for n, o in kwargs['operators'].items():
                self.assemble_inner(o[0], o[1])
        # otherwise just use defaults
        if 'operators' not in kwargs:
            
            operator = Operator.Dxy_mask
            weight =  self.interpolation_weights['dxy'] / \
                             1#(4*self.support.step_vector[0]*self.support.step_vector[1])
            self.assemble_inner(operator, weight )
            operator = Operator.Dyz_mask
            weight = self.interpolation_weights['dyz'] / \
                             1#(4*self.support.step_vector[1]*self.support.step_vector[2])
            self.assemble_inner(operator, weight)
            operator = Operator.Dxz_mask
            weight =  self.interpolation_weights['dxz'] / \
                             1#(4*self.support.step_vector[0]*self.support.step_vector[2])
            self.assemble_inner(operator, weight)
            operator = Operator.Dxx_mask
            weight = self.interpolation_weights['dxx'] \
                             / 1#self.support.step_vector[0]**2
            self.assemble_inner(operator,
                                weight)
            operator = Operator.Dyy_mask
            weight =  self.interpolation_weights['dyy'] / \
                             1#self.support.step_vector[1]**2
            self.assemble_inner(operator,weight)
            operator = Operator.Dzz_mask
            weight = self.interpolation_weights['dzz'] / \
                             1#self.support.step_vector[2]**2
            self.assemble_inner(operator,weight)
        self.add_norm_constraint(
            self.interpolation_weights['npw'])
        self.add_gradient_constraint(
             self.interpolation_weights['gpw'])
        self.add_vaue_constraint(
             self.interpolation_weights['cpw'])
        self.add_tangent_ctr_pts(
            self.interpolation_weights['tpw']
        )
        self.add_interface_ctr_pts(
            self.interpolation_weights['ipw']
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
            # a/=np.product(self.support.step_vector)
            self.add_constraints_to_least_squares(a.T * w,
                                                  points[inside, 3] * w,
                                                  idc[inside, :],
                                                  name='value')
    def add_interface_ctr_pts(self, w=1.0):  # for now weight all value points the same
        """
        Adds a constraint that defines all points with the same 'id' to be the same value
        Sets all P1-P2 = 0 for all pairs of points

        Parameters
        ----------
        w : double
            weight

        Returns
        -------

        """
        # get elements for points
        points = self.get_interface_constraints()
        if points.shape[0] > 1:
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
            a = self.support.position_to_dof_coefs(points[inside, :3]).T
            # create oversided array for storing constraints
            A = np.zeros((a.shape[0]*a.shape[0],a.shape[1]*2))
            interface_idc = np.zeros((a.shape[0]*a.shape[0],a.shape[1]*2),dtype=int)
            interface_idc[:] = -1
            c_i = 0
            
            for i in np.unique(points[np.logical_and(~np.isnan(points[:,3]),inside),3]):
                mask = points[inside,3] == i
                for p1 in range(points[inside][mask].shape[0]):
                    for p2 in range(p1+1,points[inside][mask].shape[0]):
                        A[c_i,:8] = a[mask][p1,:]
                        A[c_i,8:] -= a[mask][p2,:]
                        interface_idc[c_i,:8] = idc[inside,:][mask,:][p1,:]
                        interface_idc[c_i,8:] = idc[inside,:][mask,:][p2,:]
                        c_i+=1
            outside = ~np.any(interface_idc == -1, axis=1) 
            
            self.add_constraints_to_least_squares(A[outside,:] * w,
                                                    np.zeros(A[outside,:].shape[0]),
                                                    interface_idc[outside, :], name='interface')

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
            norm = np.linalg.norm(T,axis=2)
            T/=norm[:,:,None]
            strike_vector, dip_vector = get_vectors(points[inside, 3:6])
            A = np.einsum('ij,ijk->ik', strike_vector.T, T)

            B = np.zeros(points[inside, :].shape[0])
            self.add_constraints_to_least_squares(A * w, B, idc[inside, :], name='gradient')
            A = np.einsum('ij,ijk->ik', dip_vector.T, T)
            self.add_constraints_to_least_squares(A * w, B, idc[inside, :], name='gradient')

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
            # T*=np.product(self.support.step_vector)
            # T/=self.support.step_vector[0]
            w /= 3
            self.add_constraints_to_least_squares(T[:, 0, :] * w,
                                                  points[inside, 3] * w ,
                                                  idc[inside, :], name='norm')
            self.add_constraints_to_least_squares(T[:, 1, :] * w,
                                                  points[inside, 4] * w ,
                                                  idc[inside, :], name='norm')
            self.add_constraints_to_least_squares(T[:, 2, :] * w,
                                                  points[inside, 5] * w ,
                                                  idc[inside, :], name='norm')

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

            #normalise element vector to unit vector for dot product
            T = self.support.calcul_T(points[inside, :3])
            norm = np.linalg.norm(T,axis=1)
            T/=norm[:,None,:]
            # normalise vector to unit vector for dot product
            vector[inside,:3] /= np.linalg.norm(vector[inside,:3],axis=1)[:,None]
            # dot product of vector and element gradient 
            A = np.einsum('ij,ijk->ik', vector[inside, :3], T)
            
            B = np.zeros(points[inside, :].shape[0])
            self.add_constraints_to_least_squares(A * w, B, idc[inside, :], name='gradient orthogonal')

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
        inside = ~np.any(idc == -1, axis=1)#np.ones(a.shape[0],dtype=bool)#
        # a[idc==-1] = 0
        # idc[idc==-1] = 0
        B = np.zeros(global_indexes.shape[1])
        self.add_constraints_to_least_squares(a[inside, :] * w ,
                                              B[inside],
                                              idc[inside, :],
                                              name='regularisation'
                                              )
        return
