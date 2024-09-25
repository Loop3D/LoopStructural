"""
FiniteDifference interpolator
"""

import numpy as np

from ..utils import get_vectors
from ._discrete_interpolator import DiscreteInterpolator
from ..interpolators import InterpolatorType

from ._operator import Operator

from LoopStructural.utils import getLogger

logger = getLogger(__name__)


class FiniteDifferenceInterpolator(DiscreteInterpolator):
    def __init__(self, grid, data={}):
        """
        Finite difference interpolation on a regular cartesian grid

        Parameters
        ----------
        grid : StructuredGrid
        """
        self.shape = "rectangular"
        DiscreteInterpolator.__init__(self, grid, data=data)
        # default weights for the interpolation matrix are 1 in x,y,z and
        # 1/
        self.set_interpolation_weights(
            {
                "dxy": 0.7,
                "dyz": 0.7,
                "dxz": 0.7,
                "dxx": 1.0,
                "dyy": 1.0,
                "dzz": 1.0,
                "dx": 1.0,
                "dy": 1.0,
                "dz": 1.0,
                "cpw": 1.0,
                "gpw": 1.0,
                "npw": 1.0,
                "tpw": 1.0,
                "ipw": 1.0,
            }
        )

        # grid.step_vector[2]
        self.type = InterpolatorType.FINITE_DIFFERENCE

    def setup_interpolator(self, **kwargs):
        """

        Parameters
        ----------
        kwargs
            possible kwargs are weights for the different masks and masks.

        Notes
        -----
        Default masks are the second derivative in x,y,z direction and the second
        derivative of x wrt y and y wrt z and z wrt x. Custom masks can be used
        by specifying the operator as a 3d numpy array
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
        self.reset()
        for key in kwargs:
            self.up_to_date = False
            if "regularisation" in kwargs:
                self.interpolation_weights["dxy"] = 0.1 * kwargs["regularisation"]
                self.interpolation_weights["dyz"] = 0.1 * kwargs["regularisation"]
                self.interpolation_weights["dxz"] = 0.1 * kwargs["regularisation"]
                self.interpolation_weights["dxx"] = 0.1 * kwargs["regularisation"]
                self.interpolation_weights["dyy"] = 0.1 * kwargs["regularisation"]
                self.interpolation_weights["dzz"] = 0.1 * kwargs["regularisation"]
            self.interpolation_weights[key] = kwargs[key]
        # either use the default operators or the ones passed to the function
        operators = kwargs.get(
            "operators", self.support.get_operators(weights=self.interpolation_weights)
        )
        for k, o in operators.items():
            self.assemble_inner(o[0], o[1], name=k)

        self.add_norm_constraints(self.interpolation_weights["npw"])
        self.add_gradient_constraints(self.interpolation_weights["gpw"])
        self.add_value_constraints(self.interpolation_weights["cpw"])
        self.add_tangent_constraints(self.interpolation_weights["tpw"])
        self.add_interface_constraints(self.interpolation_weights["ipw"])
        self.add_value_inequality_constraints()
        self.add_inequality_pairs_constraints(
            kwargs.get('inequality_pairs', None),
            upper_bound=kwargs.get('inequality_pair_upper_bound', np.finfo(float).eps),
            lower_bound=kwargs.get('inequality_pair_lower_bound', -np.inf),
        )

    def copy(self):
        """
        Create a new identical interpolator

        Returns
        -------
        returns a new empy interpolator from the same support
        """
        return FiniteDifferenceInterpolator(self.support)

    def add_value_constraints(self, w=1.0):
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
                points[:, : self.support.dimension]
            )
            # print(points[inside,:].shape)
            gi = np.zeros(self.support.n_nodes, dtype=int)
            gi[:] = -1
            gi[self.region] = np.arange(0, self.nx, dtype=int)
            idc = np.zeros(node_idx.shape)
            idc[:] = -1
            idc[inside, :] = gi[node_idx[inside, :]]
            inside = np.logical_and(~np.any(idc == -1, axis=1), inside)
            a = self.support.position_to_dof_coefs(points[inside, : self.support.dimension])
            # a *= w
            # a/=self.support.enp.product(self.support.step_vector)
            self.add_constraints_to_least_squares(
                a,
                points[inside, 3],
                idc[inside, :],
                w=w * points[inside, 4],
                name="value",
            )
            if np.sum(inside) <= 0:
                logger.warning(
                    f"{np.sum(~inside)} \
                        value constraints not added: outside of model bounding box"
                )

    def add_interface_constraints(self, w=1.0):
        """
        Adds a constraint that defines all points
        with the same 'id' to be the same value
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
            vertices, c, tetras, inside = self.support.get_element_for_location(
                points[:, : self.support.dimension]
            )
            # calculate volume of tetras
            # vecs = vertices[inside, 1:, :] - vertices[inside, 0, None, :]
            # vol = np.abs(np.linalg.det(vecs)) / 6
            A = c[inside, :]
            # A *= vol[:,None]
            idc = tetras[inside, :]

            for unique_id in np.unique(points[np.logical_and(~np.isnan(points[:, 3]), inside), 3]):
                mask = points[inside, 3] == unique_id
                ij = np.array(
                    np.meshgrid(
                        np.arange(0, A[mask, :].shape[0]),
                        np.arange(0, A[mask, :].shape[0]),
                    )
                ).T.reshape(-1, 2)
                interface_A = np.hstack([A[mask, :][ij[:, 0], :], -A[mask, :][ij[:, 1], :]])
                interface_idc = np.hstack([idc[mask, :][ij[:, 0], :], idc[mask, :][ij[:, 1], :]])

                # now map the index from global to region create array size of mesh
                # initialise as np.nan, then map points inside region to 0->nx
                gi = np.zeros(self.support.n_nodes).astype(int)
                gi[:] = -1

                gi[self.region] = np.arange(0, self.nx)
                interface_idc = gi[interface_idc]
                outside = ~np.any(interface_idc == -1, axis=1)
                self.add_constraints_to_least_squares(
                    interface_A[outside, :],
                    np.zeros(interface_A[outside, :].shape[0]),
                    interface_idc[outside, :],
                    w=w,
                    name="interface_{}".format(unique_id),
                )

    def add_gradient_constraints(self, w=1.0):
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

            node_idx, inside = self.support.position_to_cell_corners(points[:, :3])
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

            (
                vertices,
                T,
                elements,
                inside_,
            ) = self.support.get_element_gradient_for_location(points[inside, :3])
            # normalise constraint vector and scale element matrix by this
            norm = np.linalg.norm(points[:, 3:6], axis=1)
            points[:, 3:6] /= norm[:, None]
            T /= norm[inside, None, None]
            # calculate two orthogonal vectors to constraint (strike and dip vector)
            strike_vector, dip_vector = get_vectors(points[inside, 3:6])
            A = np.einsum("ij,ijk->ik", strike_vector.T, T)
            B = np.zeros(points[inside, :].shape[0])
            self.add_constraints_to_least_squares(A, B, idc[inside, :], w=w, name="gradient")
            A = np.einsum("ij,ijk->ik", dip_vector.T, T)
            self.add_constraints_to_least_squares(A, B, idc[inside, :], w=w, name="gradient")
            if np.sum(inside) <= 0:
                logger.warning(
                    f" {np.sum(~inside)} \
                        norm constraints not added: outside of model bounding box"
                )

    def add_norm_constraints(self, w=1.0):
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
            node_idx, inside = self.support.position_to_cell_corners(points[:, :3])
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
            (
                vertices,
                T,
                elements,
                inside_,
            ) = self.support.get_element_gradient_for_location(points[inside, :3])
            # T*=np.product(self.support.step_vector)
            # T/=self.support.step_vector[0]
            w /= 3
            self.add_constraints_to_least_squares(
                T[:, 0, :],
                points[inside, 3],
                idc[inside, :],
                w=w,
                name="norm",
            )
            self.add_constraints_to_least_squares(
                T[:, 1, :],
                points[inside, 4],
                idc[inside, :],
                w=w,
                name="norm",
            )
            self.add_constraints_to_least_squares(
                T[:, 2, :],
                points[inside, 5],
                idc[inside, :],
                w=w,
                name="norm",
            )
            if np.sum(inside) <= 0:
                logger.warning(
                    f"{np.sum(~inside)} \
                        norm constraints not added: outside of model bounding box"
                )
            self.up_to_date = False

    def add_gradient_orthogonal_constraints(
        self,
        points: np.ndarray,
        vectors: np.ndarray,
        w: float = 1.0,
        B: float = 0,
        name="gradient orthogonal",
    ):
        """
        constraints scalar field to be orthogonal to a given vector

        Parameters
        ----------
        points : np.darray
            location to add gradient orthogonal constraint
        vector : np.darray
            vector to be orthogonal to, should be the same shape as points
        w : double
        B : np.array

        Returns
        -------

        """
        if points.shape[0] > 0:

            # calculate unit vector for orientation data
            node_idx, inside = self.support.position_to_cell_corners(points[:, :3])
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
            # normalise vector and scale element gradient matrix by norm as well
            norm = np.linalg.norm(vectors, axis=1)
            vectors[norm > 0, :] /= norm[norm > 0, None]

            # normalise element vector to unit vector for dot product
            (
                vertices,
                T,
                elements,
                inside_,
            ) = self.support.get_element_gradient_for_location(points[inside, :3])
            T[norm > 0, :, :] /= norm[norm > 0, None, None]

            # dot product of vector and element gradient = 0
            A = np.einsum("ij,ijk->ik", vectors[inside, :3], T)
            B = np.zeros(points[inside, :].shape[0]) + B
            self.add_constraints_to_least_squares(A, B, idc[inside, :], w=w, name=name)

            if np.sum(inside) <= 0:
                logger.warning(
                    f"{np.sum(~inside)} \
                        gradient constraints not added: outside of model bounding box"
                )
            self.up_to_date = False

    def add_regularisation(self, operator, w=0.1):
        """

        Parameters
        ----------
        operator
        w

        Returns
        -------

        """
        self.assemble_inner(operator, w)
        # self.assemble_borders()

    # def assemble_borders(self, operator, w):

    def assemble_inner(self, operator, w, name='regularisation'):
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
        inside = ~np.any(idc == -1, axis=1)  # np.ones(a.shape[0],dtype=bool)#
        # a[idc==-1] = 0
        # idc[idc==-1] = 0
        B = np.zeros(global_indexes.shape[1])
        self.add_constraints_to_least_squares(
            a[inside, :],
            B[inside],
            idc[inside, :],
            w=w,
            name=name,
        )
        return
