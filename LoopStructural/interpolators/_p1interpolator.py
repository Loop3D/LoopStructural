"""
Piecewise linear interpolator
"""
import logging

import numpy as np

from ._discrete_interpolator import DiscreteInterpolator
from LoopStructural.utils.helper import get_vectors

logger = logging.getLogger(__name__)


class P1Interpolator(DiscreteInterpolator):
    """ """

    def __init__(self, mesh):
        """
        Piecewise Linear Interpolator
        Approximates scalar field by finding coefficients to a piecewise linear
        equation on a tetrahedral mesh. Uses constant gradient regularisation.

        Parameters
        ----------
        mesh - TetMesh
            interpolation support
        """

        self.shape = "rectangular"
        DiscreteInterpolator.__init__(self, mesh)
        # whether to assemble a rectangular matrix or a square matrix
        self.support = mesh

        self.interpolation_weights = {
            "cgw": 0.1,
            "cpw": 1.0,
            "npw": 1.0,
            "gpw": 1.0,
            "tpw": 1.0,
            "ipw": 1.0,
        }

    def add_gradient_ctr_pts(self, w=1.0):
        pass

    def add_norm_ctr_pts(self, w=1.0):
        points = self.get_norm_constraints()
        if points.shape[0] > 0:
            grad, elements, inside = self.support.evaluate_shape_derivatives(
                points[:, :3]
            )
            size = self.support.element_size[inside]
            wt = np.ones(size.shape[0])
            wt *= w * size
            # print(grad[inside,:,:].shape)
            # print(self.support.elements[elements[inside]].shape)
            elements = np.tile(self.support.elements[elements[inside]], (3, 1, 1))

            elements = elements.swapaxes(0, 1)
            # elements = elements.swapaxes(0,2)
            grad = grad.swapaxes(1, 2)

            self.add_constraints_to_least_squares(
                grad[inside, :, :] * wt[:, None, None],
                points[inside, 3:5] * wt[:, None],
                elements,
                name="norm",
            )

        pass

    def add_ctr_pts(self, w=1.0):
        points = self.get_value_constraints()
        if points.shape[0] > 1:
            N, elements, inside = self.support.evaluate_shape(points[:, :3])
            size = self.support.element_size[elements[inside]]
            wt = np.ones(size.shape[0])
            wt *= w * size
            self.add_constraints_to_least_squares(
                N[inside, :] * wt[:, None],
                points[inside, 3] * wt,
                self.support.elements[elements[inside], :],
                name="value",
            )

    def minimize_edge_jumps(
        self, w=0.1, vector_func=None
    ):  # NOTE: imposes \phi_T1(xi)-\phi_T2(xi) dot n =0
        # iterate over all triangles
        # flag inidicate which triangles have had all their relationships added
        v1 = self.support.nodes[self.support.shared_elements][:, 0, :]
        v2 = self.support.nodes[self.support.shared_elements][:, 1, :]
        bc_t1 = self.support.barycentre[self.support.shared_element_relationships[:, 0]]
        bc_t2 = self.support.barycentre[self.support.shared_element_relationships[:, 1]]
        norm = self.support.shared_element_norm
        shared_element_size = self.support.shared_element_size

        # evaluate normal if using vector func for cp2
        if vector_func:
            norm = vector_func((v1 + v2) / 2)
        # evaluate the shape function for the edges for each neighbouring triangle
        Dt, tri1 = self.support.evaluate_shape_derivatives(
            bc_t1, elements=self.support.shared_element_relationships[:, 0]
        )
        Dn, tri2 = self.support.evaluate_shape_derivatives(
            bc_t2, elements=self.support.shared_element_relationships[:, 1]
        )

        # constraint for each cp is triangle - neighbour create a Nx12 matrix
        const_t = np.einsum("ij,ijk->ik", norm, Dt)
        const_n = -np.einsum("ij,ijk->ik", norm, Dn)
        # const_t_cp2 = np.einsum('ij,ikj->ik',normal,cp2_Dt)
        # const_n_cp2 = -np.einsum('ij,ikj->ik',normal,cp2_Dn)

        const = np.hstack([const_t, const_n])

        # get vertex indexes
        tri_cp1 = np.hstack([self.support.elements[tri1], self.support.elements[tri2]])
        # tri_cp2 = np.hstack([self.support.elements[cp2_tri1],self.support.elements[tri2]])
        # add cp1 and cp2 to the least squares system
        self.add_constraints_to_least_squares(
            const * shared_element_size[:, None] * w,
            np.zeros(const.shape[0]),
            tri_cp1,
            name="edge jump",
        )
        # p2.add_constraints_to_least_squares(const_cp2*e_len[:,None]*w,np.zeros(const_cp1.shape[0]),tri_cp2, name='edge jump cp2')
