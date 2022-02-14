"""
Piecewise linear interpolator
"""
import logging

import numpy as np

from LoopStructural.interpolators import DiscreteInterpolator
from LoopStructural.utils.helper import get_vectors

logger = logging.getLogger(__name__)


class P2Interpolator(DiscreteInterpolator):
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
        self.interpolator_type = "P2"
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
        points = self.get_gradient_constraints()
        if points.shape[0] > 0:
            raise NotImplementedError
            
            # grad, elements = self.support.evaluate_shape_derivatives(points[:, :2])
            # inside = elements > -1
            # area = self.support.element_area(elements[inside])
            # wt = np.ones(area.shape[0])
            # wt *= w * area
            # A = np.einsum("ikj,ij->ik", grad[:, :], points[:, 3:5])
            # # A = np.einsum('ij,ikj->ik', grad[:,:], points[:,3:5])
            # A *= area[idx, None]
            # # # A *= fold_orientation
            # B = np.zeros(A.shape[0])
            # idc = p2_mesh.elements[idx, :]
            # self.add_constraints_to_least_squares(A, B, elements)
        # p2.add_constraints_to_least_squares(A, B, idc,'fold ori')
        pass

    def add_norm_ctr_pts(self, w=1.0):
        points = self.get_norm_constraints()
        if points.shape[0] > 0:
            grad, elements = self.support.evaluate_shape_derivatives(points[:, :2])
            inside = elements > -1
            area = self.support.element_area(elements[inside])
            wt = np.ones(area.shape[0])
            wt *= w * area
            # print(grad[inside,:,:].shape)
            # print(self.support.elements[elements[inside]].shape)
            elements = np.tile(self.support.elements[elements[inside]], (2, 1, 1))

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

    def add_gradient_orthogonal_constraint(self, points, vector, w=1.0, B=0):
        """
        constraints scalar field to be orthogonal to a given vector

        Parameters
        ----------
        position
        normals
        w
        B

        Returns
        -------

        """
        pass
        # if points.shape[0] > 0:
        #     vertices, element_gradients, tetras, inside = self.support.get_tetra_gradient_for_location(points[:,:3])
        #     #e, inside = self.support.elements_for_array(points[:, :3])
        #     #nodes = self.support.nodes[self.support.elements[e]]
        #     vector /= np.linalg.norm(vector,axis=1)[:,None]
        #     vecs = vertices[:, 1:, :] - vertices[:, 0, None, :]
        #     vol = np.abs(np.linalg.det(vecs))  # / 6
        #     # d_t = self.support.get_elements_gradients(e)
        #     norm = np.linalg.norm(element_gradients, axis=2)
        #     element_gradients /= norm[:, :, None]

        #     A = np.einsum('ij,ijk->ik', vector, element_gradients)

        #     A *= vol[:, None]

        #     gi = np.zeros(self.support.n_nodes).astype(int)
        #     gi[:] = -1
        #     gi[self.region] = np.arange(0, self.nx).astype(int)
        #     w /= 3
        #     idc = gi[tetras]
        #     B = np.zeros(idc.shape[0])+B
        #     outside = ~np.any(idc == -1, axis=1)
        #     self.add_constraints_to_least_squares(A[outside, :] * w,
        #                                           B[outside], idc[outside, :])

    def add_ctr_pts(self, w=1.0):
        points = self.get_value_constraints()
        if points.shape[0] > 1:
            N, elements, mask = self.support.evaluate_shape(points[:, :3])
            # mask = elements > 0
            size = self.support.element_size[elements[mask]]
            wt = np.ones(size.shape[0])
            wt *= w 
            self.add_constraints_to_least_squares(
                N[mask, :] * wt[:, None],
                points[mask, 3] * wt,
                self.support.elements[elements[mask], :],
                name="value",
            )

    def minimise_grad_steepness(self, stren=0.0, w=0.1, maskall=True, wtfunc=None):
        """This constraint minimises the second derivative of the gradient
        mimimising the 2nd derivative should prevent high curvature solutions
        It is not added on the borders

        Parameters
        ----------
        sten : float, optional
            [description], by default 0.0
        w : float, optional
            [description], by default 0.1
        """
        elements = np.arange(0, len(self.support.elements))
        mask = np.ones(self.support.neighbours.shape[0], dtype=bool)
        if maskall == False:
            mask[:] = np.all(self.support.neighbours > 3, axis=1)
        # vertices = self.support.nodes[self.support.elements[elements[mask], :], :]
        # barycentre = np.mean(vertices, axis=1)
        # M_t = np.ones((vertices.shape[0], 3, 3))
        # M_t[:, :, 1:] = vertices[:, :3, :]
        # area = np.abs(np.linalg.det(M_t)) * 0.5
        d2 = self.support.evaluate_shape_d2(elements[mask])
        # d2 is [ele_idx, deriv, node]
        wt = np.ones(d2.shape[0])
        
        wt *= w #* self.support.element_size[mask]
        if wtfunc:
            wt = wtfunc(self.support.barycentre) * self.support.element_size[mask]
        idc = self.support.elements[elements[mask]]
        for i in range(d2.shape[1]):
            self.add_constraints_to_least_squares(
                d2[:, i, :] * wt[i, None, None],
                np.zeros(d2.shape[0]),
                idc[:, :],
                name=f"grad_steepness_{i}",
            )
 
    def minimize_edge_jumps(
        self, w=0.1, wtfunc=None, vector_func=None
    ):  # NOTE: imposes \phi_T1(xi)-\phi_T2(xi) dot n =0
        # iterate over all triangles
        
        cp, weight = self.support.get_quadrature_points(3)

        norm = self.support.shared_element_norm
        shared_element_size = self.support.shared_element_size
        
        # evaluate normal if using vector func for cp1
        for i in range(cp.shape[1]):
            if callable(vector_func):
                norm = vector_func(cp[:,i,:])
            # evaluate the shape function for the edges for each neighbouring triangle
            cp_Dt, cp_tri1 = self.support.evaluate_shape_derivatives(
                cp[:,i, :], elements=self.support.shared_element_relationships[:, 0]
            )
            cp_Dn, cp_tri2 = self.support.evaluate_shape_derivatives(
                cp[:,i, :], elements=self.support.shared_element_relationships[:, 1]
            )
            # constraint for each cp is triangle - neighbour create a Nx12 matrix
            const_t_cp = np.einsum("ij,ijk->ik", norm, cp_Dt)
            const_n_cp = -np.einsum("ij,ijk->ik", norm, cp_Dn)
           
            const_cp = np.hstack([const_t_cp, const_n_cp])
            tri_cp = np.hstack(
                        [self.support.elements[cp_tri1], self.support.elements[cp_tri2]]
                    )
            wt = np.zeros(tri_cp.shape[0])
            wt[:] = w
            if wtfunc:
                wt = wtfunc(tri_cp)
            self.add_constraints_to_least_squares(
                const_cp * wt[:, None],
                np.zeros(const_cp.shape[0]),
                tri_cp,
                name=f"shared element jump cp{i}",
            )
        

    

    def evaluate_d2(self, evaluation_points):
        evaluation_points = np.array(evaluation_points)
        evaluated = np.zeros(evaluation_points.shape[0])
        mask = np.any(evaluation_points == np.nan, axis=1)

        if evaluation_points[~mask, :].shape[0] > 0:
            evaluated[~mask] = self.support.evaluate_d2(
                evaluation_points[~mask], self.propertyname
            )
        return evaluated
