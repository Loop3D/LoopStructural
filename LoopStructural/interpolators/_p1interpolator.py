"""
Piecewise linear interpolator
"""

import logging

import numpy as np


from ._discrete_interpolator import DiscreteInterpolator
from . import InterpolatorType
logger = logging.getLogger(__name__)


class P1Interpolator(DiscreteInterpolator):
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
        self.type = InterpolatorType.PIECEWISE_LINEAR
    def add_gradient_constraints(self, w=1.0):
        pass

    def add_norm_constraints(self, w=1.0):
        points = self.get_norm_constraints()
        if points.shape[0] > 0:
            grad, elements, inside = self.support.evaluate_shape_derivatives(points[:, :3])
            size = self.support.element_scale[elements[inside]]
            wt = np.ones(size.shape[0])
            wt *= w  # s* size
            elements = np.tile(self.support.elements[elements[inside]], (3, 1, 1))

            elements = elements.swapaxes(0, 1)
            # elements = elements.swapaxes(0, 2)
            # grad = grad.swapaxes(1, 2)
            # elements = elements.swapaxes(1, 2)

            self.add_constraints_to_least_squares(
                grad[inside, :, :],
                points[inside, 3:6],
                elements,
                w=wt,
                name="norm",
            )
            self.up_to_date = False
        pass

    def add_value_constraints(self, w=1.0):
        points = self.get_value_constraints()
        if points.shape[0] > 1:
            N, elements, inside = self.support.evaluate_shape(points[:, :3])
            size = self.support.element_size[elements[inside]]

            wt = np.ones(size.shape[0])
            wt *= w  # * size
            self.add_constraints_to_least_squares(
                N[inside, :],
                points[inside, 3],
                self.support.elements[elements[inside], :],
                w=wt,
                name="value",
            )
            self.up_to_date = False

    def minimise_edge_jumps(self, w=0.1, vector_func=None, vector=None, name="edge jump"):
        # NOTE: imposes \phi_T1(xi)-\phi_T2(xi) dot n =0
        # iterate over all triangles
        # flag inidicate which triangles have had all their relationships added
        v1 = self.support.nodes[self.support.shared_elements][:, 0, :]
        v2 = self.support.nodes[self.support.shared_elements][:, 1, :]
        bc_t1 = self.support.barycentre[self.support.shared_element_relationships[:, 0]]
        bc_t2 = self.support.barycentre[self.support.shared_element_relationships[:, 1]]
        norm = self.support.shared_element_norm
        # shared_element_scale = self.support.shared_element_scale

        # evaluate normal if using vector func for cp2
        if vector_func:
            norm = vector_func((v1 + v2) / 2)
        if vector is not None:
            if bc_t1.shape[0] == vector.shape[0]:
                norm = vector
        # evaluate the shape function for the edges for each neighbouring triangle
        Dt, tri1, inside = self.support.evaluate_shape_derivatives(
            bc_t1, elements=self.support.shared_element_relationships[:, 0]
        )
        Dn, tri2, inside = self.support.evaluate_shape_derivatives(
            bc_t2, elements=self.support.shared_element_relationships[:, 1]
        )
        # constraint for each cp is triangle - neighbour create a Nx12 matrix
        const_t = np.einsum("ij,ijk->ik", norm, Dt)
        const_n = -np.einsum("ij,ijk->ik", norm, Dn)
        # const_t_cp2 = np.einsum('ij,ikj->ik',normal,cp2_Dt)
        # const_n_cp2 = -np.einsum('ij,ikj->ik',normal,cp2_Dn)
        # shared_element_size = self.support.shared_element_size
        # const_t /= shared_element_size[:, None]  # normalise by element size
        # const_n /= shared_element_size[:, None]  # normalise by element size
        const = np.hstack([const_t, const_n])

        # get vertex indexes
        tri_cp1 = np.hstack([self.support.elements[tri1], self.support.elements[tri2]])
        # tri_cp2 = np.hstack([self.support.elements[cp2_tri1],self.support.elements[tri2]])
        # add cp1 and cp2 to the least squares system

        self.add_constraints_to_least_squares(
            const,
            np.zeros(const.shape[0]),
            tri_cp1,
            w=w,
            name=name,
        )
        self.up_to_date = False
        # p2.add_constraints_to_least_squares(const_cp2*e_len[:,None]*w,np.zeros(const_cp1.shape[0]),tri_cp2, name='edge jump cp2')

    def setup_interpolator(self, **kwargs):
        """
        Searches through kwargs for any interpolation weights and updates
        the dictionary.
        Then adds the constraints to the linear system using the
        interpolation weights values
        Parameters
        ----------
        kwargs -
            interpolation weights

        Returns
        -------

        """
        # can't reset here, clears fold constraints
        self.reset()
        for key in kwargs:
            if "regularisation" in kwargs:
                self.interpolation_weights["cgw"] = kwargs["regularisation"]
            self.up_to_date = False
            self.interpolation_weights[key] = kwargs[key]
        if self.interpolation_weights["cgw"] > 0.0:
            self.up_to_date = False
            self.minimise_edge_jumps(self.interpolation_weights["cgw"])
            #     direction_feature=kwargs.get("direction_feature", None),
            #     direction_vector=kwargs.get("direction_vector", None),
            # )
            # self.minimise_grad_steepness(
            #     w=self.interpolation_weights.get("steepness_weight", 0.01),
            #     wtfunc=self.interpolation_weights.get("steepness_wtfunc", None),
            # )
            logger.info(
                "Using constant gradient regularisation w = %f" % self.interpolation_weights["cgw"]
            )

        logger.info(
            "Added %i gradient constraints, %i normal constraints,"
            "%i tangent constraints and %i value constraints"
            % (self.n_g, self.n_n, self.n_t, self.n_i)
        )
        self.add_gradient_constraints(self.interpolation_weights["gpw"])
        self.add_norm_constraints(self.interpolation_weights["npw"])
        self.add_value_constraints(self.interpolation_weights["cpw"])
        self.add_tangent_constraints(self.interpolation_weights["tpw"])
        self.add_value_inequality_constraints()
        self.add_inequality_pairs_constraints()
        # self.add_interface_constraints(self.interpolation_weights["ipw"])

    def add_gradient_orthogonal_constraints(
        self,
        points: np.ndarray,
        vectors: np.ndarray,
        w: float = 1.0,
        b: float = 0,
        name='undefined gradient orthogonal constraint',
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
            grad, elements, inside = self.support.evaluate_shape_derivatives(points[:, :3])
            size = self.support.element_size[elements[inside]]
            wt = np.ones(size.shape[0])
            wt *= w * size
            elements = self.support.elements[elements[inside], :]
            # elements = np.tile(self.support.elements[elements[inside]], (3, 1, 1))

            # elements = elements.swapaxes(0, 1)
            # elements = elements.swapaxes(0, 2)
            # grad = grad.swapaxes(1, 2)
            # elements = elements.swapaxes(1, 2)
            norm = np.linalg.norm(vectors, axis=1)
            vectors[norm > 0, :] /= norm[norm > 0, None]
            A = np.einsum("ij,ijk->ik", vectors[inside, :3], grad[inside, :, :])
            B = np.zeros(points[inside, :].shape[0]) + b
            self.add_constraints_to_least_squares(A, B, elements, w=wt, name=name)
            if np.sum(inside) <= 0:
                logger.warning(
                    f"{np.sum(~inside)} \
                        gradient constraints not added: outside of model bounding box"
                )
            self.up_to_date = False

    def add_interface_constraints(self, w: float = 1):
        raise NotImplementedError
