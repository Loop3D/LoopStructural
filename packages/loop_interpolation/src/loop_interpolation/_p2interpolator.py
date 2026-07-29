"""
Piecewise quadratic interpolator
"""

import logging
from typing import Optional, Callable

import numpy as np

from ._discrete_interpolator import DiscreteInterpolator
from . import InterpolatorType

logger = logging.getLogger(__name__)


class P2Interpolator(DiscreteInterpolator):
    """ """

    def __init__(self, mesh):
        """
        Piecewise Quadratic Interpolator
        Approximates scalar field by finding coefficients to a piecewise quadratic
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
        self.type = InterpolatorType.PIECEWISE_QUADRATIC

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
        self.reset()
        regularisation_config = self.resolve_regularisation_config(
            regularisation=kwargs.get("regularisation", None),
            directional_regularisation=kwargs.get("directional_regularisation", None),
        )
        self._apply_isotropic_regularisation_weight(
            regularisation_config.isotropic,
            ("cgw",),
        )
        self._apply_interpolation_weight_kwargs(
            kwargs,
            skip_keys=("regularisation", "directional_regularisation"),
        )

        if self.interpolation_weights["cgw"] > 0.0:
            self.up_to_date = False
            self.minimise_edge_jumps(self.interpolation_weights["cgw"])
            self.minimise_grad_steepness(
                w=self.interpolation_weights.get("steepness_weight", 0.01),
                wtfunc=self.interpolation_weights.get("steepness_wtfunc", None),
            )
            logger.info(
                "Using constant gradient regularisation w = %f" % self.interpolation_weights["cgw"]
            )
        self.add_directional_regularisation(regularisation_config.directional)

        logger.info(
            "Added %i gradient constraints, %i normal constraints,"
            "%i tangent constraints and %i value constraints"
            % (self.n_g, self.n_n, self.n_t, self.n_i)
        )
        self.add_gradient_constraints(self.interpolation_weights["gpw"])
        self.add_norm_constraints(self.interpolation_weights["npw"])
        self.add_value_constraints(self.interpolation_weights["cpw"])
        self.add_tangent_constraints(self.interpolation_weights["tpw"])
        # self.add_interface_constraints(self.interpolation_weights["ipw"])
        return self.finalize_setup_diagnostics_report()

    def copy(self):
        return P2Interpolator(self.support)

    def add_gradient_constraints(self, w: float = 1.0):
        points = self.get_gradient_constraints()
        if points.shape[0] > 0:
            grad, elements = self.support.evaluate_shape_derivatives(points[:, : self.dimensions])
            inside = elements > -1
            area = self.support.element_size[elements[inside]]
            wt = np.ones(area.shape[0])
            wt *= w * area
            A = np.einsum(
                "ikj,ij->ik",
                grad[inside, :],
                points[inside, self.dimensions : self.dimensions * 2],
            )
            B = np.zeros(A.shape[0])
            elements = self.support.elements[elements[inside]]
            self.add_constraints_to_least_squares(A * wt[:, None], B, elements, name="gradient")

    def add_gradient_orthogonal_constraints(
        self, points: np.ndarray, vector: np.ndarray, w=1.0, B=0
    ):
        """
        Constraints scalar field to be orthogonal to a given vector

        Parameters
        ----------
        position
        normals
        w
        B

        Returns
        -------

        """
        if points.shape[0] > 0:
            grad, elements = self.support.evaluate_shape_derivatives(points[:, :3])
            inside = elements > -1
            area = self.support.element_size[elements[inside]]
            wt = np.ones(area.shape[0])
            wt *= w * area
            A = np.einsum("ijk,ij->ik", grad[inside, :], vector[inside, :])
            B = np.zeros(A.shape[0])
            elements = self.support.elements[elements[inside]]
            self.add_constraints_to_least_squares(
                A * wt[:, None], B, elements, name="gradient orthogonal"
            )

    def add_norm_constraints(self, w: float = 1.0):
        points = self.get_norm_constraints()
        if points.shape[0] > 0:
            grad, elements = self.support.evaluate_shape_derivatives(points[:, : self.dimensions])
            inside = elements > -1
            area = self.support.element_size[elements[inside]]
            wt = np.ones(area.shape[0])
            wt *= w * area
            elements = np.tile(self.support.elements[elements[inside]], (self.dimensions, 1, 1))
            elements = elements.swapaxes(0, 1)
            self.add_constraints_to_least_squares(
                grad[inside, :, :] * wt[:, None, None],
                points[inside, self.dimensions : self.dimensions * 2] * wt[:, None],
                elements,
                name="norm",
            )

    def add_value_constraints(self, w: float = 1.0):
        points = self.get_value_constraints()
        if points.shape[0] > 0:
            N, elements, mask = self.support.evaluate_shape(points[:, : self.dimensions])
            # mask = elements > 0
            size = self.support.element_size[elements[mask]]
            wt = np.ones(size.shape[0])
            wt *= w
            self.add_constraints_to_least_squares(
                N[mask, :],
                points[mask, self.dimensions],
                self.support.elements[elements[mask], :],
                w=wt,
                name="value",
            )

    def minimise_grad_steepness(
        self,
        w: float = 0.1,
        maskall: bool = False,
        wtfunc: Optional[Callable[[np.ndarray], np.ndarray]] = None,
    ):
        """This constraint minimises the second derivative of the gradient
        mimimising the 2nd derivative should prevent high curvature solutions
        It is not added on the borders

        Parameters
        ----------
        w : float, optional
            [description], by default 0.1
        maskall : bool, default False
            whether to apply on all elements or just internal elements (default)
        wtfunc :  callable, optional
            a function that returns the weight to be applied at xyz. Called on the barycentre
            of the tetrahedron
        """
        elements = np.arange(0, len(self.support.elements))
        mask = np.ones(self.support.neighbours.shape[0], dtype=bool)
        if not maskall:
            mask[:] = np.all(self.support.neighbours > 3, axis=1)

        d2 = self.support.evaluate_shape_d2(elements[mask])
        # d2 shape is [ele_idx, deriv, node]
        wt = np.ones(d2.shape[0])
        wt *= w  # * self.support.element_size[mask]
        if callable(wtfunc):
            logger.info("Using function to weight gradient steepness")
            wt = wtfunc(self.support.barycentre) * self.support.element_size[mask]
        idc = self.support.elements[elements[mask]]
        for i in range(d2.shape[1]):
            self.add_constraints_to_least_squares(
                d2[:, i, :],
                np.zeros(d2.shape[0]),
                idc[:, :],
                w=wt,
                name=f"gradsteepness_{i}",
            )

    def minimise_edge_jumps(
        self,
        w: float = 0.1,
        wtfunc: Optional[Callable[[np.ndarray], np.ndarray]] = None,
        vector_func: Optional[Callable[[np.ndarray], np.ndarray]] = None,
        quadrature_points: Optional[int] = None,
    ):
        """_summary_

        Parameters
        ----------
        w : float, optional
            _description_, by default 0.1
        wtfunc : callable, optional
            _description_, by default None
        vector_func : callable, optional
            _description_, by default None
        """
        # NOTE: imposes \phi_T1(xi)-\phi_T2(xi) dot n =0
        # iterate over all triangles

        get_qp_kwargs = {} if quadrature_points is None else {"npts": quadrature_points}
        cp, weight = self.support.get_quadrature_points(**get_qp_kwargs)

        norm = self.support.shared_element_norm

        # evaluate normal if using vector func for cp1
        for i in range(cp.shape[1]):
            if callable(vector_func):
                norm = vector_func(cp[:, i, :])
            # evaluate the shape function for the edges for each neighbouring triangle
            cp_Dt, cp_tri1 = self.support.evaluate_shape_derivatives(
                cp[:, i, :], elements=self.support.shared_element_relationships[:, 0]
            )
            cp_Dn, cp_tri2 = self.support.evaluate_shape_derivatives(
                cp[:, i, :], elements=self.support.shared_element_relationships[:, 1]
            )
            # constraint for each cp is triangle - neighbour create a Nx12 matrix
            const_t_cp = np.einsum("ij,ijk->ik", norm, cp_Dt)
            const_n_cp = -np.einsum("ij,ijk->ik", norm, cp_Dn)

            const_cp = np.hstack([const_t_cp, const_n_cp])
            tri_cp = np.hstack([self.support.elements[cp_tri1], self.support.elements[cp_tri2]])
            wt = np.zeros(tri_cp.shape[0])
            wt[:] = w * weight[:, i]
            if wtfunc:
                wt = wtfunc(tri_cp)
            self.add_constraints_to_least_squares(
                const_cp,
                np.zeros(const_cp.shape[0]),
                tri_cp,
                w=wt,
                name=f"shared element jump cp{i}",
            )

    def evaluate_d2(self, evaluation_points: np.ndarray) -> np.ndarray:
        """Evaluate second derivatives of the interpolant at given points.

        Parameters
        ----------
        evaluation_points : np.ndarray
            Array of shape (n_points, 3) containing point coordinates

        Returns
        -------
        np.ndarray
            Array of shape (n_points, 6) containing second derivatives:
            [d2x, dxdy, d2y, dxdz, dydz, d2z]
        """
        evaluation_points = np.array(evaluation_points)
        mask = np.isnan(evaluation_points).any(axis=1)
        valid_points = evaluation_points[~mask, :]

        if valid_points.shape[0] == 0:
            # Preserve historical shape when all rows are invalid.
            return np.zeros(evaluation_points.shape[0])

        valid_d2 = np.asarray(self.support.evaluate_d2(valid_points, self.c))
        if valid_d2.ndim == 1:
            evaluated = np.zeros(evaluation_points.shape[0])
            evaluated[~mask] = valid_d2
            return evaluated

        evaluated = np.full((evaluation_points.shape[0], valid_d2.shape[1]), np.nan)
        evaluated[~mask] = valid_d2
        return evaluated

    def add_interface_constraints(self, w: float = 1):
        raise NotImplementedError
