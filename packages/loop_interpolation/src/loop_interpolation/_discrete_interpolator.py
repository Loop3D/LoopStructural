"""
Discrete interpolator base for least squares
"""

from abc import abstractmethod
from collections import defaultdict
from typing import Callable, Optional, Union
import logging
from time import perf_counter

import numpy as np
from scipy import sparse  # import sparse.coo_matrix, sparse.bmat, sparse.eye
from scipy.sparse.linalg import LinearOperator
from ._interpolatortype import InterpolatorType
from ._regularisation import (
    DirectionalRegularisation,
    RegularisationConfig,
    coerce_regularisation_config,
)

from ._diagnostics import ConstraintDiagnosticsReport, ConstraintFamilyDiagnostics
from ._geological_interpolator import GeologicalInterpolator
from . import _solver_strategy
from . import _solver_pipeline
from loop_common.logging import get_logger as getLogger

logger = getLogger(__name__)


class DiscreteInterpolator(GeologicalInterpolator):
    """ """

    def __init__(self, support, data=None, c=None, up_to_date=False):
        """
        Base class for a discrete interpolator e.g. piecewise linear or finite difference which is
        any interpolator that solves the system using least squares approximation

        Parameters
        ----------
        support
            A discrete mesh with, nodes, elements, etc
        """
        GeologicalInterpolator.__init__(self, data=data, up_to_date=up_to_date)
        self.B = []
        self.support = support
        self.dimensions = support.dimension
        self.c = (
            np.array(c)
            if c is not None and np.array(c).shape[0] == self.support.n_nodes
            else np.zeros(self.support.n_nodes)
        )
        # Region masking is no longer supported; all support nodes are active.
        self.region_function = lambda xyz: np.ones(xyz.shape[0], dtype=bool)

        self.shape = "rectangular"
        if self.shape == "square":
            self.B = np.zeros(self.dof)
        self.c_ = 0

        self.solver = "cg"

        self.eq_const_C = []
        self.eq_const_row = []
        self.eq_const_col = []
        self.eq_const_d = []

        self.equal_constraints = {}
        self.eq_const_c = 0
        self.ineq_constraints = {}
        self.ineq_const_c = 0

        self.non_linear_constraints = []
        self.constraints = {}
        self.interpolation_weights = {}
        logger.info("Creating discrete interpolator with {} degrees of freedom".format(self.dof))
        self.type = InterpolatorType.BASE_DISCRETE
        self.apply_scaling_matrix = True
        self.add_ridge_regulatisation = True
        self.ridge_factor = 1e-8
        self.solver_history = None
        self.latest_solve_timing = {}

    def get_last_solve_timing(self) -> dict:
        """Return timing metrics captured during the most recent solve."""
        return dict(self.latest_solve_timing)

    def set_nelements(self, nelements: int) -> int:
        return self.support.set_nelements(nelements)

    @property
    def n_elements(self) -> int:
        """Number of elements in the interpolator

        Returns
        -------
        int
            number of elements, positive
        """
        return self.support.n_elements

    @property
    def dof(self) -> int:
        """Number of degrees of freedom for the interpolator

        Returns
        -------
        int
            number of degrees of freedom, positve
        """
        return int(self.support.n_nodes)

    @property
    def region(self) -> np.ndarray:
        """The active region of the interpolator. A boolean
        mask for all elements that are interpolated

        Returns
        -------
        np.ndarray

        """

        return np.ones(self.support.n_nodes, dtype=bool)

    @property
    def region_map(self):
        return np.arange(self.support.n_nodes, dtype=int)

    def set_region(self, region=None):
        """
        Set the region of the support the interpolator is working on

        Parameters
        ----------
        region - function(position)
            return true when in region, false when out

        Returns
        -------

        """
        # evaluate the region function on the support to determine
        # which nodes are inside update region map and degrees of freedom
        # self.region_function = region
        logger.info(
            "Region masking has been removed; the interpolator always uses all support nodes (%s DOF).",
            self.dof,
        )

    def set_interpolation_weights(self, weights):
        """
        Set the interpolation weights dictionary

        Parameters
        ----------
        weights - dictionary
            Entry of new weights to assign to self.interpolation_weights

        Returns
        -------

        """
        for key in weights:
            self.up_to_date = False
            self.interpolation_weights[key] = weights[key]

    def _apply_isotropic_regularisation_weight(self, value: Optional[float], keys: tuple[str, ...]):
        if value is None:
            return
        for key in keys:
            self.interpolation_weights[key] = value

    def _apply_interpolation_weight_kwargs(self, kwargs: dict, skip_keys: tuple[str, ...] = ()):
        for key, value in kwargs.items():
            self.up_to_date = False
            if key in skip_keys:
                continue
            self.interpolation_weights[key] = value

    def _normalise_constraint_family(self, name: str, inequality: bool = False) -> str:
        name = name.lower()
        if name.startswith("value"):
            return "value"
        if name.startswith("gradient"):
            return "gradient"
        if name.startswith("norm"):
            return "normal"
        if name.startswith("tangent"):
            return "tangent"
        if name.startswith("interface"):
            return "interface"
        if name.startswith("inequality_value"):
            return "inequality_value"
        if name.startswith("inequality_pairs"):
            return "inequality_pairs"
        if name.startswith("fold"):
            return "fold"
        if name.startswith("d") and len(name) in (3, 4, 5, 6, 7, 8):
            return "regularisation"
        if inequality:
            return "inequality"
        return "other"

    def _build_constraint_diagnostics_report(self) -> ConstraintDiagnosticsReport:
        report = super()._build_constraint_diagnostics_report()
        outside = dict(report.outside_model_points)

        raw_points = {
            "value": int(self.get_value_constraints().shape[0]),
            "gradient": int(self.get_gradient_constraints().shape[0]),
            "normal": int(self.get_norm_constraints().shape[0]),
            "tangent": int(self.get_tangent_constraints().shape[0]),
            "interface": int(self.get_interface_constraints().shape[0]),
            "inequality_value": int(self.get_inequality_value_constraints().shape[0]),
            "inequality_pairs": int(self.get_inequality_pairs_constraints().shape[0]),
        }

        row_counts = defaultdict(int)
        weight_values = defaultdict(list)

        for name, constraint in self.constraints.items():
            family = self._normalise_constraint_family(name)
            matrix = constraint.get("matrix")
            if matrix is not None:
                row_counts[family] += int(matrix.shape[0])
            w = constraint.get("w")
            if w is not None:
                w_arr = np.asarray(w, dtype=float)
                if w_arr.size > 0:
                    weight_values[family].append(w_arr)

        for name, constraint in self.ineq_constraints.items():
            family = self._normalise_constraint_family(name, inequality=True)
            matrix = constraint.get("matrix")
            if matrix is not None:
                row_counts[family] += int(matrix.shape[0])

        families = {}
        family_names = set(raw_points.keys()) | set(row_counts.keys()) | set(report.families.keys())

        for family_name in sorted(family_names):
            rows = int(row_counts.get(family_name, 0))
            source_points = int(raw_points.get(family_name, 0))
            outside_points = int(outside.get(family_name, 0))

            dropped = None
            if family_name in ("value", "normal", "tangent", "inequality_value"):
                dropped = max(source_points - rows, 0)
            elif family_name == "gradient":
                dropped = max(source_points * 2 - rows, 0)

            if len(weight_values.get(family_name, [])) > 0:
                all_weights = np.concatenate(weight_values[family_name])
                mean_w = float(np.mean(all_weights))
                min_w = float(np.min(all_weights))
                max_w = float(np.max(all_weights))
            else:
                mean_w = None
                min_w = None
                max_w = None

            families[family_name] = ConstraintFamilyDiagnostics(
                name=family_name,
                active=rows > 0,
                row_count=rows,
                dropped_rows=dropped,
                effective_weight_mean=mean_w,
                effective_weight_min=min_w,
                effective_weight_max=max_w,
                source_point_count=source_points,
                outside_model_point_count=outside_points,
            )

        return ConstraintDiagnosticsReport(
            interpolator_type=report.interpolator_type,
            families=families,
            region_coverage=report.region_coverage,
            outside_model_points=outside,
        )

    def finalize_setup_diagnostics_report(self) -> ConstraintDiagnosticsReport:
        self.latest_diagnostics_report = self._build_constraint_diagnostics_report()
        return self.latest_diagnostics_report

    def _pre_solve(self):
        """
        Pre solve function to be run before solving the interpolation
        """
        self.c = np.zeros(self.support.n_nodes)
        self.c[:] = np.nan
        return True

    def _post_solve(self):
        """Post solve function(s) to be run after the solver has been called"""
        self.clear_constraints()
        return True

    def clear_constraints(self):
        """
        Clear the constraints from the interpolator, this makes sure we are not storing
        the constraints after the solver has been run
        """
        self.constraints = {}
        self.ineq_constraints = {}
        self.equal_constraints = {}

    def reset(self):
        """
        Reset the interpolation constraints

        """
        self.constraints = {}
        self.c_ = 0
        self.regularisation_scale = np.ones(self.dof)
        logger.info("Resetting interpolation constraints")

    def add_constraints_to_least_squares(self, A, B, idc, w=1.0, name="undefined"):
        """
        Adds constraints to the least squares system. Automatically works
        out the row
        index given the shape of the input arrays

        Parameters
        ----------
        A : numpy array / list
            RxC numpy array of constraints where C is number of columns,R rows
        B : numpy array /list
            B values array length R
        idc : numpy array/list
            RxC column index

        Returns
        -------
        list of constraint ids

        """
        A = np.array(A)
        B = np.array(B)
        idc = np.array(idc)
        n_rows = A.shape[0]
        # logger.debug('Adding constraints to interpolator: {} {} {}'.format(A.shape[0]))
        # print(A.shape,B.shape,idc.shape)
        if A.shape != idc.shape:
            logger.error(f"Cannot add constraints: A and indexes have different shape : {name}")
            return

        if len(A.shape) > 2:
            n_rows = A.shape[0] * A.shape[1]
            if isinstance(w, np.ndarray):
                w = np.tile(w, (A.shape[1]))
            A = A.reshape((A.shape[0] * A.shape[1], A.shape[2]))
            idc = idc.reshape((idc.shape[0] * idc.shape[1], idc.shape[2]))
            B = B.reshape((A.shape[0]))
            # w = w.reshape((A.shape[0]))

        # Check for nan before any row normalisation, which would otherwise
        # zero out nan rows in A and mask this check further down.
        if np.any(np.isnan(idc)) or np.any(np.isnan(A)) or np.any(np.isnan(B)):
            logger.warning("Constraints contain nan not adding constraints: {}".format(name))
            return

        # normalise by rows of A
        # Should this be done? It should make the solution more stable
        length = np.linalg.norm(A, axis=1)
        # length[length>0] = 1.
        B[length > 0] /= length[length > 0]
        A[length > 0, :] /= length[length > 0, None]
        if isinstance(w, (float, int)):
            w = np.ones(A.shape[0]) * w
        if not isinstance(w, np.ndarray):
            raise BaseException("w must be a numpy array")

        if w.shape[0] != A.shape[0]:
            raise BaseException("Weight array does not match number of constraints")
        rows = np.arange(0, n_rows).astype(int)
        base_name = name
        while name in self.constraints:
            count = 0
            if "_" in name:
                count = int(name.split("_")[1]) + 1
            name = base_name + "_{}".format(count)

        rows = np.tile(rows, (A.shape[-1], 1)).T
        self.constraints[name] = {
            "matrix": sparse.coo_matrix(
                (A.flatten(), (rows.flatten(), idc.flatten())), shape=(n_rows, self.dof)
            ).tocsc(),
            "b": B.flatten(),
            "w": w,
        }

    @abstractmethod
    def add_gradient_orthogonal_constraints(
        self, points: np.ndarray, vectors: np.ndarray, w: float = 1.0
    ):
        pass

    def get_regularisation_sample_points(self) -> np.ndarray:
        raise NotImplementedError(
            f"{self.__class__.__name__} does not define regularisation sample points"
        )

    def _add_directional_regularisation(
        self,
        weight: float,
        vectors: np.ndarray,
        name: str = "directional regularisation",
    ):
        raise NotImplementedError(
            f"{self.__class__.__name__} does not implement directional regularisation"
        )

    def resolve_regularisation_config(
        self,
        regularisation=None,
        directional_regularisation=None,
    ) -> RegularisationConfig:
        return coerce_regularisation_config(
            regularisation=regularisation,
            directional_regularisation=directional_regularisation,
        )

    def add_directional_regularisation(
        self,
        directional_regularisation,
    ) -> tuple[DirectionalRegularisation, ...]:
        directional_terms = self.resolve_regularisation_config(
            directional_regularisation=directional_regularisation
        ).directional
        if len(directional_terms) == 0:
            return directional_terms

        sample_points = np.asarray(self.get_regularisation_sample_points(), dtype=float)
        expected_shape = sample_points.shape

        for term in directional_terms:
            if term.weight == 0:
                continue

            vectors = term.direction(sample_points) if callable(term.direction) else term.direction
            vectors = np.asarray(vectors, dtype=float)
            if vectors.shape != expected_shape:
                logger.warning(
                    "%s: directional regularisation vectors must have shape %s, got %s. Skipping.",
                    term.name,
                    expected_shape,
                    vectors.shape,
                )
                continue
            self._add_directional_regularisation(term.weight, vectors, name=term.name)

        return directional_terms

    def calculate_residual_for_constraints(self):
        """Calculates Ax-B for all constraints added to the interpolator
        This could be a proxy to identify which constraints are controlling the model

        Returns
        -------
        np.ndarray
            vector of Ax-B
        """
        residuals = {}
        for constraint_name, constraint in self.constraints.items():
            residuals[constraint_name] = constraint["matrix"].dot(self.c) - constraint["b"].flatten()
        return residuals

    def add_inequality_constraints_to_matrix(
        self, A: np.ndarray, bounds: np.ndarray, idc: np.ndarray, name: str = "undefined"
    ):
        """Adds constraints for a matrix where the linear function
        l < Ax > u constrains the objective function


        Parameters
        ----------
        A : numpy array
            matrix of coefficients
        bounds : numpy array
            n*3 lower, upper, 1
        idc : numpy array
            index of constraints in the matrix
        Returns
        -------

        """
        A = np.asarray(A, dtype=float)
        idc = np.asarray(idc, dtype=int)
        bounds = np.asarray(bounds, dtype=float)

        if A.ndim != 2 or idc.ndim != 2:
            raise ValueError("A and idc must be 2D arrays")
        if A.shape != idc.shape:
            raise ValueError("A and idc must have the same shape")
        if bounds.ndim != 2 or bounds.shape[0] != A.shape[0] or bounds.shape[1] not in (2, 3):
            raise ValueError("bounds must have shape (n_rows, 2) or (n_rows, 3)")
        if np.any(idc < 0) or np.any(idc >= self.dof):
            raise ValueError("Inequality constraint indices are out of range")

        rows = np.arange(0, idc.shape[0])
        rows = np.tile(rows, (A.shape[-1], 1)).T

        self.ineq_constraints[name] = {
            "matrix": sparse.coo_matrix(
                (A.flatten(), (rows.flatten(), idc.flatten())), shape=(rows.shape[0], self.dof)
            ).tocsc(),
            "bounds": bounds,
        }

    def add_value_inequality_constraints(self, w: float = 1.0):
        points = self.get_inequality_value_constraints()
        # check that we have added some points
        if points.shape[0] > 0:
            coords = points[:, : self.support.dimension]
            vertices, a, element, inside = self.support.get_element_for_location(coords)
            a = a[inside]
            cols = self.support.elements[element[inside]]
            bounds = points[inside, self.support.dimension : self.support.dimension + 2]
            self.add_inequality_constraints_to_matrix(a, bounds, cols, "inequality_value")

    def add_inequality_pairs_constraints(
        self,
        w: float = 1.0,
        upper_bound=1.0,  # np.finfo(float).eps,
        lower_bound=-np.inf,
        pairs: Optional[list] = None,
    ):

        points = self.get_inequality_pairs_constraints()
        if points.shape[0] > 0:
            # assemble a list of pairs in the model
            # this will make pairs even across stratigraphic boundaries
            # TODO add option to only add stratigraphic pairs
            if not pairs:
                pairs = {}
                k = 0
                for i in np.unique(points[:, self.support.dimension]):
                    for j in np.unique(points[:, self.support.dimension]):
                        if i == j:
                            continue
                        if tuple(sorted([i, j])) not in pairs:
                            pairs[tuple(sorted([i, j]))] = k
                            k += 1
                pairs = list(pairs.keys())
            for pair in pairs:
                upper_points = points[points[:, self.support.dimension] == pair[0]]
                lower_points = points[points[:, self.support.dimension] == pair[1]]

                upper_coords = upper_points[:, : self.support.dimension]
                lower_coords = lower_points[:, : self.support.dimension]
                upper_interpolation = self.support.get_element_for_location(upper_coords)
                lower_interpolation = self.support.get_element_for_location(lower_coords)
                if (~upper_interpolation[3]).sum() > 0:
                    logger.warning(
                        f"Upper points not in mesh {upper_points[~upper_interpolation[3]]}"
                    )
                if (~lower_interpolation[3]).sum() > 0:
                    logger.warning(
                        f"Lower points not in mesh {lower_points[~lower_interpolation[3]]}"
                    )
                ij = np.array(
                    [
                        *np.meshgrid(
                            np.arange(0, int(upper_interpolation[3].sum()), dtype=int),
                            np.arange(0, int(lower_interpolation[3].sum()), dtype=int),
                        )
                    ],
                    dtype=int,
                )

                ij = ij.reshape(2, -1).T
                rows = np.arange(0, ij.shape[0], dtype=int)
                rows = np.tile(rows, (upper_interpolation[1].shape[-1], 1)).T
                rows = np.hstack([rows, rows])
                a = upper_interpolation[1][upper_interpolation[3]][ij[:, 0]]
                a = np.hstack([a, -lower_interpolation[1][lower_interpolation[3]][ij[:, 1]]])
                cols = np.hstack(
                    [
                        self.support.elements[
                            upper_interpolation[2][upper_interpolation[3]][ij[:, 0]]
                        ],
                        self.support.elements[
                            lower_interpolation[2][lower_interpolation[3]][ij[:, 1]]
                        ],
                    ]
                )

                bounds = np.zeros((ij.shape[0], 2))
                bounds[:, 0] = lower_bound
                bounds[:, 1] = upper_bound

                self.add_inequality_constraints_to_matrix(
                    a, bounds, cols, f"inequality_pairs_{pair[0]}_{pair[1]}"
                )

    def add_inequality_feature(
        self,
        feature: Callable[[np.ndarray], np.ndarray],
        lower: bool = True,
        mask: Optional[np.ndarray] = None,
    ):
        """Add an inequality constraint to the interpolator using an existing feature.
        This will make the interpolator greater than or less than the exising feature.
        Evaluate the feature at the interpolation nodes.
        Can provide a boolean mask to restrict to only some parts

        Parameters
        ----------
        feature : BaseFeature
            the feature that will be used to constraint the interpolator
        lower : bool, optional
            lower or upper constraint, by default True
        mask : np.ndarray, optional
            restrict the nodes to evaluate on, by default None
        """
        # add inequality value for the nodes of the mesh
        # flag lower determines whether the feature is a lower bound or upper bound
        # mask is just a boolean array determining which nodes to apply it to

        value = feature(self.support.nodes)
        if mask is None:
            mask = np.ones(value.shape[0], dtype=bool)
        l = np.zeros(value.shape[0]) - np.inf
        u = np.zeros(value.shape[0]) + np.inf
        mask = np.logical_and(mask, ~np.isnan(value))
        if lower:
            l[mask] = value[mask]
        if not lower:
            u[mask] = value[mask]
        bounds = np.column_stack([l, u])

        self.add_inequality_constraints_to_matrix(
            np.ones((value.shape[0], 1)),
            bounds,
            np.arange(0, self.dof, dtype=int)[:, None],
            name="inequality_feature",
        )

    def add_equality_constraints(self, node_idx, values, name="undefined"):
        """
        Adds hard constraints to the least squares system. For now this just
        sets
        the node values to be fixed using a lagrangian.

        Parameters
        ----------
        node_idx : numpy array/list
            int array of node indexes
        values : numpy array/list
            array of node values

        Returns
        -------

        """
        idc = np.asarray(node_idx, dtype=int)
        values = np.asarray(values)
        inside = np.logical_and(idc >= 0, idc < self.dof)

        self.equal_constraints[name] = {
            "A": np.ones(idc[inside].shape[0]),
            "B": values[inside],
            "col": idc[inside],
            # "w": w,
            "row": np.arange(self.eq_const_c, self.eq_const_c + idc[inside].shape[0]),
        }
        self.eq_const_c += idc[inside].shape[0]

    def add_tangent_constraints(self, w=1.0):
        """Adds the constraints :math:`f(X)\cdotT=0`

        Parameters
        ----------
        w : double


        Returns
        -------

        """
        points = self.get_tangent_constraints()
        if points.shape[0] > 0:
            self.add_gradient_orthogonal_constraints(points[:, :3], points[:, 3:6], w)

    def _assemble_explicit_constraint_matrix(self) -> tuple[sparse.spmatrix, np.ndarray]:
        """Stack every recorded explicit constraint (``self.constraints``,
        i.e. data constraints and - unless ``regularisation_matrix_free`` is
        active - regularisation constraints too) into a single sparse ``A, B``
        pair. Extracted from ``build_matrix`` so the matrix-free-regularisation
        + ``cg`` fast path (``_solve_with_cg_fused_regularisation``) can reuse
        the data-only assembly without duplicating this loop.
        """
        mats = []
        bs = []
        for c in self.constraints.values():
            if len(c["w"]) == 0:
                continue
            mats.append(c["matrix"].multiply(c["w"][:, None]))
            bs.append(c["b"] * c["w"])
        if len(mats) == 0:
            A = sparse.csr_matrix((0, self.dof), dtype=float)
            B = np.zeros(0, dtype=float)
        else:
            A = sparse.vstack(mats)
            B = np.hstack(bs)
        return A, B

    def build_matrix(self):
        """
        Assemble constraints into interpolation matrix. Adds equaltiy
        constraints
        using lagrange modifiers if necessary

        Parameters
        ----------
        damp: bool
            Flag whether damping should be added to the diagonal of the matrix
        Returns
        -------
        Interpolation matrix and B
        """

        A, B = self._assemble_explicit_constraint_matrix()

        # Matrix-free regularisation (FiniteDifferenceInterpolator opt-in path,
        # see `regularisation_matrix_free` / `get_regularisation_linear_operator`).
        # When active, the structured-grid regularisation stencils were recorded
        # as `matrix_free_regularisation_blocks` instead of being added to
        # `self.constraints`, so `A`/`B` above only contain the explicit data
        # constraints (value/gradient/normal/tangent/interface/... and anything
        # else, e.g. directional regularisation, that was not diverted). Combine
        # them with the matrix-free regularisation `LinearOperator` into a single
        # combined `LinearOperator` so cg/lsmr can solve against the full system
        # without ever materialising the regularisation stencils as a sparse
        # matrix. This branch is a no-op (returns the explicit sparse matrix,
        # exactly as before) unless the flag is set and blocks were recorded.
        use_matrix_free_regularisation = bool(
            getattr(self, "regularisation_matrix_free", False)
        ) and bool(getattr(self, "matrix_free_regularisation_blocks", None))
        if use_matrix_free_regularisation:
            reg_operator = self.get_regularisation_linear_operator()
            if reg_operator is not None:
                n_data = A.shape[0]
                A, B = self._combine_explicit_matrix_with_linear_operator(A, B, reg_operator)
                logger.info(
                    "Interpolation matrix-free system is %d x %d (%d explicit rows + "
                    "%d matrix-free regularisation rows)",
                    A.shape[0],
                    A.shape[1],
                    n_data,
                    reg_operator.shape[0],
                )
                return A, B

        logger.info(f"Interpolation matrix is {A.shape[0]} x {A.shape[1]}")
        return A, B

    def _combine_explicit_matrix_with_linear_operator(
        self,
        A: sparse.spmatrix,
        B: np.ndarray,
        reg_operator: LinearOperator,
    ) -> tuple[LinearOperator, np.ndarray]:
        """Combine an explicit sparse data-constraint matrix with a matrix-free
        regularisation ``LinearOperator`` into a single ``LinearOperator``.

        ``matvec``/``rmatvec`` stack the explicit data rows on top of the
        matrix-free regularisation rows, matching the row order
        ``sparse.vstack`` would otherwise use. The regularisation rows target
        zero (see ``_assemble_operator``'s ``B = np.zeros(...)``), so the
        combined right-hand side is the data ``B`` followed by zeros.
        """
        A = A.tocsr()
        n_data = A.shape[0]
        n_reg = reg_operator.shape[0]
        dof = self.dof

        def matvec(x):
            x = np.asarray(x).reshape(-1)
            return np.concatenate([A @ x, reg_operator.matvec(x)])

        def rmatvec(y):
            y = np.asarray(y).reshape(-1)
            return A.T @ y[:n_data] + reg_operator.rmatvec(y[n_data:])

        combined = LinearOperator(
            shape=(n_data + n_reg, dof), matvec=matvec, rmatvec=rmatvec, dtype=float
        )
        combined_b = np.concatenate([np.asarray(B, dtype=float).reshape(-1), np.zeros(n_reg)])
        return combined, combined_b

    def _materialise_matrix_free_regularisation_blocks(self) -> None:
        """Fallback used when a solver that cannot consume a ``LinearOperator``
        (currently: ``admm``) is selected while ``regularisation_matrix_free=True``.

        Converts every recorded matrix-free regularisation block back into an
        explicit sparse constraint (the same rows ``_assemble_operator`` would
        have built with ``regularisation_matrix_free=False``) and clears the
        matrix-free state so the ordinary explicit ``build_matrix`` path is used
        for this solve. This mirrors the existing ``apply_scaling_matrix``
        fallback pattern in ``FiniteDifferenceInterpolator.setup_interpolator``.
        """
        blocks = getattr(self, "matrix_free_regularisation_blocks", None)
        if not blocks:
            return
        for name, block in blocks.items():
            idc = block["idc"]
            values = block["values"]
            row_w = block["w"]
            a = np.tile(values, (idc.shape[0], 1))
            b = np.zeros(idc.shape[0])
            self.add_constraints_to_least_squares(a, b, idc, w=row_w, name=name)
        self.matrix_free_regularisation_blocks = {}
        self.regularisation_matrix_free = False

    def compute_column_scaling_matrix(self, A: sparse.csr_matrix) -> sparse.dia_matrix:
        """Compute column scaling matrix S for matrix A so that A @ S has columns with unit norm.

        Parameters
        ----------
        A : sparse.csr_matrix
            interpolation matrix

        Returns
        -------
        scipy.sparse.dia_matrix
            diagonal scaling matrix S
        """
        col_norms = sparse.linalg.norm(A, axis=0)
        scaling_factors = np.ones(A.shape[1])
        mask = col_norms > 0
        scaling_factors[mask] = 1.0 / col_norms[mask]
        S = sparse.diags(scaling_factors)
        return S

    def add_equality_block(self, A, B):
        if len(self.equal_constraints) > 0:
            ATA = A.T.dot(A)
            ATB = A.T.dot(B)
            logger.info(f"Equality block is {self.eq_const_c} x {self.dof}")
            # solving constrained least squares using
            # | ATA CT | |c| = b
            # | C   0  | |y|   d
            # where A is the interpoaltion matrix
            # C is the equality constraint matrix
            # b is the interpolation constraints to be honoured
            # in a least squares sense
            # and d are the equality constraints
            # c are the node values and y are the
            # lagrange multipliers#
            a = []
            rows = []
            cols = []
            b = []
            for c in self.equal_constraints.values():
                b.extend((c["B"]).tolist())
                aa = c["A"].flatten()
                mask = aa == 0
                a.extend(aa[~mask].tolist())
                rows.extend(c["row"].flatten()[~mask].tolist())
                cols.extend(c["col"].flatten()[~mask].tolist())

            C = sparse.coo_matrix(
                (np.array(a), (np.array(rows), cols)),
                shape=(self.eq_const_c, self.dof),
                dtype=float,
            ).tocsr()

            d = np.array(b)
            ATA = sparse.bmat([[ATA, C.T], [C, None]])
            ATB = np.hstack([ATB, d])

            return ATA, ATB

    def build_inequality_matrix(self):
        mats = []
        bounds = []
        for c in self.ineq_constraints.values():
            mats.append(c["matrix"])
            bounds.append(c["bounds"])
        if len(mats) == 0:
            return sparse.csr_matrix((0, self.dof), dtype=float), np.zeros((0, 3))
        Q = sparse.vstack(mats)
        bounds = np.vstack(bounds)
        return Q, bounds

    def _remove_constraints_with_prefix(self, prefix: str) -> None:
        to_remove = [name for name in self.constraints if name.startswith(prefix)]
        for name in to_remove:
            self.constraints.pop(name, None)

    def _constant_norm_gradient_data(self):
        support = self.support
        if not hasattr(support, "elements") or not hasattr(support, "barycentre"):
            return None

        try:
            element_indices = np.arange(support.elements.shape[0], dtype=int)
            _, gradient, elements, inside = support.get_element_gradient_for_location(
                support.barycentre[element_indices]
            )
        except Exception as err:
            logger.debug("Unable to build constant-norm gradient rows: %s", err)
            return None

        inside = np.asarray(inside, dtype=bool)
        gradient = np.asarray(gradient, dtype=float)
        if inside.shape[0] != gradient.shape[0]:
            inside = np.ones(gradient.shape[0], dtype=bool)
        if not np.any(inside):
            return None

        element_ids = np.asarray(elements, dtype=int)[inside]
        gradient = gradient[inside]
        idc = np.asarray(support.elements[element_ids], dtype=int)
        values = np.asarray(self.c[idc], dtype=float)
        grad_vec = np.einsum("ijk,ik->ij", gradient, values)
        grad_norm = np.linalg.norm(grad_vec, axis=1)
        valid = np.isfinite(grad_norm) & (grad_norm > 1e-10)
        if not np.any(valid):
            return None

        volume = np.ones(np.sum(valid), dtype=float)
        if hasattr(support, "element_size"):
            raw_volume = np.asarray(support.element_size, dtype=float)
            if raw_volume.ndim == 0:
                volume = np.full(np.sum(valid), float(raw_volume), dtype=float)
            else:
                volume = raw_volume[element_ids][valid]
            volume = np.maximum(volume, 1e-12)

        return {
            "gradient": gradient[valid],
            "grad_vec": grad_vec[valid],
            "grad_norm": grad_norm[valid],
            "idc": idc[valid],
            "volume": volume,
        }

    def _run_constant_norm_polish(
        self,
        solver_kwargs: dict,
        iterations: int,
        base_weight: float,
        target_norm: Optional[float],
    ) -> None:
        if iterations <= 0 or base_weight <= 0.0:
            return

        if not np.all(np.isfinite(self.c)):
            logger.warning("Skipping constant-norm polish because current solution is not finite")
            return

        prefix = "__constant_norm_polish__"
        self._remove_constraints_with_prefix(prefix)

        stable_solver_kwargs = dict(solver_kwargs or {})
        stable_solver_kwargs.pop("x0", None)

        for i in range(int(iterations)):
            grad_data = self._constant_norm_gradient_data()
            if grad_data is None:
                logger.warning("Constant-norm polish stopped: gradient rows unavailable")
                break

            if target_norm is None:
                norm_target = float(np.median(grad_data["grad_norm"]))
            else:
                norm_target = float(target_norm)

            unit_grad = grad_data["grad_vec"] / grad_data["grad_norm"][:, None]
            A = np.einsum("ij,ijk->ik", unit_grad, grad_data["gradient"])
            A = A / grad_data["volume"][:, None]
            b = np.full(A.shape[0], norm_target, dtype=float) / grad_data["volume"]
            iter_weight = float(base_weight) * float(i + 1)

            self.add_constraints_to_least_squares(
                A,
                b,
                grad_data["idc"],
                w=np.full(A.shape[0], iter_weight, dtype=float),
                name=f"{prefix}{i}",
            )

            x0_seed = np.asarray(self.c[self.region], dtype=float).copy()
            polish_solver_kwargs = dict(stable_solver_kwargs)
            polish_solver_kwargs["x0"] = lambda _support, x0=x0_seed: np.array(x0, copy=True)
            polish_solver_kwargs["constant_norm_iterations"] = 0

            solved = self.solve_system("admm", solver_kwargs=polish_solver_kwargs)
            if not solved:
                logger.warning("Constant-norm polish terminated because ADMM solve failed")
                break

            updated = self._constant_norm_gradient_data()
            if updated is not None:
                logger.info(
                    "Constant-norm iter %d: target=%.3e, grad_norm mean=%.3e, std=%.3e",
                    i + 1,
                    norm_target,
                    float(np.mean(updated["grad_norm"])),
                    float(np.std(updated["grad_norm"])),
                )

        self._remove_constraints_with_prefix(prefix)

    def _normalise_solver_choice(
        self,
        solver: Optional[Union[Callable[[sparse.csr_matrix, np.ndarray], np.ndarray], str]],
    ):
        return _solver_strategy.resolve_solver_choice(solver, logger)

    def _assemble_main_system(self, timing: dict) -> tuple[sparse.spmatrix, np.ndarray]:
        return _solver_pipeline.assemble_main_system(self.build_matrix, timing)

    def _preprocess_main_system(
        self,
        A: sparse.spmatrix,
        b: np.ndarray,
        timing: dict,
    ) -> tuple[sparse.spmatrix, np.ndarray, Optional[sparse.spmatrix]]:
        return _solver_pipeline.preprocess_main_system(
            A=A,
            b=b,
            add_ridge_regularisation=self.add_ridge_regulatisation,
            ridge_factor=self.ridge_factor,
            apply_scaling_matrix=self.apply_scaling_matrix,
            compute_column_scaling_matrix_fn=self.compute_column_scaling_matrix,
            logger=logger,
            timing=timing,
        )

    def _assemble_inequality_system(self, timing: dict) -> tuple[sparse.spmatrix, np.ndarray]:
        return _solver_pipeline.assemble_inequality_system(self.build_inequality_matrix, timing)

    def _solve_with_callable(
        self,
        solver: Callable[[sparse.csr_matrix, np.ndarray], np.ndarray],
        A: sparse.spmatrix,
        b: np.ndarray,
        timing: dict,
    ) -> bool:
        self.c, ok = _solver_strategy.solve_with_callable(solver, A, b, timing, logger)
        return ok

    def _solve_with_cg(
        self,
        A: sparse.spmatrix,
        b: np.ndarray,
        tol: Optional[float],
        solver_kwargs: dict,
        timing: dict,
    ) -> bool:
        self.c, ok = _solver_strategy.solve_with_cg(A, b, tol, solver_kwargs, timing, logger)
        return ok

    def _use_fused_cg_regularisation(self, solver_choice) -> bool:
        """Whether the fused-kernel matrix-free-regularisation CG fast path
        (``_solve_with_cg_fused_regularisation``) should be used for this
        solve, instead of the generic ``build_matrix`` + ``A.T @ A`` path.

        Only ever true for ``solver_choice == "cg"`` with
        ``regularisation_matrix_free=True`` and recorded blocks (i.e. exactly
        the case that otherwise builds the combined rectangular
        ``LinearOperator`` via ``_combine_explicit_matrix_with_linear_operator``
        and squares it generically in ``_solver_strategy.solve_with_cg``).
        Interpolators that never set `matrix_free_regularisation_blocks`
        (P1/P2/anything else) or that don't define the fused-operator builder
        (anything but ``FiniteDifferenceInterpolator``) always take the
        existing, unmodified path.
        """
        return bool(
            solver_choice == "cg"
            and getattr(self, "regularisation_matrix_free", False)
            and getattr(self, "matrix_free_regularisation_blocks", None)
            and not getattr(self, "apply_scaling_matrix", False)
            and hasattr(self, "_build_fused_cg_regularisation_operator")
        )

    def _solve_with_cg_fused_regularisation(
        self,
        tol: Optional[float],
        solver_kwargs: dict,
        timing: dict,
    ) -> bool:
        """CG fast path for ``regularisation_matrix_free=True``: assemble the
        normal-equations system ``(A_data^T A_data + R_reg^T R_reg) x =
        A_data^T b_data`` directly, using the fused single-kernel + boundary-
        corrected regularisation operator
        (``FiniteDifferenceInterpolator._build_fused_cg_regularisation_operator``)
        for the ``R_reg^T R_reg`` term, instead of assembling a combined
        rectangular ``LinearOperator`` (data rows stacked on top of
        regularisation rows) and letting ``scipy.sparse.linalg.cg`` square it
        generically via ``A.T @ A`` (which re-derives the same normal
        equations but pays for a matvec-then-rmatvec pass over the
        regularisation rows every CG iteration instead of one fused
        convolution call).

        Ridge regularisation (``add_ridge_regulatisation``) contributes exactly
        ``ridge_factor**2 * x`` to the normal-equations LHS (and nothing to
        the RHS) for any ``x``; it's added directly to the precomputed data
        Gram matrix's diagonal rather than appended as ``ridge_factor * I``
        rows to ``A_data`` first (as ``_solver_pipeline.preprocess_main_system``
        does for the non-matrix-free case) purely so the one-off
        ``A_data.T @ A_data`` product below stays a small, cheap
        (n_data-row-driven) sparse-sparse product instead of a
        ``(n_data + dof)``-row one; both are exact and mathematically
        equivalent. ``apply_scaling_matrix`` is always forced off before this
        path can be selected (see ``_use_fused_cg_regularisation``), matching
        ``FiniteDifferenceInterpolator.setup_interpolator``'s existing
        column-scaling fallback.

        The data-constraint block's normal-equations contribution
        (``A_data.T @ A_data``) is precomputed ONCE as an explicit sparse
        matrix here (data constraint rows are comparatively few - value/
        gradient/norm/etc. constraints, not the O(dof) regularisation rows -
        so this product is cheap and its fill-in bounded), so each CG
        iteration is a single fast precomputed-sparse-matrix matvec for the
        data term plus the fused regularisation term, instead of two
        separate sparse matvecs (``A_data @ x`` then ``A_data.T @ (...)``)
        every iteration.
        """
        assembly_started = perf_counter()
        A_data, b_data = self._assemble_explicit_constraint_matrix()
        A_data = A_data.tocsr()
        dof = self.dof

        AtA_data = (A_data.T @ A_data).tocsr()
        rhs = np.asarray(A_data.T @ b_data).reshape(-1)
        if self.add_ridge_regulatisation:
            AtA_data = AtA_data + sparse.eye(dof, format="csr") * (self.ridge_factor**2)

        reg_operator = self._build_fused_cg_regularisation_operator()
        timing["assembly_seconds"] = perf_counter() - assembly_started
        timing["matrix_rows"] = int(A_data.shape[0])
        timing["matrix_cols"] = int(A_data.shape[1])
        timing["matrix_nnz"] = int(A_data.nnz)

        preprocess_started = perf_counter()

        def matvec(x):
            x = np.asarray(x, dtype=float).reshape(-1)
            out = AtA_data @ x
            if reg_operator is not None:
                out = out + reg_operator.matvec(x)
            return np.asarray(out).reshape(-1)

        normal_operator = LinearOperator(
            shape=(dof, dof), matvec=matvec, rmatvec=matvec, dtype=float
        )
        timing["preprocess_seconds"] = perf_counter() - preprocess_started

        self.c, ok = _solver_strategy.solve_with_cg_normal_equations(
            normal_operator, rhs, tol, solver_kwargs, timing, logger
        )
        return ok

    def _solve_with_lsmr(
        self,
        A: sparse.spmatrix,
        b: np.ndarray,
        tol: Optional[float],
        solver_kwargs: dict,
        timing: dict,
    ) -> bool:
        self.c, ok = _solver_strategy.solve_with_lsmr(A, b, tol, solver_kwargs, timing, logger)
        return ok

    def _extract_admm_kwargs(self, solver_kwargs: dict, admm_solve) -> tuple[str, dict]:
        return _solver_strategy.extract_admm_kwargs(solver_kwargs, admm_solve)

    def _solve_with_admm(
        self,
        A: sparse.spmatrix,
        b: np.ndarray,
        Q: sparse.spmatrix,
        bounds: np.ndarray,
        solver_kwargs: dict,
        timing: dict,
        constant_norm_iterations: int,
        constant_norm_weight: float,
        constant_norm_target: Optional[float],
    ) -> bool:
        logger.info("Solving using admm")

        constant_norm_solver_kwargs = dict(solver_kwargs)
        if Q is None:
            logger.warning("No inequality constraints, using lsmr")
            return self.solve_system("lsmr", solver_kwargs=solver_kwargs)

        try:
            c, history, ok = _solver_strategy.solve_with_admm(
                A,
                b,
                Q,
                bounds,
                solver_kwargs,
                timing,
                self.support,
                logger,
            )
            if not ok:
                return False
            self.c = c
            self.solver_history = history

            if constant_norm_iterations > 0 and constant_norm_weight > 0.0:
                self._run_constant_norm_polish(
                    solver_kwargs=constant_norm_solver_kwargs,
                    iterations=constant_norm_iterations,
                    base_weight=constant_norm_weight,
                    target_norm=constant_norm_target,
                )
            return True
        except ValueError as e:
            logger.error(f"ADMM solver failed: {e}")
            return False
        except ImportError:
            logger.warning("Cannot import embedded admm solver. Use lsmr or cg")
            return False

    def solve_system(
        self,
        solver: Optional[Union[Callable[[sparse.csr_matrix, np.ndarray], np.ndarray], str]] = None,
        tol: Optional[float] = None,
        solver_kwargs: Optional[dict] = None,
    ) -> bool:
        """
        Main entry point to run the solver and update the node value
        attribute for the
        discreteinterpolator class

        Parameters
        ----------
        solver : string/callable
            solver 'cg' conjugate gradient, 'lsmr' or callable function
        solver_kwargs
            kwargs for solver check scipy documentation for more information

        Returns
        -------
        bool
            True if the interpolation is run

        """
        if not self._pre_solve():
            raise ValueError("Pre solve failed")

        solve_started = perf_counter()
        solver_choice = self._normalise_solver_choice(solver)
        solver_kwargs = dict(solver_kwargs or {})
        constant_norm_iterations, constant_norm_weight, constant_norm_target = (
            _solver_pipeline.extract_constant_norm_options(solver_kwargs, logger)
        )
        timing = _solver_pipeline.init_timing(solver_choice)

        self.solver_history = None

        if solver_choice == "admm" and getattr(self, "regularisation_matrix_free", False) and getattr(
            self, "matrix_free_regularisation_blocks", None
        ):
            logger.warning(
                "regularisation_matrix_free=True is not supported with solver='admm' "
                "(ADMM's inequality-constrained solve needs an explicit sparse system matrix); "
                "falling back to explicit regularisation assembly for this solve."
            )
            self._materialise_matrix_free_regularisation_blocks()

        Q, bounds = self._assemble_inequality_system(timing)
        # Legacy compatibility: ignore deprecated backend switch.
        solver_kwargs.pop("backend", None)

        scaling_matrix = None
        if self._use_fused_cg_regularisation(solver_choice):
            # Fast path: build the CG normal equations directly, using a
            # fused single-kernel + boundary-corrected regularisation
            # operator instead of the generic combined-LinearOperator
            # `build_matrix` + `A.T @ A` route (still used, unchanged, for
            # every other solver/case below). See
            # `_solve_with_cg_fused_regularisation` / `FiniteDifferenceInterpolator.
            # _build_fused_cg_regularisation_operator`.
            self.up_to_date = self._solve_with_cg_fused_regularisation(
                tol, solver_kwargs, timing
            )
        else:
            A, b = self._assemble_main_system(timing)
            A, b, scaling_matrix = self._preprocess_main_system(A, b, timing)

            if callable(solver_choice):
                self.up_to_date = self._solve_with_callable(solver_choice, A, b, timing)
            elif solver_choice == "cg":
                self.up_to_date = self._solve_with_cg(A, b, tol, solver_kwargs, timing)
            elif solver_choice == "lsmr":
                self.up_to_date = self._solve_with_lsmr(A, b, tol, solver_kwargs, timing)
            elif solver_choice == "admm":
                self.up_to_date = self._solve_with_admm(
                    A,
                    b,
                    Q,
                    bounds,
                    solver_kwargs,
                    timing,
                    constant_norm_iterations,
                    constant_norm_weight,
                    constant_norm_target,
                )
            else:
                logger.error(f"Unknown solver {solver_choice}")
                self.up_to_date = False

        # self._post_solve()
        # apply scaling matrix to solution
        if scaling_matrix is not None:
            self.c = scaling_matrix @ self.c
        self.latest_solve_timing = _solver_pipeline.finalize_timing(
            timing=timing,
            solve_started=solve_started,
            up_to_date=bool(self.up_to_date),
        )
        return self.up_to_date

    def update(self) -> bool:
        """
        Check if the solver is up to date, if not rerun interpolation using
        the previously used solver. If the interpolation has not been run
        before it will
        return False

        Returns
        -------
        bool

        """
        if self.solver is None:
            logging.debug("Cannot rerun interpolator")
            return False
        if not self.up_to_date:
            self.setup_interpolator()
            self.up_to_date = self.solve_system(self.solver)
            return self.up_to_date
        return bool(self.up_to_date)

    def evaluate_value(self, locations: np.ndarray) -> np.ndarray:
        """Evaluate the value of the interpolator at location

        Parameters
        ----------
        evaluation_points : np.ndarray
            location to evaluate the interpolator

        Returns
        -------
        np.ndarray
            value of the interpolator
        """
        self.update()
        evaluation_points = np.array(locations)
        evaluated = np.zeros(evaluation_points.shape[0])
        mask = np.any(np.isnan(evaluation_points), axis=1)

        if evaluation_points[~mask, :].shape[0] > 0:
            evaluated[~mask] = self.support.evaluate_value(evaluation_points[~mask], self.c)
        return evaluated

    def evaluate_gradient(self, locations: np.ndarray) -> np.ndarray:
        """
        Evaluate the gradient of the scalar field at the evaluation points
        Parameters
        ----------
        evaluation_points : np.array
            xyz locations to evaluate the gradient

        Returns
        -------

        """
        self.update()
        if locations.shape[0] > 0:
            return self.support.evaluate_gradient(locations, self.c)
        return np.zeros((0, 3))

    def to_dict(self):
        return {
            "type": self.type.name,
            "support": self.support.to_dict(),
            "c": np.asarray(self.c).tolist(),
            **super().to_dict(),
            # 'region_function':self.region_function,
        }

    def vtk(self):
        if self.up_to_date is False:
            self.update()
        return self.support.vtk({"c": self.c})

    def surfaces(self, value, **kwargs):
        from loop_common.geometry._surface import Surface

        if self.up_to_date is False:
            self.update()
        mesh = self.support.vtk({"c": self.c})
        contour = mesh.contour([value], scalars="c")
        if contour.n_points == 0:
            return []
        triangles = contour.faces.reshape(-1, 4)[:, 1:]
        return [Surface(vertices=contour.points.copy(), triangles=triangles)]
