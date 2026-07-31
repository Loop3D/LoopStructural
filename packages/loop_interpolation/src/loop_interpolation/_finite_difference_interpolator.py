"""
FiniteDifference interpolator
"""
from __future__ import annotations

import numpy as np
from loop_common.logging import get_logger as getLogger
from loop_common.math import get_vectors
from scipy import ndimage, signal
from scipy.sparse.linalg import LinearOperator
from scipy.spatial import KDTree

from ._discrete_interpolator import DiscreteInterpolator
from ._interpolatortype import InterpolatorType
from ._operator import Operator

logger = getLogger(__name__)

# The six interior second-derivative stencil families that qualify for the
# fused single-kernel CG regularisation fast path (see
# `FiniteDifferenceInterpolator._build_fused_cg_regularisation_operator`).
# Keyed by the family name `setup_interpolator`/`support.get_operators` use.
_FUSED_FAMILY_MASKS = {
    "dxx": Operator.Dxx_mask,
    "dyy": Operator.Dyy_mask,
    "dzz": Operator.Dzz_mask,
    "dxy": Operator.Dxy_mask,
    "dxz": Operator.Dxz_mask,
    "dyz": Operator.Dyz_mask,
}


def _normalised_mask_xyz(mask: np.ndarray) -> np.ndarray:
    """(3,3,3) mask in the Operator's native (z,y,x) index order -> normalised
    (unit L2 norm over nonzero taps) mask in (x,y,z) order, i.e. matching
    ``arr.reshape((nx, ny, nz), order='F')`` (axis0=x, axis1=y, axis2=z)."""
    mask = np.asarray(mask, dtype=float)
    active = mask != 0
    norm = np.linalg.norm(mask[active])
    mask_norm = mask / norm if norm > 0 else mask.copy()
    return np.ascontiguousarray(mask_norm.transpose(2, 1, 0))


def _composite_kernel(mask_xyz: np.ndarray) -> np.ndarray:
    """Autocorrelation of a (3,3,3) kernel -> (5,5,5) symmetric composite
    kernel representing "correlate then convolve with the same kernel", i.e.
    one application of ``rmatvec(matvec(x))`` for a single translation-
    invariant stencil family, away from any domain boundary."""
    return signal.correlate(mask_xyz, mask_xyz, mode="full")


def _dist_to_true_edge(nx: int, ny: int, nz: int) -> np.ndarray:
    """For every dof (F-order flat index, ``gi = i + nx*j + nx*ny*k``), the
    minimum number of grid steps to the nearest TRUE domain edge along any
    axis (0 for the outermost node layer, 1 for the next layer in, ...)."""
    xi = np.arange(nx)
    yi = np.arange(ny)
    zi = np.arange(nz)
    dx = np.minimum(xi, nx - 1 - xi)
    dy = np.minimum(yi, ny - 1 - yi)
    dz = np.minimum(zi, nz - 1 - zi)
    d = np.minimum(np.minimum(dx[:, None, None], dy[None, :, None]), dz[None, None, :])
    return d.ravel(order="F")


def compute_weighting(grid_points, gradient_constraint_points, alpha=10.0, sigma=1.0):
    """
    Compute weights for second derivative regularization based on proximity to gradient constraints.

    Parameters:
        grid_points (ndarray): (N, 3) array of 3D coordinates for grid cells.
        gradient_constraint_points (ndarray): (M, 3) array of 3D coordinates for gradient constraints.
        alpha (float): Strength of weighting increase.
        sigma (float): Decay parameter for Gaussian-like influence.

    Returns:
        weights (ndarray): (N,) array of weights for each grid point.
    """
    # Build a KDTree with the gradient constraint locations
    tree = KDTree(gradient_constraint_points)

    # Find the distance from each grid point to the nearest gradient constraint
    distances, _ = tree.query(grid_points, k=1)

    # Compute weighting function (higher weight for nearby points)
    weights = 1 + alpha * np.exp(-(distances**2) / (2 * sigma**2))

    return weights


class FiniteDifferenceInterpolator(DiscreteInterpolator):
    def __init__(self, grid, data=None):
        """
        Finite difference interpolation on a regular cartesian grid

        Parameters
        ----------
        grid : StructuredGrid
        """
        self.shape = "rectangular"
        DiscreteInterpolator.__init__(self, grid, data=data)
        self.set_interpolation_weights(
            {
                "dxy": 1.0,
                "dyz": 1.0,
                "dxz": 1.0,
                "dxx": 1.0,
                "dyy": 1.0,
                "dzz": 1.0,
                "dx": 1.0,
                "dy": 1.0,
                "dz": 1.0,
                "cpw": 1.0,
                "gpw": 1.0,
                "npw": 10.0,
                "tpw": 1.0,
                "ipw": 1.0,
            }
        )

        self.type = InterpolatorType.FINITE_DIFFERENCE
        self.use_regularisation_weight_scale = False
        self.regularisation_weight_sigma = None
        self.interface_pair_mode = "auto"
        self.interface_pair_threshold = 200
        self.interface_pair_fallback = "star"
        # Matrix-free regularisation (opt-in, StructuredGrid mask-based path only).
        # See `_assemble_operator` / `get_regularisation_linear_operator`.
        self.regularisation_matrix_free = False
        self.matrix_free_regularisation_blocks: dict = {}

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
        regularisation_config = self.resolve_regularisation_config(
            regularisation=kwargs.get("regularisation", None),
            directional_regularisation=kwargs.get("directional_regularisation", None),
        )
        self._apply_isotropic_regularisation_weight(
            regularisation_config.isotropic,
            ("dxy", "dyz", "dxz", "dxx", "dyy", "dzz"),
        )
        self._apply_interpolation_weight_kwargs(
            kwargs,
            skip_keys=("regularisation", "directional_regularisation"),
        )
        # either use the default operators or the ones passed to the function
        operators = kwargs.get(
            "operators", self.support.get_operators(weights=self.interpolation_weights)
        )

        self.use_regularisation_weight_scale = kwargs.get("use_regularisation_weight_scale", False)
        self.regularisation_weight_sigma = kwargs.get("regularisation_weight_sigma", None)
        self.matrix_free_regularisation_blocks = {}
        self.regularisation_matrix_free = bool(kwargs.get("regularisation_matrix_free", False))
        if self.regularisation_matrix_free and self.apply_scaling_matrix:
            logger.warning(
                "regularisation_matrix_free=True was requested together with "
                "apply_scaling_matrix=True. compute_column_scaling_matrix() needs an explicit "
                "sparse matrix (it computes per-column norms) and cannot operate on a matrix-free "
                "LinearOperator, so the structured-grid regularisation stencils will fall back to "
                "the explicit assembly path. Set apply_scaling_matrix=False to use the matrix-free "
                "path instead."
            )
            self.regularisation_matrix_free = False
        elif self.regularisation_matrix_free:
            logger.warning(
                "regularisation_matrix_free=True: the structured-grid regularisation stencils "
                "(dxx, dyy, dzz, dxy, dxz, dyz and the border second-derivative terms) will be "
                "assembled as matrix-free LinearOperator blocks (self.matrix_free_regularisation_blocks, "
                "get_regularisation_linear_operator()) instead of explicit sparse matrices. "
                "build_matrix()/solve_system() combine these blocks with the explicit data-constraint "
                "matrix into a single LinearOperator for the 'cg' and 'lsmr' solvers. 'admm' cannot "
                "consume a LinearOperator and will fall back to explicit regularisation assembly "
                "(with a warning) for that solve."
            )
        self.interface_pair_mode = kwargs.get("interface_pair_mode", self.interface_pair_mode)
        self.interface_pair_threshold = int(
            kwargs.get("interface_pair_threshold", self.interface_pair_threshold)
        )
        self.interface_pair_fallback = kwargs.get(
            "interface_pair_fallback", self.interface_pair_fallback
        )
        self.add_norm_constraints(self.interpolation_weights["npw"])
        self.add_gradient_constraints(self.interpolation_weights["gpw"])
        self.add_value_constraints(self.interpolation_weights["cpw"])
        self.add_tangent_constraints(self.interpolation_weights["tpw"])
        self.add_interface_constraints(self.interpolation_weights["ipw"])
        self.add_value_inequality_constraints()
        self.add_inequality_pairs_constraints(
            pairs=kwargs.get("inequality_pairs", None),
            upper_bound=kwargs.get("inequality_pair_upper_bound", np.finfo(float).eps),
            lower_bound=kwargs.get("inequality_pair_lower_bound", -np.inf),
        )
        for k, o in operators.items():
            self.assemble_inner(o[0], o[1], name=k)
        self.add_directional_regularisation(regularisation_config.directional)
        self.assemble_borders()
        return self.finalize_setup_diagnostics_report()

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
            idc = np.asarray(node_idx, dtype=int)
            a = self.support.position_to_dof_coefs(points[inside, : self.support.dimension])
            # a *= w
            # a/=self.support.enp.product(self.support.step_vector)
            self.add_constraints_to_least_squares(
                a,
                points[inside, self.support.dimension],
                idc[inside, :],
                w=w * points[inside, self.support.dimension + 1],
                name="value",
            )
            if np.sum(inside) <= 0:
                logger.warning(
                    f"{np.sum(~inside)} \
                        value constraints not added: outside of model bounding box"
                )

    def _interface_pair_indices(self, n: int) -> np.ndarray:
        """Return index pairs for interface constraints with scalable defaults."""
        if n < 2:
            return np.zeros((0, 2), dtype=int)

        mode = str(self.interface_pair_mode).lower()
        threshold = max(2, int(self.interface_pair_threshold))
        fallback = str(self.interface_pair_fallback).lower()
        if fallback not in ("star", "chain"):
            fallback = "star"

        if mode == "auto":
            mode = fallback if n > threshold else "all"

        if mode == "all":
            ii, jj = np.triu_indices(n, k=1)
            return np.column_stack([ii, jj]).astype(int)
        if mode == "star":
            return np.column_stack([np.zeros(n - 1, dtype=int), np.arange(1, n, dtype=int)])
        if mode == "chain":
            return np.column_stack([np.arange(0, n - 1, dtype=int), np.arange(1, n, dtype=int)])

        raise ValueError(
            f"Unknown interface_pair_mode '{self.interface_pair_mode}'. Use one of: auto, all, star, chain"
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
            node_idx, inside = self.support.position_to_cell_corners(
                points[:, : self.support.dimension]
            )
            idc = np.asarray(node_idx, dtype=int)[inside, :]
            A = self.support.position_to_dof_coefs(points[inside, : self.support.dimension])
            for unique_id in np.unique(
                points[
                    np.logical_and(~np.isnan(points[:, self.support.dimension]), inside),
                    self.support.dimension,
                ]
            ):
                mask = points[inside, self.support.dimension] == unique_id
                pair_idx = self._interface_pair_indices(int(np.sum(mask)))
                if pair_idx.shape[0] == 0:
                    continue
                interface_A = np.hstack([A[mask, :][pair_idx[:, 0], :], -A[mask, :][pair_idx[:, 1], :]])
                interface_idc = np.hstack(
                    [idc[mask, :][pair_idx[:, 0], :], idc[mask, :][pair_idx[:, 1], :]]
                )
                self.add_constraints_to_least_squares(
                    interface_A,
                    np.zeros(interface_A.shape[0]),
                    interface_idc,
                    w=w,
                    name=f"interface_{unique_id}",
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

            node_idx, inside = self.support.position_to_cell_corners(
                points[:, : self.support.dimension]
            )
            # calculate unit vector for node gradients
            # this means we are only constraining direction of grad not the
            # magnitude
            idc = np.asarray(node_idx, dtype=int)

            (
                _vertices,
                T,
                _elements,
                _inside,
            ) = self.support.get_element_gradient_for_location(
                points[inside, : self.support.dimension]
            )
            # normalise constraint vector and scale element matrix by this
            norm = np.linalg.norm(
                points[:, self.support.dimension : self.support.dimension + self.support.dimension],
                axis=1,
            )
            points[:, 3:6] /= norm[:, None]
            T /= norm[inside, None, None]
            # calculate two orthogonal vectors to constraint (strike and dip vector)
            strike_vector, dip_vector = get_vectors(
                points[
                    inside, self.support.dimension : self.support.dimension + self.support.dimension
                ]
            )
            A = np.einsum("ij,ijk->ik", strike_vector.T, T)
            B = np.zeros(points[inside, :].shape[0])
            self.add_constraints_to_least_squares(A, B, idc[inside, :], w=w, name="gradient")
            A = np.einsum("ij,ijk->ik", dip_vector.T, T)
            self.add_constraints_to_least_squares(A, B, idc[inside, :], w=w, name="gradient")
            # self.regularisation_scale += compute_weighting(
            #     self.support.nodes,
            #     points[inside, : self.support.dimension],
            #     sigma=self.support.nsteps[0] * 10,
            # )
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
            node_idx, inside = self.support.position_to_cell_corners(
                points[:, : self.support.dimension]
            )
            idc = np.asarray(node_idx, dtype=int)

            # calculate unit vector for node gradients and their magnitudes
            # to preserve magnitude enforcement across the split 3-component constraint
            (
                _vertices,
                T,
                _elements,
                _inside,
            ) = self.support.get_element_gradient_for_location(
                points[inside, : self.support.dimension]
            )

            sigma = self.regularisation_weight_sigma
            if sigma is None:
                sigma = self.support.nsteps[0] * 10

            self.regularisation_scale += compute_weighting(
                self.support.nodes,
                points[inside, : self.support.dimension],
                sigma=sigma,
            )
            # Apply optional per-constraint weights from the points array.
            # For normal constraints, row format is xyz|nx ny nz|w.
            point_weights = np.ones(np.sum(inside), dtype=float)
            if points.shape[1] > self.support.dimension * 2:
                point_weights = points[inside, self.support.dimension * 2]

            # Each normal constraint is split into one row per axis below, so
            # divide by the number of axes here - otherwise a single
            # orientation point ends up weighted `dimension` times more
            # strongly than the same-`w` value/gradient constraints, enough
            # to dominate the least-squares system when data is sparse.
            w = w / self.support.dimension
            if isinstance(w, np.ndarray):
                constraint_weights = w[inside] * point_weights
            else:
                constraint_weights = float(w) * point_weights

            for d in range(self.support.dimension):
                self.add_constraints_to_least_squares(
                    T[:, d, :],
                    points[inside, self.support.dimension + d],
                    idc[inside, :],
                    w=constraint_weights,
                    name=f"norm_{d}",
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
        b: float = 0,
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
            node_idx, inside = self.support.position_to_cell_corners(
                points[:, : self.support.dimension]
            )
            # calculate unit vector for node gradients
            # this means we are only constraining direction of grad not the
            # magnitude
            idc = np.asarray(node_idx, dtype=int)
            # normalise vector and scale element gradient matrix by norm as well
            norm = np.linalg.norm(vectors, axis=1)
            vectors[norm > 0, :] /= norm[norm > 0, None]

            # normalise element vector to unit vector for dot product
            (
                _vertices,
                T,
                _elements,
                _inside,
            ) = self.support.get_element_gradient_for_location(
                points[inside, : self.support.dimension]
            )
            # norm_inside indexes norm over the inside-filtered subset so
            # the boolean mask aligns with T (shape n_inside x ...).
            norm_inside = norm[inside]
            T[norm_inside > 0, :, :] /= norm_inside[norm_inside > 0, None, None]

            # dot product of vector and element gradient = 0
            A = np.einsum("ij,ijk->ik", vectors[inside, : self.support.dimension], T)
            b_ = np.zeros(np.sum(inside)) + b
            self.add_constraints_to_least_squares(A, b_, idc[inside, :], w=w, name=name)

            if np.sum(inside) <= 0:
                logger.warning(
                    f"{np.sum(~inside)} \
                        gradient constraints not added: outside of model bounding box"
                )
            self.up_to_date = False

    def _full_neighbour_mask(self):
        if self.support.dimension == 2:
            return np.array([[-1, 0, 1, -1, 0, 1, -1, 0, 1], [1, 1, 1, 0, 0, 0, -1, -1, -1]])
        return np.array(
            [
                [
                    -1,
                    0,
                    1,
                    -1,
                    0,
                    1,
                    -1,
                    0,
                    1,
                    -1,
                    0,
                    1,
                    -1,
                    0,
                    1,
                    -1,
                    0,
                    1,
                    -1,
                    0,
                    1,
                    -1,
                    0,
                    1,
                    -1,
                    0,
                    1,
                ],
                [
                    -1,
                    -1,
                    -1,
                    0,
                    0,
                    0,
                    1,
                    1,
                    1,
                    -1,
                    -1,
                    -1,
                    0,
                    0,
                    0,
                    1,
                    1,
                    1,
                    -1,
                    -1,
                    -1,
                    0,
                    0,
                    0,
                    1,
                    1,
                    1,
                ],
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                ],
            ]
        )

    def _boundary_indexes(self, axis, upper=False):
        if self.support.nsteps[axis] < 3:
            return None
        indexes = self.support.global_index_to_node_index(np.arange(self.support.n_nodes))
        boundary_index = self.support.nsteps[axis] - 2 if upper else 1
        return indexes[indexes[:, axis] == boundary_index, :].T

    def _assemble_operator(self, operator, w, name="regularisation", indexes=None):
        active = operator.flatten() != 0
        if not np.any(active):
            return

        full_mask = self._full_neighbour_mask()
        active_mask = full_mask[:, active]
        operator_values = operator.flatten()[active]

        neighbour_kwargs = {"mask": active_mask}
        if indexes is not None:
            neighbour_kwargs["indexes"] = indexes
        global_indexes = self.support.neighbour_global_indexes(**neighbour_kwargs)
        if global_indexes is None or global_indexes.size == 0:
            return

        centre_kwargs = {"mask": np.zeros((self.support.dimension, 1), dtype=int)}
        if indexes is not None:
            centre_kwargs["indexes"] = indexes
        centre_indexes = self.support.neighbour_global_indexes(**centre_kwargs)

        idc = global_indexes.T

        centre_idc = centre_indexes.T[:, 0]

        row_w = (
            self.regularisation_scale[centre_idc.astype(int)] * w
            if self.use_regularisation_weight_scale
            else w
        )

        if self.regularisation_matrix_free:
            self._store_matrix_free_regularisation_block(
                name=name,
                idc=idc,
                operator_values=operator_values,
                row_w=row_w,
                centre_idc=centre_idc,
            )
            return

        a = np.tile(operator_values, (global_indexes.shape[1], 1))
        B = np.zeros(global_indexes.shape[1])
        self.add_constraints_to_least_squares(
            a,
            B,
            idc,
            w=row_w,
            name=name,
        )

    def _store_matrix_free_regularisation_block(
        self,
        name: str,
        idc: np.ndarray,
        operator_values: np.ndarray,
        row_w,
        centre_idc: np.ndarray | None = None,
    ) -> None:
        """Record a matrix-free regularisation block for one stencil family.

        Mirrors the row-normalisation that ``add_constraints_to_least_squares``
        applies (dividing each row by its L2 norm) so that the stored block is
        numerically equivalent to the ``coo_matrix`` that would otherwise be
        built for this family, without ever materialising it.

        Parameters
        ----------
        name : str
            Family name (e.g. ``"dxx"``, ``"dx_lower"``); used as the key in
            ``self.matrix_free_regularisation_blocks``.
        idc : np.ndarray
            (n_rows, n_active_offsets) global column indices, one row per
            interior/boundary grid node.
        operator_values : np.ndarray
            (n_active_offsets,) nonzero stencil coefficients, shared by every
            row (uniform-spacing StructuredGrid).
        row_w : float or np.ndarray
            Scalar or per-row weight (as computed for ``add_constraints_to_least_squares``).
        centre_idc : np.ndarray, optional
            (n_rows,) global dof index of the stencil's centre node for each
            row (i.e. the grid node the constraint is "centred" on). Only used
            by the CG matrix-free-regularisation fused fast path (see
            ``_build_fused_cg_regularisation_operator``) to identify rows near
            the true domain boundary; harmless to omit for callers that don't
            need it (e.g. border second-derivative families).
        """
        idc = np.asarray(idc, dtype=np.int64)
        n_rows = idc.shape[0]
        if n_rows == 0:
            return

        norm = np.linalg.norm(operator_values)
        values = (operator_values / norm if norm > 0 else operator_values).astype(float)

        if isinstance(row_w, np.ndarray):
            row_w_arr = np.asarray(row_w, dtype=float)
            if row_w_arr.shape[0] != n_rows:
                row_w_arr = np.broadcast_to(row_w_arr, (n_rows,)).copy()
        else:
            row_w_arr = np.full(n_rows, float(row_w))

        base_name = name
        count = 0
        while name in self.matrix_free_regularisation_blocks:
            name = f"{base_name}_{count}"
            count += 1

        block = {
            "idc": idc,
            "values": values,
            "w": row_w_arr,
        }
        if centre_idc is not None:
            block["centre_idc"] = np.asarray(centre_idc, dtype=np.int64).reshape(-1)
        self.matrix_free_regularisation_blocks[name] = block

    def _matrix_free_operator_from_blocks(self, blocks) -> LinearOperator:
        """Build a single ``LinearOperator`` from a list of regularisation blocks.

        Each block is a dict with ``idc`` (n_rows, n_offsets), ``values``
        (n_offsets,) and ``w`` (n_rows,), as produced by
        ``_store_matrix_free_regularisation_block``. Blocks are stacked along
        the row axis (equivalent to ``sparse.vstack`` of each family's
        weighted matrix in ``build_matrix``).

        ``matvec`` keeps the original per-family dense gather: within one
        family every row has the same number of stencil offsets, so
        ``values[None, :] * x[idc]`` followed by ``.sum(axis=1)`` is a
        rectangular, fully vectorised reduction with no index collisions to
        worry about - profiling showed this was already competitive with the
        explicit CSR matvec and looping over the (small, fixed) number of
        families adds negligible overhead.

        ``rmatvec`` is different: multiple rows can scatter into the *same*
        column (neighbouring grid nodes share stencil columns), so the
        reduction is a genuine scatter-add. The original implementation used
        ``np.add.at`` once per family, which is a known-slow path in numpy (it
        cannot vectorise internally because of possible index collisions).
        Profiling this class of family (many rows, few offsets per row, uniform
        stencil) found ``np.add.at`` to be ~8-12x slower per call than the
        explicit CSR ``rmatvec`` equivalent. Precomputing one flattened
        (row, col, weighted-value) triple across *every* family up front and
        reducing with a single ``np.bincount`` call (instead of six ``np.add.at``
        calls) removes that bottleneck: ``np.bincount`` performs the same
        "sum contributions landing on the same index" reduction in one
        optimised pass, and combining every family into one call also amortises
        Python-level loop overhead across the whole regularisation system
        rather than paying it per family.
        """
        blocks = list(blocks)
        dof = self.dof
        row_counts = [b["idc"].shape[0] for b in blocks]
        n_rows_total = int(sum(row_counts))

        # Flattened (row, col, weighted-value) triples across every family,
        # used only by rmatvec's combined bincount scatter-add (see above).
        flat_rows = []
        flat_cols = []
        flat_vals = []
        row_offset = 0
        for block in blocks:
            idc = block["idc"]
            values = block["values"]
            w = block["w"]
            n_rows, n_offsets = idc.shape
            if n_rows == 0 or n_offsets == 0:
                row_offset += n_rows
                continue
            # Global row index for every (row, offset) entry in this family.
            flat_rows.append(np.repeat(np.arange(row_offset, row_offset + n_rows), n_offsets))
            flat_cols.append(idc.ravel())
            # w is per-row, values is per-offset (shared across rows in a
            # uniform-spacing family); broadcast-multiply once and flatten.
            flat_vals.append((w[:, None] * values[None, :]).ravel())
            row_offset += n_rows

        if flat_rows:
            flat_rows = np.concatenate(flat_rows).astype(np.intp, copy=False)
            flat_cols = np.concatenate(flat_cols).astype(np.intp, copy=False)
            flat_vals = np.concatenate(flat_vals)
        else:
            flat_rows = np.zeros(0, dtype=np.intp)
            flat_cols = np.zeros(0, dtype=np.intp)
            flat_vals = np.zeros(0, dtype=float)

        def matvec(x):
            x = np.asarray(x).reshape(-1)
            out = np.empty(n_rows_total, dtype=float)
            offset = 0
            for block in blocks:
                idc = block["idc"]
                values = block["values"]
                w = block["w"]
                n_rows = idc.shape[0]
                gathered = values[None, :] * x[idc]
                out[offset : offset + n_rows] = w * gathered.sum(axis=1)
                offset += n_rows
            return out

        def rmatvec(y):
            y = np.asarray(y).reshape(-1)
            contributions = flat_vals * y[flat_rows]
            return np.bincount(flat_cols, weights=contributions, minlength=dof)

        return LinearOperator(
            shape=(n_rows_total, dof), matvec=matvec, rmatvec=rmatvec, dtype=float
        )

    def get_regularisation_linear_operator(self, names=None) -> LinearOperator | None:
        """Return a matrix-free ``LinearOperator`` for the recorded regularisation blocks.

        Only populated after ``setup_interpolator(..., regularisation_matrix_free=True)``
        has run (and ``apply_scaling_matrix`` is False, otherwise the explicit
        path is used, see ``setup_interpolator``).

        Parameters
        ----------
        names : iterable of str, optional
            Subset of family names (keys of ``self.matrix_free_regularisation_blocks``)
            to combine, e.g. ``["dxx", "dyy", "dzz"]``. Defaults to all recorded
            blocks, stacked in insertion order.

        Returns
        -------
        scipy.sparse.linalg.LinearOperator or None
            ``None`` if no matrix-free regularisation blocks have been recorded.
        """
        if not self.matrix_free_regularisation_blocks:
            return None
        if names is None:
            names = list(self.matrix_free_regularisation_blocks.keys())
        blocks = [
            self.matrix_free_regularisation_blocks[n]
            for n in names
            if n in self.matrix_free_regularisation_blocks
        ]
        if not blocks:
            return None
        return self._matrix_free_operator_from_blocks(blocks)

    def _build_fused_cg_regularisation_operator(self) -> LinearOperator | None:
        """Return a ``LinearOperator`` computing ``R_reg^T @ R_reg @ x`` (the
        matrix-free regularisation block's contribution to the CG normal
        equations) using a fused single-kernel convolution for the six
        interior second-derivative families (dxx/dyy/dzz/dxy/dxz/dyz), plus an
        exact boundary-shell correction, instead of the (rectangular,
        matvec-then-rmatvec) ``get_regularisation_linear_operator`` path.

        This is CG/normal-equations-specific (the fused kernel represents
        ``A^T A`` directly, not ``A``, so it has no rectangular matvec and
        cannot be used by lsmr) and self-adjoint by construction
        (``matvec is rmatvec``).

        Correctness
        -----------
        Composing two radius-1 stencils (correlate then convolve, i.e. one
        family's ``rmatvec(matvec(x))``) is itself a translation-invariant
        radius-2 stencil ONLY away from the true domain boundary: the real
        two-pass computation discards/zeroes "rows" at grid nodes that are not
        genuinely interior (see ``support.neighbour_global_indexes``'s default
        edge mask) before scattering back, which is not translation-invariant,
        so a single fused (5,5,5) kernel applied everywhere is wrong at nodes
        within 1 cell of the true edge (see module-level prototypes this
        implementation is derived from). This method fixes that exactly: the
        fused kernel is used only for dof at distance >= 2 from the true edge
        (``_dist_to_true_edge``), and the existing exact index-array
        matvec/rmatvec machinery (``_matrix_free_operator_from_blocks``,
        reused unmodified) is used for the outer 2-cell shell, restricted to
        just the rows whose stencil footprint can reach that shell (row centre
        distance <= 2) - the two contributions are masked to disjoint,
        exhaustive column sets (>=2 vs <=1) so they sum to the exact result at
        every single dof, not merely in the deep interior.

        Families are only fused if (a) their name is one of the six canonical
        interior families, (b) their stored normalised stencil coefficients
        match the canonical ``Operator.*_mask`` for that name (guards against
        a caller-supplied custom ``operators=`` override reusing a canonical
        name for a different mask), and (c) their per-row weight is uniform
        (spatially-varying weights, e.g. ``use_regularisation_weight_scale``,
        are not translation-invariant, so fusing them would be wrong). Any
        family that doesn't qualify (including all border second-derivative
        families) is instead included, in full and unrestricted, via the
        existing exact ``_matrix_free_operator_from_blocks`` two-pass - never
        dropped, only computed the "slow" (but always correct) way.
        """
        blocks = self.matrix_free_regularisation_blocks
        if not blocks:
            return None

        dof = self.dof

        def _normal_operator_from_blocks(block_list) -> LinearOperator | None:
            """Wrap the existing rectangular (rows x dof) matvec/rmatvec
            LinearOperator into a (dof x dof) normal-equations operator
            (``rmatvec(matvec(x))``), self-adjoint by construction. This is
            the "safe, unfused, but always exact" fallback used whenever
            fusion doesn't apply (2D/non-structured dof layout, or no
            families qualify for fusion) - `_build_fused_cg_regularisation_operator`
            must always return a square operator computing R_reg^T @ R_reg @ x,
            never the rectangular forward operator itself.
            """
            if not block_list:
                return None
            rect_op = self._matrix_free_operator_from_blocks(block_list)

            def _matvec(x):
                x = np.asarray(x, dtype=float).reshape(-1)
                return rect_op.rmatvec(rect_op.matvec(x))

            return LinearOperator(shape=(dof, dof), matvec=_matvec, rmatvec=_matvec, dtype=float)

        if getattr(self.support, "dimension", 3) == 3:
            nsteps = np.asarray(self.support.nsteps, dtype=int).reshape(-1)
        else:
            nsteps = np.zeros(0, dtype=int)

        if nsteps.shape[0] != 3 or int(np.prod(nsteps)) != dof:
            # Not a (genuine) 3D structured grid dof layout - safe fallback,
            # exact but unfused, identical to the pre-existing behaviour.
            return _normal_operator_from_blocks(list(blocks.values()))

        nx, ny, nz = (int(v) for v in nsteps)
        dist = _dist_to_true_edge(nx, ny, nz)
        deep_mask = dist >= 2
        shell_rows_mask = dist <= 2  # candidate row centres near the shell

        fused_terms = []  # list of (raw mask (3,3,3) z,y,x order, scalar w)
        shell_blocks = []  # per-family row-restricted blocks (exact boundary fix)
        remainder_blocks = []  # families handled entirely by the plain two-pass

        for name, block in blocks.items():
            idc = block["idc"]
            values = block["values"]
            w = block["w"]
            centre_idc = block.get("centre_idc")

            fusable = (
                name in _FUSED_FAMILY_MASKS
                and centre_idc is not None
                and idc.shape[0] > 0
                and w.shape[0] > 0
                and np.all(w == w[0])
            )
            if fusable:
                expected_xyz = _normalised_mask_xyz(_FUSED_FAMILY_MASKS[name])
                expected_values = expected_xyz[expected_xyz != 0]
                fusable = values.shape[0] == expected_values.shape[0] and np.allclose(
                    np.sort(values), np.sort(expected_values), atol=1e-12
                )

            if not fusable:
                remainder_blocks.append(block)
                continue

            fused_terms.append((_FUSED_FAMILY_MASKS[name], float(w[0])))

            row_mask = shell_rows_mask[centre_idc]
            if np.any(row_mask):
                shell_blocks.append(
                    {
                        "idc": idc[row_mask],
                        "values": values,
                        "w": w[row_mask],
                    }
                )

        if not fused_terms:
            # Nothing qualified for fusion (e.g. spatially-varying weights,
            # non-canonical masks, or no interior families at all) - fall back
            # to the exact, unfused operator for everything.
            return _normal_operator_from_blocks(list(blocks.values()))

        combined_kernel = None
        for mask, w in fused_terms:
            k = (w * w) * _composite_kernel(_normalised_mask_xyz(mask))
            combined_kernel = k if combined_kernel is None else combined_kernel + k

        shell_operator = (
            self._matrix_free_operator_from_blocks(shell_blocks) if shell_blocks else None
        )
        remainder_operator = (
            self._matrix_free_operator_from_blocks(remainder_blocks)
            if remainder_blocks
            else None
        )

        def matvec(x):
            x = np.asarray(x, dtype=float).reshape(-1)
            arr3d = x.reshape((nx, ny, nz), order="F")
            fused_out = ndimage.correlate(
                arr3d, combined_kernel, mode="constant", cval=0.0
            ).ravel(order="F")
            fused_out = np.where(deep_mask, fused_out, 0.0)
            out = fused_out
            if shell_operator is not None:
                shell_out = shell_operator.rmatvec(shell_operator.matvec(x))
                out = out + np.where(deep_mask, 0.0, shell_out)
            if remainder_operator is not None:
                out = out + remainder_operator.rmatvec(remainder_operator.matvec(x))
            return out

        # Self-adjoint by construction: every term above is an A^T A
        # contribution (fused-deep-interior, boundary-shell two-pass, and
        # remainder two-pass are each individually a valid A_i^T A_i, and the
        # sum of self-adjoint operators is self-adjoint), so rmatvec is the
        # same function as matvec.
        return LinearOperator(shape=(dof, dof), matvec=matvec, rmatvec=matvec, dtype=float)

    def assemble_borders(self):
        """Regularise the one-node-in-from-the-edge layer along each axis.

        The centred second-derivative stencils (dxx/dyy/dzz) assembled by
        ``assemble_inner`` only ever fire at nodes that are strictly interior
        on *all three* axes (``support.neighbour_global_indexes``'s default
        edge mask), so a node that is interior along this axis but sits on a
        face/edge/corner along another axis never gets a curvature constraint
        along this axis either, even though both of its neighbours here exist.

        This uses the same centred second-derivative mask as the interior
        pass (not a one-sided first-derivative/Neumann "zero slope" mask):
        the field's true gradient is generally nonzero and non-constant right
        up to the domain edge (e.g. a planar fault or a tilted contact), so
        forcing the slope to zero there fights the correct solution and can
        dominate the sparse-data regime enough to visibly flatten/curve
        surfaces well inside the domain. Zero curvature ("natural" boundary
        condition, as in a natural cubic spline) is compatible with any
        locally-linear field and only asks that curvature not be introduced
        artificially at the edge.
        """
        operators = []
        if self.support.dimension == 2:
            operators = [
                (
                    Operator.Dxx_mask[1, :, :],
                    self.interpolation_weights["dx"],
                    "dx_lower",
                    self._boundary_indexes(0, upper=False),
                ),
                (
                    Operator.Dxx_mask[1, :, :],
                    self.interpolation_weights["dx"],
                    "dx_upper",
                    self._boundary_indexes(0, upper=True),
                ),
                (
                    Operator.Dyy_mask[1, :, :],
                    self.interpolation_weights["dy"],
                    "dy_lower",
                    self._boundary_indexes(1, upper=False),
                ),
                (
                    Operator.Dyy_mask[1, :, :],
                    self.interpolation_weights["dy"],
                    "dy_upper",
                    self._boundary_indexes(1, upper=True),
                ),
            ]
        else:
            operators = [
                (
                    Operator.Dxx_mask,
                    self.interpolation_weights["dx"],
                    "dx_lower",
                    self._boundary_indexes(0, upper=False),
                ),
                (
                    Operator.Dxx_mask,
                    self.interpolation_weights["dx"],
                    "dx_upper",
                    self._boundary_indexes(0, upper=True),
                ),
                (
                    Operator.Dyy_mask,
                    self.interpolation_weights["dy"],
                    "dy_lower",
                    self._boundary_indexes(1, upper=False),
                ),
                (
                    Operator.Dyy_mask,
                    self.interpolation_weights["dy"],
                    "dy_upper",
                    self._boundary_indexes(1, upper=True),
                ),
                (
                    Operator.Dzz_mask,
                    self.interpolation_weights["dz"],
                    "dz_lower",
                    self._boundary_indexes(2, upper=False),
                ),
                (
                    Operator.Dzz_mask,
                    self.interpolation_weights["dz"],
                    "dz_upper",
                    self._boundary_indexes(2, upper=True),
                ),
            ]

        for operator, weight, name, indexes in operators:
            if weight == 0 or indexes is None or indexes.shape[1] == 0:
                continue
            self._assemble_operator(operator, weight, name=name, indexes=indexes)

    def assemble_inner(self, operator, w, name="regularisation"):
        """

        Parameters
        ----------
        operator : Operator mask (ndarray) or None for rectilinear grids
        w : double

        Returns
        -------

        """
        if operator is None:
            # Rectilinear grid: operator rows are built from the grid itself
            # using the operator name as a hint for which derivative to assemble.
            self._assemble_rectilinear_operator(name, w)
            return
        self._assemble_operator(operator, w, name=name)
        return

    def _assemble_rectilinear_operator(self, name: str, w: float):
        """Assemble a scaled FD regularisation operator for a rectilinear grid.

        Parameters
        ----------
        name : str
            Operator name, e.g. 'dxx', 'dyy', 'dzz', 'dxy', 'dyz', 'dxz'.
        w : float
            Weight applied to every row.
        """
        axis_map = {
            "dxx": (0, -1),
            "dyy": (1, -1),
            "dzz": (2, -1),
            "dxy": (0, 1),
            "dxz": (0, 2),
            "dyz": (1, 2),
        }
        if name not in axis_map:
            logger.warning(f"Unknown rectilinear operator name '{name}', skipping.")
            return
        axis, cross = axis_map[name]
        A_values, col_global, row_global = self.support.build_scaled_operator_rows(axis, cross)

        idc = np.asarray(col_global, dtype=int)
        centre_dof = np.asarray(row_global, dtype=int)

        B = np.zeros(idc.shape[0])
        row_w = (
            self.regularisation_scale[centre_dof.astype(int)] * w
            if self.use_regularisation_weight_scale
            else w
        )
        self.add_constraints_to_least_squares(
            A_values,
            B,
            idc,
            w=row_w,
            name=name,
        )

    def minimise_directional_gradient_change(
        self,
        w: float,
        vector: np.ndarray,
        name: str = "directional regularisation",
    ):
        """
        Anisotropic regularisation that penalises the directional second
        derivative ``(v·∇)²f = 0`` at each interior grid node.

        This is the finite-difference analogue of the P1
        ``minimise_edge_jumps`` with a direction vector.  For a given
        direction field ``v = (vx, vy, vz)`` sampled at every grid node, the
        constraint at each interior node is

        .. math::

            v_x^2 f_{xx} + v_y^2 f_{yy} + v_z^2 f_{zz}
            + 2 v_x v_y f_{xy} + 2 v_x v_z f_{xz} + 2 v_y v_z f_{yz} = 0

        The six second-derivative operators are each weighted by the
        corresponding squared direction component so the regularisation is
        strong along ``v`` and weak across it.

        Parameters
        ----------
        w : float
            Base regularisation weight.
        vector : np.ndarray, shape (n_nodes, 3)
            Direction field evaluated at every grid node
            (``self.support.nodes``).  Typically the fold normal, fold axis,
            or deformed-orientation vector returned by
            ``FoldEvent.get_deformed_orientation``.
        name : str
            Label stored with these constraints.
        """
        if vector is None or vector.ndim != 2 or vector.shape != (self.support.n_nodes, 3):
            logger.warning(
                f"{name}: vector must have shape ({self.support.n_nodes}, 3), "
                f"got {None if vector is None else vector.shape}.  Skipping."
            )
            return

        has_scaled_rows = hasattr(self.support, "build_scaled_operator_rows")

        # Six second-derivative operator types and the matching direction-
        # weight formula.  For RectilinearGrid we use build_scaled_operator_rows
        # which gives per-node stencil coefficients scaled by the local spacing.
        # For StructuredGrid (uniform spacing) we use the fixed Operator masks via
        # neighbour_global_indexes, which mirrors how _assemble_operator works.
        axis_map = {
            "dxx": (0, -1),
            "dyy": (1, -1),
            "dzz": (2, -1),
            "dxy": (0, 1),
            "dxz": (0, 2),
            "dyz": (1, 2),
        }
        # Operator masks for the mask-based path (StructuredGrid).
        op_masks = {
            "dxx": Operator.Dxx_mask,
            "dyy": Operator.Dyy_mask,
            "dzz": Operator.Dzz_mask,
            "dxy": Operator.Dxy_mask,
            "dxz": Operator.Dxz_mask,
            "dyz": Operator.Dyz_mask,
        }

        for op_key, (ax, cx) in axis_map.items():
            if has_scaled_rows:
                # --- RectilinearGrid path: per-node scaled stencil rows -------
                A_values, col_global, row_nodes = self.support.build_scaled_operator_rows(ax, cx)
                idc = np.asarray(col_global, dtype=int)
                row_nodes = np.asarray(row_nodes, dtype=int)
                if idc.shape[0] == 0:
                    continue
            else:
                # --- StructuredGrid path: fixed stencil mask ------------------
                operator = op_masks[op_key]
                active = operator.flatten() != 0
                if not np.any(active):
                    continue
                full_mask = self._full_neighbour_mask()
                active_mask = full_mask[:, active]
                operator_values = operator.flatten()[active]

                global_indexes = self.support.neighbour_global_indexes(mask=active_mask)
                if global_indexes is None or global_indexes.size == 0:
                    continue
                centre_indexes = self.support.neighbour_global_indexes(
                    mask=np.zeros((self.support.dimension, 1), dtype=int)
                )
                col_global = global_indexes
                row_nodes = centre_indexes.T[:, 0]  # global index of each interior centre node
                A_values = np.tile(operator_values, (col_global.shape[1], 1))
                idc = np.asarray(col_global.T, dtype=int)
                row_nodes = np.asarray(row_nodes, dtype=int)
                if idc.shape[0] == 0:
                    continue

            # Direction-component weight for this operator type.
            vx = vector[row_nodes, 0]
            vy = vector[row_nodes, 1]
            vz = vector[row_nodes, 2]
            if op_key == "dxx":
                comp_w = vx**2
            elif op_key == "dyy":
                comp_w = vy**2
            elif op_key == "dzz":
                comp_w = vz**2
            elif op_key == "dxy":
                comp_w = 2.0 * vx * vy
            elif op_key == "dxz":
                comp_w = 2.0 * vx * vz
            else:  # dyz
                comp_w = 2.0 * vy * vz

            row_w = w * comp_w
            # Skip rows where the direction weight is effectively zero.
            nonzero = np.abs(row_w) > 0.0
            if not np.any(nonzero):
                continue

            B = np.zeros(np.sum(nonzero))
            self.add_constraints_to_least_squares(
                A_values[nonzero],
                B,
                idc[nonzero],
                w=row_w[nonzero],
                name=f"{name}_{op_key}",
            )

    def get_regularisation_sample_points(self) -> np.ndarray:
        return self.support.nodes

    def _add_directional_regularisation(
        self,
        weight: float,
        vectors: np.ndarray,
        name: str = "directional regularisation",
    ):
        self.minimise_directional_gradient_change(weight, vectors, name=name)
