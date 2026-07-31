"""
Rectilinear grid support for finite difference interpolation.

A rectilinear grid has the same topology as a structured grid, but allows
varying cell sizes along each axis. Node positions along each axis are
specified explicitly as monotonically increasing 1-D arrays.
"""

from __future__ import annotations

import numpy as np

from ..logging import get_logger as getLogger
from . import SupportType
from ._3d_structured_grid import StructuredGrid

logger = getLogger(__name__)


class RectilinearGrid(StructuredGrid):
    """A rectilinear (non-uniformly spaced) 3-D structured grid.

    Unlike :class:`StructuredGrid`, the spacing between nodes may differ from
    cell to cell.  Node positions are given as three 1-D arrays *xnodes*,
    *ynodes*, *znodes* (strictly increasing).

    Parameters
    ----------
    xnodes, ynodes, znodes : array-like
        Node positions along each axis (must be strictly increasing).
    """

    def __init__(
        self,
        xnodes: np.ndarray,
        ynodes: np.ndarray,
        znodes: np.ndarray,
    ):
        xnodes = np.asarray(xnodes, dtype=float)
        ynodes = np.asarray(ynodes, dtype=float)
        znodes = np.asarray(znodes, dtype=float)

        if xnodes.ndim != 1 or ynodes.ndim != 1 or znodes.ndim != 1:
            raise ValueError("Node arrays must be 1-D.")

        self._xnodes = xnodes
        self._ynodes = ynodes
        self._znodes = znodes

        nsteps = np.array([len(xnodes) - 1, len(ynodes) - 1, len(znodes) - 1], dtype=int)
        origin = np.array([xnodes[0], ynodes[0], znodes[0]])
        with np.errstate(divide="ignore", invalid="ignore"):
            extent = np.array(
                [xnodes[-1] - xnodes[0], ynodes[-1] - ynodes[0], znodes[-1] - znodes[0]]
            )
            step_vector = np.where(nsteps > 0, extent / nsteps, 1.0)

        from ._3d_base_structured import BaseStructuredSupport

        BaseStructuredSupport.__init__(self, origin=origin, nsteps=nsteps, step_vector=step_vector)
        self.type = SupportType.RectilinearGrid
        self.regions = {}
        self.regions["everywhere"] = np.ones(self.n_nodes).astype(bool)

    # ------------------------------------------------------------------
    # Node-array accessors
    # ------------------------------------------------------------------

    @property
    def xnodes(self) -> np.ndarray:
        return self._xnodes

    @property
    def ynodes(self) -> np.ndarray:
        return self._ynodes

    @property
    def znodes(self) -> np.ndarray:
        return self._znodes

    def onGeometryChange(self):
        if self.interpolator is not None:
            self.interpolator.reset()

    # ------------------------------------------------------------------
    # Override geometry properties
    # ------------------------------------------------------------------

    @property
    def nodes(self) -> np.ndarray:
        """Return all node positions as an (n_nodes, 3) array (Fortran order)."""
        xx, yy, zz = np.meshgrid(self._xnodes, self._ynodes, self._znodes, indexing="ij")
        return np.column_stack([xx.ravel(order="F"), yy.ravel(order="F"), zz.ravel(order="F")])

    @property
    def maximum(self) -> np.ndarray:
        return np.array([self._xnodes[-1], self._ynodes[-1], self._znodes[-1]])

    @property
    def barycentre(self) -> np.ndarray:
        return self.cell_centres(np.arange(self.n_elements))

    def cell_centres(self, global_index: np.ndarray) -> np.ndarray:
        global_index = np.asarray(global_index)
        cell_idx = self.global_index_to_cell_index(global_index)
        cx = 0.5 * (self._xnodes[cell_idx[:, 0]] + self._xnodes[cell_idx[:, 0] + 1])
        cy = 0.5 * (self._ynodes[cell_idx[:, 1]] + self._ynodes[cell_idx[:, 1] + 1])
        cz = 0.5 * (self._znodes[cell_idx[:, 2]] + self._znodes[cell_idx[:, 2] + 1])
        return np.column_stack([cx, cy, cz])

    # ------------------------------------------------------------------
    # Geometry: position to cell index
    # ------------------------------------------------------------------

    def inside(self, pos: np.ndarray) -> np.ndarray:
        return np.all((pos > self.origin[None, :]) & (pos < self.maximum[None, :]), axis=1)

    def position_to_cell_index(self, pos: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return (i,j,k) cell indices and an *inside* boolean mask."""
        pos = self.check_position(pos)
        inside = self.inside(pos)

        cell_idx = np.zeros((pos.shape[0], 3), dtype=int)
        for dim, nodes in enumerate((self._xnodes, self._ynodes, self._znodes)):
            idx = np.searchsorted(nodes, pos[:, dim], side="right") - 1
            n_cells = len(nodes) - 1
            idx = np.clip(idx, 0, n_cells - 1)
            cell_idx[:, dim] = idx

        return cell_idx, inside

    def node_indexes_to_position(self, node_indexes: np.ndarray) -> np.ndarray:
        original_shape = node_indexes.shape
        ni = node_indexes.reshape((-1, 3))
        xyz = np.zeros((ni.shape[0], 3), dtype=float)
        xyz[:, 0] = self._xnodes[np.clip(ni[:, 0], 0, len(self._xnodes) - 1)]
        xyz[:, 1] = self._ynodes[np.clip(ni[:, 1], 0, len(self._ynodes) - 1)]
        xyz[:, 2] = self._znodes[np.clip(ni[:, 2], 0, len(self._znodes) - 1)]
        return xyz.reshape(original_shape)

    # ------------------------------------------------------------------
    # Local coordinates within a cell (0 -> 1 per axis)
    # ------------------------------------------------------------------

    def position_to_local_coordinates(self, pos: np.ndarray) -> np.ndarray:
        pos = np.asarray(pos)
        cell_idx, _ = self.position_to_cell_index(pos)
        local = np.zeros(pos.shape)
        for dim, nodes in enumerate((self._xnodes, self._ynodes, self._znodes)):
            x0 = nodes[cell_idx[:, dim]]
            x1 = nodes[np.minimum(cell_idx[:, dim] + 1, len(nodes) - 1)]
            h = x1 - x0
            h = np.where(h == 0, 1.0, h)
            local[:, dim] = (pos[:, dim] - x0) / h
        return local

    # ------------------------------------------------------------------
    # Gradient of shape functions
    # ------------------------------------------------------------------

    def get_element_gradient_for_location(
        self, pos: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Return (vertices, T, elements, inside) with T of shape (N, 3, 8).

        T correctly accounts for local cell size so that
        T[:, d, :] * node_values approximates df/dx_d at *pos*.
        """
        pos = np.asarray(pos)
        T = np.zeros((pos.shape[0], 3, 8))
        local_coords = self.position_to_local_coordinates(pos)
        vertices, inside = self.position_to_cell_vertices(pos)
        elements, inside = self.position_to_cell_index(pos)
        elements_global = self.global_cell_indices(elements)

        c = local_coords

        T[:, 0, 0] = (1 - c[:, 2]) * (c[:, 1] - 1)
        T[:, 0, 1] = (1 - c[:, 1]) * (1 - c[:, 2])
        T[:, 0, 2] = -c[:, 1] * (1 - c[:, 2])
        T[:, 0, 4] = -(1 - c[:, 1]) * c[:, 2]
        T[:, 0, 5] = (1 - c[:, 1]) * c[:, 2]
        T[:, 0, 6] = -c[:, 1] * c[:, 2]
        T[:, 0, 3] = c[:, 1] * (1 - c[:, 2])
        T[:, 0, 7] = c[:, 1] * c[:, 2]

        T[:, 1, 0] = (c[:, 0] - 1) * (1 - c[:, 2])
        T[:, 1, 1] = -c[:, 0] * (1 - c[:, 2])
        T[:, 1, 2] = (1 - c[:, 0]) * (1 - c[:, 2])
        T[:, 1, 4] = -(1 - c[:, 0]) * c[:, 2]
        T[:, 1, 5] = -c[:, 0] * c[:, 2]
        T[:, 1, 6] = (1 - c[:, 0]) * c[:, 2]
        T[:, 1, 3] = c[:, 0] * (1 - c[:, 2])
        T[:, 1, 7] = c[:, 0] * c[:, 2]

        T[:, 2, 0] = -(1 - c[:, 0]) * (1 - c[:, 1])
        T[:, 2, 1] = -c[:, 0] * (1 - c[:, 1])
        T[:, 2, 2] = -(1 - c[:, 0]) * c[:, 1]
        T[:, 2, 4] = (1 - c[:, 0]) * (1 - c[:, 1])
        T[:, 2, 5] = c[:, 0] * (1 - c[:, 1])
        T[:, 2, 6] = (1 - c[:, 0]) * c[:, 1]
        T[:, 2, 3] = -c[:, 0] * c[:, 1]
        T[:, 2, 7] = c[:, 0] * c[:, 1]

        # Chain rule: dN/dx = dN/d(local) / h
        for dim, nodes in enumerate((self._xnodes, self._ynodes, self._znodes)):
            h = nodes[elements[:, dim] + 1] - nodes[elements[:, dim]]
            h = np.where(h == 0, 1.0, h)
            T[:, dim, :] /= h[:, None]

        return vertices, T, elements_global, inside

    # ------------------------------------------------------------------
    # FD regularisation operators
    # ------------------------------------------------------------------

    def get_operators(self, weights: dict[str, float]) -> dict[str, tuple]:
        """Return rectilinear FD operators.

        The mask is ``None`` to signal to the FD interpolator that it should
        call :meth:`build_scaled_operator_rows` instead of the fixed stencil.
        """
        return {
            "dxx": (None, weights["dxx"]),
            "dyy": (None, weights["dyy"]),
            "dzz": (None, weights["dzz"]),
            "dxy": (None, weights["dxy"] / 4),
            "dyz": (None, weights["dyz"] / 4),
            "dxz": (None, weights["dxz"] / 4),
        }

    def build_scaled_operator_rows(
        self, axis: int, cross_axis: int = -1
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Build per-node FD stencil rows scaled for non-uniform spacing.

        Parameters
        ----------
        axis : int
            Primary axis (0=x, 1=y, 2=z).
        cross_axis : int
            If >= 0 and != axis, build d^2/(d_axis d_cross_axis).
            Otherwise, build d^2/d_axis^2.

        Returns
        -------
        A_values : (n_interior, n_stencil_pts)
        col_global_idx : (n_interior, n_stencil_pts)
        row_node_idx : (n_interior,)
        """
        nodes_list = (self._xnodes, self._ynodes, self._znodes)
        nsteps = self.nsteps  # number of *nodes* per axis

        gi_all = np.arange(self.n_nodes)
        ijk = self.global_index_to_node_index(gi_all)  # (n_nodes, 3)

        interior_mask = np.ones(self.n_nodes, dtype=bool)
        pure = cross_axis < 0 or cross_axis == axis
        if pure:
            interior_mask &= (ijk[:, axis] > 0) & (ijk[:, axis] < nsteps[axis] - 1)
        else:
            for ax in (axis, cross_axis):
                interior_mask &= (ijk[:, ax] > 0) & (ijk[:, ax] < nsteps[ax] - 1)

        interior_ijk = ijk[interior_mask]
        n_interior = interior_ijk.shape[0]
        row_node_idx = self.global_node_indices(interior_ijk)

        if pure:
            # Non-uniform second derivative:
            # f''(x_i) ~= 2[f_{i-1}/(hL*(hL+hR)) - f_i/(hL*hR) + f_{i+1}/(hR*(hL+hR))]
            nodes_ax = nodes_list[axis]
            i = interior_ijk[:, axis]
            hL = nodes_ax[i] - nodes_ax[i - 1]
            hR = nodes_ax[i + 1] - nodes_ax[i]

            coef_m = 2.0 / (hL * (hL + hR))
            coef_p = 2.0 / (hR * (hL + hR))
            coef_c = -(coef_m + coef_p)

            A_values = np.column_stack([coef_m, coef_c, coef_p])

            col_ijk = np.zeros((n_interior, 3, 3), dtype=int)
            col_ijk[:, :, :] = interior_ijk[:, None, :]
            col_ijk[:, 0, axis] -= 1
            col_ijk[:, 2, axis] += 1
            col_global_idx = self.global_node_indices(col_ijk.reshape(-1, 3)).reshape(n_interior, 3)

        else:
            # Mixed second derivative:
            # d^2f/dxdy ~= [f(i+1,j+1) - f(i+1,j-1) - f(i-1,j+1) + f(i-1,j-1)]
            #              / ((hxL+hxR) * (hyL+hyR))
            nodes_ax = nodes_list[axis]
            nodes_cx = nodes_list[cross_axis]
            ia = interior_ijk[:, axis]
            ic = interior_ijk[:, cross_axis]
            hx = nodes_ax[ia + 1] - nodes_ax[ia - 1]
            hy = nodes_cx[ic + 1] - nodes_cx[ic - 1]

            denom = hx * hy
            signs = np.array([1.0, -1.0, -1.0, 1.0])
            A_values = signs[None, :] / denom[:, None]

            da = [1, 1, -1, -1]
            dc = [1, -1, 1, -1]
            col_ijk_4 = np.zeros((n_interior, 4, 3), dtype=int)
            for k, (da_k, dc_k) in enumerate(zip(da, dc)):
                col_ijk_4[:, k, :] = interior_ijk
                col_ijk_4[:, k, axis] += da_k
                col_ijk_4[:, k, cross_axis] += dc_k
            col_global_idx = self.global_node_indices(col_ijk_4.reshape(-1, 3)).reshape(
                n_interior, 4
            )

        return A_values, col_global_idx, row_node_idx
