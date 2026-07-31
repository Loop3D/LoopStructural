"""FoldEvent — deformed-orientation geometry for fold-constrained interpolation.

Ported from LoopStructural/modelling/features/fold/_fold.py (Laurent et al., 2016).
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from loop_common.logging import get_logger

logger = get_logger(__name__)


class FoldEvent:
    """Describes fold geometry via a fold frame and two rotation-angle functions.

    The fold frame provides two coordinate fields:
    - ``foldframe.features[0]``: axial-surface scalar field (gx).  Its gradient
      (``evaluate_gradient``) gives the across-fold direction and its value
      (``evaluate_value``) parameterises the fold limb rotation angle.
    - ``foldframe.features[1]``: along-fold scalar field (gy).  Its value is
      used to parameterise the fold axis rotation angle (only when
      ``fold_axis_rotation`` is supplied instead of a constant ``fold_axis``).

    Parameters
    ----------
    foldframe
        Object exposing a ``features`` sequence of at least one (two when
        ``fold_axis_rotation`` is set) feature objects with
        ``evaluate_value(points)`` and ``evaluate_gradient(points)`` methods.
    fold_axis_rotation : callable, optional
        Function ``f(gy_values) -> angles_deg`` that gives the rotation of the
        fold axis from the gy gradient direction.  Supply either this or
        ``fold_axis``, not both.
    fold_limb_rotation : callable
        Function ``f(gx_values) -> angles_deg`` that gives the fold limb
        rotation angle from the gx gradient direction (axial-plane rotation).
    fold_axis : array-like of shape (3,), optional
        Constant fold axis direction.  Used when the axis is known a-priori and
        ``fold_axis_rotation`` is None.
    invert_norm : bool
        When True, invert dgz for points where the fold direction and axis are
        nearly parallel (dot product < 0).  Rarely needed.
    name : str
        Label for logging.
    """

    def __init__(
        self,
        foldframe,
        fold_axis_rotation: Callable | None = None,
        fold_limb_rotation: Callable | None = None,
        fold_axis: np.ndarray | None = None,
        invert_norm: bool = False,
        name: str = "Fold",
    ):
        self.foldframe = foldframe
        self.fold_axis_rotation = fold_axis_rotation
        self.fold_limb_rotation = fold_limb_rotation
        self.fold_axis = np.asarray(fold_axis, dtype=float) if fold_axis is not None else None
        self.invert_norm = invert_norm
        self.name = name

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_fold_axis_orientation(self, points: np.ndarray) -> np.ndarray:
        """Return unit fold-axis vectors at *points* (N×3)."""
        if self.fold_axis_rotation is not None:
            logger.debug("FoldEvent: evaluating fold axis via rotation function")
            dgx = self.foldframe.features[0].evaluate_gradient(points)
            dgy = self.foldframe.features[1].evaluate_gradient(points)
            valid_x = np.all(~np.isnan(dgx), axis=1)
            valid_y = np.all(~np.isnan(dgy), axis=1)
            nrm = np.linalg.norm(dgx[valid_x], axis=1)
            dgx[valid_x] /= np.where(nrm > 1e-12, nrm, 1.0)[:, None]
            nrm = np.linalg.norm(dgy[valid_y], axis=1)
            dgy[valid_y] /= np.where(nrm > 1e-12, nrm, 1.0)[:, None]
            gy = self.foldframe.features[1].evaluate_value(points)
            R1 = self._rot_mat(-dgx, self.fold_axis_rotation(gy))
            fold_axis = np.einsum("ijk,ki->kj", R1, dgy)
            nrm = np.linalg.norm(fold_axis, axis=1, keepdims=True)
            nrm = np.where(nrm > 1e-12, nrm, 1.0)
            return fold_axis / nrm

        if self.fold_axis is not None:
            logger.debug("FoldEvent: using constant fold axis")
            return np.tile(self.fold_axis, (points.shape[0], 1))

        raise ValueError(
            f"FoldEvent '{self.name}': either fold_axis or fold_axis_rotation must be set"
        )

    def get_deformed_orientation(
        self, points: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return the three fold-geometry direction vectors at *points*.

        Returns
        -------
        fold_direction : np.ndarray (N, 3)
            Vector in the deformed foliation plane, orthogonal to the fold axis.
        fold_axis : np.ndarray (N, 3)
            Fold hinge direction.
        dgz : np.ndarray (N, 3)
            Across-fold direction (perpendicular to both fold_direction and fold_axis).
        """
        if self.fold_limb_rotation is None:
            raise ValueError(
                f"FoldEvent '{self.name}': fold_limb_rotation must be set before calling "
                "get_deformed_orientation()"
            )

        fold_axis = self.get_fold_axis_orientation(points)

        gx = self.foldframe.features[0].evaluate_value(points)
        dgx = self.foldframe.features[0].evaluate_gradient(points)
        mask = np.all(~np.isnan(dgx), axis=1)
        nrm = np.linalg.norm(dgx[mask], axis=1)
        dgx[mask] /= np.where(nrm > 1e-12, nrm, 1.0)[:, None]

        dgz = np.full_like(dgx, np.nan)
        dgz[mask] = np.cross(dgx[mask], fold_axis[mask], axisa=1, axisb=1)
        nrm = np.linalg.norm(dgz[mask], axis=1)
        nrm = np.where(nrm > 1e-12, nrm, 1.0)
        dgz[mask] /= nrm[:, None]

        R2 = self._rot_mat(fold_axis, self.fold_limb_rotation(gx))
        fold_direction = np.einsum("ijk,ki->kj", R2, dgx)
        fd_nrm = np.linalg.norm(fold_direction, axis=1, keepdims=True)
        fd_nrm = np.where(fd_nrm > 1e-12, fd_nrm, 1.0)
        fold_direction = fold_direction / fd_nrm

        if self.invert_norm:
            d = np.einsum("ij,ij->i", fold_direction, fold_axis)
            dgz[mask & (d < 0)] *= -1.0

        return fold_direction, fold_axis, dgz

    def vtk(self):
        import pyvista as pv
        points = self.foldframe.vtk().points

        fold_direction, fold_axis, dgz = self.get_deformed_orientation(points)
        points = pv.PolyData(points)

        points['fold_direction'] = fold_direction
        points['fold_axis'] = fold_axis
        points['dgz'] = dgz
        points['folded_normal'] = np.cross(fold_direction, fold_axis, axisa=1, axisb=1)
        return points
    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _rot_mat(axis: np.ndarray, angle: np.ndarray) -> np.ndarray:
        """Rodrigues rotation matrices — returns (3, 3, N) array.

        Parameters
        ----------
        axis : (N, 3) unit vectors
        angle : (N,) angles in degrees
        """
        c = np.cos(np.deg2rad(angle))
        s = np.sin(np.deg2rad(angle))
        C = 1.0 - c
        x, y, z = axis[:, 0], axis[:, 1], axis[:, 2]
        xs, ys, zs = x * s, y * s, z * s
        xC, yC, zC = x * C, y * C, z * C
        xyC, yzC, zxC = x * yC, y * zC, z * xC

        R = np.zeros((3, 3, len(angle)))
        R[0, 0] = x * xC + c
        R[0, 1] = xyC - zs
        R[0, 2] = zxC + ys
        R[1, 0] = xyC + zs
        R[1, 1] = y * yC + c
        R[1, 2] = yzC - xs
        R[2, 0] = zxC - ys
        R[2, 1] = yzC + xs
        R[2, 2] = z * zC + c
        return R
