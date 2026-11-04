"""Shared analytic data and reference functions for coordinate frame tests."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel

ANCHORS = np.array(
    [
        (0.0, 0.0, 0.0),
        (100.0, 200.0, 300.0),
        (-100.0, -200.0, -300.0),
        (500000.0, 6500000.0, 100.0),
    ]
)
ROTATIONS_DEG = (0.0, 30.0, 90.0)


def rot_z(deg: float) -> np.ndarray:
    """Return a rotation matrix for a rotation around the z axis."""
    angle = np.deg2rad(deg)
    cosine, sine = np.cos(angle), np.sin(angle)
    return np.array([[cosine, -sine, 0.0], [sine, cosine, 0.0], [0.0, 0.0, 1.0]])


def tp(length: float, width: float) -> float:
    """Return a scale-aware floating point tolerance."""
    return 1e-10 * length + 32 * np.finfo(float).eps * max(width, length)


def model_with_frame(
    origin: np.ndarray, maximum: np.ndarray, rotation_deg: float
) -> GeologicalModel:
    """Create a model with a rotated local bounding-box frame."""
    model = GeologicalModel(origin, maximum)
    model.bounding_box.set_local_transform(local_origin=origin, rotation_matrix=rot_z(rotation_deg))
    return model


def _fault_data(anchor: np.ndarray, offsets: tuple[int, ...]) -> pd.DataFrame:
    y_values, z_values = np.meshgrid(offsets, offsets, indexing="ij")
    local = np.column_stack((np.zeros(y_values.size), y_values.ravel(), z_values.ravel()))
    points = anchor + local
    count = points.shape[0]
    return pd.DataFrame(
        {
            "X": points[:, 0],
            "Y": points[:, 1],
            "Z": points[:, 2],
            "val": np.zeros(count),
            "coord": np.zeros(count),
            "feature_name": ["fault"] * count,
            "w": np.ones(count),
        }
    )


def r2_fault_data(anchor: np.ndarray) -> pd.DataFrame:
    """Return fault observations for requirement R2."""
    return _fault_data(anchor, (-8, -4, 0, 4, 8))


def r5_fault_data(anchor: np.ndarray) -> pd.DataFrame:
    """Return fault observations for requirement R5."""
    return _fault_data(anchor, (-4, -2, 0, 2, 4))


def r2_probes(anchor: np.ndarray) -> np.ndarray:
    """Return the 19-cubed probe points around an anchor."""
    offsets = np.arange(-9, 10, dtype=float)
    grid = np.stack(np.meshgrid(offsets, offsets, offsets, indexing="ij"), axis=-1)
    return anchor + grid.reshape(-1, 3)


def r2_frame(points: np.ndarray, anchor: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the R2 frame values and gradients at world points."""
    local = np.asarray(points) - anchor
    x, y, z = local.T
    values = np.column_stack((x / 20, -z / 20, y / 20))
    gradients = np.broadcast_to(
        np.array([[1 / 20, 0, 0], [0, 0, -1 / 20], [0, 1 / 20, 0]]),
        (len(local), 3, 3),
    ).copy()
    return values, gradients


def r3_value(points: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    """Evaluate the R3 analytic scalar field."""
    local = np.asarray(points) - anchor
    x, y, z = local.T
    return x + 0.05 * y**2 + z


def r3_gradient(points: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    """Evaluate the R3 analytic scalar field gradient."""
    local = np.asarray(points) - anchor
    gradient = np.zeros((len(local), 3), dtype=float)
    gradient[:, 0] = 1.0
    gradient[:, 1] = 0.1 * local[:, 1]
    gradient[:, 2] = 1.0
    return gradient


def r4_frame(
    points: np.ndarray, anchor: np.ndarray
) -> tuple[np.ndarray, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Evaluate the R4 cylindrical frame values and gradients."""
    local = np.asarray(points) - anchor
    x, y, z = local.T
    radius_origin = 40.0
    world_x = radius_origin + x
    radius = np.sqrt(world_x**2 + y**2)
    values = np.column_stack(
        (
            (radius - radius_origin) / 20,
            -z / 20,
            radius_origin * np.arctan2(y, world_x) / 20,
        )
    )
    gradients = (
        np.column_stack((world_x / (20 * radius), y / (20 * radius), np.zeros(len(local)))),
        np.broadcast_to([0.0, 0.0, -1 / 20], (len(local), 3)).copy(),
        np.column_stack(
            (
                -radius_origin * y / (20 * radius**2),
                radius_origin * world_x / (20 * radius**2),
                np.zeros(len(local)),
            )
        ),
    )
    return values, gradients


class AnalyticReference:
    """Evaluate an analytic scalar field and gradient within an optional box."""

    def __init__(
        self,
        value_fn: Callable[[np.ndarray], np.ndarray],
        gradient_fn: Callable[[np.ndarray], np.ndarray],
        name: str,
        support_box: tuple[np.ndarray, np.ndarray] | None = None,
    ) -> None:
        self.value_fn = value_fn
        self.gradient_fn = gradient_fn
        self.name = name
        self.support_box = support_box

    def _supported(self, points: np.ndarray) -> np.ndarray:
        if self.support_box is None:
            return np.ones(len(points), dtype=bool)
        lower, upper = self.support_box
        return np.all((points >= lower) & (points <= upper), axis=1)

    def evaluate_value(self, points: np.ndarray) -> np.ndarray:
        """Evaluate values, returning NaN outside the closed support box."""
        points = np.asarray(points)
        values = np.asarray(self.value_fn(points), dtype=float).copy()
        values[~self._supported(points)] = np.nan
        return values

    def evaluate_gradient(self, points: np.ndarray, ignore_regions: bool = True) -> np.ndarray:
        """Evaluate gradients, returning NaN outside the closed support box."""
        del ignore_regions
        points = np.asarray(points)
        gradients = np.asarray(self.gradient_fn(points), dtype=float).copy()
        gradients[~self._supported(points)] = np.nan
        return gradients

    def __call__(self, points: np.ndarray) -> np.ndarray:
        """Evaluate the reference values."""
        return self.evaluate_value(points)
