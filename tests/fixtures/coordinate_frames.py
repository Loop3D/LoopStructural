"""Shared analytic data and reference functions for coordinate frame tests."""

from __future__ import annotations

import copy
import itertools
from collections.abc import Callable

import numpy as np
import pandas as pd
from loop_common.supports import P2TetMesh, StructuredGrid, TetMesh
from loop_interpolation import (
    ConstantNormFDIInterpolator,
    ConstantNormP1Interpolator,
    FiniteDifferenceInterpolator,
    P1Interpolator,
    P2Interpolator,
)

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
PERSIST_ENGINES = ("fdi", "pli", "p2")
LIFECYCLE_ENGINES = PERSIST_ENGINES + ("cn_fdi", "cn_pli")
QUERY_OFFSETS = np.array([[2.2, 3.7, 4.3], [6.1, 8.2, 9.1]])


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


R1_GRADIENT = np.array([0.2, -0.3, 0.4])


def _tile(vector: np.ndarray, points: np.ndarray) -> np.ndarray:
    """Repeat one world vector for each point."""
    return np.tile(np.asarray(vector, dtype=float), (len(points), 1))


def _unit(vectors: np.ndarray) -> np.ndarray:
    """Normalise each row of an (N, 3) array to unit length."""
    return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)


def r1_value(points: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    """Evaluate the R1 analytic scalar field."""
    local = np.asarray(points) - anchor
    return 2 + local @ R1_GRADIENT


def r1_gradient(points: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    """Evaluate the R1 analytic scalar field gradient."""
    del anchor
    return _tile(R1_GRADIENT, points)


def r1_observations(anchor: np.ndarray) -> pd.DataFrame:
    """Return 27 scalar observations of the R1 field around an anchor."""
    x_values, y_values, z_values = np.meshgrid(
        (1.0, 4.0, 7.0), (1.0, 5.0, 9.0), (1.0, 6.0, 11.0), indexing="ij"
    )
    local = np.column_stack((x_values.ravel(), y_values.ravel(), z_values.ravel()))
    points = anchor + local
    count = points.shape[0]
    return pd.DataFrame(
        {
            "X": points[:, 0],
            "Y": points[:, 1],
            "Z": points[:, 2],
            "val": r1_value(points, anchor),
            "feature_name": ["r1"] * count,
            "w": np.ones(count),
        }
    )


def r1_interpolator(
    engine: str, anchor: np.ndarray, rotation_deg: float, exact: bool = True
) -> tuple[object, pd.DataFrame]:
    """Return an R1 interpolator and its observations in a rotated frame.

    Parameters
    ----------
    engine : str
        One of ``LIFECYCLE_ENGINES``.
    anchor : np.ndarray
        World origin of the model frame.
    rotation_deg : float
        Rotation of the local frame around the z axis, in degrees.
    exact : bool
        If True, set the node values to the R1 field and skip the solve.
        If False, solve the system with the conjugate gradient solver.
    """
    model = model_with_frame(anchor, anchor + (8, 10, 12), rotation_deg)
    bounding_box = model.bounding_box
    corners = np.array(list(itertools.product(*zip(bounding_box.origin, bounding_box.maximum))))
    local = bounding_box.project(corners)
    lower, upper = local.min(axis=0), local.max(axis=0)
    cells = (5, 6, 7) if exact else (6, 6, 6)
    step = (upper - lower) / cells
    if engine == "fdi":
        interp = FiniteDifferenceInterpolator(
            StructuredGrid(origin=lower, nsteps_cells=cells, step_vector=step)
        )
    elif engine == "pli":
        interp = P1Interpolator(TetMesh(origin=lower, nsteps_cells=cells, step_vector=step))
    elif engine == "p2":
        interp = P2Interpolator(P2TetMesh(origin=lower, nsteps=cells, step_vector=step))
    elif engine == "cn_fdi":
        interp = ConstantNormFDIInterpolator(
            StructuredGrid(origin=lower, nsteps_cells=cells, step_vector=step)
        )
    elif engine == "cn_pli":
        interp = ConstantNormP1Interpolator(
            TetMesh(origin=lower, nsteps_cells=cells, step_vector=step)
        )
    else:
        raise ValueError(f"Unknown engine: {engine}")
    if engine.startswith("cn_"):
        interp.norm_length = np.sqrt(0.29)
    interp.bounding_box = bounding_box
    observations = r1_observations(anchor)
    interp.set_value_constraints(observations[["X", "Y", "Z", "val", "w"]].to_numpy())
    if exact:
        interp.c = r1_value(interp.local_to_world_points(interp.support.nodes), anchor)
        interp.up_to_date = True
    else:
        interp.setup()
        assert interp.solve_system(
            solver="cg", solver_kwargs={"rtol": 1e-10, "atol": 1e-12, "maxiter": 20000}
        )
    return interp, observations


def strip_frame_keys(payload: dict) -> dict:
    """Return a deep copy of a payload without the frame keys (a legacy payload)."""
    legacy = copy.deepcopy(payload)
    for key in ("format_version", "constraint_space", "transform", "bounding_box"):
        legacy.pop(key, None)
    return legacy


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


class AnalyticFold:
    """Return analytic fold orientation vectors and record the points asked for."""

    def __init__(
        self,
        deformed_fn: Callable[[np.ndarray], np.ndarray],
        axis_fn: Callable[[np.ndarray], np.ndarray],
        dgz_fn: Callable[[np.ndarray], np.ndarray],
    ) -> None:
        self.deformed_fn = deformed_fn
        self.axis_fn = axis_fn
        self.dgz_fn = dgz_fn
        self.seen_points: list[np.ndarray] = []

    def get_deformed_orientation(
        self, points: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return the deformed, axis and dgz vectors at world points."""
        self.seen_points.append(np.array(points, copy=True))
        return self.deformed_fn(points), self.axis_fn(points), self.dgz_fn(points)

    @classmethod
    def constant(cls, deformed: np.ndarray, axis: np.ndarray, dgz: np.ndarray) -> AnalyticFold:
        """Return a fold that gives three fixed world vectors everywhere."""
        return cls(
            lambda points: _tile(deformed, points),
            lambda points: _tile(axis, points),
            lambda points: _tile(dgz, points),
        )

    @classmethod
    def r4(cls, anchor: np.ndarray) -> AnalyticFold:
        """Return a fold that uses the R4 cylindrical frame."""
        return cls(
            lambda points: _unit(r4_frame(points, anchor)[1][2]),
            lambda points: _tile((0.0, 0.0, 1.0), points),
            lambda points: _unit(r4_frame(points, anchor)[1][0]),
        )
