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
from LoopStructural.modelling.features.fault import FaultSegment

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


def world_foliation_model(
    anchor: np.ndarray, rotation_deg: float
) -> tuple[GeologicalModel, object, pd.DataFrame]:
    """Return a model, an unsolved foliation and its world data in a rotated frame.

    The data frame holds 27 value rows, 4 gradient rows and 4 normal rows.
    All locations are world coordinates. The function does not run the interpolation.
    """
    model = model_with_frame(anchor, anchor + (8, 10, 12), rotation_deg)
    columns = ["X", "Y", "Z", "feature_name", "val", "nx", "ny", "nz", "gx", "gy", "gz", "w"]
    values = r1_observations(anchor)
    offsets = np.array([(2, 3, 4), (5, 6, 7), (3, 8, 5), (6, 4, 9)], dtype=float)
    points = anchor + offsets
    unit = R1_GRADIENT / np.linalg.norm(R1_GRADIENT)

    def _direction_rows(prefix: str) -> pd.DataFrame:
        rows = pd.DataFrame(points, columns=["X", "Y", "Z"])
        rows["feature_name"] = "r1"
        rows["w"] = 1.0
        for axis, component in zip("xyz", unit):
            rows[f"{prefix}{axis}"] = component
        return rows

    data = pd.concat(
        [values, _direction_rows("g"), _direction_rows("n")], ignore_index=True
    ).reindex(columns=columns)
    model.data = data
    feature = model.create_and_add_foliation("r1", interpolatortype="FDI", nelements=1000)
    return model, feature, data


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


R2_SLIP = 3.0
R2_STEPS = 20
R2_DISPLACEMENT_TOL = 0.03
CURVED_RADIUS = 5.0
CURVED_GX = 0.2


def r2_displacement_section(anchor: np.ndarray) -> np.ndarray:
    """Return the 1369 world points of the R2 section at y' = 0.

    The local x' and z' coordinates are -9, -8.5, ..., 9 (37 values each).
    """
    offsets = np.arange(-9, 9.5, 0.5)
    x_values, z_values = np.meshgrid(offsets, offsets, indexing="ij")
    local = np.column_stack((x_values.ravel(), np.zeros(x_values.size), z_values.ravel()))
    return anchor + local


def r2_section_masks(anchor: np.ndarray, section: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the masks (smooth, plane) of the section points.

    The smooth mask selects ``abs(x') >= 0.5`` (1332 points). The plane mask
    selects the points on the fault plane (37 points).
    """
    smooth = np.abs((np.asarray(section) - anchor)[:, 0]) >= 0.5
    return smooth, ~smooth


def r2_analytic_fault(
    anchor: np.ndarray,
    steps: int = R2_STEPS,
    displacement: float = R2_SLIP,
    gradient_override: Callable[[int, np.ndarray, np.ndarray], np.ndarray] | None = None,
) -> FaultSegment:
    """Return a fault segment whose three frame coordinates are the R2 analytic frame.

    Parameters
    ----------
    anchor : np.ndarray
        World origin of the R2 frame.
    steps : int
        Number of integration steps of the fault.
    displacement : float
        Maximum displacement of the fault.
    gradient_override : callable, optional
        Function ``(coordinate, points, gradient)`` that returns a replacement for the
        gradient of one coordinate. It gets the current points of each call.
    """

    def _value(coordinate: int) -> Callable[[np.ndarray], np.ndarray]:
        return lambda points: r2_frame(points, anchor)[0][:, coordinate]

    def _gradient(coordinate: int) -> Callable[[np.ndarray], np.ndarray]:
        def _evaluate(points: np.ndarray) -> np.ndarray:
            gradient = r2_frame(points, anchor)[1][:, coordinate, :]
            if gradient_override is not None:
                gradient = gradient_override(coordinate, points, gradient)
            return gradient

        return _evaluate

    features = [AnalyticReference(_value(i), _gradient(i), f"r2_coordinate_{i}") for i in range(3)]
    return FaultSegment("r2_fault", features, steps=steps, displacement=displacement)


def _cubic_profile(t: np.ndarray) -> np.ndarray:
    """Return ``f(t) = 1 - 3 t^2 + 2 t^3``."""
    return 1 - 3 * t**2 + 2 * t**3


def r2_closed_form_displacement(
    section: np.ndarray, anchor: np.ndarray, slip: float = R2_SLIP
) -> np.ndarray:
    """Return the R2 displacement vectors by hand, with no library profile.

    For ``x' > 0`` the vector is ``(0, 0, -0.5 * slip * f(t))``, for ``x' < 0`` it is
    ``(0, 0, +0.5 * slip * f(t))`` and for ``x' = 0`` it is zero. Here ``t = abs(x') / 20``.
    """
    x = (np.asarray(section) - anchor)[:, 0]
    vectors = np.zeros((len(x), 3))
    vectors[:, 2] = -np.sign(x) * 0.5 * slip * _cubic_profile(np.abs(x) / 20)
    return vectors


def r2_curved_slip_fault(anchor: np.ndarray, steps: int) -> FaultSegment:
    """Return a fault segment with a slip direction that turns about a vertical axis.

    Coordinate 0 is the constant ``CURVED_GX`` and coordinate 2 is the constant 0.
    The gradient of coordinate 1 is ``(-(y - cy), (x - cx), 0) / (20 r)`` about the
    vertical axis through ``anchor``. The displacement is ``R2_SLIP``. By hand,
    ``d = 0.5 * 3 * f(0.2) = 1.344``.
    """

    def _constant(value: float) -> Callable[[np.ndarray], np.ndarray]:
        return lambda points: np.full(len(points), value)

    def _zero_gradient(points: np.ndarray) -> np.ndarray:
        return np.zeros((len(points), 3))

    def _turning_gradient(points: np.ndarray) -> np.ndarray:
        local = np.asarray(points) - anchor
        radius = np.hypot(local[:, 0], local[:, 1])
        gradient = np.zeros((len(local), 3))
        gradient[:, 0] = -local[:, 1] / (20 * radius)
        gradient[:, 1] = local[:, 0] / (20 * radius)
        return gradient

    features = [
        AnalyticReference(_constant(CURVED_GX), _zero_gradient, "curved_coordinate_0"),
        AnalyticReference(_constant(0.0), _turning_gradient, "curved_coordinate_1"),
        AnalyticReference(_constant(0.0), _zero_gradient, "curved_coordinate_2"),
    ]
    return FaultSegment("curved_fault", features, steps=steps, displacement=R2_SLIP)


def euler_circle_closed_form(radius: float, d: float, steps: int) -> tuple[float, float]:
    """Return the final radius and the turn angle of an Euler walk about a circle.

    Each step has the length ``s = d / steps`` and goes along the tangent. Then
    ``r_k^2 = radius^2 + k s^2`` and the angle gains ``arctan(s / r_k)`` in step k.
    """
    step_length = d / steps
    angle = 0.0
    for k in range(steps):
        angle += np.arctan(step_length / np.sqrt(radius**2 + k * step_length**2))
    return float(np.sqrt(radius**2 + steps * step_length**2)), float(angle)


def build_r2_fault(
    anchor: np.ndarray,
    rotation_deg: float,
    interpolatortype: str,
    explicit_center: bool = True,
    displacement: float = R2_SLIP,
) -> FaultSegment:
    """Build the R2 fault in a rotated frame, with a non-zero displacement.

    This is a copy of ``_fault`` in ``test_fault_frame_analytic.py``. Only the
    displacement and the number of steps (``R2_STEPS``) are different.
    """
    half_size = np.full(3, 10.0)
    model = model_with_frame(anchor - half_size, anchor + half_size, rotation_deg)
    data = r2_fault_data(anchor)
    model.data = data
    fault = model.create_and_add_fault(
        "fault",
        displacement,
        data=data,
        interpolatortype=interpolatortype,
        fault_center=anchor if explicit_center else None,
        fault_normal_vector=np.array([1.0, 0.0, 0.0]),
        fault_slip_vector=np.array([0.0, 0.0, -1.0]),
        major_axis=20,
        minor_axis=10,
        intermediate_axis=20,
        fault_buffer=0.2,
        nelements=20_000 if interpolatortype == "FDI" else 1_000,
        fault_trace_anisotropy=0,
        fault_dip_anisotropy=0,
    )
    fault.steps = R2_STEPS
    # These are the CG solve settings of SC-003. The builder ignores
    # top-level rtol, atol and maxiter, so they go in solver_kwargs.
    fault.builder.build(
        solver="cg",
        solver_kwargs={"rtol": 1e-10, "atol": 1e-12, "maxiter": 20_000},
        w1=1,
        w2=1,
        w3=1,
        step=10,
    )
    return fault


def faulted_foliation_model(
    anchor: np.ndarray, rotation_deg: float, engine: str
) -> tuple[GeologicalModel, FaultSegment, object]:
    """Return an updated model with one fault and one foliation in a rotated frame.

    The data are those of a unit cube. All locations are shifted by ``anchor``.

    Parameters
    ----------
    anchor : np.ndarray
        World origin of the model frame.
    rotation_deg : float
        Rotation of the local frame around the z axis, in degrees.
    engine : str
        Interpolator type of both features, ``"FDI"`` or ``"PLI"``.
    """
    model = model_with_frame(anchor, anchor + 1, rotation_deg)
    rows = [
        [0.5, 0.5, 0.5, 1, 0, 0, 0, "fault", 0],
        [0.5, 0.5, 0.5, 0, 0, 1, 1, "fault", 0],
        [0.5, 0.5, 0.5, 0, 1, 0, 2, "fault", 0],
    ]
    for x in (0.2, 0.8):
        for y in (0.2, 0.8):
            rows.append([x, y, 0.3, np.nan, np.nan, np.nan, np.nan, "strati", 0.0])
            rows.append([x, y, 0.7, np.nan, np.nan, np.nan, np.nan, "strati", 1.0])
    rows.append([0.3, 0.3, 0.5, 0, 0, 1, np.nan, "strati", np.nan])
    data = pd.DataFrame(
        rows, columns=["X", "Y", "Z", "nx", "ny", "nz", "coord", "feature_name", "val"]
    )
    data[["X", "Y", "Z"]] = data[["X", "Y", "Z"]].to_numpy(dtype=float) + anchor
    model.data = data
    fault = model.create_and_add_fault("fault", 0.1, interpolatortype=engine, nelements=1e3)
    strati = model.create_and_add_foliation("strati", interpolatortype=engine, nelements=1e3)
    model.update()
    return model, fault, strati


def residuals_of_field(interpolator: object, nodal_values: np.ndarray) -> dict[str, np.ndarray]:
    """Return ``matrix @ nodal_values - b`` for each constraint of an interpolator."""
    return {
        name: constraint["matrix"] @ nodal_values - constraint["b"]
        for name, constraint in interpolator.constraints.items()
    }
