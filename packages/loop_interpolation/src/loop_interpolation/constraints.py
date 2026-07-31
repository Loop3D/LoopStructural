from __future__ import annotations

import logging
from typing import ClassVar, Union

import numpy as np
from loop_common.base import NumpyArray
from pydantic import BaseModel, ConfigDict, Field, model_validator
from pydantic import ValidationError as PydanticValidationError

from ._validation import (
    DtypeError,
    ShapeError,
    ValidationError,
    VectorError,
    _check_finite,
    _drop_nan_data_rows,
    _ensure_float_array,
)

_logger = logging.getLogger(__name__)


def _construct(cls, **kwargs):
    """Construct a constraint model, re-raising the library's own exception
    types instead of the generic ``pydantic.ValidationError`` that wraps
    them when raised from inside a ``@model_validator``.
    """
    try:
        return cls(**kwargs)
    except PydanticValidationError as exc:
        msg = "; ".join(err.get("msg", "") for err in exc.errors())
        if "zero or near-zero magnitude" in msg:
            raise VectorError(msg) from exc
        if "must have shape" in msg or "must contain the same number of rows" in msg:
            raise ShapeError(msg) from exc
        if "could not be converted to float64" in msg:
            raise DtypeError(msg) from exc
        raise ValidationError(msg) from exc


def _to_float_array(points: np.ndarray, name: str) -> np.ndarray:
    try:
        return np.asarray(points, dtype=float)
    except (TypeError, ValueError) as e:
        raise DtypeError(f"{name} could not be converted to float64. Error: {e}") from e


class BaseConstraint(BaseModel):
    """Base pydantic model for interpolation constraints."""

    model_config = ConfigDict(
        arbitrary_types_allowed=True, validate_assignment=True, extra="forbid"
    )

    weight_name: ClassVar[str] = "constraint"

    def _set_field(self, name: str, value) -> None:
        object.__setattr__(self, name, value)

    def _validate_weights(self) -> None:
        if np.isscalar(self.weights):
            return

        weights = _ensure_float_array(self.weights, f"{self.weight_name} weights")
        if weights.ndim != 1 or weights.shape[0] != self.points.shape[0]:
            raise ShapeError(
                f"{self.weight_name} weights must be scalar or shape ({self.points.shape[0]},), got {weights.shape}"
            )
        if weights.shape[0] == 0:
            self._set_field("weights", weights)
            return

        nan_mask = ~np.isfinite(weights)
        if np.any(nan_mask):
            n_replaced = int(np.sum(nan_mask))
            _logger.warning(
                "%d NaN weight value(s) replaced with 1.0 (default weight). Provide explicit finite weights to suppress this warning.",
                n_replaced,
            )
            weights = weights.copy()
            weights[nan_mask] = 1.0

        _check_finite(weights, f"{self.weight_name} weights")
        self._set_field("weights", weights)

    @staticmethod
    def _drop_invalid_rows(
        points: np.ndarray,
        data_arrays: list[np.ndarray],
        name: str,
    ) -> tuple[np.ndarray, list[np.ndarray], np.ndarray]:
        normalized_arrays = [arr.reshape(-1, 1) if arr.ndim == 1 else arr for arr in data_arrays]
        combined = np.hstack([points, *normalized_arrays])
        filtered = _drop_nan_data_rows(combined, slice(0, combined.shape[1]), name)
        if filtered.shape[0] == combined.shape[0]:
            mask = np.ones(combined.shape[0], dtype=bool)
        else:
            mask = np.isfinite(combined).all(axis=1)

        if not data_arrays:
            return filtered[:, : points.shape[1]], [], mask

        new_arrays = []
        start = points.shape[1]
        for arr in data_arrays:
            width = 1 if arr.ndim == 1 else arr.shape[1]
            chunk = filtered[:, start : start + width]
            if arr.ndim == 1:
                chunk = chunk.reshape(-1)
            new_arrays.append(chunk)
            start += width
        return filtered[:, : points.shape[1]], new_arrays, mask


def _weights_to_column(weights: Union[float, np.ndarray], n_rows: int) -> np.ndarray:
    if n_rows == 0:
        return np.empty((0, 1), dtype=float)
    if np.isscalar(weights):
        return np.full((n_rows, 1), float(weights), dtype=float)

    arr = np.asarray(weights, dtype=float)
    if arr.ndim != 1 or arr.shape[0] != n_rows:
        raise ValueError(f"Weights must be scalar or shape ({n_rows},), got {arr.shape}")
    return arr.reshape(-1, 1)


class ValueConstraint(BaseConstraint):
    weight_name: ClassVar[str] = "Value constraint"
    points: NumpyArray = Field(default_factory=lambda: np.empty((0, 3), dtype=float))
    values: NumpyArray = Field(default_factory=lambda: np.empty((0,), dtype=float))
    weights: Union[float, NumpyArray] = 1.0

    @model_validator(mode="after")
    def check_shapes(self):
        self._set_field("points", _ensure_float_array(self.points, "Value constraint points"))
        self._set_field("values", _ensure_float_array(self.values, "Value constraint values"))

        if self.points.ndim != 2:
            raise ShapeError("Value constraint points must have shape (N, D)")
        if self.values.ndim != 1:
            raise ShapeError("Value constraint values must have shape (N,)")
        if self.points.shape[0] != self.values.shape[0]:
            raise ShapeError(
                "Value constraint points and values must contain the same number of rows"
            )

        self._validate_weights()
        points, values_list, keep_mask = self._drop_invalid_rows(
            self.points,
            [self.values],
            "value constraint",
        )
        self._set_field("points", points)
        self._set_field("values", values_list[0])
        if not np.isscalar(self.weights):
            self._set_field("weights", self.weights[keep_mask])
            self._validate_weights()

        if self.points.shape[0] == 0:
            return self

        _check_finite(self.points, "Position (X, Y, Z)")
        _check_finite(self.values, "Value column")
        return self

    def to_array(self) -> np.ndarray:
        n_rows = self.points.shape[0]
        weights = _weights_to_column(self.weights, n_rows)
        return np.hstack([self.points, self.values.reshape(-1, 1), weights])

    @classmethod
    def from_array(cls, points: np.ndarray, dimensions: int = 3) -> ValueConstraint:
        pts = _to_float_array(points, "Value constraint array")
        if pts.ndim != 2:
            raise ShapeError("Value constraint array must be 2D")
        if pts.shape[1] == dimensions + 1:
            return _construct(cls, points=pts[:, :dimensions], values=pts[:, dimensions], weights=1.0)
        if pts.shape[1] == dimensions + 2:
            return _construct(
                cls,
                points=pts[:, :dimensions],
                values=pts[:, dimensions],
                weights=pts[:, dimensions + 1],
            )
        raise ShapeError(
            f"Value constraint array must have {dimensions + 1} or {dimensions + 2} columns"
        )


class GradientConstraint(BaseConstraint):
    weight_name: ClassVar[str] = "Gradient constraint"
    points: NumpyArray = Field(default_factory=lambda: np.empty((0, 3), dtype=float))
    vectors: NumpyArray = Field(default_factory=lambda: np.empty((0, 3), dtype=float))
    weights: Union[float, NumpyArray] = 1.0
    is_normal: bool = False
    drop_invalid_rows: bool = True
    @model_validator(mode="after")
    def check_shapes(self):
        label = "Normal" if self.is_normal else "Gradient"
        vector_label = (
            "Normal vector (nx, ny, nz)" if self.is_normal else "Gradient vector (gx, gy, gz)"
        )

        self._set_field("points", _ensure_float_array(self.points, f"{label} constraint points"))
        self._set_field("vectors", _ensure_float_array(self.vectors, f"{label} constraint vectors"))

        if self.points.ndim != 2:
            raise ShapeError(f"{label} constraint points must have shape (N, D)")
        if self.vectors.ndim != 2:
            raise ShapeError(f"{label} constraint vectors must have shape (N, D)")
        if self.points.shape != self.vectors.shape:
            raise ShapeError(f"{label} constraint points and vectors must have matching shape")

        self._validate_weights()
        points, vectors_list, keep_mask = self._drop_invalid_rows(
            self.points,
            [self.vectors],
            f"{label.lower()} constraint",
        )
        self._set_field("points", points)
        self._set_field("vectors", vectors_list[0])
        if not np.isscalar(self.weights):
            self._set_field("weights", self.weights[keep_mask])
            self._validate_weights()

        if self.points.shape[0] == 0:
            return self

        _check_finite(self.points, "Position (X, Y, Z)")
        _check_finite(self.vectors, vector_label)

        magnitudes = np.linalg.norm(self.vectors, axis=1)
        zero_mag = magnitudes < 1e-14
        if np.any(zero_mag):
            if self.drop_invalid_rows:
                _logger.warning(
                    "Dropping %d %s constraints with zero or near-zero magnitude vectors.",
                    int(np.sum(zero_mag)),
                    label.lower(),
                )
                self._set_field("vectors", self.vectors[~zero_mag])
                self._set_field("points", self.points[~zero_mag])
                if not np.isscalar(self.weights):
                    self._set_field("weights", self.weights[~zero_mag])
                    self._validate_weights()
            else:

                zero_indices = np.where(zero_mag)[0]
                raise VectorError(
                    f"Found {int(np.sum(zero_mag))} {label.lower()} constraints with zero or near-zero magnitude. "
                    f"{label} vectors must have non-zero length. "
                    f"Zero-magnitude vectors at indices: {zero_indices}"
                )
        return self

    def to_array(self) -> np.ndarray:
        n_rows = self.points.shape[0]
        weights = _weights_to_column(self.weights, n_rows)
        return np.hstack([self.points, self.vectors, weights])

    @classmethod
    def from_array(
        cls, points: np.ndarray, dimensions: int = 3, is_normal: bool = False
    ) -> GradientConstraint:
        pts = _to_float_array(points, "Gradient constraint array")
        if pts.ndim != 2:
            raise ShapeError("Gradient constraint array must be 2D")
        if pts.shape[1] == dimensions * 2:
            return _construct(
                cls,
                points=pts[:, :dimensions],
                vectors=pts[:, dimensions : 2 * dimensions],
                weights=1.0,
                is_normal=is_normal,
            )
        if pts.shape[1] == dimensions * 2 + 1:
            return _construct(
                cls,
                points=pts[:, :dimensions],
                vectors=pts[:, dimensions : 2 * dimensions],
                weights=pts[:, 2 * dimensions],
                is_normal=is_normal,
            )
        raise ShapeError(
            f"Gradient constraint array must have {dimensions * 2} or {dimensions * 2 + 1} columns"
        )


class InequalityConstraint(BaseConstraint):
    weight_name: ClassVar[str] = "Inequality constraint"
    points: NumpyArray = Field(default_factory=lambda: np.empty((0, 3), dtype=float))
    bounds: NumpyArray = Field(default_factory=lambda: np.empty((0, 2), dtype=float))
    weights: Union[float, NumpyArray] = 1.0

    @model_validator(mode="after")
    def check_shapes(self):
        self._set_field("points", _ensure_float_array(self.points, "Inequality constraint points"))
        self._set_field("bounds", _ensure_float_array(self.bounds, "Inequality constraint bounds"))

        if self.points.ndim != 2:
            raise ShapeError("Inequality constraint points must have shape (N, D)")
        if self.bounds.ndim != 2 or self.bounds.shape[1] != 2:
            raise ShapeError("Inequality constraint bounds must have shape (N, 2)")
        if self.points.shape[0] != self.bounds.shape[0]:
            raise ShapeError(
                "Inequality constraint points and bounds must contain the same number of rows"
            )

        self._validate_weights()
        _check_finite(self.points, "Position (X, Y, Z)")
        _check_finite(self.bounds, "Bound values")

        lower_bounds = self.bounds[:, 0]
        upper_bounds = self.bounds[:, 1]
        invalid_bounds = lower_bounds >= upper_bounds
        if np.any(invalid_bounds):
            invalid_indices = np.where(invalid_bounds)[0]
            raise ValidationError(
                f"Found {int(np.sum(invalid_bounds))} inequality constraints with lower_bound >= upper_bound. "
                "For each constraint, lower_bound must be strictly less than upper_bound. "
                f"Invalid constraints at indices: {invalid_indices}"
            )
        return self

    def to_array(self) -> np.ndarray:
        n_rows = self.points.shape[0]
        weights = _weights_to_column(self.weights, n_rows)
        return np.hstack([self.points, self.bounds, weights])

    @classmethod
    def from_array(cls, points: np.ndarray, dimensions: int = 3) -> InequalityConstraint:
        pts = _to_float_array(points, "Inequality constraint array")
        if pts.ndim != 2:
            raise ShapeError("Inequality constraint array must be 2D")
        if pts.shape[1] == dimensions + 2:
            return _construct(
                cls, points=pts[:, :dimensions], bounds=pts[:, dimensions : dimensions + 2]
            )
        if pts.shape[1] == dimensions + 3:
            return _construct(
                cls,
                points=pts[:, :dimensions],
                bounds=pts[:, dimensions : dimensions + 2],
                weights=pts[:, dimensions + 2],
            )
        raise ShapeError(
            f"Inequality constraint array must have {dimensions + 2} or {dimensions + 3} columns"
        )


class InequalityPair(BaseConstraint):
    """Represents pairwise ordering rows [x, y, z, pair_id, weight?]."""

    weight_name: ClassVar[str] = "Inequality pairs constraint"
    points: NumpyArray = Field(default_factory=lambda: np.empty((0, 3), dtype=float))
    pair_ids: NumpyArray = Field(default_factory=lambda: np.empty((0,), dtype=float))
    weights: Union[float, NumpyArray] = 1.0

    @model_validator(mode="after")
    def check_shapes(self):
        self._set_field(
            "points", _ensure_float_array(self.points, "Inequality pairs constraint points")
        )
        self._set_field(
            "pair_ids", _ensure_float_array(self.pair_ids, "Inequality pairs constraint ids")
        )

        if self.points.ndim != 2:
            raise ShapeError("Inequality pairs constraint points must have shape (N, D)")
        if self.pair_ids.ndim != 1:
            raise ShapeError("Inequality pairs constraint pair_ids must have shape (N,)")
        if self.points.shape[0] != self.pair_ids.shape[0]:
            raise ShapeError(
                "Inequality pairs constraint points and pair_ids must contain the same number of rows"
            )

        self._validate_weights()
        _check_finite(self.points, "Position (X, Y, Z)")
        _check_finite(self.pair_ids, "Pair id column")
        return self

    def to_array(self) -> np.ndarray:
        n_rows = self.points.shape[0]
        weights = _weights_to_column(self.weights, n_rows)
        return np.hstack([self.points, self.pair_ids.reshape(-1, 1), weights])

    @classmethod
    def from_array(cls, points: np.ndarray, dimensions: int = 3) -> InequalityPair:
        pts = _to_float_array(points, "Inequality pair constraint array")
        if pts.ndim != 2:
            raise ShapeError("Inequality pair constraint array must be 2D")
        if pts.shape[1] == dimensions + 1:
            return _construct(
                cls, points=pts[:, :dimensions], pair_ids=pts[:, dimensions], weights=1.0
            )
        if pts.shape[1] == dimensions + 2:
            return _construct(
                cls,
                points=pts[:, :dimensions],
                pair_ids=pts[:, dimensions],
                weights=pts[:, dimensions + 1],
            )
        raise ShapeError(
            f"Inequality pair array must have {dimensions + 1} or {dimensions + 2} columns"
        )


class InterfaceConstraint(BaseConstraint):
    weight_name: ClassVar[str] = "Interface constraint"
    points: NumpyArray = Field(default_factory=lambda: np.empty((0, 3), dtype=float))
    interface_ids: NumpyArray = Field(default_factory=lambda: np.empty((0,), dtype=float))
    weights: Union[float, NumpyArray] = 1.0

    @model_validator(mode="after")
    def check_shapes(self):
        self._set_field("points", _ensure_float_array(self.points, "Interface constraint points"))
        self._set_field(
            "interface_ids", _ensure_float_array(self.interface_ids, "Interface constraint ids")
        )

        if self.points.ndim != 2:
            raise ShapeError("Interface constraint points must have shape (N, D)")
        if self.interface_ids.ndim != 1:
            raise ShapeError("Interface constraint interface_ids must have shape (N,)")
        if self.points.shape[0] != self.interface_ids.shape[0]:
            raise ShapeError(
                "Interface constraint points and interface_ids must contain the same number of rows"
            )

        self._validate_weights()
        points, interface_ids_list, keep_mask = self._drop_invalid_rows(
            self.points,
            [self.interface_ids],
            "interface constraint",
        )
        self._set_field("points", points)
        self._set_field("interface_ids", interface_ids_list[0])
        if not np.isscalar(self.weights):
            self._set_field("weights", self.weights[keep_mask])
            self._validate_weights()

        if self.points.shape[0] == 0:
            return self

        _check_finite(self.points, "Position (X, Y, Z)")
        _check_finite(self.interface_ids, "Interface id column")
        return self

    def to_array(self) -> np.ndarray:
        n_rows = self.points.shape[0]
        weights = _weights_to_column(self.weights, n_rows)
        return np.hstack([self.points, self.interface_ids.reshape(-1, 1), weights])

    @classmethod
    def from_array(cls, points: np.ndarray, dimensions: int = 3) -> InterfaceConstraint:
        pts = _to_float_array(points, "Interface constraint array")
        if pts.ndim != 2:
            raise ShapeError("Interface constraint array must be 2D")
        if pts.shape[1] == dimensions + 1:
            return _construct(
                cls, points=pts[:, :dimensions], interface_ids=pts[:, dimensions], weights=1.0
            )
        if pts.shape[1] == dimensions + 2:
            return _construct(
                cls,
                points=pts[:, :dimensions],
                interface_ids=pts[:, dimensions],
                weights=pts[:, dimensions + 1],
            )
        raise ShapeError(
            f"Interface constraint array must have {dimensions + 1} or {dimensions + 2} columns"
        )
