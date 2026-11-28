from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np

EXCLUSION_EXPLANATIONS: dict[str, str] = {
    "invalid_reference_value": "The reference value of the point is not a finite number.",
    "invalid_reference_direction": "The reference direction of the point is zero or not finite.",
    "evaluation_error": "The model could not be evaluated at the point.",
    "outside_model": "The point is outside the model bounds. The count is in points.",
    "unattributed": "The rows were dropped for an unknown reason. The count is in matrix rows.",
}
UNKNOWN_EXPLANATION = "No explanation is registered for this reason."


@dataclass(frozen=True)
class ConstraintFamilyDiagnostics:
    name: str
    active: bool
    row_count: int
    dropped_rows: int | None
    effective_weight_mean: float | None
    effective_weight_min: float | None
    effective_weight_max: float | None
    source_point_count: int = 0
    outside_model_point_count: int = 0


@dataclass(frozen=True)
class RegionCoverageDiagnostics:
    total_support_nodes: int
    active_region_nodes: int
    inactive_region_nodes: int
    active_fraction: float


@dataclass(frozen=True)
class ResidualSummary:
    """Statistics of the residuals of one constraint family.

    Attributes
    ----------
    name : str
        Name of the residual set.
    kind : str
        Either ``"physical"`` or ``"assembled"``.
    unit : str
        One of ``"field"``, ``"degree"`` or ``"assembled"``.
    count : int
        Number of finite residuals.
    nan_count : int
        Number of residuals that are not finite.
    rms : float or None
        Root mean square of the finite residuals. None when ``count`` is 0.
    max_abs : float or None
        Largest absolute finite residual. None when ``count`` is 0.
    mean_abs : float or None
        Mean absolute finite residual. None when ``count`` is 0.
    """

    name: str
    kind: str
    unit: str
    count: int
    nan_count: int
    rms: float | None
    max_abs: float | None
    mean_abs: float | None


@dataclass(frozen=True)
class ExclusionRecord:
    """One reason why rows or points were excluded.

    Attributes
    ----------
    family : str
        Name of the constraint family.
    reason : str
        Short key of the reason.
    count : int
        Number of excluded items.
    explanation : str
        One sentence that explains the reason.
    """

    family: str
    reason: str
    count: int
    explanation: str


def summarise_residuals(name: str, kind: str, unit: str, residuals: np.ndarray) -> ResidualSummary:
    """Summarise a set of residuals.

    Parameters
    ----------
    name : str
        Name of the residual set.
    kind : str
        Either ``"physical"`` or ``"assembled"``.
    unit : str
        One of ``"field"``, ``"degree"`` or ``"assembled"``.
    residuals : numpy.ndarray
        The residuals. The array is flattened and is not changed.

    Returns
    -------
    ResidualSummary
        Statistics of the finite residuals. Values that are not finite
        are only counted.
    """
    flat = np.asarray(residuals, dtype=float).ravel()
    finite = flat[np.isfinite(flat)]
    nan_count = int(flat.size - finite.size)
    if finite.size == 0:
        return ResidualSummary(name, kind, unit, 0, nan_count, None, None, None)
    absolute = np.abs(finite)
    return ResidualSummary(
        name=name,
        kind=kind,
        unit=unit,
        count=int(finite.size),
        nan_count=nan_count,
        rms=float(np.sqrt(np.mean(finite**2))),
        max_abs=float(np.max(absolute)),
        mean_abs=float(np.mean(absolute)),
    )


def direction_cosines(model_vectors: np.ndarray, data_vectors: np.ndarray) -> np.ndarray:
    """Compute the cosine of the angle between paired vectors.

    Parameters
    ----------
    model_vectors : numpy.ndarray
        Array of shape (n, d).
    data_vectors : numpy.ndarray
        Array of shape (n, d).

    Returns
    -------
    numpy.ndarray
        Array of shape (n,). A row with a zero or not finite norm gives NaN.

    Raises
    ------
    ValueError
        If the two arrays do not have the same shape.
    """
    model = np.asarray(model_vectors, dtype=float)
    data = np.asarray(data_vectors, dtype=float)
    if model.shape != data.shape:
        raise ValueError(f"Shape mismatch: model vectors {model.shape}, data vectors {data.shape}.")
    model_norm = np.linalg.norm(model, axis=-1)
    data_norm = np.linalg.norm(data, axis=-1)
    valid = np.isfinite(model_norm) & np.isfinite(data_norm) & (model_norm > 0) & (data_norm > 0)
    result = np.full(model_norm.shape, np.nan)
    dot = np.sum(model * data, axis=-1)
    result[valid] = dot[valid] / (model_norm[valid] * data_norm[valid])
    return result


def build_exclusion_records(
    excluded_rows: dict[str, dict[str, int]],
    families: dict[str, ConstraintFamilyDiagnostics],
    rows_per_point: dict[str, int] | None = None,
) -> tuple[ExclusionRecord, ...]:
    """Build the records of excluded rows and points.

    Parameters
    ----------
    excluded_rows : dict
        Maps a family name to a dict of reason and count.
    families : dict
        Maps a family name to its diagnostics.
    rows_per_point : dict, optional
        Number of matrix rows that one point makes, for each family.
        The default is 1 for each family.

    Returns
    -------
    tuple of ExclusionRecord
        Records sorted by family and reason. A count of 0 makes no record.
    """
    rows_per_point = rows_per_point or {}
    records = []

    def add(family: str, reason: str, count: int) -> None:
        """Add one record when the count is more than 0."""
        if count > 0:
            explanation = EXCLUSION_EXPLANATIONS.get(reason, UNKNOWN_EXPLANATION)
            records.append(ExclusionRecord(family, reason, int(count), explanation))

    for family, reasons in excluded_rows.items():
        for reason, count in reasons.items():
            add(family, reason, count)
    for name, diagnostics in families.items():
        outside = diagnostics.outside_model_point_count
        add(name, "outside_model", outside)
        if diagnostics.dropped_rows is not None:
            k = rows_per_point.get(name, 1)
            add(name, "unattributed", diagnostics.dropped_rows - outside * k)
    return tuple(sorted(records, key=lambda r: (r.family, r.reason)))


@dataclass(frozen=True)
class ConstraintDiagnosticsReport:
    interpolator_type: str
    families: dict[str, ConstraintFamilyDiagnostics] = field(default_factory=dict)
    region_coverage: RegionCoverageDiagnostics | None = None
    outside_model_points: dict[str, int] = field(default_factory=dict)
    physical_residuals: dict[str, ResidualSummary] = field(default_factory=dict)
    assembled_residuals: dict[str, ResidualSummary] = field(default_factory=dict)
    exclusions: tuple[ExclusionRecord, ...] = ()

    @property
    def total_rows(self) -> int:
        return int(sum(f.row_count for f in self.families.values()))

    @property
    def active_families(self) -> dict[str, ConstraintFamilyDiagnostics]:
        return {k: v for k, v in self.families.items() if v.active}

    def to_dict(self) -> dict:
        data = {
            "interpolator_type": self.interpolator_type,
            "total_rows": self.total_rows,
            "families": {
                name: {
                    "active": family.active,
                    "row_count": family.row_count,
                    "dropped_rows": family.dropped_rows,
                    "effective_weight_mean": family.effective_weight_mean,
                    "effective_weight_min": family.effective_weight_min,
                    "effective_weight_max": family.effective_weight_max,
                    "source_point_count": family.source_point_count,
                    "outside_model_point_count": family.outside_model_point_count,
                }
                for name, family in self.families.items()
            },
            "region_coverage": (
                None
                if self.region_coverage is None
                else {
                    "total_support_nodes": self.region_coverage.total_support_nodes,
                    "active_region_nodes": self.region_coverage.active_region_nodes,
                    "inactive_region_nodes": self.region_coverage.inactive_region_nodes,
                    "active_fraction": self.region_coverage.active_fraction,
                }
            ),
            "outside_model_points": dict(self.outside_model_points),
        }
        if self.physical_residuals:
            data["physical_residuals"] = {
                key: asdict(value) for key, value in self.physical_residuals.items()
            }
        if self.assembled_residuals:
            data["assembled_residuals"] = {
                key: asdict(value) for key, value in self.assembled_residuals.items()
            }
        if self.exclusions:
            data["exclusions"] = [asdict(record) for record in self.exclusions]
        return data

    def summary(self) -> str:
        lines = [
            f"Constraint diagnostics for {self.interpolator_type}",
            f"Total rows: {self.total_rows}",
            "Active families:",
        ]
        active = self.active_families
        if len(active) == 0:
            lines.append("  - none")
        for name, family in sorted(active.items()):
            dropped = "unknown" if family.dropped_rows is None else str(family.dropped_rows)
            mean_weight = (
                "n/a"
                if family.effective_weight_mean is None
                else f"{family.effective_weight_mean:.4g}"
            )
            lines.append(
                f"  - {name}: rows={family.row_count}, dropped={dropped}, mean_weight={mean_weight}"
            )

        if self.region_coverage is not None:
            lines.append(
                "Region coverage: "
                f"{self.region_coverage.active_region_nodes}/{self.region_coverage.total_support_nodes} "
                f"({self.region_coverage.active_fraction:.1%})"
            )

        if len(self.outside_model_points) > 0:
            lines.append("Outside-model points:")
            for key, value in sorted(self.outside_model_points.items()):
                lines.append(f"  - {key}: {value}")

        for title, block in (
            ("Physical residuals:", self.physical_residuals),
            ("Assembled residuals:", self.assembled_residuals),
        ):
            if block:
                lines.append(title)
                for key, item in sorted(block.items()):
                    lines.append(
                        f"  - {key}: count={item.count}, nan={item.nan_count}, "
                        f"rms={item.rms}, max_abs={item.max_abs} ({item.unit})"
                    )
        if self.exclusions:
            lines.append("Exclusions:")
            for record in self.exclusions:
                lines.append(f"  - {record.family}: {record.reason}={record.count}")

        return "\n".join(lines)
