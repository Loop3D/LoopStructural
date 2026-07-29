from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional


@dataclass(frozen=True)
class ConstraintFamilyDiagnostics:
    name: str
    active: bool
    row_count: int
    dropped_rows: Optional[int]
    effective_weight_mean: Optional[float]
    effective_weight_min: Optional[float]
    effective_weight_max: Optional[float]
    source_point_count: int = 0
    outside_model_point_count: int = 0


@dataclass(frozen=True)
class RegionCoverageDiagnostics:
    total_support_nodes: int
    active_region_nodes: int
    inactive_region_nodes: int
    active_fraction: float


@dataclass(frozen=True)
class ConstraintDiagnosticsReport:
    interpolator_type: str
    families: Dict[str, ConstraintFamilyDiagnostics] = field(default_factory=dict)
    region_coverage: Optional[RegionCoverageDiagnostics] = None
    outside_model_points: Dict[str, int] = field(default_factory=dict)

    @property
    def total_rows(self) -> int:
        return int(sum(f.row_count for f in self.families.values()))

    @property
    def active_families(self) -> Dict[str, ConstraintFamilyDiagnostics]:
        return {k: v for k, v in self.families.items() if v.active}

    def to_dict(self) -> dict:
        return {
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
            "region_coverage": None
            if self.region_coverage is None
            else {
                "total_support_nodes": self.region_coverage.total_support_nodes,
                "active_region_nodes": self.region_coverage.active_region_nodes,
                "inactive_region_nodes": self.region_coverage.inactive_region_nodes,
                "active_fraction": self.region_coverage.active_fraction,
            },
            "outside_model_points": dict(self.outside_model_points),
        }

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

        return "\n".join(lines)
