"""Summarise a completed Bermagui comparison and render review figures."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
import numpy as np


def decision(rows):
    """Require comparable fold events and a lower median with no more unknowns."""
    result = {}
    for target in ("S2", "S0"):
        candidates = [r for r in rows if r["target"] == target]
        if not candidates:
            raise ValueError(f"Missing {target} results")
        if len({r["expected_normal_sha256"] for r in candidates}) != 1:
            raise ValueError(f"{target} candidates used different fold events")
        if len({r["buffer"] for r in candidates}) != 1:
            raise ValueError(f"{target} candidates used different buffers")
        (reference,) = [r for r in candidates if r["engine"] == "DFI"]
        eligible = [
            r
            for r in candidates
            if r["engine"] == "P2Fold"
            and r["rule"] == "gauss4"
            and r["unknowns"] <= reference["unknowns"]
        ]
        if not eligible:
            raise ValueError(f"No {target} P2Fold candidate within the DFI nodal budget")
        best = min(eligible, key=lambda r: r["median_angle_deg"])
        result[target] = {
            "passed": best["median_angle_deg"] < reference["median_angle_deg"],
            "dfi_unknowns": reference["unknowns"],
            "p2_unknowns": best["unknowns"],
            "dfi_median_deg": reference["median_angle_deg"],
            "p2_median_deg": best["median_angle_deg"],
        }
    result["passed"] = all(result[t]["passed"] for t in ("S2", "S0"))
    return result


def render(runtime, output, rows):
    output.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), layout="constrained")
    for i, target in enumerate(("S2", "S0")):
        for engine, rule, color, marker in (
            ("DFI", "gauss4", "black", "s"),
            ("P2Fold", "centroid", "tab:orange", "o"),
            ("P2Fold", "gauss4", "tab:blue", "o"),
        ):
            group = sorted(
                (
                    r
                    for r in rows
                    if r["target"] == target and r["engine"] == engine and r["rule"] == rule
                ),
                key=lambda r: r["unknowns"],
            )
            label = "DFI" if engine == "DFI" else "P2Fold " + rule
            for j, (key, title) in enumerate(
                (
                    ("median_angle_deg", "Median axial direction error (deg)"),
                    ("fraction_above_30_deg", "Share above 30 deg (%)"),
                    ("solve_seconds", "Solve time (s)"),
                )
            ):
                factor = 100 if j == 1 else 1
                axes[i, j].plot(
                    [r["unknowns"] for r in group],
                    [r[key] * factor for r in group],
                    color=color,
                    marker=marker,
                    label=label,
                )
                axes[i, j].set(xscale="log", xlabel="Unknowns", ylabel=title, title=target)
                axes[i, j].xaxis.set_major_locator(LogLocator(base=10, subs=(1, 3)))
                axes[i, j].xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x / 1000:g}k"))
                axes[i, j].xaxis.set_minor_formatter(NullFormatter())
                axes[i, j].grid(alpha=0.2)
        axes[i, 0].legend(fontsize=8)
    fig.suptitle("Bermagui: corrected library fold weights, fixed upstream frames, buffer 0")
    fig.savefig(output / "comparison.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(13, 7), layout="constrained")
    for i, target in enumerate(("S2", "S0")):
        baseline = next(r for r in rows if r["target"] == target and r["engine"] == "DFI")
        for j, (engine, rule) in enumerate(
            (("DFI", "gauss4"), ("P2Fold", "centroid"), ("P2Fold", "gauss4"))
        ):
            row = next(
                r
                for r in rows
                if r["target"] == target
                and r["engine"] == engine
                and r["rule"] == rule
                and r["requested_nelements"] == baseline["requested_nelements"]
            )
            label = f"{target}_{engine}_{rule}_{row['requested_nelements']}"
            with np.load(runtime / (label + ".npz")) as z:
                p, a = z["points"], z["angle_deg"]
                nearest_y = np.unique(p[:, 1])[np.argmin(np.abs(np.unique(p[:, 1]) - 5965125))]
                m = p[:, 1] == nearest_y
                ax = axes[i, j]
                artist = ax.pcolormesh(
                    np.unique(p[m, 0]) - 238275,
                    np.unique(p[m, 2]),
                    a[m].reshape(len(np.unique(p[m, 0])), -1).T,
                    shading="nearest",
                    vmin=0,
                    vmax=60,
                    cmap="magma_r",
                )
                ax.set(
                    title=f"{target} {engine} {rule if engine != 'DFI' else ''}\n"
                    f"{row['unknowns']:,} unknowns; median {row['median_angle_deg']:.2f} deg",
                    xlabel="Easting minus 238275 (m)",
                    ylabel="Elevation (m)",
                )
    fig.colorbar(artist, ax=axes, label="Axial direction error (deg)", shrink=0.8)
    fig.savefig(output / "sections.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout="constrained")
    for ax, target in zip(axes, ("S2", "S0")):
        row = next(r for r in rows if r["target"] == target and r["engine"] == "DFI")
        label = f"{target}_DFI_gauss4_{row['requested_nelements']}"
        with np.load(runtime / (label + ".npz")) as z:
            ax.plot(z["curve_x"], z["curve_y"], label="Fitted limb profile")
            ax.scatter(z["profile_x"], z["profile_y"], color="tab:red", label="Observations")
        ax.set(
            title=f"{target}: shared profile; RMS {row['profile_rms_mod180_deg']:.2f} deg",
            xlabel="Driver field coordinate (arbitrary scale)",
            ylabel="Rotation angle (deg)",
            ylim=(-100, 100),
        )
        ax.legend(fontsize=8)
    fig.savefig(output / "profiles.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runtime", type=Path)
    parser.add_argument("--output", type=Path, default=Path("validation/p2_fold"))
    args = parser.parse_args()
    rows = [json.loads(p.read_text()) for p in sorted(args.runtime.glob("S*.json"))]
    outcome = decision(rows)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "results.json").write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")
    (args.output / "decision.json").write_text(json.dumps(outcome, indent=2, sort_keys=True) + "\n")
    render(args.runtime, args.output, rows)
    print(json.dumps(outcome, indent=2))


if __name__ == "__main__":
    main()
