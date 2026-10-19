"""Run a (case, engine) pair and reduce it to numbers that can be pinned.

Three layers are pinned:
  * the probe-grid field and gradient  (npz)  - catches any change in the answer
  * the metrics against the analytic truth (json) - tells you if the change is better or worse
  * a fingerprint of every constraint family (json) - tells you which family changed
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import engines, geomodels

HERE = Path(__file__).parent
GOLDEN_DIR = HERE / "golden_files"
OUT_DIR = HERE / "out"

PROBE_N = 15  # probe grid: PROBE_N^3 points over the unit cube, boundary included


def probe_points() -> np.ndarray:
    ax = np.linspace(0.0, 1.0, PROBE_N)
    # order="F"-like: x fastest
    zz, yy, xx = np.meshgrid(ax, ax, ax, indexing="ij")
    return np.c_[xx.ravel(), yy.ravel(), zz.ravel()]


def eval_mask(points: np.ndarray, case) -> np.ndarray:
    lo, hi = np.array(case.eval_box[0]), np.array(case.eval_box[1])
    return np.all((points >= lo - 1e-12) & (points <= hi + 1e-12), axis=1)


def angle_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Signed-direction angle: 180 means the gradient is reversed."""
    na, nb = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.sum(a * b, axis=1) / (na * nb)
    return np.rad2deg(np.arccos(np.clip(c, -1.0, 1.0)))


def metrics(interp, case, value, grad, pts) -> dict:
    m = eval_mask(pts, case)
    truth, tgrad = case.truth(pts), case.truth_grad(pts)
    err = value - truth
    ang = angle_deg(grad, tgrad)
    nrm_ratio = np.linalg.norm(grad, axis=1) / np.linalg.norm(tgrad, axis=1)
    dv = interp.evaluate_value(case.values[:, :3]) - case.values[:, 3]
    dn = angle_deg(interp.evaluate_gradient(case.normals[:, :3]), case.normals[:, 3:6])
    return {
        "nan_fraction": float(np.mean(~np.isfinite(value[m]))),
        "value_rms": float(np.sqrt(np.nanmean(err[m] ** 2))),
        "value_max": float(np.nanmax(np.abs(err[m]))),
        "angle_median": float(np.nanmedian(ang[m])),
        "angle_p90": float(np.nanpercentile(ang[m], 90)),
        "angle_max": float(np.nanmax(ang[m])),
        "flipped_fraction": float(np.nanmean(ang[m] > 90.0)),
        "grad_norm_ratio_median": float(np.nanmedian(nrm_ratio[m])),
        "data_value_rms": float(np.sqrt(np.nanmean(dv**2))),
        "data_normal_angle_max": float(np.nanmax(dn)),
    }


def compute(case_name: str, engine: str, solver_kwargs: dict | None = None) -> dict:
    case = geomodels.build(case_name)
    interp, info = engines.run(engine, case, solver_kwargs)
    pts = probe_points()
    value = np.asarray(interp.evaluate_value(pts), float)
    grad = np.asarray(interp.evaluate_gradient(pts), float)
    return {
        "case": case,
        "interp": interp,
        "info": info,
        "points": pts,
        "value": value,
        "gradient": grad,
        "metrics": metrics(interp, case, value, grad, pts),
        "fingerprint": engines.constraint_fingerprint(interp),
    }


def _stem(case: str, engine: str) -> str:
    return f"{case}__{engine}"


def save_golden(case_name: str, engine: str, result: dict) -> None:
    GOLDEN_DIR.mkdir(exist_ok=True)
    np.savez_compressed(
        GOLDEN_DIR / (_stem(case_name, engine) + ".npz"),
        value=result["value"],
        gradient=result["gradient"],
    )
    path = GOLDEN_DIR / (_stem(case_name, engine) + ".json")
    path.write_text(
        json.dumps(
            {
                "case": case_name,
                "engine": engine,
                "dof": result["info"]["dof"],
                "metrics": result["metrics"],
                "constraints": result["fingerprint"],
            },
            indent=1,
            sort_keys=True,
        )
        + "\n"
    )


def load_golden(case_name: str, engine: str):
    npz = GOLDEN_DIR / (_stem(case_name, engine) + ".npz")
    js = GOLDEN_DIR / (_stem(case_name, engine) + ".json")
    if not npz.exists() or not js.exists():
        return None, None
    with np.load(npz) as z:
        arrays = {k: z[k] for k in z.files}
    return arrays, json.loads(js.read_text())


def all_pairs():
    for name, factory in geomodels.ALL.items():
        for engine in factory().engines:
            yield name, engine
