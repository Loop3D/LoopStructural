"""Validate library P2Fold against Bermagui data held outside this repository.

Run from the repository root. Each candidate gets a fresh subprocess and the
same saved upstream DFI frames, so peak RSS and profile fits are comparable.
No prototype class, interpolator replacement, or private data is copied here.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import time

import numpy as np

ORIGIN = np.array([238275.0, 5965000.0, -50.0])
MAXIMUM = np.array([238500.0, 5965250.0, 50.0])
AXIS = [-0.31882007, -0.46674375, -0.06258398]
SEED = 12345
SOLVER = {"maxiter": 20000, "rtol": 1e-8}
FILES = {s: f"TrimData{s}.csv" for s in ("S4", "S3", "S2")}
FILES["S0"] = "TrimDataS0_LA2.csv"


def reseed():
    from loop_common.math import rng

    random.seed(SEED)
    np.random.seed(SEED)
    rng.bit_generator.state = np.random.default_rng(SEED).bit_generator.state


def peak_rss_mib():
    if os.name != "nt":
        import resource

        value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return value / (1024**2 if sys.platform == "darwin" else 1024)

    class Counters(ctypes.Structure):
        _fields_ = [("cb", ctypes.c_ulong), ("PageFaultCount", ctypes.c_ulong)] + [
            (name, ctypes.c_size_t)
            for name in (
                "PeakWorkingSetSize",
                "WorkingSetSize",
                "QuotaPeakPagedPoolUsage",
                "QuotaPagedPoolUsage",
                "QuotaPeakNonPagedPoolUsage",
                "QuotaNonPagedPoolUsage",
                "PagefileUsage",
                "PeakPagefileUsage",
            )
        ]

    counts = Counters()
    counts.cb = ctypes.sizeof(counts)
    handle = ctypes.windll.kernel32.GetCurrentProcess()
    ctypes.windll.psapi.GetProcessMemoryInfo.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(Counters),
        ctypes.c_ulong,
    ]
    if not ctypes.windll.psapi.GetProcessMemoryInfo(handle, ctypes.byref(counts), counts.cb):
        raise ctypes.WinError()
    return counts.PeakWorkingSetSize / 1024**2


def data(root, name):
    import pandas as pd

    return pd.read_csv(root / "Data" / "Processed" / FILES[name])


def save_json(path, payload):
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")


def prepare(args):
    import dill
    from LoopStructural import BoundingBox, GeologicalModel

    reseed()
    model = GeologicalModel(BoundingBox(origin=ORIGIN, maximum=MAXIMUM))
    model.set_model_data(data(args.project, "S4"))
    s4 = model.create_and_add_fold_frame(
        "S4",
        nelements=10000,
        buffer=0.5,
        solver="cg",
        damp=True,
        solver_kwargs=SOLVER,
    )
    s3 = model.create_and_add_folded_fold_frame(
        "S3",
        limb_wl=100,
        data=data(args.project, "S3"),
        fold_frame=s4,
        av_fold_axis=True,
        nelements=10000,
        buffer=0.3,
        lag=0.05,
        solver="cg",
        solver_kwargs=SOLVER,
        interpolatortype="DFI",
    )
    s3[0].evaluate_value((ORIGIN + 1)[None, :])
    reseed()
    s2 = model.create_and_add_folded_fold_frame(
        "S2",
        limb_wl=80,
        data=data(args.project, "S2"),
        fold_frame=s3,
        av_fold_axis=True,
        nelements=10000,
        buffer=args.buffer,
        lag=0.05,
        solver="cg",
        solver_kwargs=SOLVER,
        interpolatortype="DFI",
    )
    # Build all coordinates used by a downstream fold before serialising.
    for frame in (s3, s2):
        for feature in frame:
            feature.evaluate_value((ORIGIN + 1)[None, :])
    with (args.out / "drivers.pkl").open("wb") as stream:
        dill.dump(model, stream)
    save_json(
        args.out / "inputs.json",
        {
            "seed": SEED,
            "origin": ORIGIN.tolist(),
            "maximum": MAXIMUM.tolist(),
            "target_buffer": args.buffer,
            "solver": "cg",
            "solver_kwargs": SOLVER,
            "input_sha256": {
                name: hashlib.sha256(
                    (args.project / "Data" / "Processed" / file).read_bytes()
                ).hexdigest()
                for name, file in FILES.items()
            },
            "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        },
    )


def unit(vectors):
    with np.errstate(divide="ignore", invalid="ignore"):
        return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)


def worker(args):
    import dill
    from loop_interpolation import DiscreteFoldInterpolator, P2FoldInterpolator

    with (args.out / "drivers.pkl").open("rb") as stream:
        model = dill.load(stream)
    reseed()
    driver = model["S3"] if args.target == "S2" else model["S2"]
    options = (
        {"limb_wl": 80, "av_fold_axis": True}
        if args.target == "S2"
        else {
            "limb_wl": None,
            "fold_axis": AXIS,
        }
    )
    folded = model.create_and_add_folded_fold_frame(
        args.target + "_validation",
        data=data(args.project, args.target),
        fold_frame=driver,
        interpolatortype=args.engine,
        nelements=args.nelements,
        buffer=args.buffer,
        lag=0.05,
        solver="cg",
        solver_kwargs=dict(SOLVER),
        fold_weights={"sample_rule": args.rule} if args.engine == "P2Fold" else {},
        **options,
    )
    feature = folded[0]
    interp = feature.builder.interpolator
    assert type(interp) is (
        P2FoldInterpolator if args.engine == "P2Fold" else DiscreteFoldInterpolator
    )
    solve_seconds = []
    original = interp.solve_system

    def timed_solve(*pos, **kw):
        start = time.perf_counter()
        ok = original(*pos, **kw)
        solve_seconds.append(time.perf_counter() - start)
        if not ok:
            raise RuntimeError("Candidate CG did not converge")
        return ok

    interp.solve_system = timed_solve
    start = time.perf_counter()
    feature.evaluate_value((ORIGIN + 1)[None, :])
    build_seconds = time.perf_counter() - start
    assert solve_seconds, "Lazy build did not solve the candidate"
    # Same 60 x 60 x 30 grid and inner 70% box as the prototype.
    axes = [np.linspace(a, b, n) for a, b, n in zip(ORIGIN, MAXIMUM, (60, 60, 30))]
    points = np.array(np.meshgrid(*axes, indexing="ij")).reshape(3, -1).T
    fraction = (points - ORIGIN) / (MAXIMUM - ORIGIN)
    points = points[np.all((fraction > 0.15) & (fraction < 0.85), axis=1)]
    gradient = feature.evaluate_gradient(points)
    fdir, axis, _ = folded.fold.get_deformed_orientation(points)
    expected = unit(np.cross(unit(fdir), unit(axis)))
    angles = np.degrees(
        np.arccos(np.clip(np.abs(np.sum(unit(gradient) * expected, axis=1)), 0.0, 1.0))
    )
    finite = np.isfinite(angles)
    if not finite.all():
        raise RuntimeError(f"Nonfinite inner probes: {np.sum(~finite)}/{len(points)}")
    profile = folded.fold.fold_limb_rotation
    xs = np.asarray(profile.fold_frame_coordinate)
    ys = np.asarray(profile.rotation_angle)
    used = np.isfinite(xs) & np.isfinite(ys) & (np.abs(ys) < 89)
    residual = (profile(xs[used]) - ys[used] + 90) % 180 - 90
    row = {
        "target": args.target,
        "engine": args.engine,
        "rule": args.rule,
        "requested_nelements": args.nelements,
        "buffer": args.buffer,
        "unknowns": int(interp.dof),
        "tets": int(interp.support.n_elements),
        "median_angle_deg": float(np.median(angles)),
        "fraction_above_30_deg": float(np.mean(angles > 30)),
        "solve_seconds": float(sum(solve_seconds)),
        "build_seconds": build_seconds,
        "peak_rss_mib": peak_rss_mib(),
        "finite_probe_fraction": float(finite.mean()),
        "mean_gradient_norm": float(np.linalg.norm(gradient, axis=1).mean()),
        "profile_rms_mod180_deg": float(np.sqrt(np.mean(residual**2))),
        "expected_normal_sha256": hashlib.sha256(np.round(expected, 8).tobytes()).hexdigest(),
        "profile_params": {k: float(v) for k, v in profile.params.items()},
    }
    label = f"{args.target}_{args.engine}_{args.rule}_{args.nelements}"
    np.savez_compressed(
        args.out / (label + ".npz"),
        points=points,
        angle_deg=angles,
        gradient=gradient,
        expected=expected,
        profile_x=xs,
        profile_y=ys,
        curve_x=np.linspace(xs.min(), xs.max(), 2000),
        curve_y=profile(np.linspace(xs.min(), xs.max(), 2000)),
    )
    save_json(args.out / (label + ".json"), row)
    print("RESULT", json.dumps(row), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("validation/p2_fold/runtime"))
    parser.add_argument("--buffer", type=float, default=0.0)
    parser.add_argument("--mode", choices=("prepare", "worker", "all"), default="all")
    parser.add_argument("--target", choices=("S2", "S0"), default="S2")
    parser.add_argument("--engine", choices=("DFI", "P2Fold"), default="DFI")
    parser.add_argument("--rule", choices=("centroid", "gauss4"), default="gauss4")
    parser.add_argument("--nelements", type=int, default=10000)
    args = parser.parse_args()
    args.project = args.project.resolve()
    args.out = args.out.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.mode == "prepare":
        prepare(args)
    elif args.mode == "worker":
        worker(args)
    else:
        base = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--project",
            str(args.project),
            "--out",
            str(args.out),
            "--buffer",
            str(args.buffer),
        ]
        jobs = [("prepare", [])]
        for target, nelements in (("S2", 10000), ("S0", 100000)):
            configs = [("DFI", "gauss4", nelements)]
            configs += [
                ("P2Fold", rule, nel)
                for nel in (nelements, nelements // 10)
                for rule in ("centroid", "gauss4")
            ]
            for engine, rule, nel in configs:
                jobs.append(
                    (
                        f"{target}_{engine}_{rule}_{nel}",
                        [
                            "--target",
                            target,
                            "--engine",
                            engine,
                            "--rule",
                            rule,
                            "--nelements",
                            str(nel),
                        ],
                    )
                )
        for label, extra in jobs:
            print("START", label, flush=True)
            with (args.out / (label + ".log")).open("w") as log:
                subprocess.run(
                    base + ["--mode", "prepare" if label == "prepare" else "worker"] + extra,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
            print("DONE", label, flush=True)
        rows = [json.loads(p.read_text()) for p in sorted(args.out.glob("S*.json"))]
        with (args.out / "results.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
