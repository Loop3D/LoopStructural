"""Command line for the golden suite.

  python -m tests.golden.run                      run everything, write figures + report.html to out/
  python -m tests.golden.run --update             also (re)write the golden files from the current code
  python -m tests.golden.run --case dome          one case
  python -m tests.golden.run --engine P2          one engine
  python -m tests.golden.run --vti                also write one .vti per run (open in ParaView / pyvista)
"""
from __future__ import annotations

import argparse
import html
import json

import numpy as np

from . import geomodels, harness, render

METRIC_COLS = ["value_rms", "angle_median", "angle_p90", "flipped_fraction", "data_value_rms", "data_normal_angle_max"]


def write_vti(case_name: str, engine: str, result: dict, n: int = 33):
    import pyvista as pv

    case = result["case"]
    u = np.linspace(0, 1, n)
    zz, yy, xx = np.meshgrid(u, u, u, indexing="ij")
    pts = np.c_[xx.ravel(), yy.ravel(), zz.ravel()]
    v = np.asarray(result["interp"].evaluate_value(pts), float)
    g = np.asarray(result["interp"].evaluate_gradient(pts), float)
    grid = pv.ImageData(dimensions=(n, n, n), spacing=(1 / (n - 1),) * 3, origin=(0, 0, 0))
    grid["solved"] = v
    grid["truth"] = case.truth(pts)
    grid["error"] = v - case.truth(pts)
    grid["gradient"] = g
    grid["truth_gradient"] = case.truth_grad(pts)
    grid["angle_error_deg"] = harness.angle_deg(g, case.truth_grad(pts))
    d = harness.OUT_DIR / "vti"
    d.mkdir(parents=True, exist_ok=True)
    grid.save(d / f"{case_name}__{engine}.vti")


def report(all_results: dict, golden_status: dict):
    rows = []
    for case_name, per_engine in all_results.items():
        case = geomodels.build(case_name)
        rows.append(f"<h2>{html.escape(case.title)}</h2><p>{html.escape(case.why)}</p>")
        rows.append(f'<img src="figures/{case_name}.png" style="max-width:100%">')
        head = "".join(f"<th>{c}</th>" for c in METRIC_COLS)
        body = ""
        for eng, r in per_engine.items():
            st = golden_status.get((case_name, eng), "")
            cells = "".join(f"<td>{r['metrics'][c]:.4g}</td>" for c in METRIC_COLS)
            body += f"<tr><td>{eng}</td><td>{r['info']['dof']}</td>{cells}<td>{st}</td></tr>"
        rows.append(f"<table><tr><th>engine</th><th>dof</th>{head}<th>vs golden</th></tr>{body}</table>")
    page = (
        "<!doctype html><meta charset=utf-8><title>Interpolation golden report</title>"
        "<style>body{font:14px system-ui;margin:20px;max-width:1500px}table{border-collapse:collapse;margin:8px 0 28px}"
        "td,th{border:1px solid #ccc;padding:3px 8px;text-align:right}td:first-child{text-align:left}"
        "h2{margin-top:36px}</style><h1>Interpolation golden report</h1>"
        "<p>Dashed black = truth horizons. Orange = solved horizons. See README.md for how to read each panel.</p>"
        + "".join(rows)
    )
    path = harness.OUT_DIR / "report.html"
    path.write_text(page, encoding="utf-8")
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--update", action="store_true", help="rewrite golden files")
    ap.add_argument("--case")
    ap.add_argument("--engine")
    ap.add_argument("--vti", action="store_true")
    a = ap.parse_args()

    all_results: dict[str, dict] = {}
    status: dict = {}
    for case_name, engine in harness.all_pairs():
        if (a.case and a.case != case_name) or (a.engine and a.engine != engine):
            continue
        r = harness.compute(case_name, engine)
        all_results.setdefault(case_name, {})[engine] = r
        arrays, meta = harness.load_golden(case_name, engine)
        if arrays is None:
            status[(case_name, engine)] = "no golden"
        else:
            dv = np.nanmax(np.abs(arrays["value"] - r["value"]))
            status[(case_name, engine)] = "same" if dv < 1e-6 else f"CHANGED (max dv {dv:.2g})"
        if a.update:
            harness.save_golden(case_name, engine, r)
            status[(case_name, engine)] = "updated"
        if a.vti:
            write_vti(case_name, engine, r)
        m = r["metrics"]
        print(f"{case_name:22s} {engine:7s} rms {m['value_rms']:.4f}  ang med {m['angle_median']:6.2f}  "
              f"p90 {m['angle_p90']:6.2f}  solve {r['info']['t_solve']:.1f}s  {status[(case_name, engine)]}")

    for case_name, per_engine in all_results.items():
        render.render_case(case_name, per_engine)
    path = report(all_results, status)
    print(f"\nreport: {path}")


if __name__ == "__main__":
    main()
