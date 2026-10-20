

def pytest_addoption(parser):
    parser.addoption(
        "--update-golden",
        action="store_true",
        default=False,
        help="rewrite golden_files/ from the current code instead of comparing",
    )
    parser.addoption(
        "--report",
        action="store_true",
        default=False,
        help="write out/test_report.html with the resulting surfaces for every test",
    )


# --------------------------------------------------------------------------- HTML test report
import base64  # noqa: E402
import html as _html  # noqa: E402
import time  # noqa: E402

_RECORDS: list[dict] = []

_FOLD_ENGINES = [("fold_sparse", e, False) for e in ("FD", "P1", "DFI", "FDFold", "P2Fold")]
_BY_NAME = {
    "test_fold_constraints_help": _FOLD_ENGINES,
    "test_refolded_fold_constraints_help": [("refolded_sparse", e, False) for e in ("FD", "P1", "DFI", "FDFold", "P2Fold")],
    "test_p2_not_worse_than_p1_on_smooth_dome": [("dome", "P1", False), ("dome", "P2", False)],
    "test_known_p2_normal_weight_has_effect": [("planar_noisy", "P2", False)],
    "test_dfi_is_repeatable_with_reseed": [("fold_sparse", "DFI", False)],
    "test_dfi_is_repeatable_without_reseed": [("fold_sparse", "DFI", False)],
}


def _pairs_for(item):
    name = item.originalname or item.name
    params = getattr(getattr(item, "callspec", None), "params", {})
    if name in _BY_NAME:
        return _BY_NAME[name]
    if "case" in params and "engine" in params:
        return [(params["case"], params["engine"], False)]
    if name == "test_known_default_solve_is_converged":
        return [("planar_dip", params["engine"], True), ("planar_dip", params["engine"], False)]
    if name in ("test_known_update_keeps_fold_weights", "test_known_points_on_lower_boundary_are_finite"):
        return [("fold_sparse", params["engine"], False)]
    return []


def pytest_runtest_makereport(item, call):
    if call.when != "call" and not (call.when == "setup" and call.excinfo is not None):
        return
    if not item.config.getoption("--report"):
        return
    outcome = "passed"
    detail = ""
    if call.excinfo is not None:
        xfail = item.get_closest_marker("xfail")
        if xfail is not None and call.excinfo.errisinstance(Exception) and call.when == "call":
            outcome, detail = "xfailed (known bug)", xfail.kwargs.get("reason", "")
        elif call.excinfo.errisinstance(__import__("_pytest.outcomes", fromlist=["Skipped"]).Skipped):
            outcome, detail = "skipped", str(call.excinfo.value)
        else:
            outcome, detail = "FAILED", str(call.excinfo.value)[:1500]
    elif item.get_closest_marker("xfail") is not None and item.get_closest_marker("xfail").kwargs.get("strict"):
        outcome, detail = "XPASS (bug fixed? remove the marker)", item.get_closest_marker("xfail").kwargs.get("reason", "")
    _RECORDS.append({"item": item, "outcome": outcome, "detail": detail, "t": call.duration,
                     "doc": (item.function.__doc__ or "").strip()})


def _b64(data: bytes) -> str:
    return "data:image/png;base64," + base64.b64encode(data).decode()


def pytest_sessionfinish(session, exitstatus):
    if not session.config.getoption("--report") or not _RECORDS:
        return
    import tempfile
    from pathlib import Path

    from . import harness, render, test_golden

    img_cache: dict = {}
    store: dict = {}  # each image is embedded once, however many tests show it

    def images(case, engine, default):
        key = (case, engine, default)
        if key not in img_cache:
            r = test_golden.result(case, engine, default)
            surf = render.render_surfaces(case, engine, r["interp"])
            with tempfile.TemporaryDirectory() as d:
                sec = render.render_case(case, {engine: r}, out_dir=d).read_bytes()
            k = f"{case}|{engine}|{int(default)}"
            store[k + "|s"], store[k + "|x"] = _b64(surf), _b64(sec)
            img_cache[key] = (k + "|s", k + "|x", r["metrics"])
        return img_cache[key]

    badge = {"passed": "#2e7d32", "skipped": "#757575", "FAILED": "#c62828"}
    blocks = []
    counts: dict[str, int] = {}
    for rec in _RECORDS:
        item, out = rec["item"], rec["outcome"]
        counts[out.split(" ")[0]] = counts.get(out.split(" ")[0], 0) + 1
        colour = badge.get(out, "#ef6c00")
        figs = ""
        for case, engine, default in _pairs_for(item):
            surf, sec, m = images(case, engine, default)
            label = f"{case} / {engine}" + (" / package-default solver" if default else " / converged solver")
            figs += (
                f'<div class="pair"><div class="lab">{_html.escape(label)} &mdash; value rms {m["value_rms"]:.3g}, '
                f'gradient error median {m["angle_median"]:.1f} deg</div>'
                f'<img class="s" data-k="{surf}"><img class="x" data-k="{sec}"></div>'
            )
        blocks.append(
            f'<section><h3><span class="b" style="background:{colour}">{_html.escape(out)}</span> '
            f'{_html.escape(item.name)} <small>{rec["t"]:.2f}s</small></h3>'
            + (f"<p>{_html.escape(rec['doc'])}</p>" if rec["doc"] else "")
            + (f'<pre>{_html.escape(rec["detail"])}</pre>' if rec["detail"] else "")
            + figs + "</section>"
        )
    summary = ", ".join(f"{v} {k}" for k, v in sorted(counts.items()))
    legend = (
        "<p><b>Left:</b> solved horizons (green / orange / purple = s 0.3 / 0.5 / 0.7) over true horizons (black mesh). "
        "White spheres = value data, red arrows = normals. <b>Right:</b> XZ section at y = 0.5 "
        "(field, gradient error, gradient arrows) and plan-view error. See README.md for the panels.</p>"
    )
    page = (
        "<!doctype html><meta charset=utf-8><title>Interpolation golden test report</title>"
        "<style>body{font:14px system-ui;margin:20px auto;max-width:1700px;padding:0 16px}"
        "section{border-top:1px solid #ccc;padding:6px 0 14px}h3{margin:10px 0 4px;font-size:15px}"
        ".b{color:#fff;padding:1px 8px;border-radius:10px;font-size:12px}small{color:#777}"
        "pre{background:#f5f5f5;padding:8px;white-space:pre-wrap;font-size:12px}"
        ".pair{margin:6px 0}.lab{font-weight:600;margin:4px 0}"
        "img.s{width:520px;vertical-align:top}img.x{width:min(1100px,calc(100% - 540px));vertical-align:top}"
        "@media(max-width:1200px){img.s,img.x{width:100%}}</style>"
        f"<h1>Interpolation golden test report</h1><p>{time.strftime('%Y-%m-%d %H:%M')} &mdash; {summary}</p>"
        + legend + "".join(blocks)
        + "<script>const IMG=" + __import__("json").dumps(store)
        + ";document.querySelectorAll('img[data-k]').forEach(i=>i.src=IMG[i.dataset.k]);</script>"
    )
    path = harness.OUT_DIR / "test_report.html"
    path.parent.mkdir(exist_ok=True)
    path.write_text(page, encoding="utf-8")
    print(f"\ntest report: {path}")
