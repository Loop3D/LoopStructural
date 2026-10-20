"""Golden tests for the interpolation code.

Run:   python -m pytest tests/golden -q
Update golden files after an intended change:   python -m pytest tests/golden --update-golden
Figures for debugging:   python -m tests.golden.run   (see README.md; run from packages/loop_interpolation)

Layers, in the order you should read a failure:
  test_physical_bounds / test_data_honoured   the answer is geologically wrong
  test_matches_golden_field                   the answer changed (a diff figure is written to out/failures)
  test_matches_golden_constraints             names the constraint family whose matrix changed
  test_known_*                                bugs from agent/REVIEW_interpolation_hardening.md, strict xfail:
                                              the test FAILS when the bug is fixed - then remove the
                                              marker and update the golden files.
"""
from __future__ import annotations

import functools

import numpy as np
import pytest

from . import engines, geomodels, harness

pytestmark = pytest.mark.filterwarnings("ignore::DeprecationWarning")

PAIRS = list(harness.all_pairs())
IDS = [f"{c}-{e}" for c, e in PAIRS]


@functools.lru_cache(maxsize=None)
def result(case: str, engine: str, default_solver: bool = False):
    return harness.compute(case, engine, solver_kwargs={} if default_solver else None)


# --------------------------------------------------------------------------- physical
# (value_rms, angle_median_deg) upper bounds. Planar truth must be exact. The others are
# loose limits set about 1.5x above what the code gives today. Tighten them as the code improves.
BOUNDS = {
    ("planar_dip", "*"): (1e-4, 0.1),
    ("planar_extrapolation", "*"): (1e-4, 0.1),
    ("planar_noisy", "*"): (0.03, 3.0),
    ("dome", "*"): (0.09, 17.0),
    ("fold_section", "*"): (0.12, 32.0),
    ("fold_sparse", "DFI"): (0.15, 6.0),
    ("fold_sparse", "FDFold"): (0.20, 15.0),
    ("refolded_sparse", "DFI"): (0.20, 8.0),
    ("refolded_sparse", "FDFold"): (0.13, 12.0),
}


def _bound(case, engine):
    return BOUNDS.get((case, engine)) or BOUNDS.get((case, "*"))


@pytest.mark.parametrize("case,engine", PAIRS, ids=IDS)
def test_solves_and_is_finite_inside(case, engine):
    r = result(case, engine)
    assert r["info"]["solved"]
    assert r["metrics"]["nan_fraction"] == 0.0


@pytest.mark.parametrize("case,engine", PAIRS, ids=IDS)
def test_physical_bounds(case, engine):
    bound = _bound(case, engine)
    if bound is None:
        pytest.skip("reference only: no fold constraint, the error is expected to be large")
    m = result(case, engine)["metrics"]
    assert m["value_rms"] <= bound[0], f"value rms {m['value_rms']:.4g} > {bound[0]}"
    assert m["angle_median"] <= bound[1], f"median gradient error {m['angle_median']:.2f} deg > {bound[1]}"
    assert m["flipped_fraction"] == 0.0, "part of the gradient field points the wrong way (stratigraphic inversion)"


@pytest.mark.parametrize("case,engine", PAIRS, ids=IDS)
def test_data_honoured(case, engine):
    if case == "planar_noisy":
        pytest.skip("noisy data are smoothed on purpose")
    m = result(case, engine)["metrics"]
    # Refolded: 9 sparse data against fold constraints on every element. The fold constraints win (FDFold ends 40 deg from one normal; that is pinned, not endorsed).
    value_tol, angle_tol = (0.05, 45.0) if case == "refolded_sparse" else (0.03, 8.0)
    assert m["data_value_rms"] < value_tol
    assert m["data_normal_angle_max"] < angle_tol


def test_fold_constraints_help():
    """The reason the fold interpolators exist."""
    plain = {e: result("fold_sparse", e)["metrics"]["angle_median"] for e in ("P1", "FD")}
    assert result("fold_sparse", "DFI")["metrics"]["angle_median"] < 0.25 * plain["P1"]
    assert result("fold_sparse", "FDFold")["metrics"]["angle_median"] < 0.5 * plain["FD"]


def test_refolded_fold_constraints_help():
    """Same, with a hinge plunge that changes along the axis (fold_axis_rotation path)."""
    plain = {e: result("refolded_sparse", e)["metrics"]["angle_median"] for e in ("P1", "FD")}
    assert result("refolded_sparse", "DFI")["metrics"]["angle_median"] < 0.25 * plain["P1"]
    assert result("refolded_sparse", "FDFold")["metrics"]["angle_median"] < 0.4 * plain["FD"]


def test_refolded_fold_event_matches_truth():
    """The fold event itself (no solve): expected normal vs true gradient. The floor of any solve.

    About 2 deg median. The limb angle is a function of gx only, the exact one also depends on the plunge.
    """
    case = geomodels.build("refolded_sparse")
    P = harness.probe_points()
    P = P[harness.eval_mask(P, case)]
    fd, axis, dgz = engines.make_fold_event(case).get_deformed_orientation(P)
    ang = harness.angle_deg(np.cross(fd, axis), case.truth_grad(P))
    ang = np.minimum(ang, 180 - ang)  # the sign of the normal is not set by the fold event
    assert np.median(ang) < 3.0 and ang.max() < 6.0
    assert np.all(np.sum(dgz * case.truth_grad(P), axis=1) < 0), "dgz must point against the normals"


def test_p2_not_worse_than_p1_on_smooth_dome():
    assert result("dome", "P2")["metrics"]["value_rms"] <= result("dome", "P1")["metrics"]["value_rms"]


# --------------------------------------------------------------------------- golden
@pytest.mark.parametrize("case,engine", PAIRS, ids=IDS)
def test_matches_golden_field(case, engine, request):
    r = result(case, engine)
    if request.config.getoption("--update-golden"):
        harness.save_golden(case, engine, r)
        pytest.skip("golden updated")
    golden, _ = harness.load_golden(case, engine)
    if golden is None:
        pytest.fail("no golden file: run with --update-golden and review the figures first")
    ok = np.allclose(r["value"], golden["value"], atol=1e-6, rtol=0, equal_nan=True) and np.allclose(
        r["gradient"], golden["gradient"], atol=1e-5, rtol=0, equal_nan=True
    )
    if not ok:
        from . import render  # needs matplotlib

        fig = render.render_diff(case, engine, r["points"], golden, r)
        dv = np.nanmax(np.abs(r["value"] - golden["value"]))
        pytest.fail(f"field changed, max |dv| = {dv:.3g}. Diff figure: {fig}")


def _compare_family(name, now, ref):
    problems = []
    for key in ("shape", "nnz"):
        if now[key] != ref[key]:
            problems.append(f"{name}.{key}: {ref[key]} -> {now[key]}")
    for key in ("A_sum", "A_fro", "b_sum", "b_norm", "w_min", "w_max", "WA_fro", "A_rowmoment"):
        # A_rowmoment is a signed sum that cancels to about 1e-9, so it carries
        # platform noise (numpy/BLAS ordering). A real row reorder moves it by far more.
        atol = 1e-6 if key == "A_rowmoment" else 1e-10
        if not np.isclose(now[key], ref[key], rtol=1e-8, atol=atol):
            problems.append(f"{name}.{key}: {ref[key]:.10g} -> {now[key]:.10g}")
    return problems


@pytest.mark.parametrize("case,engine", PAIRS, ids=IDS)
def test_matches_golden_constraints(case, engine, request):
    if request.config.getoption("--update-golden"):
        pytest.skip("golden updated by the field test")
    _, meta = harness.load_golden(case, engine)
    if meta is None:
        pytest.fail("no golden file: run with --update-golden")
    now, ref = result(case, engine)["fingerprint"], meta["constraints"]
    problems = []
    for name in sorted(set(now) | set(ref)):
        if name not in ref:
            problems.append(f"new constraint family: {name}")
        elif name not in now:
            problems.append(f"constraint family disappeared: {name}")
        else:
            problems += _compare_family(name, now[name], ref[name])
    assert not problems, "constraint matrices changed:\n  " + "\n  ".join(problems)


# --------------------------------------------------------------------------- known bugs
# Each test states the CORRECT behaviour. It fails today and is marked xfail(strict=True).
# When a fix lands the test XPASSes, strict turns that into a failure: remove the marker,
# then run --update-golden and review the diff figures.


@pytest.mark.parametrize("engine", ["P1", "P2"])
def test_default_solve_is_converged(engine):
    default, converged = result("planar_dip", engine, True), result("planar_dip", engine)
    assert np.nanmax(np.abs(default["value"] - converged["value"])) < 1e-3


def test_p2_normal_weight_has_effect():
    case = geomodels.build("planar_noisy")
    out = []
    for npw in (1.0, 100.0):
        engines.reseed()
        it = engines.make_interpolator("P2", case)
        engines.add_data(it, case)
        it.setup_interpolator(npw=npw)
        it.solve_system(engines.SOLVER, solver_kwargs=engines.SOLVER_KWARGS)
        out.append(it.evaluate_value(harness.probe_points()))
    assert np.nanmax(np.abs(out[0] - out[1])) > 1e-4


@pytest.mark.parametrize("engine", ["DFI", "FDFold"])
def test_update_keeps_fold_weights(engine):
    case = geomodels.build("fold_sparse")
    engines.reseed()
    it = engines.make_interpolator(engine, case)
    engines.add_data(it, case)
    it.setup_interpolator(fold_weights={"fold_orientation": 50.0})
    before = {k: float(np.mean(v["w"])) for k, v in it.constraints.items() if "orientation" in k}
    it.up_to_date = False
    it.update()
    after = {k: float(np.mean(v["w"])) for k, v in it.constraints.items() if "orientation" in k}
    assert before == after


def test_dfi_is_repeatable_without_reseed():
    """Fixed in the dev tree: DFI no longer shuffles rows (review section 3, item 6)."""
    case = geomodels.build("fold_sparse")
    P = harness.probe_points()
    out = []
    for _ in range(2):  # no reseed on purpose
        it = engines.make_interpolator("DFI", case)
        engines.add_data(it, case)
        it.setup_interpolator()
        it.solve_system(engines.SOLVER, solver_kwargs=engines.SOLVER_KWARGS)
        out.append(it.evaluate_value(P))
    assert np.nanmax(np.abs(out[0] - out[1])) < 1e-9


def test_dfi_is_repeatable_with_reseed():
    a, b = harness.compute("fold_sparse", "DFI"), harness.compute("fold_sparse", "DFI")
    assert np.allclose(a["value"], b["value"], atol=1e-9, equal_nan=True)


def test_p2_mesh_stores_its_geometry():
    from loop_common.supports import P2UnstructuredTetMesh

    mesh = P2UnstructuredTetMesh(origin=np.zeros(3), step_vector=np.ones(3) / 4, nsteps_cells=np.full(3, 4))
    assert hasattr(mesh, "origin")  # inside() and the setup diagnostics need it


@pytest.mark.xfail(strict=True, reason="new finding F2: P1/DFI return NaN for points on the lower x/y/z faces")
@pytest.mark.parametrize("engine", ["P1", "DFI"])
def test_known_points_on_lower_boundary_are_finite(engine):
    case = "fold_sparse"
    r = result(case, engine)
    on_lower = np.any(np.isclose(r["points"], 0.0), axis=1)
    assert np.all(np.isfinite(r["value"][on_lower]))
