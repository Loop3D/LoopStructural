# Golden tests for the interpolation code

Run from `packages/loop_interpolation`, against the source tree (`pip install -e .` or `PYTHONPATH=src:../loop_common/src`).
These are step 1 of the interpolation hardening review (pin behaviour before fixing bugs).
The figures and `--report` need matplotlib and pyvista; the tests themselves need only numpy, scipy and pytest.
Open `--update-golden` diffs in review: the golden files are the Rust parity oracle.
Each model has an analytic truth, so every run is scored against the right answer, not only
against the last run.

## Use

```
python -m pytest tests/golden -q               # about 20 s
python -m pytest tests/golden --update-golden  # accept the current answers
python -m pytest tests/golden --report         # out/test_report.html: every test with its surfaces (about 70 s)
python -m tests.golden.run                     # figures + out/report.html
python -m tests.golden.run --case fold_sparse --engine DFI --vti   # one run + a 3D .vti
```

## Models (unit cube, fixed seeds)

| Case | Truth | Data | What it catches |
|---|---|---|---|
| `planar_dip` | plane, dip 30 | 3 horizons, 6 boreholes, 3 normals | sign, scale, weight and placement bugs. Must be exact. |
| `planar_extrapolation` | same plane | data only for x < 0.35 | regularisation as an extrapolator |
| `planar_noisy` | same plane | z +-0.01, normals +-3 deg | smoothing versus over-fitting |
| `dome` | Gaussian dome, relief 0.2 | 3 horizons, 9 boreholes, 5 normals | crest flattening, curved horizons |
| `fold_section` | cylindrical fold (wavelength 0.6, amplitude 0.12) | one cross-section at y = 0.5 | plain interpolation along the fold axis |
| `fold_sparse` | same fold | 2 value points + 1 normal | fold constraints (`DFI`, `FDFold`, `P2Fold`) against none (`FD`, `P1`, `P2`) |
| `refolded_sparse` | F2 fold (axis y) refolded by F3 (axis x): dome and basin, `s = z - A sin(k2 x) - B sin(k3 y)` | 3 contacts + 6 normals | the `fold_axis_rotation` path: the hinge plunge changes along the axis (as at Bermagui) |

Engines: `FD` (17^3 nodes), `P1` (17^3 nodes), `P2` (8^3 cells, 4.4k unknowns), `DFI` (P1 + fold), `FDFold` (FD + fold), `P2Fold` (P2 + fold, 8^3 cells / 4401 unknowns).
Fold frame is analytic (gx = x, gy = y), fold axis (0, -1, 0). With that axis `dgz` points against the
normals, which is the library convention (`fold_norm = -1`).

### Refolded fold: limits of the model

The fold event takes the limb angle as a function of `gx` only. For this surface the exact angle is
`tan(theta) = -h'(x) cos(plunge(y))`, so the constraint direction is about 2 deg off (median, max 4.7)
even before any solve. `test_refolded_fold_event_matches_truth` pins that floor. Results: median gradient
error is 5 deg (DFI) and 8 deg (FDFold) against 33 to 40 deg without fold constraints.

In this case the fold constraints override the sparse orientation data. FDFold ends 40 deg from one of the six
normals, and DFI 14 deg. `test_data_honoured` allows this for `refolded_sparse` only. Worth a look: a data
normal that is ignored by 40 deg is a weighting question (review item 9).

## What is pinned

`golden_files/<case>__<engine>.npz` value and gradient on a 15^3 probe grid.
`golden_files/<case>__<engine>.json` metrics against truth, and a fingerprint of every constraint
family (shape, nnz, sums, norms, row-order moment, weights). These are the matrices the Rust port must match.

Solves are pinned to CG with `maxiter=20000, rtol=1e-10`. The package default stops early (see K1).

## Reading the figures (`out/figures/<case>.png`, one row per engine)

1. Field on the XZ section at y = 0.5. Dashed black = true horizons, orange = solved horizons,
   white dots = value data near the section, black arrows = normals near the section.
2. Gradient direction error in degrees. 0 is good, 180 is reversed (stratigraphic inversion).
3. Unit gradient arrows. Solved (blue) over true (grey).
4. Plan view at z = 0.5. Signed value error. Blue = too low, red = too high.

A failing field test also writes `out/failures/<case>__<engine>__diff.png`: where the field moved
and by how much, against the golden file. A failing constraint test names the family whose matrix changed.

## Findings from building these tests (all confirmed by running code)

Bugs from the review that the suite demonstrated. The tests are now normal tests (`test_default_solve_is_converged`,
`test_p2_normal_weight_has_effect`, `test_update_keeps_fold_weights`, `test_dfi_is_repeatable_*`).

| Test | Review item | Status |
|---|---|---|
| `test_default_solve_is_converged` | bug 3 | Fixed (Phase 0.4). CG returns `False` at the iteration limit. Default `rtol` is 1e-8: P2 `planar_dip` error 0.65 -> 2.5e-4, P1 1.8e-2 -> 5e-7. |
| `test_p2_normal_weight_has_effect` | bug 4 | Fixed (Phase 0.3). `npw` 1 -> 100 moves the P2 field by 7.4e-4 (was 1e-7). |
| `test_update_keeps_fold_weights` | bug 11 | Fixed (Phase 0.5). `update()` replays the kwargs of the last `setup_interpolator` call. |
| `test_dfi_is_repeatable_without_reseed` | section 3, item 6 | Fixed. DFI no longer shuffles rows. |
| (no golden test; see `test_p0_nan_constraints_skipped.py`) | bug 6 | Fixed (Phase 0.6). NaN rows are dropped, not the whole family. |
| `TestP2GradSteepnessMask` (`test_p2_interpolator.py`) | bug 5 | Fixed (Phase 0.2). Mask is `neighbours >= 0`. P2 fields moved by at most 2.3e-3. |

New findings, not in the review:

- **F1.** (Fixed, Phase 0.1.) `P2UnstructuredTetMesh(origin=..., step_vector=..., nsteps_cells=...)` did not store those
  attributes, so `inside()` and the diagnostics report raised `AttributeError`.
- **F2.** `TetMesh` (P1 and DFI) returns NaN for points on the lower x, y or z faces (strict `>` in `inside`).
  631 of 3375 probe points. FD and P2 return values. The golden files keep the NaN. Metrics use the inner box.
- **F3.** Fold amplitude is damped. The median ratio of solved to true gradient magnitude is 0.66 (`fold_section`),
  0.57 (DFI) and 0.35 (FDFold) on `fold_sparse`. The direction is good, the scale is not.
  Not tested because it may be intended; `grad_norm_ratio_median` is in the metrics. Relevant to the
  "decide the intended energy first" step.
- In `fold_sparse`, DFI and FDFold get the horizon shape right but not the height. The plan-view value error
  is up to about 0.05 (DFI) and 0.08 (FDFold), and it follows the fold wavelength. I did not find the cause.

## Adding a case

Add a function to `geomodels.py` returning a `Case` and add it to `ALL`. Run
`python -m tests.golden.run --case <name>`, read the figure, then `pytest --update-golden`.
Choose bounds in `test_golden.py::BOUNDS` about 1.5x above what the code gives.

## Not covered yet

Faults (equality constraints, review bug 1), inequality constraints (bug 2), tangent and interface
constraints, regularisation variants, additional `fold_weights` options,
`mask_fn`, and meshes that are not unit cubes. DFI on a cube mesh only: no graded or unstructured mesh.


## P2 fold coverage (Phase 3)

P2Fold is pinned on `fold_sparse` and `refolded_sparse`, including fields,
gradients and all constraint-family fingerprints. Only the four new P2Fold
baseline files were generated; existing golden files were preserved.
The figures were reviewed: fold shape and gradient directions agree well, with
a remaining vertical offset in the refolded section and spatially varying value
error. Neither gauss4 result has inverted gradients in the inner evaluation box.
Bounds are approximately 1.5 times the observed RMS and median direction errors;
the refolded bound stays above the fold event's approximately 2 degree model floor.

| Case | Rule | Value RMS | Median gradient error | P90 gradient error |
|---|---|---:|---:|---:|
| fold_sparse | centroid | 0.00329 | 1.364 deg | 5.303 deg |
| fold_sparse | gauss4 | 0.01652 | 0.856 deg | 2.124 deg |
| refolded_sparse | centroid | 0.03068 | 4.422 deg | 12.571 deg |
| refolded_sparse | gauss4 | 0.05507 | 2.179 deg | 5.795 deg |

These comparisons use the same fixed seed, explicit unit-cube mesh bounds
(zero buffer), data, solver tolerance and 4401 unknowns. Gauss4 improves gradient
direction in both cases, now checked by a regression test. It does **not** improve
value RMS or fit to the sparse normals; no universal accuracy claim is made.
Centroid also has one inverted inner probe on the refolded case.
The DFI comparison uses a separate budget-matched mesh with 15 x 15 x 16
cells (4352 unknowns) versus P2Fold's 4401, a 1.13% difference. These mesh
families do not have identical node counts. Both supports have explicit unit-cube
bounds and no buffer, asserted in the test, and use the same fixed seed. The
existing, larger DFI golden mesh and its baseline are unchanged.

Unit tests pin total weighted operator energy (squared Frobenius norm) for each
fold family between quadrature rules, weight sensitivity (10 to 50), fold norm
sign against DFI and `resolve_fold_norm_target`, update persistence, constant-angle
folds without orientation data, and dropped NaN row counts. The energy test does
not assert equal energy for every field: normalising the shape-derivative rows
changes the residual at different sample locations. There is no LinearOperator
path here, so the optional adjoint test does not apply.
Model integration checks folded-fold-frame construction without DFI fallback and
preservation of the P2Fold interpolator type after save/load.


Validation on Windows / Python 3.13 (shared environment, source paths explicitly
pointed at this worktree): interpolation package **577 passed, 34 skipped,
2 xfailed in 56.82 s**; golden suite **143 passed, 9 skipped, 2 xfailed in
27.92 s**; P2Fold model integration **5 passed in 11.55 s**. The final golden run was standalone; the package and integration runs overlapped
with the root suite. These are wall-clock observations rather than a benchmark.
For sandboxed runs, pass `--basetemp=<worktree-local-directory>` when the default
pytest temporary directory is inaccessible. Existing lower-boundary NaNs remain
strict expected failures. The new P2Fold cases add two small 8^3-cell solves,
plus weight-sensitivity, quadrature-comparison and budget-matched DFI checks.

The full root suite also passed: **727 passed, 3 skipped, 21 warnings in
407.56 s (6 min 48 s)**, using a worktree-local pytest temporary directory.
No test failures remain.
