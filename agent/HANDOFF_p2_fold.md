# P2 fold handoff: phase 4 validation

## Current status - 20 October 2026

Phases 0-3 are implemented on `master`. The original prototype narrative below
is retained as historical context and is **superseded** for current results and
weighting rules by [the corrected validation report](../validation/p2_fold/README.md).
The library class and factory registration are in use; an engine swap is no
longer required.

The phase 4 decision gate passed on both targets with the same fixed upstream
frames, seed 12345, explicit target buffer 0, pinned CG settings, and S0
`limb_wl=None`:

| Target | Method | Unknowns | Median axial error (deg) | Above 30 deg (%) |
|---|---|---:|---:|---:|
| S2 | DFI | 12,992 | 6.15 | 0.009 |
| S2 | P2Fold gauss4 | 9,711 | 4.69 | 0.000 |
| S0 | DFI | 112,560 | 20.63 | 37.282 |
| S0 | P2Fold gauss4 | 85,673 | 11.14 | 19.379 |

The report also records the requested full P2 meshes (S2 1e4, S0 1e5), both
sampling rules, tets, solve time, peak resident memory, input hashes and figures.
The full S0 gauss4 run has 765,347 unknowns, median 8.69 deg, 16.1% above 30 deg,
92.5 s solve time and about 6.6 GiB peak RSS. The remaining S0 error tail matters.
Unregularised centroid is unstable (about 60 deg median), even with converged CG.

Corrected rules: use `w*sqrt(q_i)` after row normalisation, not `w/n_samples`.
The common volume factor cancels and is removed. The gauss4 squared residual
weight is four times the prototype's at the same named family weight. Actual
residuals can still differ between sampling rules. `nelements` is currently a
structured-grid cell budget; realised P2 tets and unknowns exceed that request.
The older statement that it is a literal tet count is inaccurate.

These results are target-interpolator comparisons on a fixed **DFI** upstream
chain. The old all-P2 S2/S0 table changed the S0 driver and profile, so its
6.20 deg / 0.94% S0 result is not a direct before/after weight comparison.
Field scale remains arbitrary. Buffer, F2, F3 and profile-fit issues remain open.
Documentation and a bundled-data example are prepared; merge and post-merge
release checks precede Phase 5. Root dependency ranges currently exclude the
pending 0.1.0 workspace-package releases and need alignment before publication.

---

## Historical prototype handoff (superseded)


## Goal
Make the quadratic (P2) tetra interpolator usable for fold modelling (`create_and_add_folded_fold_frame`).
Today P2 has no fold constraints, and the model factory forces `"DFI"` (linear) for folded features.

## Why do it
Tests on the Bermagui S2 and S0 folded foliations (`debug_s0/`). The metric is the angle between the unit
gradient of the solved field and the direction expected from the fold constraints, `cross(fold_direction, fold_axis)`.

| Mesh | Unknowns | Median error |
|---|---|---|
| Linear DFI | 12k | 6 to 12 deg |
| P2, centroid constraint only | 9k to 28k | 9 to 17 deg (does not improve with refinement) |
| P2, 4-point quadrature | 9k to 28k | 4 to 8 deg |

P2 with quadrature is better than linear DFI at the same or fewer unknowns. The hinge-zone error nearly disappears.
All values are from single runs. The model chain varies by about 3 deg between runs, so compare only inside one run.

## What the prototype does
File: `debug_s0/scripts/p2_fold.py`, class `P2FoldInterpolator(P2Interpolator)`.
It uses the loop_interpolation package from `~/dev/LoopStructural` (version 1.8.x).

1. Put sample points inside every tet. The best rule is the degree-2 rule, 4 points with barycentric coordinates
   (a, a, a, b), a = 0.1381966011250105, b = 0.5854101966249685.
   The centroid alone is not enough: the P2 gradient is linear inside a tet, so one point does not control it.
   Using all 10 nodes is slightly worse (shared vertices get too much weight).
2. Get the shape-function derivatives with `support.evaluate_shape_derivatives(points, elements=idx)`. Shape: `(n, 3, 10)`.
3. Add three constraint families, as in `DiscreteFoldInterpolator`:
   - gradient . fold_direction = 0 (w = 10)
   - gradient . fold_axis = 0 (w = 10)
   - gradient . dgz = fold_norm (w = 1, default fold_norm = -1 in the new version)
4. Row = `einsum(dN, vector) * tet_volume`, right side = `target * tet_volume`.
   Give each sample weight `w / n_samples`, so the total weight per tet does not depend on the rule.
5. Regularisation: leave it off or use the P2 default. Fold-frame versions did not help (see below).

The prototype is attached to the builder by swapping `builder.interpolator` and `feature.interpolator` after
`create_and_add_folded_fold_frame`. That is a workaround. The real change should go in the library.

## Suggested library work
- Add a class (for example `P2FoldInterpolator`) in `packages/loop_interpolation/src/loop_interpolation/`.
  Use `_discrete_fold_interpolator.py` (shared setup in `_fold_setup.py`) and `_fd_fold_interpolator.py` as templates.
- Register a new `InterpolatorType`, string key, `interpolator_map` entry and `support_interpolator_map` entry
  (support = `P2UnstructuredTetMesh`, built from `origin`, `step_vector`, `nsteps_cells`) in `loop_interpolation/__init__.py`.
- Stop forcing `"DFI"`: `LoopStructural/modelling/core/_model_feature_factory.py` (about line 248),
  `geological_model.py` (about line 784), `modelling/features/builders/_folded_feature_builder.py` (about line 75).
- Make `nelements` clear for P2. 1e4 gives 56k tets and 86k unknowns. 1e5 gives 526k tets and 765k unknowns (solve about 2.5 min).
- Expose the sample rule (`centroid` or `gauss4`) in `fold_weights`. Default `gauss4`.
- Return NaN, not `LinAlgError ("Singular matrix")`, for points outside the mesh in `evaluate_shape_derivatives`.
- Add tests (see below) and a short doc note.

## Tests to add
- The quadrature rule integrates a quadratic exactly on a reference tet.
- Constant fold angle, no orientation data: median angle to the expected direction below a fixed limit.
- At equal unknowns, P2 error is not worse than linear DFI on a synthetic fold.
- Fold norm sign: use `resolve_fold_norm_target` as `DiscreteFoldInterpolator` does.

## Findings to keep in mind
- Constraint rows are normalised by `add_constraints_to_least_squares`. The size of `fold_direction` does not matter
  (the "sum, not norm" scaling in old versions is harmless). Weights `w` are applied after the normalisation.
- The three fold constraints fix the gradient completely. Regularisation weights 0.1 and 0.01 change nothing.
  Large weights (about 100 to 1000) make the fold fit worse. Fold-frame face-jump and curvature terms gave
  at most a small gain over the isotropic terms. Do not add them for now.
- The norm constraint is almost never satisfied, and the field scale is arbitrary (mean gradient 0.1 to 0.5).
  The direction metric is not affected, but a wavelength in S2 units is not portable between runs.
- `fold_norm` default changed from +1 (old) to -1 (new). With +1, the old version gave a field with
  opposite polarity to the data normals (80% of S2 data flipped).
- `buffer` is applied twice (model box, then support). `buffer=0.2` gives a mesh about 12 times the box volume,
  and only 7% of the elements are inside the box. `buffer=0.0` roughly halves the error at equal element count.
  With `buffer=0.0` the field is valid only inside the box and boundary points return NaN.
- Builds are lazy. Evaluate a feature before you read `fold` attributes.
- Some steps use random shuffles. Do not compare numbers between separate runs.

## Separate issue: rotation-angle profile fit (not P2)
- The fit changes between runs. 2 of 8 rebuilds gave a degenerate profile (angle 0 everywhere).
- `limb_wl` is only the optimiser start. Values 2 to 5 stay in a poor local minimum
  (misfit 40 to 48 deg). Values 8 to 50 and auto converge to the same fit (misfit 11 deg).
- Full P2 run (S2 1e4, S0 1e5, 765k unknowns), S0 `limb_wl=5` against `limb_wl=None` (auto):
  profile misfit 36 deg against 12.7 deg; S0 median error 11.8 against 6.2 deg; share above 30 deg 19.5% against 0.9%;
  roughness 44 against 30 deg per 10 m; data match 0.77 (2 of 7 flipped) against 0.92 (none flipped).
  Plots: `debug_s0/images/p2_s2_s0_rotation_profiles.png` and `..._wlauto.png`.
- S2 observations contradict each other (9 deg and 69 deg only 0.7 S3 units apart), so its profile fits poorly.
- Always plot the profile with its observations: `debug_s0/scripts/p2_s2_s0.py` does this.

## Reproduce
- Environment: a venv with numpy, pandas, scipy, scikit-image, scikit-learn, pydantic, pyyaml, pyvista, dill, matplotlib
  (and trame, trame-vtk, trame-vuetify, trame-pyvista, nest_asyncio2 for HTML export).
- Set `PYTHONPATH` to `~/dev/LoopStructural;~/dev/LoopStructural/packages/loop_common/src;~/dev/LoopStructural/packages/loop_interpolation/src`.
- Scripts in `debug_s0/scripts/`:
  - `p2_fold.py` is the prototype class and the comparison runs.
  - `s0_p2.py` compares linear DFI and P2 on S0.
  - `p2_s2_s0.py` is S2 (P2, 1e4) and S0 (P2, 1e5) with profile plots (`--wl=auto` for the automatic wavelength).
  - `s0_wavelength_sweep.py` is the `limb_wl` sweep.
  - `compat_chain.py` runs the same chain in the old (1.6.27) and new (1.8.x) versions.
- Results: `debug_s0/newversion/*.csv` and `*.log`, images in `debug_s0/images/`, interactive pages in `debug_s0/html/`.
- Data: `Data/Processed/TrimDataS4.csv`, `S3`, `S2`, `S0_LA2`. The model box is origin (238275, 5965000, -50), maximum (238500, 5965250, 50).

## Caveats
- Single runs. The metric measures agreement with the fold constraints, not geological quality.
- The prototype has no anisotropic fold regularisation, and it was tested only on this one dataset (5 S2 and 7 S0 points).
