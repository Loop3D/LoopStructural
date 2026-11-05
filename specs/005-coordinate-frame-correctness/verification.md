# Verification: 005-coordinate-frame-correctness

This file records the evidence for FR-029, FR-031 and FR-032 of spec 005.

## Pre-fix failures

Method. A throw-away worktree at `ebc7c699` (the commit before the spec) held
only the 4 new test files and `tests/fixtures/coordinate_frames.py`. The
`tests/__init__.py` and `tests/fixtures/__init__.py` files already exist at
`ebc7c699`, so no `__init__.py` was copied. The run used
`PYTHONPATH=<tmp>;<tmp>/packages/loop_common/src;<tmp>/packages/loop_interpolation/src`.
Before each run, `python -c "import loop_common, LoopStructural"` showed that
both imports came from the throw-away worktree. The worktree is removed
(`git worktree remove --force`, then `git worktree prune`).

| Test file | Result at `ebc7c699` | Time |
|---|---|---|
| `packages/loop_common/tests/test_rigid_transform_contract.py` | Collection error: `ImportError: cannot import name 'RigidRotationError' from 'loop_common.geometry'`. No test runs. | 2 s |
| `tests/unit/modelling/test_reference_constraint_frames.py` | 278 failed, 34 passed | 7 s |
| `tests/unit/modelling/test_fault_frame_analytic.py` | 53 failed, 0 passed | 88 s |
| `tests/unit/modelling/test_fault_buffer_bounds.py` | 84 failed, 0 passed | 12 s |

### Rigid-transform contract (no R1 fixture; FR-003, FR-013, FR-014, FR-021; SC-001, SC-011, SC-012)
The file cannot import pre-fix, because `RigidRotationError` does not exist.
All 42 contract checks are "not runnable" pre-fix. None is "reproduced" or
"not reproduced" by a test result. The missing class is itself the evidence.

### Reference constraints (R3; FR-003 to FR-007, FR-009, FR-013, FR-024; SC-007, SC-008)
| Test function | Failed / passed | Representative failing id and value |
|---|---|---|
| `test_world_local_conversion_methods` | 24 / 0 | `[frame_case0-FDI]`: `AttributeError: 'FiniteDifferenceInterpolator' object has no attribute 'world_to_local_points'` (same for `P1Interpolator`) |
| `test_reference_gradient_samples_use_world_points_and_local_vectors` | 66 / 6 | first mismatch at `[0, 0]`: -3.807 (actual) vs -2.770 (expected), atol 1e-10 |
| `test_orthogonal_rows_affine_residual` | 48 / 24 | `[frame_case1-None-FDI]`: 210 of 210 residuals about -1.687 vs 0.0 |
| `test_inequality_feature_bounds_world` | 44 / 4 | `[frame_case1-True-FDI]`: bound -4.9999 vs -4.0176 (336 of 336 differ) |
| `test_equality_region_selects_world_nodes` | 24 / 0 | `'reference_equality'` is not in the constraint dictionary (rows are under `'undefined'`) |
| `test_restricted_reference_exclusion_counts` | 24 / 0 | `'restricted_equality'` is not in the constraint dictionary |
| `test_set_interpolation_geometry_encloses_rotated_box` | 24 / 0 | `[frame_case0-FDI]`: grid maximum [6.4, 8.33, 10.29] does not enclose the box corner [0, 10, 12]; also `ValueError: Cannot resize the grid` in some cells |
| `test_support_expansion_contains_restored_points` | 24 / 0 | `[frame_case0-FDI]`: grid maximum [27.2, 8.67, 13.71] does not enclose the restored points |

### Fault frame analytic (R2, R4; FR-001, FR-002, FR-010, FR-011, FR-026; SC-003, SC-013)
All 5 single tests fail. The 48 planar cells all fail (12 for each of
FDI/PLI x centre mode).

| Test function | Failed / passed | Representative failing id and value |
|---|---|---|
| `test_inferred_centre_keeps_elevation` | 1 / 0 | inferred centre [100, 200, 0] vs expected [100, 200, 300] |
| `test_explicit_centre_rows_not_projected` | 1 / 0 | no explicit centre row equals the centre [100, 200, 300]; the rows were projected |
| `test_slip_vector_made_orthogonal` | 1 / 0 | slip [0.7071, 0, -0.7071] vs expected [0, 0, -1] |
| `test_parallel_slip_raises` | 1 / 0 | `DID NOT RAISE ValueError` |
| `test_extent_direction_row_at_centre` | 1 / 0 | the count of rows at the centre is 0, expected 1 |
| `test_planar_fault_frame_matches_analytic` | 48 / 0 | `[FDI-True-0.0-anchor0]`: maximum field error 0.830 (limit 1e-6). Error range over cells: 1.7e-3 to 0.84. 12 cells show exactly 0.45. |

### Fault buffer bounds (R5; FR-012, FR-013; SC-005)
| Test function | Failed / passed | Representative failing id and value |
|---|---|---|
| `test_buffer_bounds_match_request` | 36 / 0 | `[0.0-True-0.0-FDI-20000]`: lower bound 0 vs the allowed -4 (bounds do not cover the buffer) |
| `test_buffers_give_distinct_bounds` | 12 / 0 | `[True-0.0-FDI-20000]`: same bound failure |
| `test_property_change_keeps_buffer_and_forced_geometry` | 36 / 0 | `ValueError: Cannot resize the grid` ("proposed number of steps is [-2 13 13]") and bound assertion failures |

### Findings that are "not reproduced"
None in the three runnable files: each test function has at least one
failing cell. The partial passes (for example 24 of 72 in
`test_orthogonal_rows_affine_residual`) are cells where the case does not
exercise the defect.

## S2 pre-fix failures

Method. A throw-away worktree at `9e7b05f7` (the S2 start commit) held the 3 S2
test files and `tests/fixtures/coordinate_frames.py`. The
`test_reference_constraint_frames.py` file is the S1 file with the P2 cells; it
ran with `-k P2`. The run used
`PYTHONPATH=<tmp>;<tmp>/packages/loop_common/src;<tmp>/packages/loop_interpolation/src`.
Before the runs, `loop_interpolation`, `loop_common` and `LoopStructural` were
imported from the throw-away worktree. The worktree is removed
(`git worktree remove --force`, then `git worktree prune`).

| Test file | Result at `9e7b05f7` | Time |
|---|---|---|
| `test_fold_constraint_frames.py` | 274 failed, 98 passed | 16.7 s |
| `test_constant_norm_frames.py` | 13 failed, 13 passed | 5.2 s |
| `test_reference_constraint_frames.py -k P2` | 156 passed, 0 failed | 9.1 s |

### Fold constraints (FR-008; SC-007, SC-008, SC-013)
| Test function | Failed / total | Representative failing id and value |
|---|---|---|
| `test_fold_geometry_sampled_at_world_points` | 44 / 48 | `[frame_case1-fdi_fold]`: 672 of 1008 differ; `[0, 0]`: -5.0 vs -4.330 (atol 5e-9) |
| `test_fold_rows_affine_residual` | 32 / 48 | `[frame_case1-fdi_fold]`: 336 of 336 residuals differ; first -0.843 vs 0.0 (atol 1e-10) |
| `test_fold_rows_follow_curved_directions_fdi` | 11 / 12 | `[frame_case1-fdi_fold]`: 336 of 336 differ; first -0.0 vs -0.559 |
| `test_fold_rows_follow_curved_directions_pli_p2` | 33 / 36 | `[frame_case1-pli_fold]`: 3360 of 4200 differ; first 0.664 vs 0.223 |
| `test_fold_mask_receives_world_points` | 44 / 48 | `[frame_case1-fdi_fold]`: 672 of 1008 differ; first -5.0 vs -4.330 |
| `test_fold_norm_alignment_uses_world_normals` | 16 / 48 | `[frame_case2-fdi_fold]`: `assert [1.0] == [-1.0]` (the sign flips in the wrong frame) |
| `test_directional_regularisation_world_contract` | 38 / 72 | `[frame_case1-fdi-array]`: 672 of 1008 differ; first 0.828 vs 0.532 |
| `test_fold_regularisation_world_contract` | 44 / 48 | `[frame_case1-fdi_fold-nomask]`: 336 of 336 differ; first 1.0 vs 0.829 |
| `test_p2_steepness_wtfunc_receives_world_barycentres` | 12 / 12 | all frame cases: shape (1050, 3) vs expected (690, 3) |

The passing cells (98 of 372) are cells of the unrotated anchor, or cells that
do not exercise the defect. Each test function has at least one failing cell.

### Constant norm (FR-008; SC-002)
| Test function | Failed / total | Representative failing id and value |
|---|---|---|
| `test_unrotated_anchor_zero_control` | 1 / 2 | `[cn_fdi]`: `IndexError: invalid index to scalar variable.` |
| `test_constant_norm_fit_recovers_r1_in_rotated_frame` | 12 / 24 | all 12 `cn_fdi` cells: the same `IndexError`; all 12 `cn_pli` cells pass |

This is an existing defect, not a frame defect. `ConstantNormFDIInterpolator`
raised `IndexError` for every frame, because `StructuredGrid.element_size` is a
scalar. It also fails on the unrotated control. Task T032 fixes it. The PLI
constant-norm control meets SC-002 pre-fix (value 1.8e-12, gradient 1.3e-10).
So `cn_pli` is "coverage, no defect".

### Reference constraints on P2 (FR-008; SC-007, SC-008)
All 156 P2 cells pass pre-fix. The S1 changes already give the P2 reference
path the world contract. These cells are "coverage, no defect": they show that
P2 has the same behaviour as FDI and PLI, and that the S2 fixtures work for P2.

### Findings that are "not reproduced"
None in the fold file and the `cn_fdi` cells. The `cn_pli` cells and the 156 P2
reference cells are "coverage, no defect". They are not defects that the
fix removed.

## Post-fix results

Spec worktree, HEAD `10fd4dea`, worktree `PYTHONPATH`. A check of
`loop_common.__file__` and `LoopStructural.__file__` showed that both imports
come from `C:\Users\lachl\dev\LoopStructural-005`.

| Command | Result |
|---|---|
| `python -m pytest packages/loop_common/tests packages/loop_interpolation/tests -p no:cacheprovider -q` | 912 passed, 47 skipped, 2 xfailed, 7 warnings in 142 s |
| The 4 S1 test files in one run | 491 passed, 7 warnings in 291 s |
| `test_rigid_transform_contract.py` | 42 passed |
| `test_reference_constraint_frames.py` | 312 passed |
| `test_fault_frame_analytic.py` | 53 passed |
| `test_fault_buffer_bounds.py` | 84 passed |
| `ruff check LoopStructural packages` (orchestrator) | All checks passed |
| Root `pytest` (orchestrator, full) | 1204 passed, 5 skipped, 21 warnings in 1082.66 s. Baseline at `3881ef32`: 755 passed, 5 skipped, 21 warnings in 1205.98 s. |

Test-setup corrections made during the run (facts for the record). No
tolerance was relaxed.
- T003: the round-trip test passes `dimensions=dimensions`.
- T011: the slip-row expectation is unit slip / intermediate_axis (spec R2 table).
- T007: the affine residual computes the exact field without large-coordinate
  cancellation. The atol of 1e-10 is kept.
- T011: the helper passes the SC-003 CG settings through `solver_kwargs`.
  Top-level values were ignored, so the solver used atol 2e-9.

## S2 post-fix results

Spec worktree `feat/005-s2-fold-frames`, HEAD `662d2713`, worktree `PYTHONPATH`.
A check showed that `loop_interpolation.__file__` and `LoopStructural.__file__`
come from `C:\Users\lachl\dev\LoopStructural-005-s2`.

| Command | Result |
|---|---|
| `test_fold_constraint_frames.py` | 372 passed in 12.2 s |
| `test_constant_norm_frames.py` | 26 passed in 8.8 s |
| `test_reference_constraint_frames.py` (S1 file, now with P2 cells) | 468 passed (312 S1 cells and 156 P2 cells) in 18.9 s |
| The 3 S2 files in one run | 866 passed in 35.6 s |
| Existing package tests `test_fd_fold_interpolator.py` (27), `test_discrete_fold_interpolator.py` (21), `test_p2_fold_interpolator.py` (20), `test_regularisation_api.py` (7) | pass with no change; `test_native_parity.py` is skipped (see Backend matrix) |
| `ruff check LoopStructural packages` (orchestrator, HEAD `662d2713`) | All checks passed |
| `pytest packages/loop_common/tests packages/loop_interpolation/tests` (orchestrator) | 912 passed, 47 skipped, 2 xfailed, 7 warnings in 126.56 s |
| Root `pytest` (orchestrator, full) | 1760 passed, 5 skipped, 13 warnings in 1052.40 s. Baseline at `9e7b05f7`: 1206 passed, 5 skipped, 21 warnings in 937.41 s. The 554 new passes are the S2 cells. |

Test groups and the requirements that they check:
- Fold rows, fold mask, fold norm alignment, fold geometry, fold and
  directional regularisation, P2 steepness weights: FR-008; SC-007, SC-008,
  SC-013.
- Constant norm: FR-008; SC-002.
- P2 cells of the reference file: FR-008; SC-007, SC-008.
- Tolerances: SC-015. No tolerance was relaxed in S2.

Test-setup correction (fact for the record). The first version of
`_set_world_normal` in the fold test file projected the normal row twice. It was
fixed before T025 to T028 were accepted. The pre-fix run above uses the
corrected file.

## Golden review (C-004)

No golden change in S1 for `tests/unit/modelling/test__fault_builder.py`.
No assertion in that file encodes an S1 defect. The file checks the
`update_geometry` bounds, and the fault normal and slip vectors. It does not
check a fault centre or a default buffer. The commented centre assertions
stay unchanged.

One golden change is made in S1 (task T006, commit d2909aa5).

| Item | Old value | New value | Defect | Expected physical effect |
|---|---|---|---|---|
| `packages/loop_common/tests/data/structured_support_reference.npz`: 27 values `<cls>/<grid>/after_origin.nsteps`, `after_maximum.nsteps`, `after_maximum.maximum` for cls StructuredGrid, TetMesh, P2TetMesh and grid unit, aniso, small | unit: origin nsteps [4 4 4], maximum nsteps [6 6 6], maximum [6.142857 x3]. aniso: [3 6 4], [5 7 5], [14.6 6.909091 7.857143]. small: [5 5 5], [8 8 8], [5.0 x3] | unit: [5 5 5], [8 8 8], [9.0 x3]. aniso: [4 7 5], [7 9 7], [17.4 10.545455 11.285714]. small: [6 6 6], [10 10 10], [6.714286 x3] | The origin and maximum setters stored a cell count as the node count. The support lost one cell and did not reach the requested maximum (FR-014). Old values did not enclose the requested bounds (for example unit: requested maximum [7.714286 x3]). | The support now encloses the requested bounds, within one cell, outward. |

The other keys of the file keep their values from before S1 (the 33 keys that `--update` re-recorded with round-off changes were restored after the spec review). No tolerance was relaxed.

## S2 golden review (C-004)

No golden change in S2. `git diff --stat 9e7b05f7..HEAD -- tests packages/*/tests`
shows only these files:
- `tests/fixtures/coordinate_frames.py` (new helpers, +85 lines).
- `tests/unit/modelling/test_fold_constraint_frames.py` (new, 411 lines).
- `tests/unit/modelling/test_constant_norm_frames.py` (new, 69 lines).
- `tests/unit/modelling/test_reference_constraint_frames.py` (S1 file, +31 -13). The change adds the `P2` engine to the builder helper and to the `engine` parametrisation. No assertion and no tolerance changed.

No existing test file and no data file in `packages/*/tests` changed.

## Backend matrix

S1 tests use Python backends only.

| Requirement | Backends | Cells |
|---|---|---|
| R2 fault frame analytic | Python FDI, Python PLI (P1) | 4 anchors x 3 rotations x 2 centre modes x 2 backends = 48 planar cells, plus 5 single tests |
| R3 reference constraints | FDI, PLI (P1) | reference file, sample rules None / centroid / gauss4, all frame cases (312 cells) |
| R5 fault buffers | FDI, PLI | buffers 0, 0.2, 0.5 x 3 rotations x 2 centre modes (84 cells) |

No native, browser or RBF backend is in S1 scope. FR-030 states that such
backends are not a prerequisite.

## S2 backend matrix

| Requirement | Backends | Cells |
|---|---|---|
| Fold rows, mask, norm alignment, geometry, regularisation (FR-008; SC-007, SC-008, SC-013) | Python FDI-fold, PLI discrete fold, P2-fold (`fdi_fold`, `pli_fold`, `p2_fold`, with the sample rules in the test ids) | `test_fold_constraint_frames.py`: 372 cells over the frame cases |
| P2 steepness weight function (SC-013) | P2 | 12 cells |
| Reference constraints on P2 (FR-008; SC-007, SC-008) | P2 (`P2Interpolator`, `P2TetMesh`) | 156 cells (`-k P2`) |
| Constant norm (FR-008; SC-002) | `cn_fdi` (`ConstantNormFDIInterpolator`), `cn_pli` (`ConstantNormInterpolator`) | 26 cells: 4 anchors x 3 rotations x 2 backends, plus 1 unrotated zero control for each backend |

Native. `packages/loop_interpolation/tests/test_native_parity.py` is skipped,
because the native extension `loop_interpolation_core` is not installed. The
adapter `_native.py` has the same world contract (T024), but no test covers it
under a rotated frame. The native path rejects rotated supports. This is an
open coverage gap. FR-030 states that native is not a prerequisite.

## Rollback checks

On a throw-away worktree at `ebc7c699` (imports checked as above):
`python -m pytest tests/unit/modelling/test__fault_builder.py packages/loop_common/tests/test_bounding_box.py -p no:cacheprovider -q`
gives 15 passed in 0.5 s. The two files pass on the pre-spec code, so the
rollback target is sound.

## S2 rollback checks

At the S2 start commit `9e7b05f7` the full root suite passed:
`pytest` gives 1206 passed, 5 skipped, 21 warnings in 937.41 s (run by the
orchestrator on 2026-11-04; collected before the wave 8 edits). This run includes all S1 tests. The S1 tests
therefore pass on the code from before S2, and the rollback target is sound.
The full suite was not run again at the start commit. The S1 test files were
not run alone in a throw-away worktree at `9e7b05f7`.

## Test run time (OP2)

| Item | Time |
|---|---|
| All 4 S1 files, one run | 291 s test time (4 min 57 s wall) |
| Slowest file: `test_fault_frame_analytic.py` | 264 s test time (4 min 31 s wall) |
| `test_fault_buffer_bounds.py` | 20 s (26 s wall) |
| `test_reference_constraint_frames.py` | 5.6 s (12 s wall) |
| `test_rigid_transform_contract.py` | 0.1 s (4 s wall) |

## S2 test run time (OP2)

| Item | Time |
|---|---|
| All 3 S2 files, 866 cells, one run | 35.6 s test time (about 45 s wall, measured while the orchestrator ran the gates in parallel) |
| `test_fold_constraint_frames.py` (372 cells) | 12.2 s |
| `test_reference_constraint_frames.py` (468 cells) | 18.9 s |
| `test_constant_norm_frames.py` (26 cells) | 8.8 s |

The times come from runs under load, so they are an upper bound. The slowest
single cell is 0.65 s.
