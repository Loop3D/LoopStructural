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

## S3 pre-fix failures

Method. HEAD was `caf0a0e7`. The two new S3 test files ran one at a time with
`python -m pytest -p no:cacheprovider <file> -q --tb=line`, with the worktree
`PYTHONPATH`. Before the first run, `loop_interpolation.__file__` and
`LoopStructural.__file__` were both checked: they come from
`LoopStructural-005-s3`. No S3 library code exists at this commit. Only the stubs
and the two error classes exist. The persistence file gives 166 failed and
3 passed. The copy and lifecycle file gives 43 failed and 6 passed.

| Test function | Failed / total | Representative id and first error |
|---|---|---|
| `test_payload_declares_frame` (persistence) | 36 / 36 | `[fdi-a0-r0]`: `KeyError: 'format_version'` (missing name) |
| `test_round_trip_world_values` | 33 / 36 | `[fdi-a0-r30]`: world error 1.9011 against a limit of 4.4e-12. 3 identity-frame cells pass (coverage, no defect) |
| `test_local_arrays_unchanged_after_two_cycles` | 36 / 36 | `[fdi-a0-r0]`: `KeyError: 'format_version'` (missing name) |
| `test_json_and_yaml_carry_frame` | 3 / 3 | `[fdi]`: `KeyError: 'transform'` (missing name) |
| `test_legacy_payload_requires_choice` | 3 / 3 | `[fdi]`: `DID NOT RAISE InterpolatorFormatError` |
| `test_legacy_with_supplied_transform` | 36 / 36 | `[fdi-a0-r0]`: `TypeError: from_dict() got an unexpected keyword argument 'local_origin'` (missing name) |
| `test_legacy_with_identity_declaration` | 3 / 3 | `[fdi]`: `TypeError: from_dict() got an unexpected keyword argument 'declare_identity'` (missing name) |
| `test_conflicting_inputs` | 4 / 4 | `[identity-and-origin]`: `TypeError: from_dict() got an unexpected keyword argument 'declare_identity'` (missing name) |
| `test_malformed_metadata` | 8 / 8 | `[missing-transform]`: `TypeError: create_interpolator() got an unexpected keyword argument 'format_version'` (missing name) |
| `test_unsupported_future_version` | 2 / 2 | `[3]`: `TypeError: create_interpolator() got an unexpected keyword argument 'format_version'` (missing name) |
| `test_creation_call_unchanged` | 1 / 1 | `TypeError: InterpolatorFactory.from_dict() got an unexpected keyword argument 'declare_identity'` (missing name) |
| `test_old_reader_would_name_format_version` | 1 / 1 | `StopIteration`: the payload has no `format_version` key (missing name) |
| `test_copy_keeps_frame` (copy and lifecycle) | 36 / 36 | `[fdi-anchor0-rot0]`: `assert None is BoundingBox(...)`, the copy has no bounding box |
| `test_stale_frame_rejected_without_state_change` | 5 / 5 | `[fdi]`: `DID NOT RAISE StaleFrameError` |
| `test_restored_transform_is_usable` | 0 / 1 | Passes now (coverage, no defect), 1 cell |
| `test_rebuild_in_new_frame_preserves_world_evaluation` | 0 / 5 | Pass now (coverage, no defect), 5 cells |
| `test_fault_restoration_uses_current_frame` | 2 / 2 | `[FDI]`: `DID NOT RAISE StaleFrameError` |

Behavioural records from the production path (a short script, not saved in the
repository):

| Record | Measured value |
|---|---|
| (1) `fdi`, `ANCHORS[3]`, rotation 30, `exact=True`: `from_dict(to_dict())` against the original, at `ANCHORS[3] + QUERY_OFFSETS` | The original gives `[3.05, 4.4]`. The loaded interpolator gives `[nan, nan]`. E = NaN (2 of 2 values are NaN) |
| (2) `interp.copy().bounding_box is None` for `fdi` at the same cell | `True` |
| (3) `fdi`, `ANCHORS[1]`, rotation 0, `exact=False`: `evaluate_value` before and after `bb.set_local_transform(local_origin=ANCHORS[1] + (3, -2, 1), rotation_matrix=rot_z(30.0))` | No error is raised. The value before is `[3.05, 4.4]`. The value after is `[nan, nan]`. D = NaN (the value changes from finite to NaN) |

The identity-frame round-trip cells and the rebuild cells pass before S3. They are coverage, not defect evidence.

Since C-018, `test_legacy_payload_requires_choice` is replaced by `test_legacy_payload_without_choice_warns_and_loads_identity`, `test_legacy_identity_values_match_source_for_identity_frame` and `test_strict_legacy_switch_raises`. The pre-fix record above keeps the old name.

## S3 rollback checks

Method. The rollback target is the S3 start commit `91ffca4a` (the S2 head, with no S3 code). A throw-away worktree was made at this commit with `git worktree add --detach ../LoopStructural-005-s3-rollback 91ffca4a`. The file `tests/fixtures/coordinate_frames.py` was copied from the S3 worktree into it, and nothing else, because the old commit has no `r1_interpolator`. Each run used its own `PYTHONPATH`. Before each run, the script printed `loop_interpolation.__file__` and `LoopStructural.__file__`. For the old code, both paths were in `LoopStructural-005-s3-rollback`. For the S3 code, both paths were in `LoopStructural-005-s3`. The scripts were short and were not saved in the repository. The throw-away worktree was removed at the end.

Check 1. The old `to_dict` saved a legacy payload (no `format_version`) for `fdi`, `pli` and `p2` with `r1_interpolator(engine, anchor, rotation, exact=True)`. The S3 code loaded each payload with `GeologicalInterpolator.from_json`. It used the `local_origin` and `local_rotation` of the S3 interpolator of the same cell, or `declare_identity=True`. The errors are the maximum absolute difference from the fresh S3 interpolator at `anchor + QUERY_OFFSETS`. The limits (SC-006) are `1e-12 * F` = 4.4e-12 for values (F = 4.4) and `1e-12 * max(1, abs(R1_GRADIENT).max())` = 1.0e-12 for gradients.

| Engine | Cell | Load mode | Value error | Gradient error |
|---|---|---|---|---|
| `fdi` | `ANCHORS[3]`, rotation 30 | supplied transform | 0.0 | 0.0 |
| `pli` | `ANCHORS[3]`, rotation 30 | supplied transform | 0.0 | 0.0 |
| `p2` | `ANCHORS[3]`, rotation 30 | supplied transform | 0.0 | 0.0 |
| `fdi` | `ANCHORS[0]`, rotation 0 | supplied transform | 0.0 | 0.0 |
| `pli` | `ANCHORS[0]`, rotation 0 | supplied transform | 0.0 | 0.0 |
| `p2` | `ANCHORS[0]`, rotation 0 | supplied transform | 0.0 | 0.0 |
| `fdi` | `ANCHORS[0]`, rotation 0 | `declare_identity=True` | 0.0 | 0.0 |
| `pli` | `ANCHORS[0]`, rotation 0 | `declare_identity=True` | 0.0 | 0.0 |
| `p2` | `ANCHORS[0]`, rotation 0 | `declare_identity=True` | 0.0 | 0.0 |

All 9 errors are inside the limits. With C-018, a legacy payload loaded with no choice (no transform and no `declare_identity`) loads with the identity frame and emits one `FutureWarning` for `fdi`, `pli` and `p2` (3 of 3). The message starts: "legacy interpolator payload has no transform; assuming the identity frame". The loaded interpolator has no bounding box. For the rotated cell (`ANCHORS[3]`, rotation 30), the world values are NaN. This is the behaviour of releases before S3, and it is the reason for the warning. With `_STRICT_LEGACY` set to true, the call raises `legacy_no_transform` (the behaviour of the next major release).

Check 2. The S3 code wrote `r1_interpolator("fdi", ANCHORS[3], 30.0, exact=True)[0].to_json()` (it has `format_version`). The old reader (`GeologicalInterpolator.from_json(text)` at `91ffca4a`) raised `TypeError` with this exact message: `InterpolatorFactory.create_interpolator() got an unexpected keyword argument 'format_version'`. The call returned no object: the result variable stayed unset (OP5). The error names `format_version`, as expected.

## S3 post-fix results

Spec worktree `feat/005-s3-frame-persistence`, HEAD `67b1e48a`, worktree
`PYTHONPATH`. A check showed that `loop_interpolation.__file__` comes from
`C:\Users\lachl\dev\LoopStructural-005-s3`. The memory on the machine was low.
The orchestrator ran the full gates. The implementer did not run them again.

| Command | Result |
|---|---|
| `ruff check LoopStructural packages` (run by the orchestrator) | All checks passed |
| `pytest packages/loop_common/tests packages/loop_interpolation/tests -p no:cacheprovider -q` (run by the orchestrator) | 912 passed, 47 skipped, 2 xfailed, 7 warnings in 108.69 s. The S3 tests are in the root suite. The result equals the S2 result (912/47/2). |
| Root `pytest -p no:cacheprovider -q`, full (run by the orchestrator) | 1978 passed, 5 skipped, 13 warnings in 797.81 s. Before S3 (S2 head `662d2713`): 1760 passed, 5 skipped, 13 warnings in 1052.40 s. The 218 new passes equal the 218 new S3 cells (169 in `test_frame_persistence.py`, 49 in `test_frame_copy_and_lifecycle.py`). |
| `test_frame_persistence.py` (run by the implementer) | 169 passed in 5.20 s |
| `test_frame_copy_and_lifecycle.py` (run by the implementer) | 49 passed in 31.27 s |

No `StaleFrameError` appeared in the full root run. No library path raised it.
The check is not weakened.

Requirements and the tests that check them:

| Requirement | Tests |
|---|---|
| FR-016 (frame and format version in the payload) | `test_payload_declares_frame`, `test_json_and_yaml_carry_frame`, `test_old_reader_would_name_format_version` |
| FR-017 (no second world conversion on load) | `test_round_trip_world_values`, `test_local_arrays_unchanged_after_two_cycles` |
| FR-018 (copy keeps the frame) | `test_copy_keeps_frame` |
| FR-019 (legacy file needs a choice; malformed, unsupported and conflicting input fail) | `test_legacy_payload_requires_choice`, `test_legacy_with_supplied_transform`, `test_legacy_with_identity_declaration`, `test_conflicting_inputs`, `test_malformed_metadata`, `test_unsupported_future_version` |
| FR-020 (transform change: reject or rebuild) | `test_stale_frame_rejected_without_state_change`, `test_restored_transform_is_usable`, `test_rebuild_in_new_frame_preserves_world_evaluation`, `test_fault_restoration_uses_current_frame` |
| FR-021 (rotation check) | Checked in S1 (`test_rigid_transform_contract.py`). The S3 files test the saved and supplied rotations in `test_conflicting_inputs` and `test_malformed_metadata`. |
| FR-028 (public API kept) | `test_creation_call_unchanged`; the existing tests pass with no change (see the S3 golden review) |
| SC-006 (save and load accuracy; arrays unchanged after 2 cycles) | `test_round_trip_world_values`, `test_local_arrays_unchanged_after_two_cycles` |
| SC-012 (1 failure, 0 state change; legacy needs a choice) | `test_legacy_payload_requires_choice`, `test_conflicting_inputs`, `test_malformed_metadata`, `test_stale_frame_rejected_without_state_change` |
| SC-016 (origin and rotation change: reject or rebuild) | `test_stale_frame_rejected_without_state_change`, `test_rebuild_in_new_frame_preserves_world_evaluation` |
| SC-017 (old-format input keeps its load path; new format fails in the old reader with a clear error) | `test_legacy_with_supplied_transform`, `test_legacy_with_identity_declaration`, `test_unsupported_future_version`, `test_old_reader_would_name_format_version` |

Cost of the frame check. A script outside the repository built
`r1_interpolator("fdi", ANCHORS[3], 30.0, exact=True)` and 100 points
(`ANCHORS[3] + uniform(1, 7)`, seed 0). It timed 1000 calls, 5 repeats, and kept
the minimum. The "check removed" runs set
`GeologicalInterpolator._check_frame = lambda self: None` in the same process.

| Call | With the check | Check removed | Difference |
|---|---|---|---|
| `evaluate_value` | 166.3 us per call | 160.3 us per call | 6.0 us per call (about 4 %) |
| `evaluate_gradient` | 346.5 us per call | 343.4 us per call | 3.1 us per call (about 1 %) |

These numbers are for 100 points. The measure runs in one process and has some
noise. The check does not loop over the points. Its cost is therefore a fixed
cost for each call.

003 key-time check. The 003 plan has a key-time gate (T042, the time of
`fault_set_fingerprint`). It measures the key of the fault evaluation cache. It
is a different measure from the frame check of 005. The implementer did not run
it, because of the memory limit. There is no 003 measure of the frame check.

## S3 golden review (C-004)

Command: `git diff --stat 91ffca4a269252c32cefe077261edd27f10a4b14..HEAD -- tests packages/*/tests`.
Result: 3 files, 598 insertions, 0 deletions.

| File | Change |
|---|---|
| `tests/fixtures/coordinate_frames.py` | 79 lines added (helpers for the S3 tests) |
| `tests/unit/modelling/test_frame_copy_and_lifecycle.py` | new file (188 lines) |
| `tests/unit/modelling/test_frame_persistence.py` | new file (331 lines) |

No golden change in S3. No existing test file changed.

Command: `git grep -n "from_dict\|from_json\|from_yaml" -- tests packages/*/tests`.
Existing tests that load a saved interpolator payload:
- `packages/loop_interpolation/tests/test_geological_interpolator.py`: `test_interpolator_json_yaml_round_trip` loads text from `to_json()` and `to_yaml()` (new format). `test_geological_interpolator_from_dict_delegates_to_factory` patches `InterpolatorFactory.from_dict` and loads nothing.
- `packages/loop_interpolation/tests/test_p2_fold_registration.py`: `test_to_dict_from_dict_round_trip` loads a payload from `to_dict()` (new format).
- `tests/unit/interpolator/test_interpolator_factory.py`: `test_from_dict_builds_interpolator` and `test_from_dict_missing_type_raises` call `InterpolatorFactory.from_dict` with a hand-written creation dict (`type`, `boundingbox`, `nelements`). This is the factory creation call, not a saved interpolator payload. The creation call is unchanged (`test_creation_call_unchanged`).

The other matches load other objects (`LoopEntity`, `PointSet`, `FaultTopology`, the stratigraphic column types, the model recipe). They are not interpolator payloads.

No existing test loads a payload without `format_version`.

## S3 backend matrix

Engines: `fdi`, `pli`, `p2`, `cn_fdi`, `cn_pli`, `p2_fold`. Native engines and
Surfe are in S6. PLI has no `copy()`. Fold families do not save the fold object.
The `cn_*` engines have lifecycle cells only.

| Engine | persistence | copy | stale-frame lifecycle |
|---|---|---|---|
| `fdi` | tested | tested | tested |
| `pli` | tested | not applicable (no `copy()`) | tested |
| `p2` | tested | tested | tested |
| `cn_fdi` | not applicable (lifecycle cells only) | not applicable | tested |
| `cn_pli` | not applicable (lifecycle cells only) | not applicable | tested |
| `p2_fold` | not applicable (the fold object is not saved) | tested (the copy shares the fold object) | not applicable (copy cells only) |

OP6 (T037). `git grep -n 'interpolator.*copy(' -- LoopStructural/modelling/features`
finds only two numpy `.copy()` calls, on `support.origin` and `support.maximum`
in `_geological_feature_builder.py` (lines 541 and 542). No call of
`interpolator.copy()` exists in the feature code.

OP7 note. The P2 value-only fit of R1 misses SC-002 in every frame. In the
unrotated control the value error is 7.8e-8 against a limit of 7.5e-8, and the
gradient error is 3.9e-7 against a limit of 1e-8. At 30 degrees the value
error is 2e-7 to 7.6e-7. FDI and PLI stay below 3.2e-9 (value) and 1.7e-9
(gradient). Issue #399 records it. The S3 rebuild test uses the exact nodal
control for `p2` (maintainer decision 2026-11-05).

## S3 test run time (OP2)

| Item | Time |
|---|---|
| `test_frame_persistence.py` (169 cells) | 5.20 s |
| `test_frame_copy_and_lifecycle.py` (49 cells) | 31.27 s (the fault restoration test takes about 25 s) |
| Both S3 files | about 36 s test time |
| Full root suite (orchestrator) | 797.81 s (1052.40 s at the S2 head; the two runs had different machine load, so the times do not compare) |
| Package suites (orchestrator) | 108.69 s |

## S3 change C-018 (backward-compatible legacy load)

Date: 2026-11-05.

Behaviour. The legacy load stays backward compatible on the 1.x line. A legacy
payload (no `format_version`) with no choice of the caller loads with the
identity frame, as older releases did. It emits one `FutureWarning`. The warning
names the choices and says that the next major release raises an error. The
switch `_STRICT_LEGACY` in `loop_interpolation._frame_persistence` makes the
call raise `legacy_no_transform`. A partial transform, a conflict, a malformed
payload and an unsupported version still raise `InterpolatorFormatError`.

Tests (focused runs of the orchestrator):

| Command | Result |
|---|---|
| `tests/unit/modelling/test_frame_persistence.py` | 184 passed (after the review fix below) |
| `test_frame_copy_and_lifecycle.py -k "not fault_restoration"` | 48 passed |
| packages `test_geological_interpolator.py` and `test_p2_fold_registration.py` | 43 passed |
| `ruff` | pass |

The full suites are run again after this change (see the next section).

## S3 final gates after C-018

The orchestrator ran the full gates on the working tree after the C-018 change (HEAD `dddaff50` plus the review fix).

| Command | Result |
|---|---|
| `ruff check LoopStructural packages` | All checks passed |
| `pytest packages/loop_common/tests packages/loop_interpolation/tests -p no:cacheprovider -q` | 912 passed, 47 skipped, 2 xfailed, 7 warnings in 103.85 s |
| Root `pytest -p no:cacheprovider -q`, full | 1994 passed, 5 skipped, 13 warnings in 817.74 s |

Before C-018 the root suite gave 1981 passed. The difference of 13 comes from `test_frame_persistence.py`: the 3 cells of `test_legacy_payload_requires_choice` became 9 cells (the warn test, the identity-values test and the strict-switch test, 3 engines each), and 7 cells came with the review fix below. The file now has 184 cells.

Review fix. The first review of the C-018 change found that the warning pointed at a library file for `from_json` and `from_yaml`, and one frame too high for a direct `InterpolatorFactory.from_dict` call. It also found that the warning came before a later failure of the same call. The fix moves the warning to `warn_legacy_identity()`, called as the last step of a successful load. The function computes the `stacklevel` as the first frame outside the `loop_interpolation` package. The tests `test_legacy_warning_points_at_caller` (4 entry points) and `test_legacy_warning_not_emitted_when_load_fails` (3 engines) failed before the fix (6 of 7 cells) and pass now. The gates above ran after the fix.

No `StaleFrameError` and no `FutureWarning` from library code appeared in the full run. The only `FutureWarning` for a legacy load comes from the tests that load a legacy payload with no choice.

## S4 pre-fix failures

Method. The S4 consumer test file ran on the unfixed library code. No S4 library change exists at this commit. Another worker adds the deprecation tests (T056) to the same file at the same time, so the run used a filter and not the whole file. The exact command was `python -m pytest -p no:cacheprovider tests/unit/modelling/test_world_location_consumers.py -k "feature_data_locations or example_observation or exported_data or intrusion_classification or fault_ellipse" -q --tb=line`, with the worktree `PYTHONPATH` (`LoopStructural` came from `LoopStructural-005-s4`). The result is 30 failed and 5 passed. The `rescale` and `scale` tests did not run.

| Test function | Failed / total | Representative id and first error |
|---|---|---|
| `test_feature_data_locations_are_world` | 11 / 12 | `[0-30.0]`: values differ from the input rows. `[0, 0]` is 1.3660 against 1.0 (limit 4.7e-08). The 11 failing cells share this error type (the locations are converted a second time). `[0-0.0]` passes (anchor 0, no rotation: the extra conversion is the identity; coverage, no defect) |
| `test_spec_example_observation_is_not_shifted` | 3 / 3 | `[0.0]`: `assert 374.1657 <= 4.7e-08`. The cells `[30.0]` and `[90.0]` fail in the same way (374.32 and 376.07) |
| `test_exported_data_keeps_world_locations` | 4 / 4 | `[0]`: values differ from the input rows. `[0, 0]` is 3.2321 against 2.0. The cells `[1]`, `[2]` and `[3]` fail in the same way (the CSV files hold converted locations) |
| `test_intrusion_classification_uses_world_points` | 12 / 12 | `[0-0.0]`: `AttributeError: 'IntrusionFrameBuilder' object has no attribute '_network_model_values'` (missing name). All 12 cells share this error |
| `test_2d_viewer_fault_ellipse_centre_is_world` | 0 / 4 | Passes now (coverage, no defect), 4 cells |

Totals by group: the first three tests give 18 failed and 1 passed. The intrusion test gives 12 failed. The fault ellipse test gives 4 passed.

Behavioural measure from the production path (a short script, not saved in the repository). It used `world_foliation_model` with `ANCHORS[1]` = (100, 200, 300). It measured the distance between the observation `anchor + (2, 3, 4)` = (102, 203, 304) and the nearest location in the gradient object of `feature.get_data()`.

| Rotation (deg) | Distance | Nearest location |
|---|---|---|
| 0 | 374.17 | (202, 403, 604) |
| 30 | 374.32 | (289.83, 324.80, 604) |
| 90 | 376.07 | (303, 98, 604) |

At rotation 0, the nearest location is the observation plus the anchor (102 + 100, 203 + 200, 304 + 300). The expected distance is 0 (limit 4.7e-08). The shift is about 374 length units in each cell.

## S4 post-fix results

Spec worktree `LoopStructural-005-s4`, HEAD `2e84bfcc` plus the S4 working tree
changes. The worktree `PYTHONPATH` was set. A check showed that
`LoopStructural.__file__` comes from `LoopStructural-005-s4`. The memory on the
machine was low, so each command ran alone and no full suite ran in this task.
The orchestrator runs the full gates in T068, and T068 adds the gate counts to
this section.

Pre-fix result of T056. Before the fixes, 25 cells of the deprecation tests
failed. `test_rescale_warns_and_keeps_behaviour` (12 cells) and
`test_scale_warns_and_keeps_behaviour` (12 cells) failed with "DID NOT WARN".
`test_internal_callers_do_not_warn` failed with `AttributeError:
'IntrusionFrameBuilder' object has no attribute '_network_model_values'`.

| Command | Result |
|---|---|
| `pytest -p no:cacheprovider tests/unit/modelling/test_world_location_consumers.py -q` | 59 passed, 1 skipped in 2.55 s. The skip is `test_internal_callers_do_not_warn`: `pyevtk` is not in the venv. |
| Same file, with `pyevtk` on `PYTHONPATH` (see below) | 60 passed in 2.28 s |
| `pytest -p no:cacheprovider tests/unit/modelling/intrusions/test_intrusions.py -q` | 9 passed, 2 warnings in 19.81 s (the warnings are an old `setup()` deprecation in the test, not from S4) |
| `pytest -p no:cacheprovider tests/unit/io/test_exporters.py -q`, with `pyevtk` on `PYTHONPATH` | 12 passed, 1 failed in 0.12 s (pre-existing, see below) |
| `pytest -p no:cacheprovider tests/unit/test_public_api_contract.py -q` | 2 passed in 0.02 s |

The new test file has 60 cells: `test_feature_data_locations_are_world` 12,
`test_spec_example_observation_is_not_shifted` 3,
`test_exported_data_keeps_world_locations` 4,
`test_intrusion_classification_uses_world_points` 12,
`test_2d_viewer_fault_ellipse_centre_is_world` 4 (coverage, no defect),
`test_rescale_warns_and_keeps_behaviour` 12, `test_scale_warns_and_keeps_behaviour` 12,
`test_internal_callers_do_not_warn` 1. The fault ellipse test passed before and
after the fixes: coverage, no defect.

`pyevtk` fact. `pyevtk` is not in the venv. For `test_exporters.py` and for
`test_internal_callers_do_not_warn`, the run added
`C:/Users/lachl/AppData/Local/Temp/claude/C--Users-lachl-dev-dev-tools/f206cbcd-e5d0-45be-a4c8-5f5c89bb26f1/scratchpad/pyevtk_pkg`
to `PYTHONPATH`. Without it, the file skips the `pyevtk` tests.

Failure in `test_exporters.py`: `test_write_vol_gocad_is_broken_for_real_bounding_box`.
Pre-existing, not S4. The test expects `pytest.raises(LoopValueError)` with the
class from `LoopStructural.utils.exceptions`. The error comes from
`loop_common.geometry._bounding_box`, which raises its own, different
`LoopValueError` class (it derives from `ValueError`, not from the
`LoopStructural` class). The classes differ, so `pytest.raises` does not match.
The message is "slice(None, None, None) is not a valid bounding box name". The
S4 changes do not touch this path. The run with the S4 changes stashed was not
done (stash is not allowed in this task); the reason comes from reading the
test and the failure text. No existing test was changed.

Deviations.
- T059 uses `stacklevel=3`, because the `public_api` wrapper adds one frame.
- The T056 `scale` check for `inplace` records the existing behaviour: `scale`
  copies its input, and the array of the caller does not change. FR-033 says
  keep the behaviour.

S4 gate counts (T068). At HEAD `1ba71f17` (S4 code final), the root `pytest`
gives 2066 passed, 1 failed, 4 skipped. The one failure is the pre-existing
`test_write_vol_gocad_is_broken_for_real_bounding_box` above. The package tests
give 912 passed, 47 skipped, 2 xfailed. See `## S4 final gates`.

## S4 golden review (C-004)

Command: `git grep -l "get_data\|\.rescale(\|\.scale(" -- 'tests/*.py' 'packages/*.py'`.
The test files in the result were run one at a time with the worktree
`PYTHONPATH`:

| File | Result |
|---|---|
| `tests/unit/modelling/test_world_location_consumers.py` | 60 passed (new file) |
| `tests/unit/interpolator/test_geological_interpolator.py` | 24 passed in 1.08 s |
| `tests/unit/io/test_geoh5.py` | 1 skipped in 0.05 s |
| `tests/unit/utils/test_helper.py` | 9 passed in 0.08 s |
| `tests/unit/modelling/test_geological_feature_builder.py` | 18 passed in 0.14 s |
| `packages/loop_interpolation/tests/test_geological_interpolator.py` | 40 passed in 6.27 s |
| `packages/loopstructural_visualisation/tests/test_plotting.py` | 16 passed in 5.07 s |
| `tests/unit/modelling/intrusions/test_intrusions.py` | 9 passed in 19.81 s |
| `tests/unit/test_public_api_contract.py` | 2 passed in 0.02 s |
| `tests/unit/io/test_exporters.py` (with `pyevtk`) | 12 passed, 1 failed (pre-existing, not S4) |

No existing test needs a new expected value because of S4. The only failure,
`test_write_vol_gocad_is_broken_for_real_bounding_box`, has a cause that S4 does
not touch (two `LoopValueError` classes).

No golden change in S4.

## S4 rollback checks

Method. The rollback target is the S4 start commit `9b3507a9` (`start_commit` of stage S4 in `run-state.json`; old code, no S4 code). A throw-away worktree was made at this commit with `git worktree add --detach ../LoopStructural-005-s4-rollback 9b3507a9a4d4c7d7a903a08c7d33619a7546bf61`. Nothing was copied into it. The S4 version of `tests/fixtures/coordinate_frames.py` was copied to a scratch directory outside the repository and imported from there (the old commit has no `world_foliation_model`). Each tree had its own `PYTHONPATH` (tree, `packages/loop_common/src`, `packages/loop_interpolation/src`, `packages/loopstructural_visualisation/src`). Before each run the script printed `LoopStructural.__file__` and checked that it is inside the tree. For the old code the path was in `LoopStructural-005-s4-rollback`. For the S4 code (HEAD `686322e9`) the path was in `LoopStructural-005-s4`. The script printed `loopstructuralvisualisation.__file__` as well; it was inside the same tree in each run.

| Check | Old code (`9b3507a9`) | S4 code | Result |
|---|---|---|---|
| (a) `from loopstructuralvisualisation import Loop2DView, Loop3DView` | works | works | pass |
| (a) `inspect.signature(Loop2DView.add_data)` | `(self, feature, val=True, grad=True, unfault=False, dip=True, **kwargs)` | same text | pass (equal) |
| (a) `Loop3DView.add_data` | does not exist | does not exist | see note |
| (b) first gradient location of `world_foliation_model(ANCHORS[1], 0)`, `get_data()` | `(202, 403, 604)` | `(102, 203, 304)` | pass (spec example) |
| (c) `get_data()` with identity frame (`ANCHORS[0]`, rotation 0): 3 objects, locations and vectors | reference | max abs difference 0.0 for `r1+_gradient` (4 rows), `r1_norm` (4 rows), `r1_value` (27 rows) | pass (equal) |
| (d) `model.rescale(p)` warnings | 0 | 1: `GeologicalModel.rescale is deprecated ...` | pass |
| (d) `model.scale(p)` warnings | 0 | 1: `GeologicalModel.scale is deprecated ...` | pass |
| (d) values of `rescale` and `scale` | `rescale(p)` and `scale(p)` for `p = [[102, 203, 304], [110, 215, 320]]` | equal to the old values (exact) | pass |

Note on (a). The task text names `Loop3DView.add_data`. This method does not exist in either tree (`hasattr` is false). `add_data` exists only in `Loop2DView`. In `Loop3DView` the only `add_*` method defined in the class is `add_mesh(self, *args, **kwargs)`, and it is equal in both trees. The check used the signatures of `Loop2DView.add_data` and `Loop3DView.add_mesh`. The task text must name `Loop2DView.add_data`.

Check (b) used `ANCHORS[1]` = (100, 200, 300), rotation 0, and the gradient constraint at `anchor + (2, 3, 4)`. The old code shows the point shifted by the origin; the S4 code shows it at the world location.

Cleanup. The throw-away worktree was removed with `git worktree remove --force C:/Users/lachl/dev/LoopStructural-005-s4-rollback`, then `git worktree prune` ran. `git worktree list` shows no `-rollback` entry.

## S4 backend matrix

Engines and consumers that S4 touches. Native engines, Surfe and the browser
viewer are not in S4.

| Item | Level | Evidence |
|---|---|---|
| Python FDI foliation data (`get_data` locations, exported data) | mandatory | `test_feature_data_locations_are_world`, `test_spec_example_observation_is_not_shifted`, `test_exported_data_keeps_world_locations` |
| Fault centre in the 2D viewer | coverage (no defect found) | `test_2d_viewer_fault_ellipse_centre_is_world` passed before and after the fixes |
| Intrusion classification | mandatory | `test_intrusion_classification_uses_world_points` (12 cells) |
| Native engines | not in S4 | S6 |
| Surfe | not in S4 | S6 |
| Browser viewer | not in S4 | S6 |

## S4 test run time (OP2)

| Item | Time |
|---|---|
| Full root suite (orchestrator, HEAD `1ba71f17`) | 805.50 s (13 min 25 s). The S3 end run took 817.74 s. |
| Package suites (orchestrator) | 150.24 s. The S3 end run took 103.85 s; the two runs had different machine load, so the times do not compare. |

The root time is about the same as the S3 time. The root suite has 72 more
passed tests than at the S3 end: the 60 cells of the consumer file and the
exporter tests that now run, because `pyevtk` is present.

## S4 final gates

The orchestrator ran the full gates at HEAD `1ba71f17` (S4 code final). The
`PYTHONPATH` held the worktree, `packages/loop_common/src`,
`packages/loop_interpolation/src`, `packages/loopstructural_visualisation/src`
and a `pyevtk` 1.7.0 directory (the venv has no `pyevtk`).

| Command | Result |
|---|---|
| `ruff check LoopStructural packages` | All checks passed |
| `pytest packages/loop_common/tests packages/loop_interpolation/tests` | 912 passed, 47 skipped, 2 xfailed, 7 warnings in 150.24 s (same counts as the S3 end) |
| Root `pytest` | 1 failed, 2066 passed, 4 skipped, 12 warnings in 805.50 s |

The one failure is `tests/unit/io/test_exporters.py::test_write_vol_gocad_is_broken_for_real_bounding_box`.
It is pre-existing and not from S4 (two different `LoopValueError` classes; see
`## S4 post-fix results`). S3 skipped it because `pyevtk` was missing. The S3
end had 1994 passed and 5 skipped. The 72 new passes are the 60 consumer cells
and the exporter tests that now run.

Because of this one failure, the "Done when" line of T068 (`pytest` passes) is
not fully met. The failure is not in the scope of S4.

## S4 deprecated-call grep (T068)

`git grep -n "\.rescale("` shows these hits:

- `LoopStructural/modelling/features/_geological_feature.py:366`: the commented block of `get_data`.
- `tests/unit/modelling/test_world_location_consumers.py:243` and `:247`: the deprecation tests.
- `COMPAT.md:17`: the deprecation table row (text only, not a call).
- `specs/005-coordinate-frame-correctness/plan.md` and `tasks.md`: spec text, not code.
- `specs/005-coordinate-frame-correctness/verification.md:568` and `:602`: this record.

`git grep -n "\.scale("` shows these hits:

- `COMPAT.md:18`: the deprecation table row (text only).
- `specs/005-coordinate-frame-correctness/plan.md`, `tasks.md` and `verification.md` (lines 568 and 603): spec text and this record.

No code outside the commented block of `get_data` and the deprecation tests
calls `.rescale(` or `.scale(`. No other hit exists.

## S4 library deprecation-warning check (T068)

Command (with the S4 `PYTHONPATH` and the `pyevtk` directory):
`python -m pytest tests/unit/io/test_exporters.py tests/unit/modelling/intrusions tests/unit/modelling/test_world_location_consumers.py -q -p no:cacheprovider -W "error:.*is deprecated.*:DeprecationWarning"`
with the two warn-tests of `rescale` and `scale` deselected (they expect the warning).

Result: 1 failed, 57 passed, 24 deselected, 2 warnings in 23.43 s. The failure
is the pre-existing `test_write_vol_gocad_is_broken_for_real_bounding_box`. The
summary has two warnings with the text "is deprecated":
`StructuralFrameBuilder.setup() is deprecated, use .build() instead`, from
`tests/unit/modelling/intrusions/test_intrusions.py` lines 55 and 100. The tests
call the old `setup()` method. This warning is old and not from S4. No warning
of `rescale` or `scale` appears, so no S4 code path calls them.

## S4 follow-up: gocad test and viewer name

- `tests/unit/io/test_exporters.py::test_write_vol_gocad_is_broken_for_real_bounding_box` failed only because two classes are named `LoopValueError`. `BoundingBox.__getitem__` raises the class of `loop_common.geometry._bounding_box` (a `ValueError`), and the test expected the class of `LoopStructural.utils.exceptions`. The test now accepts `(LoopValueError, ValueError)`. The file passes: 13 passed (run with `pyevtk` on `PYTHONPATH`). The full root suite was not run again after this change; the failure was the only one in the S4 root run.
- The task text and the S4 notice in `COMPAT.md` named `Loop3DView.add_data`. The method is `Loop2DView.add_data`. `plan.md`, `tasks.md` and `COMPAT.md` now say so.

## S5 pre-fix failures

Method. Both `loop_common` test files ran on the unfixed code, in one call: `python -m pytest packages/loop_common/tests/test_resize_paths.py packages/loop_common/tests/test_cell_centre_sampling.py -p no:cacheprovider -q -rf` (PowerShell, the worktree `LoopStructural-005-s5` on `PYTHONPATH`). No S5 library change exists at this commit. The result is 24 failed and 76 passed. All 24 failures are in the resize file. All 60 cells of the cell-centre file pass.

Resize file (`test_resize_paths.py`): 24 failed, 76 passed.

| Test function | Failed / total | Representative id and first error |
|---|---|---|
| `test_nsteps_setter_keeps_extent` | 3 / 3 | `[G1]`: the `maximum` is [10.5, 11, 10.5] and not [10, 10, 9] after `nsteps` is set to [22, 12, 8]. `[G2]` gives 102.25, -17.75, 7.25 (start 102, -18, 7). `[G3]` gives -9795, 2250, 26.25 (start -9800, 2200, 22.5) |
| `test_step_vector_setter_encloses_extent` | 2 / 3 | `[G1]`: the `maximum` is [10.5, 9, 9] and is below the start [10, 10, 9], so the extent is lost. `[G3]` fails in the same way (the `maximum` is below the start). `[G2]` passes |
| `test_nsteps_below_two_is_rejected_without_change` | 1 / 1 | `DID NOT RAISE ValueError`. The setter accepts nsteps [1, 6, 4] |
| `test_nsteps_not_whole_is_rejected_without_change` | 1 / 1 | `DID NOT RAISE ValueError`. The setter accepts [5.5, 5, 5] and [nan, 5, 5] |
| `test_set_nelements_keeps_extent` | 3 / 3 | `[G1]`: the `maximum` is not the start `maximum`. `set_nelements(200)` gives [9.43, 10.29, 10] and not [10, 10, 9]. `[G2]` gives 102.08, -17.92, 7.08. `[G3]` gives -9811, 2270, 10 |
| `test_set_nelements_small_request_is_clamped` | 1 / 1 | `assert np.all(nsteps >= 2)` fails: `nsteps` is [2, 2, 1], and the `maximum` is [5.5, 6, 0] |
| `test_2d_counts_are_derived` | 10 / 10 | All 5 checks fail for `[H1]` and for `[H2]`: `cells`, `n_cell`, `n_elements`, `dtype`, `maximum`. See the measured values below |
| `test_rectilinear_resize_is_rejected` | 3 / 3 | `[nsteps]`: `DID NOT RAISE NotImplementedError`. `[step_vector]` and `[origin]` fail in the same way |
| `test_reapplying_same_value_10_times_changes_nothing_3d` | 0 / 12 | Passes now (coverage, no defect) |
| `test_step_vector_a_b_a_extent_drift` | 0 / 3 | Passes now (coverage, no defect). See the measured drift below |

Measured values (a short script, not saved in the repository).

| Measure | Value |
|---|---|
| `nsteps` setter, nodes 11 to 21, step 1, origin 0 | `maximum` is 10.476 (10.48) in each axis, and the expected value is 10 |
| `step_vector` setter, step x 1.5 on G1 | `maximum` is [10.5, 9, 9], the start is [10, 10, 9]. The `nsteps` become [8, 4, 3] |
| `step_vector` setter, step x 1.5 on G3 | `maximum` is [-9805, 2150, 26.25], the start is [-9800, 2200, 22.5]. The `nsteps` become [14, 2, 6] |
| `step_vector` setter, step x 1.5 on G2 | `maximum` is [102.25, -17.75, 7.25], the start is [102, -18, 7]. The extent is enclosed, so the test passes |
| A-B-A `step_vector` drift (A, then 1.5 A, then A) | G1 [1, 0, 3]. G2 [0.5, 0.5, 0.5]. G3 [0, 0, 7.5]. All are inside the limit (A + B), so the test passes |
| 2D `nsteps` set to [nx + 3.0, ny + 2.0], H1 (5, 5) | `nsteps_cells` stays [4, 4]. `n_cell_x` and `n_cell_y` stay 4. `n_elements` stays 16. The dtype is float64. The `maximum` stays [5, 5] and the expected value is [8, 7] |
| Same set on H2 (8, 4) | `nsteps_cells` stays [7, 3]. `n_elements` stays 21. The dtype is float64. The `maximum` stays [14, 3] and the expected value is [15.5, 7] |

The 2D `maximum` check fails only because the counts are stale after the set. The legacy formula (origin + nsteps x step) is the expected value in the test, so the fix does not change the formula.

Cell-centre file (`test_cell_centre_sampling.py`): 60 passed, 0 failed. Result: coverage, no defect. OP12 is not reproduced.

Note on the one-to-one test. The first version of the test used `structured_grid(local_coordinates=True)` and failed in 8 of 12 cells (rotation 30 and 90). That grid is the axis-aligned enclosure of the projected corners, and it has a different step. The premise in the task text was wrong. The test now uses the midpoints of the exported world grid and the own cell grid of the box. `_bounding_box.py` was not changed. Also, `set_local_transform` takes the argument `rotation_matrix=`.

## S5 golden review (C-004)

Order (C-004): (1) run the 3D golden test and list the failing keys; (2) check each key against the rule; (3) write the keys here; (4) run `--update`; (5) run the test again. SC-015: no tolerance was changed (`rtol=1e-9`, `atol=1e-12` stay).

Step 1. `test_structured_support_reference.py` failed in 3 cells (StructuredGrid, TetMesh, P2TetMesh). A key-by-key diff (script in the scratchpad, not in the repository) found 78 changed keys of 611. The key sets are equal. Each of the 3 classes has the same 26 changed keys: 3 grids (`unit`, `aniso`, `small`) x (`step_vector`, `to_dict.step_vector`, and `after_set_nelements`, `after_nsteps`, `after_step_vector`, `after_origin`, `after_maximum` `.step_vector`) plus `after_maximum.maximum` and `after_maximum.nsteps` where they differ (`unit` has no `nsteps` change, `aniso` and `small` do). The 2D golden (`test_structured_grid_2d_reference.py`) passes with no change (4 passed): the 2D `maximum` is the same, so no 2D key changed and `--update` was not run for it.

Step 2. Rule: each setter keeps the extent or encloses it by less than one cell. Extent = `(nsteps - 1) x step_vector`. Checked on the `StructuredGrid` keys (the 3 classes give the same numbers):

| Grid | Stage | Old extent | New extent | Result |
|---|---|---|---|---|
| unit | `set_nelements(200)` | 3.571 | 4.000 (start 4) | kept |
| unit | `nsteps` + 1 | 4.286 | 4.000 | kept |
| unit | `step_vector` x 2 | 4.286 | 4.000 | kept |
| unit | origin - 1 | max 4.714 | max 4.333 (1 cell = 1.333; enclosed by 0.333) | enclosed by less than one cell |
| unit | `maximum` + 3 (request 7.333) | 9.000 | 8.333 (over by 1.0, 1 cell = 1.333) | enclosed by less than one cell |
| aniso | `set_nelements(200)` | [12.1, 5.18, 6.29] | [13, 5, 6.5] (start [13, 5, 6.5]) | kept |
| aniso | `nsteps` + 1, `step_vector` x 2 | [12.8, 6.09, 7.14] | [13, 5, 6.5] | kept |
| aniso | origin - 1 | max [13.2, 6.91, 7.86] | max [13.5, 5.6, 7.0] (cells [1.5, 1.6, 1.5]) | enclosed by less than one cell |
| aniso | `maximum` + 3 (request [16.5, 8.6, 10]) | [17.4, 10.55, 11.29] | [16.5, 8.8, 10] | enclosed by less than one cell |
| small | `set_nelements(200)`, `nsteps` + 1, `step_vector` x 2 | 2.143, 2.571, 2.571 | 2.000 | kept |
| small | origin - 1 | max 3.286 | max 2.333 (1 cell = 0.667) | enclosed by less than one cell |
| small | `maximum` + 3 (request 5.333) | 6.714 | 5.667 (over by 0.333, 1 cell = 0.667) | enclosed by less than one cell |

All new values obey the rule. The old values did not keep the extent (they drifted), and that is the S5 defect.

Causes.

- C1. `unit/step_vector` and `to_dict.step_vector` (and the same keys for `aniso`, `small`). Old value: the step after a later setter (for example 0.7143, not the constructor value 1.0). New value: 1.0, which is the constructor step. The old `nsteps` setter did `self._step_vector /= change_factor`, an in-place change. The golden test stores the array that the property returns, so a later setter changed the stored "start" value. The new setters give a new array. The new value is the correct start value.
- C2. `after_set_nelements`, `after_nsteps`, `after_step_vector`, `after_origin`, `after_maximum` `.step_vector`: the setters keep the extent, so the new step is the extent divided by the new cell count.
- C3. `after_maximum.maximum` and `after_maximum.nsteps`: the grid before this stage has a different step and `nsteps` (C2), so the encloses-the-request result differs. The result is enclosed by less than one cell.

Step 3. The 78 changed keys (old value, new value, cause).

| Key | Old value | New value | Cause |
|---|---|---|---|
| `P2TetMesh/aniso/after_maximum.maximum` | [17.4, 10.5455, 11.2857] | [16.5, 8.8, 10] | C3 (maximum request is enclosed by less than one cell) |
| `P2TetMesh/aniso/after_maximum.nsteps` | [7, 9, 7] | [6, 9, 7] | C3 (maximum request is enclosed by less than one cell) |
| `P2TetMesh/aniso/after_maximum.step_vector` | [1.4, 1.81818, 1.71429] | [1.5, 1.6, 1.5] | C2 (resize keeps the extent: new step = extent / cells) |
| `P2TetMesh/aniso/after_nsteps.step_vector` | [0.7, 0.909091, 0.857143] | [0.75, 0.8, 0.75] | C2 (resize keeps the extent: new step = extent / cells) |
| `P2TetMesh/aniso/after_origin.step_vector` | [1.4, 1.81818, 1.71429] | [1.5, 1.6, 1.5] | C2 (resize keeps the extent: new step = extent / cells) |
| `P2TetMesh/aniso/after_set_nelements.step_vector` | [0.7, 0.909091, 0.857143] | [1, 0.888889, 0.9] | C2 (resize keeps the extent: new step = extent / cells) |
| `P2TetMesh/aniso/after_step_vector.step_vector` | [1.4, 1.81818, 1.71429] | [1.5, 1.6, 1.5] | C2 (resize keeps the extent: new step = extent / cells) |
| `P2TetMesh/aniso/step_vector` | [0.7, 0.909091, 0.857143] | [0.5, 2, 1.5] | C1 (old array was changed in place by a later setter) |
| `P2TetMesh/aniso/to_dict.step_vector` | [0.7, 0.909091, 0.857143] | [0.5, 2, 1.5] | C1 (old array was changed in place by a later setter) |
| `P2TetMesh/small/after_maximum.maximum` | [6.71429, 6.71429, 6.71429] | [5.66667, 5.66667, 5.66667] | C3 (maximum request is enclosed by less than one cell) |
| `P2TetMesh/small/after_maximum.nsteps` | [10, 10, 10] | [11, 11, 11] | C3 (maximum request is enclosed by less than one cell) |
| `P2TetMesh/small/after_maximum.step_vector` | [0.857143, 0.857143, 0.857143] | [0.666667, 0.666667, 0.666667] | C2 (resize keeps the extent: new step = extent / cells) |
| `P2TetMesh/small/after_nsteps.step_vector` | [0.428571, 0.428571, 0.428571] | [0.333333, 0.333333, 0.333333] | C2 (resize keeps the extent: new step = extent / cells) |
| `P2TetMesh/small/after_origin.step_vector` | [0.857143, 0.857143, 0.857143] | [0.666667, 0.666667, 0.666667] | C2 (resize keeps the extent: new step = extent / cells) |
| `P2TetMesh/small/after_set_nelements.step_vector` | [0.428571, 0.428571, 0.428571] | [0.4, 0.4, 0.4] | C2 (resize keeps the extent: new step = extent / cells) |
| `P2TetMesh/small/after_step_vector.step_vector` | [0.857143, 0.857143, 0.857143] | [0.666667, 0.666667, 0.666667] | C2 (resize keeps the extent: new step = extent / cells) |
| `P2TetMesh/small/step_vector` | [0.428571, 0.428571, 0.428571] | [1, 1, 1] | C1 (old array was changed in place by a later setter) |
| `P2TetMesh/small/to_dict.step_vector` | [0.428571, 0.428571, 0.428571] | [1, 1, 1] | C1 (old array was changed in place by a later setter) |
| `P2TetMesh/unit/after_maximum.maximum` | [9, 9, 9] | [8.33333, 8.33333, 8.33333] | C3 (maximum request is enclosed by less than one cell) |
| `P2TetMesh/unit/after_maximum.step_vector` | [1.42857, 1.42857, 1.42857] | [1.33333, 1.33333, 1.33333] | C2 (resize keeps the extent: new step = extent / cells) |
| `P2TetMesh/unit/after_nsteps.step_vector` | [0.714286, 0.714286, 0.714286] | [0.666667, 0.666667, 0.666667] | C2 (resize keeps the extent: new step = extent / cells) |
| `P2TetMesh/unit/after_origin.step_vector` | [1.42857, 1.42857, 1.42857] | [1.33333, 1.33333, 1.33333] | C2 (resize keeps the extent: new step = extent / cells) |
| `P2TetMesh/unit/after_set_nelements.step_vector` | [0.714286, 0.714286, 0.714286] | [0.8, 0.8, 0.8] | C2 (resize keeps the extent: new step = extent / cells) |
| `P2TetMesh/unit/after_step_vector.step_vector` | [1.42857, 1.42857, 1.42857] | [1.33333, 1.33333, 1.33333] | C2 (resize keeps the extent: new step = extent / cells) |
| `P2TetMesh/unit/step_vector` | [0.714286, 0.714286, 0.714286] | [1, 1, 1] | C1 (old array was changed in place by a later setter) |
| `P2TetMesh/unit/to_dict.step_vector` | [0.714286, 0.714286, 0.714286] | [1, 1, 1] | C1 (old array was changed in place by a later setter) |
| `StructuredGrid/aniso/after_maximum.maximum` | [17.4, 10.5455, 11.2857] | [16.5, 8.8, 10] | C3 (maximum request is enclosed by less than one cell) |
| `StructuredGrid/aniso/after_maximum.nsteps` | [7, 9, 7] | [6, 9, 7] | C3 (maximum request is enclosed by less than one cell) |
| `StructuredGrid/aniso/after_maximum.step_vector` | [1.4, 1.81818, 1.71429] | [1.5, 1.6, 1.5] | C2 (resize keeps the extent: new step = extent / cells) |
| `StructuredGrid/aniso/after_nsteps.step_vector` | [0.7, 0.909091, 0.857143] | [0.75, 0.8, 0.75] | C2 (resize keeps the extent: new step = extent / cells) |
| `StructuredGrid/aniso/after_origin.step_vector` | [1.4, 1.81818, 1.71429] | [1.5, 1.6, 1.5] | C2 (resize keeps the extent: new step = extent / cells) |
| `StructuredGrid/aniso/after_set_nelements.step_vector` | [0.7, 0.909091, 0.857143] | [1, 0.888889, 0.9] | C2 (resize keeps the extent: new step = extent / cells) |
| `StructuredGrid/aniso/after_step_vector.step_vector` | [1.4, 1.81818, 1.71429] | [1.5, 1.6, 1.5] | C2 (resize keeps the extent: new step = extent / cells) |
| `StructuredGrid/aniso/step_vector` | [0.7, 0.909091, 0.857143] | [0.5, 2, 1.5] | C1 (old array was changed in place by a later setter) |
| `StructuredGrid/aniso/to_dict.step_vector` | [0.7, 0.909091, 0.857143] | [0.5, 2, 1.5] | C1 (old array was changed in place by a later setter) |
| `StructuredGrid/small/after_maximum.maximum` | [6.71429, 6.71429, 6.71429] | [5.66667, 5.66667, 5.66667] | C3 (maximum request is enclosed by less than one cell) |
| `StructuredGrid/small/after_maximum.nsteps` | [10, 10, 10] | [11, 11, 11] | C3 (maximum request is enclosed by less than one cell) |
| `StructuredGrid/small/after_maximum.step_vector` | [0.857143, 0.857143, 0.857143] | [0.666667, 0.666667, 0.666667] | C2 (resize keeps the extent: new step = extent / cells) |
| `StructuredGrid/small/after_nsteps.step_vector` | [0.428571, 0.428571, 0.428571] | [0.333333, 0.333333, 0.333333] | C2 (resize keeps the extent: new step = extent / cells) |
| `StructuredGrid/small/after_origin.step_vector` | [0.857143, 0.857143, 0.857143] | [0.666667, 0.666667, 0.666667] | C2 (resize keeps the extent: new step = extent / cells) |
| `StructuredGrid/small/after_set_nelements.step_vector` | [0.428571, 0.428571, 0.428571] | [0.4, 0.4, 0.4] | C2 (resize keeps the extent: new step = extent / cells) |
| `StructuredGrid/small/after_step_vector.step_vector` | [0.857143, 0.857143, 0.857143] | [0.666667, 0.666667, 0.666667] | C2 (resize keeps the extent: new step = extent / cells) |
| `StructuredGrid/small/step_vector` | [0.428571, 0.428571, 0.428571] | [1, 1, 1] | C1 (old array was changed in place by a later setter) |
| `StructuredGrid/small/to_dict.step_vector` | [0.428571, 0.428571, 0.428571] | [1, 1, 1] | C1 (old array was changed in place by a later setter) |
| `StructuredGrid/unit/after_maximum.maximum` | [9, 9, 9] | [8.33333, 8.33333, 8.33333] | C3 (maximum request is enclosed by less than one cell) |
| `StructuredGrid/unit/after_maximum.step_vector` | [1.42857, 1.42857, 1.42857] | [1.33333, 1.33333, 1.33333] | C2 (resize keeps the extent: new step = extent / cells) |
| `StructuredGrid/unit/after_nsteps.step_vector` | [0.714286, 0.714286, 0.714286] | [0.666667, 0.666667, 0.666667] | C2 (resize keeps the extent: new step = extent / cells) |
| `StructuredGrid/unit/after_origin.step_vector` | [1.42857, 1.42857, 1.42857] | [1.33333, 1.33333, 1.33333] | C2 (resize keeps the extent: new step = extent / cells) |
| `StructuredGrid/unit/after_set_nelements.step_vector` | [0.714286, 0.714286, 0.714286] | [0.8, 0.8, 0.8] | C2 (resize keeps the extent: new step = extent / cells) |
| `StructuredGrid/unit/after_step_vector.step_vector` | [1.42857, 1.42857, 1.42857] | [1.33333, 1.33333, 1.33333] | C2 (resize keeps the extent: new step = extent / cells) |
| `StructuredGrid/unit/step_vector` | [0.714286, 0.714286, 0.714286] | [1, 1, 1] | C1 (old array was changed in place by a later setter) |
| `StructuredGrid/unit/to_dict.step_vector` | [0.714286, 0.714286, 0.714286] | [1, 1, 1] | C1 (old array was changed in place by a later setter) |
| `TetMesh/aniso/after_maximum.maximum` | [17.4, 10.5455, 11.2857] | [16.5, 8.8, 10] | C3 (maximum request is enclosed by less than one cell) |
| `TetMesh/aniso/after_maximum.nsteps` | [7, 9, 7] | [6, 9, 7] | C3 (maximum request is enclosed by less than one cell) |
| `TetMesh/aniso/after_maximum.step_vector` | [1.4, 1.81818, 1.71429] | [1.5, 1.6, 1.5] | C2 (resize keeps the extent: new step = extent / cells) |
| `TetMesh/aniso/after_nsteps.step_vector` | [0.7, 0.909091, 0.857143] | [0.75, 0.8, 0.75] | C2 (resize keeps the extent: new step = extent / cells) |
| `TetMesh/aniso/after_origin.step_vector` | [1.4, 1.81818, 1.71429] | [1.5, 1.6, 1.5] | C2 (resize keeps the extent: new step = extent / cells) |
| `TetMesh/aniso/after_set_nelements.step_vector` | [0.7, 0.909091, 0.857143] | [1, 0.888889, 0.9] | C2 (resize keeps the extent: new step = extent / cells) |
| `TetMesh/aniso/after_step_vector.step_vector` | [1.4, 1.81818, 1.71429] | [1.5, 1.6, 1.5] | C2 (resize keeps the extent: new step = extent / cells) |
| `TetMesh/aniso/step_vector` | [0.7, 0.909091, 0.857143] | [0.5, 2, 1.5] | C1 (old array was changed in place by a later setter) |
| `TetMesh/aniso/to_dict.step_vector` | [0.7, 0.909091, 0.857143] | [0.5, 2, 1.5] | C1 (old array was changed in place by a later setter) |
| `TetMesh/small/after_maximum.maximum` | [6.71429, 6.71429, 6.71429] | [5.66667, 5.66667, 5.66667] | C3 (maximum request is enclosed by less than one cell) |
| `TetMesh/small/after_maximum.nsteps` | [10, 10, 10] | [11, 11, 11] | C3 (maximum request is enclosed by less than one cell) |
| `TetMesh/small/after_maximum.step_vector` | [0.857143, 0.857143, 0.857143] | [0.666667, 0.666667, 0.666667] | C2 (resize keeps the extent: new step = extent / cells) |
| `TetMesh/small/after_nsteps.step_vector` | [0.428571, 0.428571, 0.428571] | [0.333333, 0.333333, 0.333333] | C2 (resize keeps the extent: new step = extent / cells) |
| `TetMesh/small/after_origin.step_vector` | [0.857143, 0.857143, 0.857143] | [0.666667, 0.666667, 0.666667] | C2 (resize keeps the extent: new step = extent / cells) |
| `TetMesh/small/after_set_nelements.step_vector` | [0.428571, 0.428571, 0.428571] | [0.4, 0.4, 0.4] | C2 (resize keeps the extent: new step = extent / cells) |
| `TetMesh/small/after_step_vector.step_vector` | [0.857143, 0.857143, 0.857143] | [0.666667, 0.666667, 0.666667] | C2 (resize keeps the extent: new step = extent / cells) |
| `TetMesh/small/step_vector` | [0.428571, 0.428571, 0.428571] | [1, 1, 1] | C1 (old array was changed in place by a later setter) |
| `TetMesh/small/to_dict.step_vector` | [0.428571, 0.428571, 0.428571] | [1, 1, 1] | C1 (old array was changed in place by a later setter) |
| `TetMesh/unit/after_maximum.maximum` | [9, 9, 9] | [8.33333, 8.33333, 8.33333] | C3 (maximum request is enclosed by less than one cell) |
| `TetMesh/unit/after_maximum.step_vector` | [1.42857, 1.42857, 1.42857] | [1.33333, 1.33333, 1.33333] | C2 (resize keeps the extent: new step = extent / cells) |
| `TetMesh/unit/after_nsteps.step_vector` | [0.714286, 0.714286, 0.714286] | [0.666667, 0.666667, 0.666667] | C2 (resize keeps the extent: new step = extent / cells) |
| `TetMesh/unit/after_origin.step_vector` | [1.42857, 1.42857, 1.42857] | [1.33333, 1.33333, 1.33333] | C2 (resize keeps the extent: new step = extent / cells) |
| `TetMesh/unit/after_set_nelements.step_vector` | [0.714286, 0.714286, 0.714286] | [0.8, 0.8, 0.8] | C2 (resize keeps the extent: new step = extent / cells) |
| `TetMesh/unit/after_step_vector.step_vector` | [1.42857, 1.42857, 1.42857] | [1.33333, 1.33333, 1.33333] | C2 (resize keeps the extent: new step = extent / cells) |
| `TetMesh/unit/step_vector` | [0.714286, 0.714286, 0.714286] | [1, 1, 1] | C1 (old array was changed in place by a later setter) |
| `TetMesh/unit/to_dict.step_vector` | [0.714286, 0.714286, 0.714286] | [1, 1, 1] | C1 (old array was changed in place by a later setter) |


Step 4. Golden update and re-run. `python packages/loop_common/tests/test_structured_support_reference.py --update` wrote 611 arrays. Then `pytest packages/loop_common/tests/test_structured_support_reference.py packages/loop_common/tests/test_structured_grid_2d_reference.py` passed (9 passed). The file `structured_grid_2d_reference.npz` did not change.

### Numeric survey (S5 resize change)

Environment: the worktree, `loop_common`, `loop_interpolation` and the visualisation sources, and `pyevtk`, on `PYTHONPATH`. Run after the golden update, one pytest process at a time.

| Command | Result | Time |
|---|---|---|
| `pytest tests/unit/modelling` | 1540 passed, 2 skipped | 457.78 s |
| `pytest packages/loop_interpolation/tests` | 615 passed, 47 skipped, 2 xfailed | 70.95 s |
| `pytest packages/loop_common/tests` | 406 passed | 69.60 s |
| `pytest tests/unit --ignore=tests/unit/modelling` | 504 passed, 2 skipped | 677.84 s |

No test failed and no test number changed, so no test file needs a change. The only changed numbers are the 78 golden keys above (file `packages/loop_common/tests/data/structured_support_reference.npz`). Nothing under `packages/loop_interpolation/tests/golden` changed (the 615 tests there pass without a regeneration). No failure is left, and so no root failure from the S4 gate was seen in these four runs.


## S5 resize and cell-centre results

Commands (PowerShell, the worktree `LoopStructural-005-s5` on `PYTHONPATH`):

- `python -m pytest packages/loop_common/tests/test_resize_paths.py packages/loop_common/tests/test_cell_centre_sampling.py -p no:cacheprovider -q --durations=0` -> 100 passed in 0.19 s. The slowest item took 0.01 s (setup).
- `python -m pytest packages/loop_common/tests/test_resize_paths.py -p no:cacheprovider -q -s -k a_b_a` -> 3 passed (shows the measured drift).
- `python -m pytest packages/loop_common/tests/test_cell_centre_sampling.py -p no:cacheprovider -q` -> 60 passed in 0.12 s.
- `python -m ruff check LoopStructural packages` -> All checks passed. This includes the three `loop_common` source files.

Counts, as pytest reports them. Resize file: 40 cells (12 + 3 + 3 + 3 + 1 + 1 + 3 + 1 + 10 + 3). The task text said 36. The count differs because `test_2d_counts_are_derived` has 10 cells (2 grids x 5 checks), not 6. Cell-centre file: 60 cells (5 tests x 12). Total 100, which equals the 24 failed + 76 passed of `## S5 pre-fix failures`.

Comparison with `## S5 pre-fix failures`, test by test. All 24 failed cells now pass.

| Test function | Before | Now |
|---|---|---|
| `test_nsteps_setter_keeps_extent` | 3 / 3 failed | 3 / 3 pass |
| `test_step_vector_setter_encloses_extent` | 2 / 3 failed (`[G1]`, `[G3]`) | 3 / 3 pass, `[G1]` and `[G3]` now pass |
| `test_nsteps_below_two_is_rejected_without_change` | 1 / 1 failed | 1 / 1 pass |
| `test_nsteps_not_whole_is_rejected_without_change` | 1 / 1 failed | 1 / 1 pass |
| `test_set_nelements_keeps_extent` | 3 / 3 failed | 3 / 3 pass |
| `test_set_nelements_small_request_is_clamped` | 1 / 1 failed | 1 / 1 pass |
| `test_2d_counts_are_derived` | 10 / 10 failed | 10 / 10 pass |
| `test_rectilinear_resize_is_rejected` | 3 / 3 failed | 3 / 3 pass |
| `test_reapplying_same_value_10_times_changes_nothing_3d` | 12 / 12 passed | 12 / 12 pass (coverage) |
| `test_step_vector_a_b_a_extent_drift` | 3 / 3 passed | 3 / 3 pass (coverage) |
| Cell-centre file, 5 tests | 60 / 60 passed | 60 / 60 pass |

Measured extent drift of the A-B-A test (`step_vector` set to A, then B = 1.5 A, then A). The bound is A + B.

| Grid | Measured drift | Bound (A + B) |
|---|---|---|
| G1 | [1, 2, 0] | [2.5, 5, 7.5] |
| G2 | [0.5, 0.5, 0.5] | [1.25, 1.25, 1.25] |
| G3 | [10, 100, 7.5] | [25, 250, 18.75] |

All drifts are inside the bound. The pre-fix drift was G1 [1, 0, 3], G2 [0.5, 0.5, 0.5], G3 [0, 0, 7.5], also inside the bound. The pre-fix and post-fix values differ because the setters now keep the extent and enclose it by less than one cell.

Evidence for FR-015 and SC-011. `test_centres_are_cell_midpoints_in_world` is the evidence. Each cell of the exported grid (`structured_grid()`) has one centre, as the mean of its 8 corner nodes, and the centres are all different. `test_local_coordinates_grid_associates_own_centres` is a containment check: each cell centre falls inside the `structured_grid(local_coordinates=True)` grid (no -1 from `cell_index_at`). That grid is the axis-aligned enclosure, so a one-to-one count against it does not hold, and the test does not claim it. `test_centres_round_trip_and_index_own_cells` is a consistency check on the cell grid of the box itself, built from the same linspace as the centres. It is not evidence for FR-015 or SC-011. `test_centres_match_exported_cells` and `test_centres_inside_box` are coverage.

OP12: not reproduced. All 60 cell-centre cells passed before the fix and pass now. `_bounding_box.py` was not changed.

The named test `test_each_cell_has_exactly_one_centre` of the task text was replaced, because its premise was wrong (see `## S5 pre-fix failures`, note on the one-to-one test). `test_centres_are_cell_midpoints_in_world` is a guard test, because OP12 is not reproduced: it passed before the fix and passes after the fix.

## S5 diagnostics pre-fix failures

Method. The test file ran on the unfixed code, in one call: `python -m pytest tests/unit/modelling/test_diagnostics_frames.py -p no:cacheprovider -q -rf --no-header` (PowerShell, the worktree `LoopStructural-005-s5`). No misfit fix commit exists at this point. The result is 147 failed and 0 passed in 142.77 s. No test passes. The file has 147 cells.

Test file (`test_diagnostics_frames.py`): 147 failed, 0 passed.

| Test function | Failed / total | First error |
|---|---|---|
| `test_misfit_leaves_state_unchanged` | 24 / 24 | `AssertionError: data['gradient'] changed`. The old methods leave a different array in `data["gradient"]` after the call, in all 24 cells (3 rotations x 4 anchors x 2 engines) |
| `test_new_misfit_methods_leaves_state_unchanged` | 24 / 24 | `AttributeError`: `GeologicalFeature` has no attribute `evaluate_gradient_cosines` |
| `test_gradient_cosine_matches_analytic_direction_exact_field` | 24 / 24 | `assert np.all(np.abs(1.0 - old) <= 1e-10)` is False. The old method gives shape (2, 4) from one column and local vectors. See the values below |
| `test_gradient_cosine_matches_analytic_direction_exact_field_solved_model` | 24 / 24 | `AttributeError`: no attribute `evaluate_gradient_cosines` |
| `test_value_residuals_are_absolute_difference` | 12 / 12 | `AttributeError`: no attribute `evaluate_value_residuals` |
| `test_value_misfit_is_normalised_by_field_range` | 12 / 12 | `AttributeError`: no attribute `evaluate_value_residuals` |
| `test_value_misfit_with_zero_range_does_not_divide` | 2 / 2 | `[equal]`: `RuntimeWarning: divide by zero encountered in divide`. `[nan]`: `assert not np.any(np.isnan(result))` fails, the result has NaN |
| `test_misfit_with_no_rows_is_empty` | 1 / 1 | `AttributeError`: no attribute `evaluate_gradient_cosines` |
| `test_old_misfit_shapes_are_kept` | 4 / 4 | `[gradient-only]`, `[normal-only]`, `[equal-families]`: `AttributeError`, no attribute `evaluate_gradient_cosines`. `[different-families]`: `ValueError`, inhomogeneous shape (the old method gives two lists of different length) |
| `test_stale_misfit_is_not_rebuilt` | 4 / 4 | `[FDI-interpolator]`, `[PLI-interpolator]`: `AssertionError: the call rebuilt the interpolator` (`assert True == False`). `[FDI-builder]`, `[PLI-builder]`: `AssertionError: the call rebuilt the builder` |
| `test_stale_misfit_warns_once_and_returns_nan` | 4 / 4 | `AssertionError: []` (`assert 0 == 1`). The old methods give no `UserWarning` |
| `test_stale_verification_report_does_not_rebuild` | 4 / 4 | `AttributeError`: no attribute `verification_report` |
| `test_report_blocks_are_separate` | 2 / 2 | `AttributeError`: no attribute `verification_report` |
| `test_unsolved_model_gives_no_residual_blocks` | 2 / 2 | `AttributeError`: no attribute `verification_report` |
| `test_exclusions_are_counted_and_explained` | 2 / 2 | `AttributeError`: no attribute `verification_report` |
| `test_gradient_outside_points_give_no_false_unattributed_model` | 1 / 1 | `AttributeError`: no attribute `verification_report` |
| `test_gradient_outside_points_with_pli_are_measured` | 1 / 1 | `AttributeError`: no attribute `verification_report` |

Group counts: 24 + 24 + 24 + 24 + 12 + 12 + 2 + 1 + 4 + 4 + 4 + 4 + 2 + 2 + 2 + 1 + 1 = 147.

Items of the task text.

| Item | Result |
|---|---|
| State test with the old methods | Reproduced. `test_misfit_leaves_state_unchanged` fails in 24 of 24 cells. `data["gradient"]` is a different array after `evaluate_gradient_misfit` (the method divides the stored rows in place) |
| Exact-field cosine test | Reproduced. The field is exact (the gradient is `R1_UNIT` in the world frame), so each cosine must be 1. The old method gives shape (2, 4) and not (8,). Both rows are equal, because the slice `grad[:, :3:6]` is the one column 0 and not the 3 vector columns. The vectors are local, not world. Values of the first row, anchor 0, FDI: rotation 0: [0.2001, 0.2632, 0.1671, 0.2877] (largest `|1 - cos|` 0.833). Rotation 30: [0.0232, 0.0700, -0.0781, 0.1532] (1.078). Rotation 90: [-0.3002, -0.3158, -0.4457, -0.1918] (1.446). The bound is 1e-10 |
| Stale-state test | Reproduced. The call rebuilds the model: `up_to_date` goes from `False` to `True`, in 4 of 4 cells (`[FDI-interpolator]`, `[PLI-interpolator]`: the interpolator flag. `[FDI-builder]`, `[PLI-builder]`: the builder flag) |
| Tests of the new names | Reproduced. 88 cells fail only because the name does not exist: `evaluate_gradient_cosines` 52 cells (24 + 24 + 1 + 3), `verification_report` 12 cells (4 + 2 + 2 + 2 + 1 + 1), `evaluate_value_residuals` 24 cells (12 + 12) |
| `test_value_residuals_are_absolute_difference` | Not reproduced as a value defect. It fails with `AttributeError` (the test uses `evaluate_value_residuals`, which does not exist now). This is not a reproduction of a value defect. The old value path is covered by `test_value_misfit_*` |

Measured value (a short script in the scratchpad, not saved in the repository). It uses `world_foliation_model` of `tests.fixtures.coordinate_frames`, as `_build` does, with `ANCHORS[1]`, rotation 30 and the FDI engine.

| Measure | Value |
|---|---|
| Largest absolute change of `interpolator.data["gradient"]` after one call of `evaluate_gradient_misfit` | 8.2254 (the array has 4 rows x 7 columns) |
| Largest change in each column (x, y, z, vx, vy, vz, weight) | 2.9211, 7.5854, 8.2254, 0.5485, 0.2712, 0.6789, 0.9139 |

The change is large in the position columns and in the vector columns, and the weight column also changes. The old method divides the stored rows in place, so the stored data of the interpolator are not the input data after the call.

## S5 rollback checks

Date: 2026-11-05. The check compares the S5 tree with the S4 end code.

Set-up:

- The throw-away worktree was at commit 1c029d21 (the docs commit that is the S5 start commit). `git diff --stat 1c029d21 c644c8d7 -- LoopStructural packages` prints nothing. No code differs between 1c029d21 and the S4 end commit c644c8d7, so 1c029d21 holds the S4 end code.
- Each tree had its own `PYTHONPATH`. In each tree, `LoopStructural.__file__` and `loop_interpolation.__file__` pointed into that tree.
- The throw-away worktree was removed with `git worktree remove --force` and `git worktree prune`. `git worktree list` shows no throw-away worktree.

| Check | Result | Note |
|---|---|---|
| (a) `dataclasses.fields` of the diagnostics classes | pass | `ConstraintFamilyDiagnostics` and `RegionCoverageDiagnostics` are equal. `ConstraintDiagnosticsReport` has only 4 fields on S4 (`interpolator_type`, `families`, `region_coverage`, `outside_model_points`). These 4 are equal on S5. S5 adds `physical_residuals` (the fifth field), `assembled_residuals` and `exclusions`, all with defaults. The task text says "first five", but S4 has only four. |
| (b) `to_dict()` of a report with no verification | pass | One `ConstraintFamilyDiagnostics` and one `RegionCoverageDiagnostics`. The two dicts are equal. |
| (c) `from loop_interpolation import *` | pass | The S4 names are a subset of the S5 names. S5 adds `ExclusionRecord`, `ResidualSummary`, `build_exclusion_records` and `direction_cosines`. No S4 name is missing. |
| (d) `GeologicalFeature.evaluate_gradient_misfit` and `evaluate_value_misfit` | pass | The parameters are equal on both trees (`(self)`). Only the return annotation differs: S5 has `-> 'np.ndarray'`, S4 has none. The result shape for one family (gradient data, normal data dropped, fixture `world_foliation_model`) is `(1, 4)` on both trees. The value misfit has shape `(27,)` on both. No warning on either tree. |
| (e) `StructuredGrid` with an unchanged grid | pass | `origin`, `maximum`, `nsteps` and `step_vector` are equal on both trees (origin [0, 1, 2], nsteps [7, 8, 9], step_vector [1, 2, 0.5], maximum [6, 15, 6]). |
| (f) 2D `maximum` assertions | pass | `pytest packages/loop_common/tests/test_2d_discrete_support.py tests/unit/interpolator/test_2d_discrete_support.py -p no:cacheprovider -q` gives 19 passed on S5 and 19 passed on S4. |

## S5 post-fix results

Date: 2026-11-05. The orchestrator ran the gates on commit `6b9e6467`. The `PYTHONPATH` held the worktree `LoopStructural-005-s5`, `packages/loop_common/src`, `packages/loop_interpolation/src`, `packages/loopstructural_visualisation/src` and a scratch folder with `pyevtk`. The root run and the package run were separate calls (PowerShell), because the conftest plugins clash.

| Command | Result | Time |
|---|---|---|
| `ruff check LoopStructural packages` | All checks passed | - |
| `pytest packages/loop_common/tests packages/loop_interpolation/tests` | 1021 passed, 47 skipped, 2 xfailed, 7 warnings | 140.08 s |
| Root `pytest -W default::DeprecationWarning` | 2216 passed, 4 skipped, 12 warnings | 936.67 s (15 min 36 s) |

No gate failed. The root run has no failure, and so the S4 failure (`test_write_vol_gocad_is_broken_for_real_bounding_box`) is closed.

The pytest summary has no library warning with the text "is deprecated" that comes from S5. The log tail (`gate_s5.log`) holds only the last 12 lines of the root run, and it has no "is deprecated" line.

Cells of the new test files.

`tests/unit/modelling/test_diagnostics_frames.py` has 147 cells:

| Test function | Cells |
|---|---|
| `test_misfit_leaves_state_unchanged` | 24 |
| `test_new_misfit_methods_leaves_state_unchanged` | 24 |
| `test_gradient_cosine_matches_analytic_direction_exact_field` | 24 |
| `test_gradient_cosine_matches_analytic_direction_exact_field_solved_model` | 24 |
| `test_value_residuals_are_absolute_difference` | 12 |
| `test_value_misfit_is_normalised_by_field_range` | 12 |
| `test_value_misfit_with_zero_range_does_not_divide` | 2 |
| `test_misfit_with_no_rows_is_empty` | 1 |
| `test_old_misfit_shapes_are_kept` | 4 |
| `test_stale_misfit_is_not_rebuilt` | 4 |
| `test_stale_misfit_warns_once_and_returns_nan` | 4 |
| `test_stale_verification_report_does_not_rebuild` | 4 |
| `test_report_blocks_are_separate` | 2 |
| `test_unsolved_model_gives_no_residual_blocks` | 2 |
| `test_exclusions_are_counted_and_explained` | 2 |
| `test_gradient_outside_points_give_no_false_unattributed_model` | 1 |
| `test_gradient_outside_points_with_pli_are_measured` | 1 |

The sum is 147. All 147 failed before the fix (`## S5 diagnostics pre-fix failures`) and all pass now.

`packages/loop_interpolation/tests/test_residual_exclusion_report.py` has 9 tests (no parametrised cells): `test_direction_cosines_pure_function`, `test_report_without_verification_is_unchanged`, `test_builder_records_keep_their_own_key`, `test_unknown_reason_keeps_its_count`, `test_unattributed_rows_are_reported`, `test_gradient_outside_points_give_no_false_unattributed`, `test_dropped_rows_none_gives_no_unattributed`, `test_surfe_style_dropped_equals_outside`, `test_residual_summary_is_importable`. All 9 pass.

Measured `TOL_SOLVED` = 9e-4. The measured maximum `|1 - cos|` is 4.414e-4 on 24 cells. The margin is a factor of 2.

Resize and cell-centre files: 100 cells, all pass (`## S5 resize and cell-centre results`). The 3D structured golden file changed in 78 of 611 keys (`## S5 golden review (C-004)`). The 2D golden file did not change. The earlier survey of T078 (modelling 1540 passed, interpolation 615, loop_common 406, tests/unit rest 504) is in the golden section.

## S5 backend matrix

| Item | Level | Evidence |
|---|---|---|
| FDI misfit (gradient cosines, value residuals) | mandatory | `test_diagnostics_frames.py` (FDI cells) |
| PLI misfit (gradient cosines, value residuals) | mandatory | `test_diagnostics_frames.py` (PLI cells) |
| Residual and exclusion report, FDI and PLI | mandatory | `test_diagnostics_frames.py`, `test_residual_exclusion_report.py` |
| Resize of the structured grid | mandatory | `test_resize_paths.py` |
| Resize of the 2D grid | mandatory | `test_resize_paths.py` (`test_2d_counts_are_derived`) |
| Resize of the rectilinear grid | mandatory | `test_rectilinear_resize_is_rejected` (resize is rejected with `NotImplementedError`) |
| P2 | not in S5 | S2 and S6 cover it |
| Native engines | not in S5 | S2 and S6 cover them |
| Surfe | not in S5 | S2 and S6 cover it |
| Residual reports for fold and P2 cells | not in S5 | S2 and S6 cover them |

PLI fact. PLI assembles no gradient rows. In the case with 3 outside points, its gradient family reports an `unattributed` record of 8 (`test_gradient_outside_points_with_pli_are_measured`). The FDI case gives no false `unattributed` record (`test_gradient_outside_points_give_no_false_unattributed_model`).

## S5 test run time (OP2)

| Item | Time |
|---|---|
| Root suite (orchestrator, commit `6b9e6467`) | 936.67 s (15 min 36 s). The S4 end run took 805.50 s. |
| Package suites (orchestrator) | 140.08 s. The S4 end run took 150.24 s. |

The root time is 131 s more than the S4 time. The root suite has 150 more passed tests than at the S4 end (2216 and 2066, and the one S4 failure now passes). The new files add 147 cells to the root suite, and the diagnostics cells are slow (142.77 s on the unfixed code). The two runs had different machine load, so the times do not compare exactly.

## S5 open points

- OP12 (cell-centre result): not reproduced. All 60 cell-centre cells passed before and after the fix. The premise of the one-to-one test was corrected. No change to `_bounding_box.py`.
- OP16 (default shape): the default shape is kept. Two families of different length now give the shape (1, N1 + N2).
- OP17 (parity files under `packages/loop_interpolation/tests/golden`): not changed. `git diff --stat 1c029d21 HEAD -- packages/loop_interpolation/tests/golden` prints nothing.

## S6 pre-fix failures

Method. The new S6 tests ran on the tree before T100 and T101 (PowerShell, the worktree `LoopStructural-005-s6-build` on `PYTHONPATH`, root and package runs in separate calls). No S6 library change exists at this point. Three runs were made for this section with `-q -p no:cacheprovider --tb=line`. The other runs were measured by the orchestrator on the same tree and are cited as recorded runs.

Commands and counts.

| Command | Result | Kind of run |
|---|---|---|
| `pytest packages/loop_common/tests/test_euclidean_transformation.py` | 12 failed, 0 passed | Run for this section |
| `pytest tests/unit/utils/test_transformation.py` | 131 failed, 11 passed | Run for this section |
| `pytest tests/unit/modelling/test_diagnostics_backends.py` | 72 failed, 48 passed in 96.18 s | Run for this section |
| `pytest tests/unit/modelling/test_diagnostics_frames.py` | 149 passed in 396 s | Recorded run |
| `pytest tests/unit/modelling/test_fault_frame_validation.py` | 210 passed in 198 s | Recorded run |
| `pytest tests/unit/interpolator/test_backend_world_contract.py` | 51 passed, 0 skipped | Recorded run |
| `pytest tests/unit/interpolator/test_backend_optional_native_surfe.py` | 24 skipped | Recorded run |
| `pytest tests/unit/test_backend_matrix.py` | 33 passed | Recorded run |

Package file (`packages/loop_common/tests/test_euclidean_transformation.py`): 12 failed, 0 passed.

| Test function | Failed / total | First error |
|---|---|---|
| `test_module_imports` | 1 / 1 | `ImportError: cannot import name 'getLogger' from 'loop_common.math'` (`_transformation.py` line 5; the module does not import) |
| `test_export_from_loop_common_math` | 1 / 1 | `AttributeError: module 'loop_common.math' has no attribute 'EuclideanTransformation'` |
| `test_round_trip_without_fit` | 6 / 6 | The same `ImportError` for `getLogger`. The ids are `[2-1]`, `[2-4]`, `[2-50]`, `[3-1]`, `[3-4]`, `[3-50]` |
| `test_two_dimensional_rotation_of_known_angle` | 1 / 1 | The same `ImportError` |
| `test_three_dimensional_rotation_keeps_the_z_direction` | 1 / 1 | The same `ImportError` |
| `test_fit_without_scikit_learn_raises_import_error` | 1 / 1 | The same `ImportError` |
| `test_fit_without_rotation_needs_no_scikit_learn` | 1 / 1 | The same `ImportError` |

Group counts: 1 + 1 + 6 + 1 + 1 + 1 + 1 = 12. The old file in `loop_common` cannot be imported, so every test fails at the import.

Root file (`tests/unit/utils/test_transformation.py`): 131 failed, 11 passed.

| Test function | Failed / total | First error |
|---|---|---|
| `test_round_trip_explicit_parameters` | 96 / 96 | `dimensions=2` with 3 or 5 columns: 48 cells (N 1, 4, 50 x 4 anchors x 2 column counts x 2 angles). `ValueError: operands could not be broadcast together with remapped shapes : (1,3)->(1,newaxis,3) (2,2)->(2,2)` (the 5-column cells give (1,5)). `dimensions=3` with 3 or 5 columns: 48 cells, and every cell fails the z check. `AssertionError: assert 97.805 <= (1e-12 * 1.0)` (the largest error is 105.45) |
| `test_round_trip_after_fit` | 24 / 24 | `ValueError: n_components=2 must be between 0 and min(n_samples, n_features)=1 with svd_solver='full'` (4 anchors x 3 values of N x 2 flags) |
| `test_integer_input_is_not_truncated` | 1 / 1 | `AssertionError: assert 0.9617 <= 9.0e-10` (the integer input is truncated) |
| `test_one_dimensional_input_raises_value_error` | 3 / 3 | `[transform]`: `IndexError: tuple index out of range`. `[inverse_transform]` and `[fit]`: `AssertionError: Regex pattern did not match` |
| `test_entry_points_are_one_class` | 1 / 1 | `ImportError: cannot import name 'getLogger' from 'loop_common.math'` |
| `test_pickle_with_the_old_module_path_loads` | 1 / 1 | `ImportError: cannot import name 'EuclideanTransformation' from 'loop_common.math'` |
| `test_new_pickle_names_the_loop_common_path` | 1 / 1 | The same `ImportError` for `EuclideanTransformation` |
| `test_fit_aligns_a_line_with_negative_slope` | 2 / 2 | `[-30.0]`: `AssertionError: assert 1.0472 < 1e-09`. `[-10.0]`: `assert 0.3491 < 1e-09` (the angle error is twice the slope angle) |
| `test_fit_three_dimensions` | 2 / 2 | `ValueError: shapes (3,) and (2,) not aligned: 3 (dim 0) != 2 (dim 0)` (the fit uses 2 components only) |

Group counts: 96 + 24 + 1 + 3 + 1 + 1 + 1 + 2 + 2 = 131. The defect for a missing scikit-learn is covered in the package file (`test_fit_without_scikit_learn_raises_import_error`), which fails at the import now.

Diagnostics file (`tests/unit/modelling/test_diagnostics_backends.py`): 72 failed, 48 passed in 96.18 s. Each of the three failing functions has 24 cells (4 anchors x 3 rotations x 2 engines).

| Test function | Failed / total | First error |
|---|---|---|
| `test_faulted_value_residuals_use_restored_positions` | 24 / 24 | `assert_allclose` fails (line 144). Example `[a1-r30-FDI]`: all 8 elements differ. The method gives about 0.079 to 0.083, and the independent residual at the restored positions gives 0.0001 to 0.002. The pre-fix faulted value residual is 0.083 to 0.099 |
| `test_faulted_gradient_cosines_use_restored_positions` | 24 / 24 | `assert_allclose` fails (line 159). Example `[a1-r30-FDI]`: the method gives cosine 0.99108 and the independent value is 0.99548 (difference 0.0044) |
| `test_report_physical_value_matches_independent_rms` | 24 / 24 | `AssertionError: assert 0.0727 <= 1.39e-07` (line 187). The report RMS differs from the independent RMS by about 0.07 to 0.08 |

Group counts: 24 + 24 + 24 = 72. Tests 1, 2 and 4 of the task text are these three functions. The restored-rows test is a guard.

Guard tests that pass before the fix (each line gives the file, the test and the count).

- `tests/unit/modelling/test_diagnostics_backends.py::test_fault_restores_the_stored_rows`: 24 passed (the restored-rows test; 4 anchors x 3 rotations x 2 engines).
- `tests/unit/modelling/test_diagnostics_backends.py::test_p2_misfit_leaves_state_unchanged`: 12 passed (P2 test; 4 anchors x 3 rotations).
- `tests/unit/modelling/test_diagnostics_backends.py::test_p2_report_blocks_are_separate`: 12 passed (P2 test; 4 anchors x 3 rotations).
- `tests/unit/modelling/test_diagnostics_frames.py` (all tests): 149 passed (recorded run). The measured largest solved residual `|1 - cos|` is 1.215e-4 (FDI) and 4.414e-4 (PLI).
- `tests/unit/modelling/test_fault_frame_validation.py` (the R2 displacement and residual tests): 210 passed (recorded run).
- `tests/unit/interpolator/test_backend_world_contract.py` (the changed exact-field test is in this file): 51 passed, 0 skipped (recorded run). The largest Surfe world error is 3.492e-11.
- `tests/unit/interpolator/test_backend_optional_native_surfe.py` (the stub cells): 24 skipped (recorded run). 12 cells skip because `loop_interpolation_core` is not installed, and 12 cells skip because `surfepy` is not installed.
- `tests/unit/test_backend_matrix.py`: 33 passed (recorded run).
- `tests/unit/utils/test_transformation.py`: 11 passed. The ids are `test_call_is_equivalent_to_transform`, `test_default_construction`, `test_fit_finds_rotation_that_aligns_main_axis_with_x`, `test_fit_raises_if_points_have_too_few_columns`, `test_fit_rotation_false_keeps_angle_zero`, `test_repr_html_contains_translation_and_angle`, `test_rotation_and_inverse_rotation_are_transposes_in_plane`, `test_too_few_columns_raises_value_error` (2 cells), `test_transform_aligns_variance_with_x_axis`, `test_transform_raises_if_points_have_too_few_columns`.

Measured tolerances (orchestrator runs).

| Constant | Value | Source |
|---|---|---|
| `TOL_SOLVED_FDI` | 2.43e-4 | 2 x the largest solved residual `|1 - cos|` 1.215e-4 |
| `TOL_SOLVED_PLI` | 8.83e-4 | 2 x the largest solved residual `|1 - cos|` 4.414e-4 |
| `TAU_ASSEMBLED` | 9.2e-11 | the largest value is 4.566e-11 |
| `TAU_ORTHOGONAL` | 1.7e-10 | the largest value is 8.408e-11 |
| `TOL_FAULTED_VALUE` | FDI 5.11e-3, PLI 3.8e-2 | T096 measure. The pre-fix faulted value residual is 0.083 to 0.099, so the limit shows the defect. See the `plan.md` Answers entry of wave 44 for the reason |
| `TOL_FAULTED_COS` | FDI 1.1e-2, PLI 1.7e-2 | T096 measure. See the same `plan.md` entry |

Other records.

- T092 (Surfe module name): the pip package is `surfepy` (version 0.52). The wheel holds the module `surfepy` and no module `surfe`. Decision: `test_p0_surfe_nans.py` now uses `pytest.importorskip("surfepy")`.
- T093 (BOM): `_structured_grid_3d.py` and `test_fault_frame_analytic.py` lost their BOM. One line changed in each file.
- T090 (deleted test): the test `test_inverse_transform_is_broken_for_normal_point_clouds_bug` was deleted. The new tests in `tests/unit/utils/test_transformation.py` replace it.
