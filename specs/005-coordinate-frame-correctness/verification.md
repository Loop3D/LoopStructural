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
