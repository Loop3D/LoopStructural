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

### Rigid-transform contract (R1, FR-001 to FR-004)
The file cannot import pre-fix, because `RigidRotationError` does not exist.
All 42 contract checks are "not runnable" pre-fix. None is "reproduced" or
"not reproduced" by a test result. The missing class is itself the evidence.

### Reference constraints (R3, FR-011 to FR-019)
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

### Fault frame analytic (R2, FR-005 to FR-010)
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

### Fault buffer bounds (R5, FR-024 to FR-028)
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

The regenerated file also re-records 33 other keys (for example evaluate_value, evaluate_gradient, rotate) that differ by at most 5.7e-14 (round-off). No tolerance was relaxed.


## Backend matrix

S1 tests use Python backends only.

| Requirement | Backends | Cells |
|---|---|---|
| R2 fault frame analytic | Python FDI, Python PLI (P1) | 4 anchors x 3 rotations x 2 centre modes x 2 backends = 48 planar cells, plus 5 single tests |
| R3 reference constraints | FDI, PLI (P1) | reference file, sample rules None / centroid / gauss4, all frame cases (312 cells) |
| R5 fault buffers | FDI, PLI | buffers 0, 0.2, 0.5 x 3 rotations x 2 centre modes (84 cells) |

No native, browser or RBF backend is in S1 scope. FR-030 states that such
backends are not a prerequisite.

## Rollback checks

On a throw-away worktree at `ebc7c699` (imports checked as above):
`python -m pytest tests/unit/modelling/test__fault_builder.py packages/loop_common/tests/test_bounding_box.py -p no:cacheprovider -q`
gives 15 passed in 0.5 s. The two files pass on the pre-spec code, so the
rollback target is sound.

## Test run time (OP2)

| Item | Time |
|---|---|
| All 4 S1 files, one run | 291 s test time (4 min 57 s wall) |
| Slowest file: `test_fault_frame_analytic.py` | 264 s test time (4 min 31 s wall) |
| `test_fault_buffer_bounds.py` | 20 s (26 s wall) |
| `test_reference_constraint_frames.py` | 5.6 s (12 s wall) |
| `test_rigid_transform_contract.py` | 0.1 s (4 s wall) |
