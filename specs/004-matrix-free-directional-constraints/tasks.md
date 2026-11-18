# Tasks: Matrix-free directional constraints

**Input**: `specs/004-matrix-free-directional-constraints/` (`spec.md`, `plan.md`)
**Prerequisites**: `plan.md` (required), `spec.md` (required for the user stories)

<!-- Written in ASD-STE100 Simplified Technical English. -->

## Orchestrator notes

- **Repository note**: The accepted constitution is at
  `.specify/memory/constitution.md`. Apply its merge gates: `pytest` and
  `ruff check LoopStructural packages`. Use the repository-fit check in
  `plan.md`.
- **Worktree**: implement this spec in its own worktree. Do not switch the
  branch of the main checkout.
- **Rules for each sub agent**: read `plan.md`; change only the files in the
  task `Files` list; for a test task, run the test and see it fail; write
  comments and docs in ASD-STE100 Simplified Technical English.
- **Tiers**: `[S]` = small, fast model. `[M]` = standard model. If an `[S]`
  task fails twice, send it to `[M]`.
- **Waves**: 16 waves. Start a wave when all `Needs` are done. Run tasks in
  one wave in parallel only when their `Files` lists do not overlap.
- **After each Rust wave**: run
  `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml`.
  Fix a failure before the next wave.
- **After each Python wave**: from `packages/loop_interpolation`, run the
  focused pytest command in the task. Fix a failure before the next wave.
- **Checkpoints**:
  - After T009, the core explicit-row tests, regularisation tests and memory
    tests pass.
  - After T013, Python can add and solve one cell family through the native
    adapter.
  - After T017, the folded feature builder selects `FDFoldRust`.
  - After T018 and T019, review golden results and record both quadrature
    measurements. T018 waits for T020 so the captured baseline uses the
    selected default.
  - T020 waits for the AC-003 and AC-006 results before choosing the default
    quadrature. This is OP1 in `plan.md`.
- **Blocked tasks**: do not start T020 until T005, T014 and T019 are done and
  their results are recorded.
- **Before done**: run Rust formatting, clippy and tests:
  `cargo fmt --manifest-path packages/loop_interpolation_core/Cargo.toml -- --check`,
  `cargo clippy --manifest-path packages/loop_interpolation_core/Cargo.toml --all-targets -- -D warnings`,
  and `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --features parallel`.
  Run the golden tests for `fold_sparse` and `refolded_sparse`, and run
  `packages/loop_interpolation/tests/test_native_parity.py`.
- **Golden files**: Keep unaffected files fixed. Source 005 C-004 allows
  reviewed corrections only for affected old results with analytic evidence.
  Keep bounds and tolerances fixed. A failed new native engine does not
  permit old-result regeneration. Add the four `FDFoldRust` golden files
  after the selected quadrature default is fixed.
- **Commits**: Conventional Commits, `type(scope): description`. Do not
  commit unless the user asks.

## Phase 1: Core test setup

- [ ] T001 [Foundation] [M] Add public API stubs for a directional cell family.
  - Do: Add public `CellQuadrature` and internal `DirectionalCellFamily`
    types in `interpolator.rs`. Add the
    `Problem.add_directional_cell_constraint` signature with `vectors`,
    optional `targets`, `weight`, `quadrature` and optional `mask`. Add the
    `Problem.add_hessian_vector_regularisation` signature with `directions`
    and `weight`, plus the public exports. The stubs must fail explicitly;
    they must not return a valid result.
  - Files: `packages/loop_interpolation_core/src/interpolator.rs`, `packages/loop_interpolation_core/src/problem.rs`, `packages/loop_interpolation_core/src/lib.rs`
  - Needs: none
  - Wave: 1
  - Done when: `cargo check --manifest-path packages/loop_interpolation_core/Cargo.toml` exits 0.

## Phase 2: Core behavior tests

- [ ] T002 [P] [US3] [M] Test explicit cell rows and Hessian-vector regularisation.
  - Do: Add unit tests in `interpolator.rs` that compare a deterministic
    centre-rule family with explicit rows for `A^T A x`, the right-hand side
    and the diagonal (AC-001). Check constant-vector centre rows against
    tangent constraints (AC-002), and test both quadrature rules. Also compare
    the three `H v` regularisation rows with explicit rows and confirm the
    existing `(v^T H v)^2` API remains available. Name the tests with the
    `cell_directional_` and `hessian_vector_` prefixes.
  - Files: `packages/loop_interpolation_core/src/interpolator.rs`
  - Needs: T001
  - Wave: 2
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --lib` fails because the family and new regularisation are stubs.
- [ ] T003 [P] [US3] [S] Test family input validation and zero vectors.
  - Do: Test wrong vector, target and mask lengths; NaN and infinite input;
    zero node vectors; a zero interpolated vector; and invalid Hessian-vector
    regularisation input. Check the named `CoreError` and confirm zero
    directions add no residual or target.
  - Files: `packages/loop_interpolation_core/tests/directional_cell_validation.rs`
  - Needs: T001
  - Wave: 2
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test directional_cell_validation` fails on the stub.
  - Context: The family accepts flat node vectors and optional targets and cell mask. Wrong lengths and nonfinite values return named `CoreError`; zero vectors are valid and add no row.
- [ ] T004 [P] [US2] [M] Test family memory estimates and budgets.
  - Do: Add one homogeneous family, one targeted family and masks to a
    problem. Check that `estimate_memory` includes their retained storage and
    that a budget below the estimate is rejected before solve allocation.
  - Files: `packages/loop_interpolation_core/tests/directional_cell_memory.rs`
  - Needs: T001
  - Wave: 2
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test directional_cell_memory` fails because family storage is not counted.
- [ ] T005 [P] [US1] [M] Test solve paths, precision, threads and the synthetic fold.
  - Do: Test `Problem.solve` and stepped `SolveSession` with a cell family.
    Compare `f32` and `f64`, and one-thread and multi-thread results. Add the
    exact AC-003 field `-z + 0.1 sin(2 pi x)` on a unit cube with 40 x 8 x 8
    cells, all three cell families, the three `|H v|^2` terms, `gauss2`, the
    specified vectors, FR-012 default weights, and unit-weight analytic value
    constraints on both z-boundary planes. Use the cell-centre angle metric.
    Check convergence with the existing CG settings and report the median
    angle.
  - Files: `packages/loop_interpolation_core/tests/directional_cell_solve.rs`
  - Needs: T001
  - Wave: 2
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --features parallel --test directional_cell_solve` fails because the family is not applied.

## Phase 3: Streamed Rust family

- [ ] T006 [US3] [M] Implement matrix-free cell rows and both quadrature rules.
  - Do: Store the family without per-row records. Apply each active cell row
    to the normal operator, right-hand side and diagonal. Reuse the corner
    gather across families. Skip zero-vector rows. Use eight cell parity
    colors for conflict-free parallel scatter. Keep the new unit tests in
    `interpolator.rs`.
  - Files: `packages/loop_interpolation_core/src/interpolator.rs`
  - Needs: T002
  - Wave: 3
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --lib cell_directional` passes the cell-row and quadrature tests.

## Phase 4: Problem API and validation

- [ ] T007 [US3] [M] Implement the `Problem` API for cell families.
  - Do: Validate the exact node vector count, optional target count, cell
    mask count, finite values and positive finite weight. Accept zero vectors.
    Store the family and return named `CoreError` values for invalid input.
  - Files: `packages/loop_interpolation_core/src/problem.rs`
  - Needs: T003, T006
  - Wave: 4
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test directional_cell_validation` passes.

## Phase 5: Memory and regularisation

- [ ] T008 [US2] [M] Count cell-family storage in memory estimates.
  - Do: Add counts for direction vectors, optional targets, masks and family
    metadata to `ProblemCounts`. Include these bytes in `estimate_memory` and
    in the pre-allocation budget check. Do not count temporary rows that are
    recomputed.
  - Files: `packages/loop_interpolation_core/src/memory.rs`, `packages/loop_interpolation_core/src/problem.rs`
  - Needs: T004, T007
  - Wave: 5
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test directional_cell_memory` passes.
- [ ] T009 [US3] [M] Add the `|H v|^2` regularisation API and term.
  - Do: Add three matrix-free residual rows per eligible node from the six
    second-derivative stencils. Add
    `Problem.add_hessian_vector_regularisation(directions, weight)` with
    input validation. Accumulate the weighted normal operator and diagonal.
    Keep `add_directional_regularisation` unchanged.
  - Files: `packages/loop_interpolation_core/src/interpolator.rs`, `packages/loop_interpolation_core/src/problem.rs`
  - Needs: T002, T006, T008
  - Wave: 6
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --lib hessian_vector` passes the reference tests.

## Phase 6: Python binding

- [ ] T010 [US3] [M] Test the NumPy cell-family and regularisation bindings.
  - Do: Test contiguous `N x 3` vectors, optional `N` targets, weight,
    quadrature and a cell mask. Also test the Hessian-vector regularisation
    method. Test a malformed shape and wrong node or mask length. Test that
    omitted targets mean zero.
  - Files: `packages/loop_interpolation/tests/test_directional_cell_binding.py`
  - Needs: T009
  - Wave: 7
  - Done when: `python -m pytest tests/test_directional_cell_binding.py -q` fails because the binding method is missing or incomplete.
- [ ] T011 [US3] [M] Expose cell families and Hessian-vector regularisation in PyO3.
  - Do: Add the binding method to `PyProblem`. Accept a two-dimensional
    contiguous vector array, optional one-dimensional targets, a quadrature
    name and an optional cell mask. Map Rust errors to `ValueError` with the
    constraint name. Expose `add_hessian_vector_regularisation` with the same
    finite-value and length checks as the Rust API.
  - Files: `packages/loop_interpolation_core/src/python.rs`
  - Needs: T010
  - Wave: 8
  - Done when: `python -m pytest tests/test_directional_cell_binding.py -q` passes with the core package built.

## Phase 7: Native adapter

- [ ] T012 [US1] [M] Test native adapter support for fold families and regularisation.
  - Do: Pass precomputed arrays in the computational basis to a fake core
    `Problem`. Check that `_native.py` sends the vectors, targets, weights,
    quadrature, mask and Hessian-vector directions unchanged. Confirm the
    nonzero scalar target remains unchanged and that the submission layer does
    not resample or project. Verify unsupported regions and settings still
    give named errors.
  - Files: `packages/loop_interpolation/tests/test_native_directional_cell.py`
  - Needs: T011
  - Wave: 9
  - Done when: `python -m pytest tests/test_native_directional_cell.py -q` fails because the adapter does not send cell families.
- [ ] T013 [US1] [M] Add cell families and regularisation to the native adapter.
  - Do: Extend `_native.py` to submit the stored family descriptions to the
    core before solving. These descriptions already contain vectors projected
    once by `FDFoldRust` setup into the computational grid basis. Do not sample
    or project geometry again in this submission layer. Send the three fold regularisation directions with
    `add_hessian_vector_regularisation`. Keep core imports lazy and preserve
    the install hint when the optional package is missing.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_native.py`
  - Needs: T012
  - Wave: 10
  - Done when: `python -m pytest tests/test_native_directional_cell.py -q` passes.

## Phase 8: Fold regularisation tests

- [ ] T014 [US1] [M] Test the Python edge-jump comparison.
  - Do: Test translated and rigidly rotated fold inputs. Confirm `FDFoldRust`
    samples geometry and masks at world positions and projects each direction
    once into the computational grid basis before submission. Check that the
    nonzero scalar target is unchanged. Compare each `|H v|^2` direction and the combined result with the
    Python edge-jump result on the AC-003 fold. Evaluate gradients at the same
    AC-003 cell centres. Report the same median angles for
    `(v^T H v)^2`. Do not assert exact equivalence for varying vectors. Also
    test that `step`, `mask_fn` and a fixed seed produce the expected cell
    mask, and that repeating a seed gives the same mask.
  - Files: `packages/loop_interpolation/tests/test_fold_directional_regularisation.py`
  - Needs: T013
  - Wave: 11
  - Done when: `python -m pytest tests/test_fold_directional_regularisation.py -q` fails because the Rust fold adapter is not implemented.

## Phase 9: Fold adapter and golden engine

- [ ] T015 [US1] [M] Implement `FDFoldRust` and add it to the golden runner.
  - Do: Create a native-only finite-difference fold interpolator. Store the
    three cell-family descriptions during setup. Sample fold geometry and
    masks at world positions. Project each fold direction exactly once into
    the computational grid basis here; `_native.py` submits these prepared
    arrays without another conversion. Add the three `|H v|^2` terms and use
    the FR-012 defaults and alignment policy. Convert `step` and `mask_fn` to
    one cell mask. Use the
    supplied seed for repeatable subsampling. Register the engine only for `fold_sparse` and
    `refolded_sparse`. Extend fingerprints with row counts, norms and weights.
    Keep unaffected existing results fixed. Any C-004 correction to an old
    result is owned by its coordinate-fix task, requires analytic evidence
    and explicit review, and does not follow from failure of this engine.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_fd_fold_rust_interpolator.py`, `packages/loop_interpolation/src/loop_interpolation/__init__.py`, `packages/loop_interpolation/tests/golden/engines.py`, `packages/loop_interpolation/tests/golden/conftest.py`, `packages/loop_interpolation/tests/golden/geomodels.py`, `packages/loop_interpolation/tests/golden/test_golden.py`
  - Needs: T014
  - Wave: 12
  - Done when: `python -m pytest tests/test_fold_directional_regularisation.py tests/golden/test_golden.py -q -k "fold_directional_regularisation or FDFoldRust"` passes and unaffected existing results and fingerprints remain fixed. Any reviewed old-result correction belongs to its coordinate-fix task; T018 creates the new baselines.

## Phase 10: Folded feature builder

- [ ] T016 [US1] [S] Test folded feature selection for `FDFoldRust`.
  - Do: Test that a folded foliation requested with
    `interpolatortype="FDFoldRust"` keeps that type and creates the Rust fold
    interpolator. Test the conversion from an existing builder with the Rust
    fold type.
  - Files: `tests/unit/modelling/features/builders/test_folded_feature_builder_rust.py`
  - Needs: T015
  - Wave: 13
  - Done when: `python -m pytest tests/unit/modelling/features/builders/test_folded_feature_builder_rust.py -q` fails because the model factory rejects or replaces the type.
  - Context: `FoldedFeatureBuilder` must preserve `interpolatortype="FDFoldRust"` and create the Rust fold interpolator. Conversion from another builder must retain the source builder's intended transform semantics and must not convert a direction twice.
- [ ] T017 [US1] [M] Register `FDFoldRust` with the model factory.
  - Do: Add the Rust fold type to the package interpolator map and to
    `FOLD_INTERPOLATOR_TYPES`. Preserve it in
    `FoldedFeatureBuilder.from_feature_builder`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_interpolatortype.py`, `packages/loop_interpolation/src/loop_interpolation/__init__.py`, `LoopStructural/modelling/core/_model_feature_factory.py`, `LoopStructural/modelling/features/builders/_folded_feature_builder.py`
  - Needs: T016
  - Wave: 14
  - Done when: `python -m pytest tests/unit/modelling/features/builders/test_folded_feature_builder_rust.py -q` passes.

## Phase 11: Golden baselines and measurements

- [ ] T018 [P] [US1] [M] Add only the two `FDFoldRust` golden baselines.
  - Do: Run the fold golden cases with CG, `maxiter=20000` and `rtol=1e-10`.
    Add the `.npz` and `.json` files for the two cases. Review the values,
    metrics and generic-family fingerprints. If a golden check fails, fix the
    implementation; if a bound cannot be met, report the failure to the owner.
    Keep unaffected files and fingerprints fixed. Source 005 C-004 permits
    reviewed corrections only for affected old results with analytic evidence.
    Keep bounds and tolerances fixed. A failed new engine does not permit old
    result regeneration. Any reviewed old-result correction belongs to its
    coordinate-fix task; this task writes only the four new files listed below.
  - Files: `packages/loop_interpolation/tests/golden/golden_files/fold_sparse__FDFoldRust.npz`, `packages/loop_interpolation/tests/golden/golden_files/fold_sparse__FDFoldRust.json`, `packages/loop_interpolation/tests/golden/golden_files/refolded_sparse__FDFoldRust.npz`, `packages/loop_interpolation/tests/golden/golden_files/refolded_sparse__FDFoldRust.json`
  - Needs: T015, T017, T020
  - Wave: 15
  - Done when: `python -m pytest tests/golden/test_golden.py -q -k "FDFoldRust"` passes and this task changes only the four new baseline files listed above. Any reviewed old-result correction is owned by its coordinate-fix task.
- [ ] T019 [P] [US2] [M] Measure cell-family memory and iteration time.
  - Do: Measure retained storage against the explicit three-row reference on
    a 100 x 100 x 100-cell grid and check that it is at least five times
    smaller. Measure peak solve memory on the AC-003 grid and check that
    `estimate_memory` is within 10 percent. Measure one CG operator iteration
    with three families for both quadrature rules and compare it with the
    existing regularisation-only path; check the 1.5x centre and 3x gauss2
    limits. Record environment, command, bytes and timings. If a target is
    missed, stop and report the failure; do not weaken the limit.
  - Files: `packages/loop_interpolation_core/tests/directional_cell_performance.rs`, `specs/004-matrix-free-directional-constraints/measurements.md`
  - Needs: T008, T009
  - Wave: 7
  - Done when: the performance test or measurement runner reports both quadrature results and `measurements.md` contains the exact measurements.

## Phase 12: Select default and final gates

- [ ] T020 [US3] [M] Select and document the default quadrature.
  - Do: Use the AC-003 and AC-006 results to select `centre` or `gauss2`.
    Set the same default in the Rust and Python APIs. Record why the choice
    meets the accuracy and iteration-time targets.
  - Files: `packages/loop_interpolation_core/src/interpolator.rs`, `packages/loop_interpolation_core/src/problem.rs`, `packages/loop_interpolation_core/src/python.rs`, `packages/loop_interpolation/src/loop_interpolation/_fd_fold_rust_interpolator.py`, `specs/004-matrix-free-directional-constraints/measurements.md`
  - Needs: T005, T014, T015, T019
  - Wave: 13
  - Blocked by: OP1
  - Done when: Rust and Python default calls select the documented rule, and the explicit rule tests still pass.
- [ ] T021 [Finish] [M] Run the final test and regression gates.
  - Do: Run Rust formatting, clippy and all core tests with and without the
    `parallel` feature. Run the package directional-family tests,
    `test_native_parity.py`, the two golden fold cases, and the folded
    feature builder test. Run constitution merge gates `pytest` and
    `ruff check LoopStructural packages`. Confirm unaffected golden files
    and fingerprints remain unchanged; any C-004 correction must be limited
    to analytically affected results.
  - Files: none
  - Needs: T001, T002, T003, T004, T005, T006, T007, T008, T009, T010, T011, T012, T013, T014, T015, T016, T017, T018, T019, T020
  - Wave: 16
  - Done when: every command in `plan.md` and the orchestrator's final checks exits 0.

## Wave table

| Wave | Tasks | Tier [S] | Tier [M] |
|---|---|---|---|
| 1 | T001 | 0 | 1 |
| 2 | T002-T005 | 1 | 3 |
| 3 | T006 | 0 | 1 |
| 4 | T007 | 0 | 1 |
| 5 | T008 | 0 | 1 |
| 6 | T009 | 0 | 1 |
| 7 | T010, T019 | 0 | 2 |
| 8 | T011 | 0 | 1 |
| 9 | T012 | 0 | 1 |
| 10 | T013 | 0 | 1 |
| 11 | T014 | 0 | 1 |
| 12 | T015 | 0 | 1 |
| 13 | T016, T020 | 1 | 1 |
| 14 | T017 | 0 | 1 |
| 15 | T018 | 0 | 1 |
| 16 | T021 | 0 | 1 |

## Critical path

T001 -> T002 -> T006 -> T007 -> T008 -> T009 -> T010 -> T011 -> T012 ->
T013 -> T014 -> T015 -> T016 -> T017 -> T018 -> T021

T019 is the parallel measurement gate for T020. The quadrature default cannot
be selected before the AC-003 test and AC-006 measurements are complete.

## Dependencies

T001 defines the Rust API surface. T002-T005 add independent failing tests.
T006 implements the streamed rows. T007 adds input validation and family
registration. T008 adds memory estimates. T009 adds Hessian-vector
regularisation. Python binding and adapter work follows those core APIs.
Fold tests precede the `FDFoldRust` class. Model registration and golden
baselines follow that class. Measurements and the synthetic-fold results
determine the default quadrature. T020 updates the default in the adapter
created by T015. T018 waits for that default and builder registration before
it captures golden results. The final task runs all gates.
