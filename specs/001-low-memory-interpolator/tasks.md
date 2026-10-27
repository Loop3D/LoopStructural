# Tasks: Low-memory interpolator for the browser

**Input**: `specs/001-low-memory-interpolator/` (`spec.md`, `plan.md`, `research.md`)
**Prerequisites**: `plan.md` (required), `spec.md` (required for the user stories)

<!-- Written in ASD-STE100 Simplified Technical English. -->

## Orchestrator notes

- **Repository note**: this repository has no `.specify/` folder, no constitution, no `specs/README.md` and no Cargo workspace. The repo-fit check in `plan.md` replaces the constitution check. Build and test commands name the package with `--manifest-path packages/loop_interpolation_core/Cargo.toml`. Do not set `CARGO_TARGET_DIR`.
- **Worktree**: implement this spec in its own worktree. Do not switch the branch of the main checkout.
- **Rules for each sub agent**: read `plan.md` (sections D1 to D8 and the open points); change only the files in the task `Files` list; for a test task, run the test and see it fail (or be skipped, where the task says so); write comments and docs in ASD-STE100 Simplified Technical English. Use `loop-rs` (`~/dev/loop-rs/crates/interpolation-fd`, commit `bff51342da8ebff7ca3c87e07c1d4fbd31279c8b`) only as the source named in the task.
- **Tiers**: `[S]` = send to a small, fast model. `[M]` = send to the standard model. If an `[S]` task fails twice, send it to `[M]`. The role names in `.agents/model-policy.json` are `worker-S` and `worker-M`.
- **Waves**: 21 waves. Start a wave when all `Needs` of its tasks are done. Run the tasks of a wave in parallel. A task with `Files: none` runs alone.
- **After each wave**: run `cargo check --manifest-path packages/loop_interpolation_core/Cargo.toml` (once the crate exists) and `python -m pytest packages/loop_interpolation/tests -q -x` after a wave that changed Python files. Fix a failure before the next wave.
- **Checkpoints**:
  - After task T011 (end of Phase 2): the ported unit tests pass with and without `--features parallel`, and `measurements.md` has `Baseline` and `After table cuts`.
  - After task T020: `bytes_per_dof` is at most 75. If it is above 95, stop and ask the user (OP1). Between 75 and 95, ask the user before Phase 7.
  - After task T023 (end of Phase 4): `cargo test` passes with all core tests.
  - After task T038: the Python suite passes in the venv with the core installed, and without it.
  - After task T045: browser tests pass and the size is at most 2 MB gzip. If the size is above, stop (OP2).
- **Blocked tasks**: do not start a task with `Blocked by` before the user confirms the open point.
- **Before done**: `cargo fmt --check`, `cargo clippy -- -D warnings`, `cargo test` (all with `--manifest-path packages/loop_interpolation_core/Cargo.toml`), and `python -m pytest packages/loop_interpolation/tests -q`.
- **Commits**: Conventional Commits, `type(scope): description`, for example `feat(loop-interpolation-core): add solve session`. Scope the commits to the package path so release-please counts them for the right package. Do not commit unless the user asks.

## Phase 1: Setup

Package scaffold, release and workspace settings. Rust is new in this repository.

- [x] T001 [P] [Setup] [S] Add the new package to the release-please files.
  - Do: Add package path `packages/loop_interpolation_core`, `release-type` `python`, `component` `loop-interpolation-core`.
  - Do: Add the same path to `exclude-paths` of the root package `.` (fix of 2026-09-27).
  - Do: Add `"packages/loop_interpolation_core": "0.0.0"` to the manifest.
  - Files: `.release-please-config.json`, `.release-please-manifest.json`
  - Needs: none
  - Wave: 1
  - Done when: `python -c "import json; c=json.load(open('.release-please-config.json')); assert 'packages/loop_interpolation_core' in c['packages']; assert 'packages/loop_interpolation_core' in c['packages']['.']['exclude-paths']; m=json.load(open('.release-please-manifest.json')); assert 'packages/loop_interpolation_core' in m"` exits 0.
- [x] T002 [P] [Setup] [S] Keep the Rust package out of the uv workspace. Ignore build output.
  - Do: In `[tool.uv.workspace]` keep `members = ["packages/*"]` and add `exclude = ["packages/loop_interpolation_core"]`.
  - Do: In `.gitignore` add `packages/loop_interpolation_core/pkg/` and `packages/loop_interpolation_core/tests_web/node_modules/`. `target/` is already ignored.
  - Files: `pyproject.toml`, `.gitignore`
  - Needs: none
  - Wave: 1
  - Done when: `python -c "import tomllib; d=tomllib.load(open('pyproject.toml','rb')); assert 'packages/loop_interpolation_core' in d['tool']['uv']['workspace']['exclude']"` exits 0.
- [x] T003 [P] [Setup] [S] Create `Cargo.toml` and an empty `lib.rs` for the core crate.
  - Do: Package name `loop-interpolation-core`, version `0.0.0`, edition `2024`.
  - Do: `[lib] name = "loop_interpolation_core"`, `crate-type = ["cdylib", "rlib"]`.
  - Do: Optional dependencies: `rayon = "1.10"`, `pyo3 = { version = "0.29", features = ["extension-module", "abi3-py39"] }`, `numpy = "0.29"`, `wasm-bindgen = "0.2"`.
  - Do: Features: `parallel = ["dep:rayon"]`, `python = ["dep:pyo3", "dep:numpy"]`, `wasm = ["dep:wasm-bindgen"]`. `default = []`.
  - Do: No `serde`, no `regex` (SC-010).
  - Do: Add `[profile.wasm-size]` with `inherits = "release"`, `opt-level = "z"`, `lto = true`, `codegen-units = 1`, `panic = "abort"`.
  - Do: `src/lib.rs` has only a crate doc comment.
  - Files: `packages/loop_interpolation_core/Cargo.toml`, `packages/loop_interpolation_core/src/lib.rs`
  - Needs: none
  - Wave: 1
  - Done when: `cargo check --manifest-path packages/loop_interpolation_core/Cargo.toml` exits 0.
- [x] T004 [P] [Setup] [S] Create `pyproject.toml` (maturin) and a README skeleton for the core package.
  - Do: Build backend `maturin>=1,<2`. Project name `loop-interpolation-core`, `requires-python = ">=3.9"`, dependency `numpy`.
  - Do: `[tool.maturin]`: `features = ["python"]`, `module-name = "loop_interpolation_core"`.
  - Do: README (ASD-STE100): purpose, link to `specs/001-low-memory-interpolator/spec.md`, and the line `Source: loop-rs crates/interpolation-fd, commit bff51342da8ebff7ca3c87e07c1d4fbd31279c8b, branch lachlan/drillhole`.
  - Files: `packages/loop_interpolation_core/pyproject.toml`, `packages/loop_interpolation_core/README.md`
  - Needs: none
  - Wave: 1
  - Done when: `python -c "import tomllib; d=tomllib.load(open('packages/loop_interpolation_core/pyproject.toml','rb')); assert d['tool']['maturin']['module-name']=='loop_interpolation_core'"` exits 0.
- [x] T005 [P] [Found] [M] Write the fixture generator and the 20x20x20 fixture file (reference solution from the matrix solver).
  - Do: Read `plan.md` section D8 (fixture format) first.
  - Do: Build a 20x20x20 `FiniteDifferenceInterpolator` problem with value, gradient, normal, tangent and interface data (fixed seed). Solve it with `solver="cg"` and `tol=1e-10`. This field is the reference solution.
  - Do: Write the flat arrays in the D8 order. Write regularisation and interpolation weights in core meaning: the cross-derivative weights `dxy`, `dxz`, `dyz` are divided by 4 (see `StructuredGrid.get_operators`). Read `loop-rs/crates/interpolation-fd/tests_python/test_drop_in_replacement.py` for the mapping of gradient, normal and tangent data.
  - Do: Record the final mapping in the module doc string. T035 copies it.
  - Files: `packages/loop_interpolation_core/fixtures/make_fixture.py`, `packages/loop_interpolation_core/fixtures/fixture_20.bin`
  - Needs: none
  - Wave: 1
  - Done when: `python -I packages/loop_interpolation_core/fixtures/make_fixture.py --check` exits 0. It reads `fixture_20.bin`, re-solves with the Python matrix solver (`solver="cg"`) and the max difference is at most 1e-6 of the field range.

## Phase 2: Foundation (port the loop-rs code, measure memory)

Phase A of the plan. Source: `~/dev/loop-rs/crates/interpolation-fd` at commit bff51342da8ebff7ca3c87e07c1d4fbd31279c8b. All tasks use the ported unit tests as the safety net.

- [x] T006 [P] [Found] [M] Copy `cg.rs`, `grid.rs`, `stencil.rs`, `interpolator.rs` and `tests/edge_cases.rs` from loop-rs. Declare the modules in `lib.rs`.
  - Do: Copy the files unchanged. Do not copy `main.rs`, `python.rs`, `tests_python/`, `Inequalities.md`, or build output (`*.pyd`, `__pycache__`).
  - Do: In `lib.rs` declare `pub mod cg; pub mod grid; pub mod interpolator; pub mod stencil;`, copy the `PAR_MIN_LEN` constant and the `pub use` lines. Use the crate name `loop_interpolation_core` in `tests/edge_cases.rs`.
  - Do: Keep `rayon` in these files for now (feature `parallel`). T008 removes it.
  - Files: `packages/loop_interpolation_core/src/cg.rs`, `packages/loop_interpolation_core/src/grid.rs`, `packages/loop_interpolation_core/src/stencil.rs`, `packages/loop_interpolation_core/src/interpolator.rs`, `packages/loop_interpolation_core/tests/edge_cases.rs`, `packages/loop_interpolation_core/src/lib.rs`
  - Needs: T003
  - Wave: 2
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --features parallel` passes (all ported unit tests and `edge_cases`).
- [x] T007 [P] [Found] [M] Remove grid rotation and the `serde` attribute (FR-020).
  - Do: Delete field `rotation`, `with_rotation_matrix`, `with_rotation_xy_degrees`, `to_local`, `direction_to_local`, the `serde` `cfg_attr` line, and every test or call site that uses them.
  - Do: Do not change any other behaviour.
  - Files: `packages/loop_interpolation_core/src/grid.rs`, `packages/loop_interpolation_core/src/interpolator.rs`, `packages/loop_interpolation_core/tests/edge_cases.rs`
  - Needs: T006
  - Wave: 3
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --features parallel` passes and `grep -rnE "rotation|to_local|direction_to_local|serde" packages/loop_interpolation_core/src packages/loop_interpolation_core/tests` prints nothing.
- [x] T008 [Found] [M] Make `rayon` optional with a one-thread shim (FR-006).
  - Do: Create `src/par.rs`. With feature `parallel` it re-exports `rayon::prelude::*`.
  - Do: Without the feature it defines traits that give the same method names on std iterators and slices: `par_iter`, `par_iter_mut`, `par_chunks_mut`, `into_par_iter`, `with_min_len` (returns `self`), and a `fold` that takes the rayon form `(identity_fn, op)`. These are the only rayon methods used (checked by `grep`).
  - Do: Replace each `use rayon::prelude::*;` with `use crate::par::*;`. Add `pub(crate) mod par;` to `lib.rs`.
  - Do: The result of every function must not change.
  - Files: `packages/loop_interpolation_core/src/par.rs`, `packages/loop_interpolation_core/src/cg.rs`, `packages/loop_interpolation_core/src/stencil.rs`, `packages/loop_interpolation_core/src/interpolator.rs`, `packages/loop_interpolation_core/src/lib.rs`
  - Needs: T007
  - Wave: 4
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml` (no features) passes, `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --features parallel` passes, and `cargo check --manifest-path packages/loop_interpolation_core/Cargo.toml --target wasm32-unknown-unknown` exits 0.
- [x] T009 [P] [Found] [M] Add the counting allocator and the baseline memory test. Record the baseline.
  - Do: `tests/common/mod.rs`: a `#[global_allocator]` wrapper of `std::alloc::System` with atomics `CURRENT` and `PEAK`, and functions `reset_peak()`, `current()`, `peak()`.
  - Do: `tests/memory.rs`: ignored test `baseline_r1` builds problem R1 (100x100x100 grid, 2000 value and 500 normal constraints, default regularisation) with the ported API, calls `solve(5000, 1e-8)`, and prints `bytes_per_dof = peak / 1e6`. Also a normal test `baseline_small` on 30^3 that only prints.
  - Do: Create `measurements.md` with the table heading `Baseline` and the value.
  - Files: `packages/loop_interpolation_core/tests/common/mod.rs`, `packages/loop_interpolation_core/tests/memory.rs`, `specs/001-low-memory-interpolator/measurements.md`
  - Needs: T008
  - Wave: 5
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --release --test memory -- --ignored --nocapture baseline_r1` prints `bytes_per_dof=` and the same value is in `specs/001-low-memory-interpolator/measurements.md` under the heading `Baseline`.
- [x] T010 [P] [Found] [M] Apply the table cuts of D4 inside `interpolator.rs`. Add `table_bytes()` (FR-001, FR-023).
  - Do: Store `shell_row_lookup`, `shell_row_indices` and `shell_output_indices` as `u32` (the top value is the empty marker).
  - Do: Add `pub fn table_bytes(&self) -> usize`: the bytes held by the interpolator's own tables (lookup, indices, directional weights, combined kernel). Add a test in `tests/memory.rs` that asserts `table_bytes() <= 12 * dof` for a 40^3 grid.
  - Do: Add the heading `After table cuts` to `measurements.md` with the new value. Do not change the CG code here (T020 rewrites it).
  - Files: `packages/loop_interpolation_core/src/interpolator.rs`, `packages/loop_interpolation_core/tests/memory.rs`, `specs/001-low-memory-interpolator/measurements.md`
  - Needs: T009
  - Wave: 6
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --release --test memory -- --ignored --nocapture baseline_r1` prints a value at least 4 bytes per dof lower than `Baseline`, and `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml` passes.
- [x] T011 [P] [Found] [M] Add the module files and the type and function stubs of D8.
  - Do: Read `plan.md` section D8. Create every type there with the exact fields and derive `Debug, Clone` (and `PartialEq` where D8 says). Implement `Display` and `std::error::Error` for `CoreError`.
  - Do: Every function in D8 returns `Err(CoreError::Unsupported { setting: "not implemented".into() })` or, if it cannot return `Result`, calls `unimplemented!()`. `Settings::default()` has the real default values of D8.
  - Do: `scalar.rs`, `python.rs`, `wasm.rs` hold a doc comment only. In `lib.rs` add `pub mod` lines for the new modules, `#[cfg(feature = "python")] mod python;`, `#[cfg(feature = "wasm")] mod wasm;` and `pub use` lines for the public types.
  - Files: `packages/loop_interpolation_core/src/error.rs`, `packages/loop_interpolation_core/src/settings.rs`, `packages/loop_interpolation_core/src/report.rs`, `packages/loop_interpolation_core/src/field.rs`, `packages/loop_interpolation_core/src/memory.rs`, `packages/loop_interpolation_core/src/session.rs`, `packages/loop_interpolation_core/src/problem.rs`, `packages/loop_interpolation_core/src/admm.rs`, `packages/loop_interpolation_core/src/scalar.rs`, `packages/loop_interpolation_core/src/python.rs`, `packages/loop_interpolation_core/src/wasm.rs`, `packages/loop_interpolation_core/src/lib.rs`
  - Needs: T008
  - Wave: 5
  - Done when: `cargo check --manifest-path packages/loop_interpolation_core/Cargo.toml` exits 0 and `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml` still passes.

## Phase 3: User story 1 core - flat input, memory estimate, sessions

Phase B of the plan, first part. Tests before code.

- [x] T012 [P] [US1] [M] Test: bad input gives an error that names the constraint (FR-018, story 1 scenario 10, edge cases).
  - Do: One test for each case. Each case builds a `Problem` and calls one `add_*` method. Check the error variant and that `constraint` is like `value[3]` (family name and row).
  - Do: Cases: NaN in a point; infinite value; array of the wrong length; zero-length vector for `add_gradient`, `add_normal`, `add_tangent`; weight 0 and weight -1; `lower > upper` for `add_value_bound` and `add_pair_bound`; grid with 2 nodes on one axis (`GridTooSmall`); grid whose node count is above `MAX_NODES` (`GridTooLarge`, no allocation); `start` with zero constraints (`NoConstraints`).
  - Do: Warnings (not errors): a point outside the grid is dropped and `warnings()` has one text with the count; two value constraints at one point with different values give a warning with the count; a step ratio above 100 gives a warning.
  - Do: A rotation cannot be given: there is no API for it. Add a doc comment that says so.
  - Files: `packages/loop_interpolation_core/tests/input_validation.rs`
  - Needs: T011
  - Wave: 6
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test input_validation` compiles and fails (the stubs return errors).
- [x] T013 [P] [US1] [M] Implement `Problem::new` and the flat `add_*` methods with validation.
  - Do: Use the ported `FiniteDifferenceInterpolator` methods for the maths. Validate first, then call the ported method. Do not add large allocations.
  - Do: Conflict count: compare each new value constraint with earlier ones at the same point with a hash of the point bits. Use a `HashMap` only if the count is below 1 000 000; above this skip the check and add a warning that says it was skipped.
  - Do: `MAX_NODES` is `isize::MAX as u64 / 8`. Use `u128` for the product of the counts.
  - Files: `packages/loop_interpolation_core/src/problem.rs`
  - Needs: T012
  - Wave: 7
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test input_validation` passes.
- [x] T014 [P] [Found] [S] Test: `MemTracker` counts bytes, tracks the peak, and refuses over budget.
  - Do: Add `#[cfg(test)] mod tests` to `memory.rs` (the stub body stays).
  - Do: Tests: `try_zeroed::<f64>(10)` raises `current` by 80 and `peak` to 80; `release(80)` sets `current` to 0 and keeps `peak` 80; with `budget = Some(79)` a request for 80 bytes gives `OverBudget` and nothing is allocated; `try_zeroed::<f64>(usize::MAX / 16)` gives `OutOfMemory`, not an abort.
  - Files: `packages/loop_interpolation_core/src/memory.rs`
  - Needs: T011
  - Wave: 6
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --lib memory::` compiles and fails.
- [x] T015 [P] [Found] [S] Implement `MemTracker` (`try_zeroed`, `release`, `peak`, `current`) with `Vec::try_reserve_exact` (FR-012).
  - Do: `try_zeroed<T: Copy + Default>(&mut self, n: usize) -> Result<Vec<T>, CoreError>`: check the budget with `u64` maths, call `try_reserve_exact`, then fill with `T::default()`.
  - Files: `packages/loop_interpolation_core/src/memory.rs`
  - Needs: T014
  - Wave: 7
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --lib memory::` passes.
- [x] T016 [P] [US4] [M] Test: memory estimate, budget refusal, and no big allocation before the check (story 4, FR-011, FR-012).
  - Do: Use `tests/common/mod.rs` (counting allocator). Add `mod common;`.
  - Do: `estimate_memory` for 100^3 with R1 data gives `total` between 40 MB and 100 MB, and all parts are above 0. It scales in a straight line: the ratio `total(60^3) / total(30^3)` is between 7.5 and 8.5. The call allocates less than 1 MB (use `reset_peak` and `peak`).
  - Do: `Problem::new` on a 200^3 grid allocates less than 1 MB.
  - Do: `start` with `budget_bytes = Some(1_000_000)` on 100^3 returns `OverBudget`; `max_nodes_per_axis` is at least 1 and a cubic grid of that size has an estimate at most the budget; the call allocates less than the budget.
  - Do: Constraint bytes are counted: adding 100 000 pair bounds raises `constraints` by at least 100 000 * 200 bytes.
  - Files: `packages/loop_interpolation_core/tests/memory_estimate.rs`
  - Needs: T013, T009, T015
  - Wave: 8
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test memory_estimate` compiles and fails.
- [x] T017 [P] [US4] [M] Implement the memory estimate, the budget check and `max_grid_for_budget`.
  - Do: `memory.rs`: `pub fn estimate(counts: &ProblemCounts, s: &Settings) -> MemoryEstimate` with arithmetic from the table in `plan.md` D4. Add 10 percent margin. `ProblemCounts` holds `dof`, counts of each constraint family, `uses_directional`.
  - Do: `pub fn max_grid_for_budget(counts: &ProblemCounts, s: &Settings, budget: u64) -> usize`: largest `n` such that an `n x n x n` grid fits. Binary search. No allocation.
  - Do: `Problem::estimate_memory` fills `ProblemCounts` and calls `estimate`. `Problem::start` calls it first and returns `OverBudget` when `budget_bytes` is set and the total is above it.
  - Files: `packages/loop_interpolation_core/src/memory.rs`, `packages/loop_interpolation_core/src/problem.rs`
  - Needs: T016
  - Wave: 9
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test memory_estimate` passes.
- [x] T018 [P] [US4] [M] Build the grid tables lazily and through `MemTracker`, not in the constructor.
  - Do: `FiniteDifferenceInterpolator::new` must not allocate any `dof`-sized vector. Move `build_shell_indices` and the combined kernel to `pub(crate) fn prepare(&mut self, tracker: &mut MemTracker) -> Result<(), CoreError>`. Call `prepare` at the start of `solve_with` in the ported code so the ported tests still pass.
  - Do: Change every `vec![0.0; dof]`-type allocation in `rhs`, `diagonal`, `scratch` and the border buffers to `tracker.try_zeroed`. Keep the old infallible methods as thin wrappers for the ported tests.
  - Files: `packages/loop_interpolation_core/src/interpolator.rs`
  - Needs: T015, T010, T016
  - Wave: 9
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml` passes (all ported tests) and `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test memory_estimate -- problem_new` passes.
- [x] T019 [P] [US1] [M] Test: the CG session equals the one-shot solve, can be cancelled, reports progress, and frees memory (story 3 core, FR-009, FR-010, FR-014).
  - Do: Problem: 30^3, value and normal constraints (fixed data in the test).
  - Do: (1) Repeat `step(7)` until `Done`: the field differs from the ported `solve(5000, 1e-10)` field by at most 1e-8 of the range.
  - Do: (2) Each `Running(progress)` has `stage == Solve`, a growing `iteration` and a falling `residual` (over 5 steps).
  - Do: (3) Drop the session after 3 steps: `common::current()` returns to within 64 KB of the value before `start`.
  - Do: (4) Drop right after `start` (cancel during setup): same memory check.
  - Do: (5) `finish()` after 3 steps returns a field and a report with `converged == false`. `best_so_far()` gives the same.
  - Do: (6) `max_iter = 3` gives `Done` with `converged == false` and the final residual in the report.
  - Do: (7) The report has `peak_bytes > 0` and `iterations` equal to the count of CG iterations.
  - Files: `packages/loop_interpolation_core/tests/session_cg.rs`
  - Needs: T013, T009, T015
  - Wave: 8
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test session_cg` compiles and fails.
- [x] T020 [US1] [M] Implement the CG state machine and `SolveSession` (`start`, `step`, `best_so_far`, `finish`). Keep 4 full vectors plus the inverse diagonal.
  - Do: Read `plan.md` D3, D4 and D8.
  - Do: `cg.rs`: `pub struct CgState { x, r, p, ap, inv_diag }` and `pub fn step(&mut self, op: &mut dyn FnMut(&[f64], &mut [f64]), iters: usize) -> CgStatus`. Do not store `z`: use `r.z = sum(r*r*inv_diag)`. Do not store `rhs`: build `r = rhs - N x0` straight into `r`. Dot products add in `f64`.
  - Do: `session.rs`: `SolveSession` owns the `Problem` (so it needs no lifetime). `Problem::start(self, s)` checks the budget, calls `prepare`, allocates through `MemTracker`, returns the session in stage `Setup`.
  - Do: `step` runs at most `max_iters` iterations and returns `Step::Running` or `Step::Done`. Dropping the session frees all memory.
  - Do: Add the line with the new value to `measurements.md` under `After CG state machine`.
  - Do: If the value is above 95, stop and report to the orchestrator (OP1).
  - Files: `packages/loop_interpolation_core/src/cg.rs`, `packages/loop_interpolation_core/src/session.rs`, `packages/loop_interpolation_core/src/problem.rs`
  - Needs: T019, T017, T018
  - Wave: 10
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test session_cg` passes and `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --release --test memory -- --ignored --nocapture baseline_r1` prints a value of at most 75 for `bytes_per_dof` (target of D4).

## Phase 4: User story 2 - inequality constraints (ADMM)

Phase B, second part. The active-set polish is not in this plan (see plan.md D5).

- [x] T021 [P] [US2] [M] Test: inequality constraints hold, inactive constraints change nothing, infeasible problems report 20 rows (story 2, FR-004, FR-015).
  - Do: Problem: 20^3 with data that makes the free answer break the order.
  - Do: (1) Value bound: violation at most 1e-4 of the range (SC-004). (2) Pair bound `A - B >= d`: same. (5) Lower and upper bound on one row: both hold.
  - Do: (4) A bound that is not active: the field differs from the field without it by at most 1e-8 of the range.
  - Do: (6) Two bounds that cannot hold together: `solve` returns `Infeasible` with at most 20 rows, sorted by falling violation.
  - Do: Every `Running` progress has `max_violation`. The final report has `max_violation`.
  - Do: Check that `Hybrid` or `ActiveSet` methods do not exist in `Settings`.
  - Files: `packages/loop_interpolation_core/tests/session_admm.rs`
  - Needs: T020
  - Wave: 11
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test session_admm` compiles and fails.
- [x] T022 [P] [US2] [M] Implement ADMM as a state machine over the CG state. Add stall detection for infeasible problems.
  - Do: Port `solve_admm` of `interpolator.rs` into `admm.rs`: outer loop `z = clip(Q x + u, lower, upper)`, `u += Q x - z`, inner CG on `N + rho Q^T Q`. The state holds `z` and `u` (one entry per inequality row) and one full-size vector for the right-hand side. No dense matrix.
  - Do: Make the helper methods the state needs `pub(crate)` in `interpolator.rs` (`accumulate_inequality_normal`, `accumulate_inequality_transpose`, `inequality_diagonal`, `inequality_bounds`).
  - Do: `session.rs`: if the problem has inequality rows, the session runs ADMM; else plain CG. Stage `CheckConstraints` is reported at each outer iteration end.
  - Do: Stall test: if the primal residual does not fall by 1 percent in 20 outer iterations and rows still violate by more than the tolerance, return `Infeasible` with the 20 largest violations.
  - Files: `packages/loop_interpolation_core/src/admm.rs`, `packages/loop_interpolation_core/src/session.rs`, `packages/loop_interpolation_core/src/interpolator.rs`
  - Needs: T021
  - Wave: 12
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test session_admm` passes.
- [ ] T023 [P] [US2] [M] Test: memory with many active constraints stays near R1 (story 2 scenario 3, SC-002, SC-011, FR-002).
  - Do: Make active pair bounds by solving once without bounds, then bounding pairs at their free difference.
  - Do: Add the ignored tests `r2` (100^3, 200 pairs) and `r3` (100^3, 20 000 pairs). Print peak bytes. Add the numbers to `measurements.md` only if the user runs them.
  - Files: `packages/loop_interpolation_core/tests/scale_inequality.rs`
  - Needs: T022
  - Wave: 13
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --release --test scale_inequality` passes (50^3, 2 500 active pair bounds, peak at most 1.2 times the peak without them). The ignored tests `r2` and `r3` run with `--ignored` and print the peaks.

## Phase 5: User story 1 - constraint families, evaluation, fixture

These tasks edit files that ADMM also edits, so they follow Phase 4.

- [x] T024 [P] [US1] [M] Test: tangent and interface properties, point evaluation, large origin, thread count (story 1 scenarios 5, 7, 8, 9; FR-007, FR-016, FR-024).
  - Do: (5) Tangent vectors: the dot product of each tangent with the field gradient at its point, divided by the gradient length, is below 1e-3. Interface pairs: the two field values differ by at most 1e-4 of the range.
  - Do: (8) `Field::evaluate` at the nodes equals the field. `evaluate_gradient` matches central differences of `evaluate` within 1e-3 of the largest gradient. A point outside the grid gives `PointOutsideGrid`.
  - Do: (9) Shift the grid origin and all points by 5 000 000 in each axis: the field differs from the unshifted field by at most 1e-6 of the range.
  - Do: (7) With `#[cfg(feature = "parallel")]`: run in a rayon pool with 1 thread and with 4 threads (`ThreadPoolBuilder`). The fields differ by at most 1e-6 of the range.
  - Files: `packages/loop_interpolation_core/tests/constraint_families.rs`
  - Needs: T020
  - Wave: 11
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test constraint_families` compiles and fails.
- [ ] T025 [P] [US1] [M] Implement `Field` and `evaluate`, `evaluate_gradient`. Fix any failure of the family tests.
  - Do: `Field` evaluate uses the ported `evaluate` and `evaluate_gradient` and the grid. Check each point is inside the grid first.
  - Do: For the large-origin case compute `(point - origin) / step` in `f64` before any cast.
  - Files: `packages/loop_interpolation_core/src/field.rs`, `packages/loop_interpolation_core/src/grid.rs`, `packages/loop_interpolation_core/src/interpolator.rs`, `packages/loop_interpolation_core/src/session.rs`
  - Needs: T024, T022
  - Wave: 13
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test constraint_families` passes (with `--features parallel` also).
- [x] T026 [P] [US1] [M] Test: the Rust solver agrees with the matrix-solver fixture (story 1 scenario 4, SC-003).
  - Do: Read `fixtures/fixture_20.bin` (format in D8, little endian, no crates). Build the `Problem`, solve with `tol = 1e-10`, compare with `expected_field`.
  - Files: `packages/loop_interpolation_core/tests/parity_fixture.rs`
  - Needs: T005, T020
  - Wave: 11
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test parity_fixture` passes. If the max difference is above 1e-6 of the range, the task is not done: report the mapping error to the orchestrator.

## Phase 6: User story 6 - Python drop-in

Phase D of the plan. It needs only phases 3 to 5. The Python tests skip when the core is not installed.

- [x] T027 [P] [US6] [M] Test: parity with the matrix solver for each constraint family, each regularisation, and node order (story 6 scenario 1, FR-025, SC-003).
  - Do: Start the file with `pytest.importorskip("loop_interpolation_core")`.
  - Do: On a 30^3 `StructuredGrid` problem, for value, gradient, normal, tangent, interface, `dxx`, `dyy`, `dzz`, `dxy`, `dxz`, `dyz`, border, and directional regularisation, one family at a time: solve with `solver="cg"` and with `solver="native"`. The fields differ by at most 1e-6 of the range.
  - Do: A node-order test: a field `f(x,y,z) = x + 2y + 3z` fixed with value constraints gives the same `interpolator.c` order as the matrix solver.
  - Do: Use the existing `InterpolatorBuilder` (`use_solver("native")`, then `setup_interpolator`, `solve`). Use `tests/fixtures/` helpers where they exist.
  - Files: `packages/loop_interpolation/tests/test_native_parity.py`
  - Needs: T004
  - Wave: 2
  - Done when: `python -m pytest packages/loop_interpolation/tests/test_native_parity.py -q` collects the tests and skips all (reason: core not installed). After T038 it passes with 0 skipped.
- [x] T028 [P] [US6] [M] Test: inequality problems agree with the Python ADMM solver (story 6 scenario 1, SC-004).
  - Do: Start with `pytest.importorskip`. Compare `solver="native"` with `solver="admm"` on a 20^3 problem with value bounds and with pair bounds (use `inequality_pairs`). The constraint violation is at most 1e-4 of the range in both, and the two fields differ by at most 1e-3 of the range.
  - Files: `packages/loop_interpolation/tests/test_native_inequality.py`
  - Needs: T004
  - Wave: 2
  - Done when: `python -m pytest packages/loop_interpolation/tests/test_native_inequality.py -q` collects the tests and skips all.
- [x] T029 [P] [US6] [S] Test: unsupported settings fail with their name; a missing core gives an install hint (story 6 scenario 3, FR-026).
  - Do: Table test over: region mask, column scaling (`apply_scaling_matrix=True`), `use_regularisation_weight_scale=True`, constant-norm iterations (`constant_norm_iterations`), `interface_pair_mode="all"`, `regularisation_matrix_free=True`. Each with `solver="native"` raises `NotImplementedError` and the message contains the setting name.
  - Do: A test that monkeypatches the import of `loop_interpolation_core` to fail: choosing `"native"` raises `ImportError` with the text `pip install loop-interpolation[native]`. This test does not need the core.
  - Files: `packages/loop_interpolation/tests/test_native_errors.py`
  - Needs: T004
  - Wave: 2
  - Done when: `python -m pytest packages/loop_interpolation/tests/test_native_errors.py -q` collects the tests.
- [x] T030 [P] [US6] [M] Test: with `native` set before setup no sparse rows are built; set after setup it still works with a warning (FR-019, plan D7.2).
  - Do: (a) Create the interpolator with `solver="native"`, set constraints, call `setup()`: `interpolator.constraints` is empty (no rows) and the diagnostics report exists.
  - Do: (b) Call `setup()` with `solver="cg"`, then `solve_system("native")`: rows exist, the field is correct, and one warning about memory is logged (`caplog`).
  - Do: (c) `update()` after (a) does not raise and the field is the same.
  - Do: (d) With no `loop_interpolation_core` installed, existing `solver="cg"`, `"lsmr"` and `"admm"` calls give the same results as before (import of `loop_interpolation` never imports the core).
  - Files: `packages/loop_interpolation/tests/test_native_no_rows.py`
  - Needs: T004
  - Wave: 2
  - Done when: `python -m pytest packages/loop_interpolation/tests/test_native_no_rows.py -q` collects the tests.
- [x] T031 [P] [US6] [S] Test: `resolve_solver_choice("native")` returns `"native"`.
  - Do: Add one test: `resolve_solver_choice("native", logger)` returns `"native"` and logs no warning.
  - Files: `packages/loop_interpolation/tests/test_solver_strategy.py`
  - Needs: none
  - Wave: 1
  - Done when: `python -m pytest packages/loop_interpolation/tests/test_solver_strategy.py -q -k native` fails.
- [x] T032 [P] [US6] [S] Accept `"native"` in `resolve_solver_choice`.
  - Do: Add `"native"` to the list `["cg", "lsmr", "admm", None]`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_solver_strategy.py`
  - Needs: T031
  - Wave: 2
  - Done when: `python -m pytest packages/loop_interpolation/tests/test_solver_strategy.py -q` passes.
- [x] T033 [P] [US6] [S] Add the optional extra `native` to the `loop-interpolation` package.
  - Do: Add `native = ["loop-interpolation-core"]` under `[project.optional-dependencies]`. Do not change `dependencies`.
  - Files: `packages/loop_interpolation/pyproject.toml`
  - Needs: none
  - Wave: 1
  - Done when: `python -c "import tomllib; d=tomllib.load(open('packages/loop_interpolation/pyproject.toml','rb')); assert d['project']['optional-dependencies']['native']==['loop-interpolation-core']"` exits 0.
- [ ] T034 [P] [US6] [M] Write the PyO3 bindings (feature `python`).
  - Do: Read `plan.md` D8 (Python API). Class `Problem` with the methods of D8. Arrays come in as `PyReadonlyArray1<f64>` and are checked for contiguity. Settings come as a `dict` with the keys of D8; unknown keys raise `ValueError`.
  - Do: `solve` returns `(numpy field, dict report)`. In `F32` mode the field is `float32`.
  - Do: `CoreError` becomes `ValueError` (input), `MemoryError` (OutOfMemory, OverBudget), `RuntimeError` (Infeasible, with `worst` in the message), `NotImplementedError` (Unsupported). Keep the text of the error.
  - Do: Take the loop-rs `python.rs` only as a guide for PyO3 0.29 syntax.
  - Files: `packages/loop_interpolation_core/src/python.rs`
  - Needs: T025, T022
  - Wave: 14
  - Done when: `cargo check --manifest-path packages/loop_interpolation_core/Cargo.toml --features python` exits 0.
- [ ] T035 [P] [US6] [M] Write the adapter `_native.py`: map the stored constraints to core calls and solve.
  - Do: Function `solve(interpolator, tol, solver_kwargs) -> bool` and function `check_supported(interpolator, setup_kwargs, solver_kwargs)`.
  - Do: Read the data with `get_value_constraints()`, `get_gradient_constraints()`, `get_norm_constraints()`, `get_tangent_constraints()`, `get_interface_constraints()`, `get_inequality_value_constraints()`, `get_inequality_pairs_constraints()`. Do not read sparse rows.
  - Do: Use the weight and gradient mapping that `fixtures/make_fixture.py` documents (cross-derivative weights divided by 4).
  - Do: Interface constraints: the core takes one weight. If the weights differ, raise `NotImplementedError` that names `interface per-point weights`.
  - Do: `solver_kwargs`: `tol`, `maxiter`, `precision`, `budget_bytes`, `rho`, `inner_max_iter`, `inner_tol`, `admm_abs_tol`, `admm_rel_tol`. Any other key raises `NotImplementedError` with the key name. Do not accept `inequality_method`.
  - Do: Return the field as `interpolator.c` (length `n_nodes`), store the report as `interpolator.latest_native_report`, and set `latest_solve_timing` keys `total` (measured with `perf_counter`).
  - Do: ImportError text: `pip install loop-interpolation[native]`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_native.py`
  - Needs: T027, T028, T029, T034, T005
  - Wave: 15
  - Done when: `python -c "import ast,sys; ast.parse(open('packages/loop_interpolation/src/loop_interpolation/_native.py').read())"` exits 0 and the module imports with no core installed (the core import is inside the function `solve`).
- [ ] T036 [P] [US6] [M] Add the `native` branch to `solve_system` and the late-solver warning.
  - Do: In `solve_system`, after `_normalise_solver_choice`: if `solver_choice == "native"`, call `_native.check_supported`, then `_native.solve`, set `self.up_to_date`, set `latest_solve_timing` through `_solver_pipeline.finalize_timing`, and return. Do not call `_assemble_main_system`.
  - Do: If rows were already built (`self.constraints` is not empty), log one warning: `solver="native" saves memory only when it is set before setup`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py`
  - Needs: T035, T032, T030
  - Wave: 16
  - Done when: `python -m pytest packages/loop_interpolation/tests/test_native_no_rows.py packages/loop_interpolation/tests/test_solver_strategy.py -q` shows no failure (skips allowed until T038).
- [x] T037 [P] [US6] [M] Skip row assembly in `setup_interpolator` when `self.solver == "native"`.
  - Do: After the resolve of weights and the regularisation config, if `self.solver == "native"`: store the weights and the setup kwargs (keep `self._setup_kwargs`), skip the `add_*` calls, `assemble_inner`, `add_directional_regularisation` and `assemble_borders`, and return `finalize_setup_diagnostics_report()`.
  - Do: Check that `finalize_setup_diagnostics_report` does not need the rows. If it does, change only the report to use the stored data when the solver is native.
  - Do: The branch must be safe to run twice (`update()`).
  - Files: `packages/loop_interpolation/src/loop_interpolation/_finite_difference_interpolator.py`
  - Needs: T030, T032
  - Wave: 3
  - Done when: `python -m pytest packages/loop_interpolation/tests/test_native_no_rows.py -q` shows no failure (skips allowed).
- [ ] T038 [US6] [M] Build the wheel into a clean venv and run all Python tests (story 6, FR-017, FR-019).
  - Do: Install `maturin` with `uv pip install maturin` in a venv in the scratchpad. Do not put the venv in the repository.
  - Do: Fix failures only in the four files of the Files list. If a fix needs another file, stop and report it.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_native.py`, `packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py`, `packages/loop_interpolation/src/loop_interpolation/_finite_difference_interpolator.py`, `packages/loop_interpolation_core/src/python.rs`
  - Needs: T036, T037, T033, T030, T027, T028, T029
  - Wave: 17
  - Done when: In a venv with `maturin develop --release -m packages/loop_interpolation_core/Cargo.toml` done: `python -m pytest packages/loop_interpolation/tests -q` passes with 0 failures and the native tests have 0 skips; with the core NOT installed `python -m pytest packages/loop_interpolation/tests -q` has no failure.

## Phase 7: User stories 1, 3, 4 - browser build

Phase C of the plan. Do not start before OP1 is answered (memory per dof from T020).

- [ ] T039 [P] [US1] [M] Test: the Worker solves the 20^3 fixture in Chromium (story 1 scenarios 1 to 4, FR-017).
  - Do: Playwright (Chromium only). `harness.html` starts `web/worker.js` as a module Worker and exposes `window.runFixture(bytes)`. The spec fetches `fixtures/fixture_20.bin`, parses it in JavaScript (D8 format), builds the problem through worker messages (D6), and compares the field with `expected_field`: max difference at most 1e-6 of the range.
  - Do: The page is served with no cross-origin isolation headers (FR-006).
  - Files: `packages/loop_interpolation_core/tests_web/package.json`, `packages/loop_interpolation_core/tests_web/playwright.config.js`, `packages/loop_interpolation_core/tests_web/harness.html`, `packages/loop_interpolation_core/tests_web/fixture.spec.js`
  - Needs: T005, T020
  - Wave: 11
  - Done when: `npm --prefix packages/loop_interpolation_core/tests_web test` fails because the package `pkg/` does not exist yet.
  - Blocked by: OP1
- [ ] T040 [P] [US3] [M] Test: page responsiveness, cancel time, budget error, second solve refused, memory reading (story 3, story 4, SC-006, FR-010).
  - Do: Problem 60^3 generated in the page. (a) A main-thread timer every 10 ms: the largest gap during the solve is below 100 ms. (b) `cancel` after 1 s: the worker answers `cancelled` within 1 s, and a new solve then works. (c) A budget of 1 MB gives an `OverBudget` message with `estimate`, `budget`, `maxNodesPerAxis`. (d) A second `solve` message while one runs gives `SolveRunning` and the first solve ends normally. (e) The worker reports `memoryBytes = WebAssembly.Memory.buffer.byteLength` with the result. (f) Progress messages arrive at least once each second.
  - Files: `packages/loop_interpolation_core/tests_web/responsiveness.spec.js`
  - Needs: T039
  - Wave: 12
  - Done when: `npm --prefix packages/loop_interpolation_core/tests_web test -- responsiveness` fails (no package).
  - Blocked by: OP1
- [ ] T041 [P] [US1] [M] Write the `wasm-bindgen` classes (feature `wasm`).
  - Do: Read `plan.md` D6 and D8. Classes `WasmProblem` and `WasmSession` with the JavaScript names of D6 (`Problem`, `Session`) set by `js_name`.
  - Do: `CoreError` becomes a thrown `JsError` whose message is JSON text: `{"kind":"OverBudget","estimate":..,"budget":..,"maxNodesPerAxis":..}`. No `serde`: write the JSON by hand.
  - Do: `Session.field()` copies out once (`Float64Array` or `Float32Array`). `Session.free()` is the wasm-bindgen free: it drops the session and frees memory.
  - Do: Add a `memory_bytes()` function that returns the wasm memory size in bytes.
  - Do: `Problem.start` consumes the problem (the Rust `start` takes `self`). After `step` returns done, `WasmSession` keeps the `Field`. `evaluate` and `evaluateGradient` call `Field::evaluate` and `Field::evaluate_gradient` on it.
  - Files: `packages/loop_interpolation_core/src/wasm.rs`
  - Needs: T025, T022, T040
  - Wave: 14
  - Done when: `cargo check --manifest-path packages/loop_interpolation_core/Cargo.toml --target wasm32-unknown-unknown --features wasm` exits 0.
  - Blocked by: OP1
- [ ] T042 [P] [US3] [M] Write the Worker wrapper `web/worker.js` (progress, cancel, one solve at a time).
  - Do: Messages in: `init`, `estimate`, `solve`, `cancel`. Messages out: `progress`, `result`, `error`, `cancelled`.
  - Do: `solve` builds the problem from transferable buffers (no copy), calls `start`, then loops `step(maxIters)` with `await new Promise(r => setTimeout(r, 0))` between steps. Pick `maxIters` so one step takes about 100 ms (start with 5, scale with the measured time of the last step).
  - Do: A `cancel` message sets a flag that the loop checks before each step. On cancel call `session.free()`.
  - Do: While a solve runs a second `solve` message answers `error` with kind `SolveRunning`.
  - Do: Post `progress` at least once each second with `stage`, `iteration`, `residual`, `maxViolation`. Post `memoryBytes` with `result`.
  - Files: `packages/loop_interpolation_core/web/worker.js`
  - Needs: T041
  - Wave: 15
  - Done when: `node --check packages/loop_interpolation_core/web/worker.js` exits 0. The browser tests that use it run in the task T045.
  - Blocked by: OP1
- [ ] T043 [P] [US1] [S] Write the TypeScript declaration file and the web README.
  - Do: Declare the message types of `worker.js` and the classes of D6. README in ASD-STE100: how to load the worker, the message list, and the note that the page must transform coordinates (FR-020).
  - Files: `packages/loop_interpolation_core/web/index.d.ts`, `packages/loop_interpolation_core/web/README.md`
  - Needs: T042
  - Wave: 16
  - Done when: `npx --prefix packages/loop_interpolation_core/tests_web tsc --noEmit packages/loop_interpolation_core/web/index.d.ts` exits 0.
  - Blocked by: OP1
- [ ] T044 [P] [US1] [M] Write the browser build script and the size check (SC-010).
  - Do: Steps: `cargo build --profile wasm-size --target wasm32-unknown-unknown --features wasm`; `wasm-bindgen --target web --out-dir pkg`; `wasm-opt -Oz` if present (skip with a message if not).
  - Do: `check_wasm_size.py` reads `pkg/*.wasm`, gzips it in memory, prints the size and exits 1 above 2 000 000 bytes. If the size is above, stop and report (OP2).
  - Files: `packages/loop_interpolation_core/scripts/build_wasm.sh`, `packages/loop_interpolation_core/scripts/check_wasm_size.py`
  - Needs: T041
  - Wave: 15
  - Done when: `bash packages/loop_interpolation_core/scripts/build_wasm.sh` creates `packages/loop_interpolation_core/pkg/` and `python -I packages/loop_interpolation_core/scripts/check_wasm_size.py` prints the gzip size and exits 0 when the size is at most 2 000 000 bytes.
  - Blocked by: OP1
- [ ] T045 [US3] [M] Run the browser tests against the built package and fix failures.
  - Do: Fix failures only in the two files of the Files list. If a fix needs another file, stop and report it.
  - Files: `packages/loop_interpolation_core/src/wasm.rs`, `packages/loop_interpolation_core/web/worker.js`
  - Needs: T044, T042, T039, T040, T043
  - Wave: 17
  - Done when: `bash packages/loop_interpolation_core/scripts/build_wasm.sh && npm --prefix packages/loop_interpolation_core/tests_web ci && npx --prefix packages/loop_interpolation_core/tests_web playwright install chromium && npm --prefix packages/loop_interpolation_core/tests_web test` passes.
  - Blocked by: OP1

## Phase 8: User story 5 - single precision

Phase E of the plan (was F). It changes the vector type, not the algorithm.

- [ ] T046 [P] [US5] [M] Test: single precision uses at most 60 percent of the memory, agrees with double precision, and reports honestly (story 5, FR-013, FR-022, SC-005).
  - Do: (1) 50^3 R1-type problem: peak of `F32` is at most 0.6 times the peak of `F64` (counting allocator). (2) The `F32` field differs from the `F64` field by at most 1e-3 of the range (SC-005). (3) For both modes `report.residual` differs from `Problem::residual_f64(&field)` (a new function that recomputes the residual in `f64`) by at most 1e-6 relative. (4) `F32` with `tol = 1e-12` gives `converged == false` and the report says so. (5) `Field` in `F32` mode holds a `Vec<f32>`.
  - Files: `packages/loop_interpolation_core/tests/precision_f32.rs`
  - Needs: T022, T025
  - Wave: 14
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test precision_f32` compiles and fails.
- [ ] T047 [P] [US5] [M] Add the `Scalar` trait and make `stencil.rs` generic over it.
  - Do: `trait Scalar: Copy + Send + Sync + Default + 'static` with `fn from_f64(f64) -> Self`, `fn to_f64(self) -> f64`, and `Add`, `Mul`, `Sub`, `AddAssign` bounds. Implement for `f32` and `f64`. No new dependency.
  - Do: Make `Stencil::forward_into`, `adjoint_accumulate_into` and the fused kernel taps generic in the vector type. Taps stay `f64`; the sum runs in `f64`.
  - Files: `packages/loop_interpolation_core/src/scalar.rs`, `packages/loop_interpolation_core/src/stencil.rs`
  - Needs: T046
  - Wave: 15
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml` passes (all `f64` results unchanged).
- [ ] T048 [P] [US5] [M] Make the regularisation operator of `interpolator.rs` generic over `Scalar`.
  - Do: Generic: `forward_all_families_into`, `adjoint_all_families_accumulate_into`, `fused_interior_regularisation_into`, `shell_regularisation_correction_into`, the border and directional passes. Keep the `f64` wrappers used by the ported tests.
  - Files: `packages/loop_interpolation_core/src/interpolator.rs`
  - Needs: T047
  - Wave: 16
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml` passes (all `f64` results unchanged).
- [ ] T049 [P] [US5] [M] Make the data and inequality operator parts generic over `Scalar`. Add `residual_f64`.
  - Do: Generic: `accumulate_data_normal`, `data_rhs`, `diagonal`, `apply_normal_operator_into`, the inequality accumulators. Add `pub fn residual_f64(&self, field: &Field) -> f64` that recomputes `|N x - rhs|` in `f64` with a streaming pass.
  - Files: `packages/loop_interpolation_core/src/interpolator.rs`
  - Needs: T048
  - Wave: 17
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml` passes.
- [ ] T050 [P] [US5] [M] Make `CgState`, `AdmmState` and the session generic. Select by `Settings.precision`. Dot products add in `f64`.
  - Do: `SolveSession` holds an enum `Inner { F64(State<f64>), F32(State<f32>) }`. `Field` uses `FieldData::F32` in single precision. The bound `F32` mode refuses `NaN` residual and reports `converged == false`.
  - Files: `packages/loop_interpolation_core/src/cg.rs`, `packages/loop_interpolation_core/src/admm.rs`, `packages/loop_interpolation_core/src/session.rs`, `packages/loop_interpolation_core/src/field.rs`
  - Needs: T049
  - Wave: 18
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --test precision_f32` passes.

## Phase 9: Finish - measurement, CI, docs, checks

- [ ] T051 [Finish] [M] Add the reference problem tests (SC-001, SC-008, SC-009, SC-012). Record the numbers.
  - Do: CI test: for 10 grid sizes from 10^3 to 100^3 the estimate is at least the peak and at most 1.25 times the peak (SC-009). Ignored test extends to 200^3.
  - Do: CI test: peak against dof for 50^3, 100^3, 150^3 fits a line with R^2 at least 0.99 (SC-012; 150^3 is `#[ignore]`d, CI uses 30^3, 45^3, 60^3).
  - Do: Ignored: R1 peak at most 100 MB (SC-001); 250^3 `F32` with `budget_bytes = 2_000_000_000` peak at most 1.5 GB (SC-008); time of R1 is printed (SC-007, report only).
  - Do: Write the field of R1 to `target/r1_field.f64` for the browser comparison.
  - Files: `packages/loop_interpolation_core/tests/reference_problems.rs`, `specs/001-low-memory-interpolator/measurements.md`
  - Needs: T023, T050, T017
  - Wave: 19
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --release --test reference_problems` passes (CI sizes). `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --release --test reference_problems -- --ignored --nocapture` prints the table that is copied to `measurements.md`.
- [ ] T052 [P] [Finish] [M] Add the browser reference run: R1 and R2 in Chromium (SC-001, SC-003, SC-006, SC-007).
  - Do: Generate the R1 data in JavaScript with the same formulas as `reference_problems.rs` (copy the formulas into both files and name them in a comment).
  - Files: `packages/loop_interpolation_core/tests_web/reference.spec.js`
  - Needs: T045, T051
  - Wave: 20
  - Done when: `npm --prefix packages/loop_interpolation_core/tests_web test -- reference` prints peak `memoryBytes` (at most 100 000 000 for R1), solve time, and the max difference from `target/r1_field.f64` (at most 1e-6 of the range). The spec is skipped unless the environment variable `RUN_REFERENCE=1` is set.
  - Blocked by: OP1
- [ ] T053 [P] [Finish] [M] Write the CI workflow with path filters.
  - Do: Trigger on `packages/loop_interpolation_core/**`, `packages/loop_interpolation/src/loop_interpolation/_native.py` and the workflow file. Also `workflow_dispatch`.
  - Do: `rust`: `cargo fmt --check`, `cargo clippy -- -D warnings` (no features, `parallel`, `wasm` on target `wasm32-unknown-unknown`, `python`), `cargo test` (with and without `parallel`).
  - Do: `wheel`: maturin build, install, run `pytest packages/loop_interpolation/tests` on ubuntu, windows and macos.
  - Do: `wasm`: `build_wasm.sh`, `check_wasm_size.py`, Playwright Chromium.
  - Do: `reference`: only on `workflow_dispatch`. Runs the ignored tests and the reference spec and writes the numbers to the job summary.
  - Files: `.github/workflows/native.yml`
  - Needs: T045, T038, T051, T044
  - Wave: 20
  - Done when: `python -c "import yaml; yaml.safe_load(open('.github/workflows/native.yml'))"` exits 0 and the file has jobs `rust`, `wheel`, `wasm`, `reference`.
- [ ] T054 [P] [Finish] [S] Write the user guide page for `solver="native"`.
  - Do: ASD-STE100 text: what the solver is for, how to install (`pip install loop-interpolation[native]`), how to select it (`solver="native"` before `setup`), the solver options, the unsupported settings list, and the memory estimate. Add `native_solver` to the toctree after `interpolation_options`.
  - Files: `docs/source/user_guide/native_solver.rst`, `docs/source/user_guide/index.rst`
  - Needs: T038
  - Wave: 18
  - Done when: `grep -q native_solver docs/source/user_guide/index.rst` exits 0 and the file `docs/source/user_guide/native_solver.rst` exists.
- [ ] T055 [P] [Finish] [S] Complete the core package README.
  - Do: Add: build commands (maturin, wasm), the three forms, the feature flags, the contract (flat arrays, grid coordinates only), and the list of differences from loop-rs (rotation removed, sessions, memory tracking).
  - Files: `packages/loop_interpolation_core/README.md`
  - Needs: T044, T038
  - Wave: 18
  - Done when: `grep -c "bff51342da8ebff7ca3c87e07c1d4fbd31279c8b" packages/loop_interpolation_core/README.md` prints 1 or more.
- [ ] T056 [Finish] [M] Run all checks and fix format and lint problems.
  - Do: Fix only format and lint problems. Name the changed files in the report. Do not change behaviour.
  - Do: Update the plan and tasks status in the index of `specs/` if an index exists (none exists today).
  - Files: none (runs alone; fixes only in files named in the task text)
  - Needs: T052, T053, T054, T055, T038, T001, T002, T026
  - Wave: 21
  - Done when: `cargo fmt --manifest-path packages/loop_interpolation_core/Cargo.toml --check`, `cargo clippy --manifest-path packages/loop_interpolation_core/Cargo.toml --all-targets -- -D warnings` (also with `--features parallel`), `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml`, `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml --features parallel`, `python -m pytest packages/loop_interpolation/tests -q`, and `python -m ruff check packages/loop_interpolation/src` all exit 0.

## Dependencies

- Phase 1 tasks have no dependencies, except that the crate files come before the port.
- The port chain is in one line: copy, remove rotation, rayon shim. The memory harness and the stubs both follow the rayon shim and run at the same time.
- The core API chain is: stubs, then tests for input and the memory tracker, then the code. Sessions need the estimate and the lazy tables. ADMM and the family checks follow the CG session. They edit the same files, so they run one after the other.
- The Python tests need only the package files, so they start early and skip until the wheel is built. The adapter needs the bindings. The two Python edits and the Python check follow.
- The browser tests start after the fixture exists. All browser code waits for OP1.
- Single precision follows ADMM and the family checks and does not need the browser.
- Phase 8 needs all earlier phases.

## Wave table
| Wave | Tasks | Tier [S] | Tier [M] |
|---|---|---|---|
| 1 | T001, T002, T003, T004, T005, T031, T033 | 6 | 1 |
| 2 | T006, T027, T028, T029, T030, T032 | 2 | 4 |
| 3 | T007, T037 | 0 | 2 |
| 4 | T008 | 0 | 1 |
| 5 | T009, T011 | 0 | 2 |
| 6 | T010, T012, T014 | 1 | 2 |
| 7 | T013, T015 | 1 | 1 |
| 8 | T016, T019 | 0 | 2 |
| 9 | T017, T018 | 0 | 2 |
| 10 | T020 | 0 | 1 |
| 11 | T021, T024, T026, T039 | 0 | 4 |
| 12 | T022, T040 | 0 | 2 |
| 13 | T023, T025 | 0 | 2 |
| 14 | T034, T041, T046 | 0 | 3 |
| 15 | T035, T042, T044, T047 | 0 | 4 |
| 16 | T036, T043, T048 | 1 | 2 |
| 17 | T038, T045, T049 | 0 | 3 |
| 18 | T050, T054, T055 | 2 | 1 |
| 19 | T051 | 0 | 1 |
| 20 | T052, T053 | 0 | 2 |
| 21 | T056 | 0 | 1 |

## Critical path

T003 -> T006 -> T007 -> T008 -> T011 -> T012 -> T013 -> T016 -> T017 -> T020 -> T021 -> T022 -> T025 -> T046 -> T047 -> T048 -> T049 -> T050 -> T051 -> T052 -> T056

(21 tasks long. It sets the minimum time.)

## Coverage of the spec

| Requirement | Tasks |
|---|---|
| FR-001, FR-023 (no matrix) | T010, T020, T051 |
| FR-002, SC-012 (linear memory) | T023, T051 |
| FR-003 (constraint kinds) | T012, T013, T024, T026, T027 |
| FR-004 (inequality kinds) | T021, T022, T028 |
| FR-005 (regularisation) | T026, T027 |
| FR-006 (no isolation, one thread) | T008, T039, T041 |
| FR-007 (thread count) | T024 |
| FR-008, FR-009, FR-010 (responsive, progress, cancel) | T019, T020, T040, T042 |
| FR-011, FR-012 (estimate, budget, no abort) | T014, T015, T016, T017, T018 |
| FR-013, FR-022, SC-005 (precision) | T046, T047, T048, T049, T050 |
| FR-014, FR-015 (report, infeasible) | T019, T021, T022 |
| FR-016 (evaluate) | T024, T025 |
| FR-017 (Python = browser) | T026, T027, T039 |
| FR-018 (bad input) | T012, T013 |
| FR-019 (no change without the solver) | T030, T038 |
| FR-020 (no coordinate change) | T007, T012, T043 |
| FR-021, FR-025 (select solver, same settings) | T032, T036, T037, T035, T033 |
| FR-024 (large origin) | T024 |
| FR-026 (unsupported setting) | T029, T036 |
| SC-001 (100 MB) | T051, T052 |
| SC-002, SC-011 | T023, T051 |
| SC-003 (agreement) | T005, T026, T027, T039 |
| SC-004 (1e-4) | T021, T028 |
| SC-006 (100 ms) | T040, T042 |
| SC-007 (60 s) | T051, T052 (report only; the release check on the reference machine is manual, OP3) |
| SC-008 (250^3 single) | T051 |
| SC-009 (estimate) | T016, T051 |
| SC-010 (2 MB) | T044, T053 |
| Story 1 scenario 3 (Firefox, Safari smoke) | manual, OP3 |
