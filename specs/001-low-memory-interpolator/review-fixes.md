# Review fixes: PR 354 (low-memory native solver)

**Input**: code review of PR 354 (10 findings), `spec.md`, `plan.md`
**Branch**: `feat/001-low-memory-interpolator` (fork `lachlangrose/LoopStructural`)

<!-- Written in ASD-STE100 Simplified Technical English. -->

## Decisions (confirmed by the user)

- **D-R1 (finding 5)**: Random sampling is accepted for the inequality pair bounds. A cap limits the rows for each pair. Use a fixed seed. Write a warning to the log when the cap applies.
- **D-R2 (finding 9)**: The version of `loop_interpolation_core` is `0.1.0` in the release manifest, `pyproject.toml` and `Cargo.toml`.

## Rules

- Work in a separate worktree of the PR branch. Use `PYTHONPATH` or a separate venv for the tests (the repo `.venv` is editable on the main checkout).
- Write a test first. See it fail. Then fix the code.
- One commit for each group (A1, A2, B, C). Conventional Commits. Scope: `loop-interpolation` for Python, `loop-interpolation-core` for Rust. Do not push before the user approves.
- Comments and docs: ASD-STE100.

## Findings and fixes

| Finding | Group | Fix | Files |
|---|---|---|---|
| 1 | A1 | Reject interpolators that the core cannot represent | `_native.py`, `_discrete_interpolator.py` |
| 2 | A1 | Make setup and solve agree on the solver | `_finite_difference_interpolator.py`, `_discrete_interpolator.py` |
| 10 | A1 | Reset history and timing on the native path | `_discrete_interpolator.py` |
| 8 | A1 | Keep the old field when the solve does not converge | `_native.py` |
| 6 | A1 | Cast the field to float64 and check its length | `_native.py` |
| 7 | A2 | Make a failed `finalize` final | `session.rs` |
| 3 | B | Filter bad rows like the matrix path | `_native.py` |
| 4 | B | Guard the interface family | `_native.py` |
| 5 | B | Cap and chunk the pair bounds | `_native.py` |
| 9 | C | Align the versions | `.release-please-manifest.json`, `.release-please-config.json` |

## Phase R1: Group A1 (Python, wrong results)

- [x] R001 [M] Test and fix finding 1: unsupported interpolators.
  - Test: a fold interpolator (`DFI`) with solver `native` raises `NotImplementedError`. A P1 or P2 interpolator raises. A plain `FiniteDifferenceInterpolator` on a uniform grid works.
  - Test: an interpolator with sparse-only constraints (from `add_constraints_to_least_squares` or `add_inequality_feature`) raises `NotImplementedError`.
  - Do: In `check_supported`, accept only the exact class `FiniteDifferenceInterpolator` (or a subclass that sets a flag `supports_native = True`). Do not use `hasattr` on the support as the only test.
  - Do: In `check_supported`, raise when `len(interpolator.constraints) > 0` and the rows are not from the raw data. Replace the warning in `solve_system` with this error.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_native.py`, `packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py`, `packages/loop_interpolation/tests/test_native.py`
  - Needs: none
  - Wave: 1
- [x] R002 [M] Test and fix finding 2: setup and solve disagree.
  - Test: setup with `native`, then `solve_system(solver="cg")` raises `ValueError` with a message to run setup again.
  - Test: setup with `cg`, then `solve_system(solver="native")` still works and warns (existing behaviour).
  - Test: setup with `native` runs `add_inequality_pairs_constraints` and the other setup steps (move the early return after them, or call them in the native branch).
  - Do: Set `self._native_setup = True` when setup skips the sparse rows. Clear it in `reset()` and at the start of each `setup_interpolator` call.
  - Do: In `solve_system`, raise when `_native_setup` is true and the solver is not `native`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_finite_difference_interpolator.py`, `packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py`, `packages/loop_interpolation/tests/test_native.py`
  - Needs: R001
  - Wave: 2
- [x] R003 [S] Test and fix finding 10: stale state.
  - Test: run `cg` (history set), then `native`. `solver_history` is `None` after the native solve. `latest_solve_timing` has no key from the `cg` run.
  - Do: In the native branch of `solve_system`, set `self.solver_history = None` and call `_solver_pipeline.init_timing("native")`. Then merge the native timing.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py`, `packages/loop_interpolation/tests/test_native.py`
  - Needs: R002
  - Wave: 3
- [x] R004 [M] Test and fix finding 8: no change of `c` on failure.
  - Test: a solve with `maxiter=1` returns `False`, and `interpolator.c` has the old value. `latest_native_report` has the partial field.
  - Do: In `_native.solve`, write `interpolator.c` only when `converged` is true. Store the partial field in the report (key `field`) when not converged.
  - Do: Make the default `maxiter` and `tol` the same as the matrix path, or write the difference in the docstring of `solve`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_native.py`, `packages/loop_interpolation/tests/test_native.py`
  - Needs: R001
  - Wave: 2
- [x] R005 [S] Test and fix finding 6: dtype and length of `c`.
  - Test: `solver_kwargs={"precision": "f32"}` gives `interpolator.c.dtype == np.float64`. A wrong length from a mocked core raises `RuntimeError`.
  - Do: `field = np.asarray(field, dtype=np.float64)`. Compare `field.shape[0]` with `support.n_nodes`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_native.py`, `packages/loop_interpolation/tests/test_native.py`
  - Needs: R004
  - Wave: 3

**Checkpoint A1**: `python -m pytest packages/loop_interpolation/tests -q` passes, with the core installed and without it.

## Phase R2: Group A2 (Rust)

- [x] R006 [M] Test and fix finding 7: `finalize` error in `SolveSession::step`.
  - Test (Rust): make `try_zeroed` fail in a budgeted f32 solve (use a small `budget_bytes`). Call `step` two times. Both calls return an `Err`. No panic. No zero-length field is reported as converged.
  - Do: Add a field `failed: Option<CoreError>` to `SolveSession`. When `finalize` returns `Err`, store the error. At the start of `step`, return the stored error.
  - Do: Reserve the `back` buffer before `x32` is dropped, so that the failure point is removed. Keep the stored error as a second guard.
  - Files: `packages/loop_interpolation_core/src/session.rs`
  - Needs: none
  - Wave: 1
  - Done when: `cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml` and `cargo clippy --manifest-path packages/loop_interpolation_core/Cargo.toml -- -D warnings` exit 0.

## Phase R3: Group B (checks in `_native.py`)

- [x] R007 [M] Test and fix finding 3: bad rows.
  - First read the matrix path. Find how it treats a negative weight, a NaN vector and a zero vector. Write the result in a comment.
  - Test: rows with a negative weight, a NaN vector and a zero vector do not raise. The native field is the same as the `cg` field on the same data (small tolerance).
  - Do: In `_add_weighted`, keep a row only if the weight is finite and more than 0, and (for vectors) all values are finite and the norm is more than 0. Log the number of dropped rows at debug level.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_native.py`, `packages/loop_interpolation/tests/test_native.py`
  - Needs: R004
  - Wave: 3
- [x] R008 [M] Test and fix finding 4: interface family.
  - Test: negative `ipw` skips the family (no error). A NaN id inside a group is dropped. A group with one valid point is dropped. Mixed `ipw` values still raise `NotImplementedError`.
  - Do: Skip the family when `ipw <= 0`. Drop the NaN ids first. Then drop groups with fewer than two points. Do this for each id.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_native.py`, `packages/loop_interpolation/tests/test_native.py`
  - Needs: R007
  - Wave: 4
- [x] R009 [M] Test and fix finding 5: pair-bound cross product.
  - First read `add_inequality_pairs_constraints`. Use the same pairing rule.
  - Test: 5 000 points for each of two rock ids. The peak Python memory (`tracemalloc`) stays below a set limit (for example 200 MB). A warning is in the log.
  - Test: a small case (below the cap) gives the same field as the matrix path.
  - Test: two runs with the same input give the same field (fixed seed).
  - Do: Add the setup keyword `inequality_pair_max_rows` (default 100 000 for each pair). If `n_a * n_b` is more than the cap, pick the pairs at random with `np.random.default_rng(seed)`. Pick the indices, do not build the full product.
  - Do: Send the rows in chunks of 100 000. Build each chunk with index arithmetic. Do not call `np.repeat` or `np.tile` on the full set.
  - Do: Write the keyword and the sampling rule in the docstring of the module.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_native.py`, `packages/loop_interpolation/tests/test_native.py`
  - Needs: R007
  - Wave: 4

**Checkpoint B**: `python -m pytest packages/loop_interpolation/tests -q` passes. The native-versus-`cg` comparison tests pass.

## Phase R4: Group C (release files)

- [x] R010 [S] Fix finding 9: align the versions to `0.1.0`.
  - Do: In `.release-please-manifest.json`, set `packages/loop_interpolation_core` to `0.1.0`.
  - Do: In `packages/loop_interpolation_core/Cargo.toml`, set `version = "0.1.0"`.
  - Do: In `.release-please-config.json`, add `extra-files` to the package `packages/loop_interpolation_core`: `{"type": "toml", "path": "Cargo.toml", "jsonpath": "$.package.version"}`. Keep the `exclude-paths` entry in the root package.
  - Do: Confirm that `pyproject.toml` of the package has `version = "0.1.0"`.
  - Files: `.release-please-manifest.json`, `.release-please-config.json`, `packages/loop_interpolation_core/Cargo.toml`
  - Needs: none
  - Wave: 1
  - Done when: `python -c "import json,tomllib; m=json.load(open('.release-please-manifest.json')); c=tomllib.load(open('packages/loop_interpolation_core/Cargo.toml','rb')); p=tomllib.load(open('packages/loop_interpolation_core/pyproject.toml','rb')); assert m['packages/loop_interpolation_core']==c['package']['version']==p['project']['version']=='0.1.0'"` exits 0.

## Waves

| Wave | Tasks | Notes |
|---|---|---|
| 1 | R001, R006, R010 | Different files. Run in parallel. |
| 2 | R002, R004 | Both need R001. They change different parts of `_native.py` and `_discrete_interpolator.py`, so run them one after the other if a merge conflict occurs. |
| 3 | R003, R005, R007 | |
| 4 | R008, R009 | Same file. Run one after the other. |

## Before done

- `python -m pytest packages/loop_interpolation/tests -q` (with and without the core installed)
- `cargo fmt --check`, `cargo clippy -- -D warnings`, `cargo test` (all with `--manifest-path packages/loop_interpolation_core/Cargo.toml`)
- Update `measurements.md` only if R009 changes a memory number.
- Ask the user before the push to the PR branch.
