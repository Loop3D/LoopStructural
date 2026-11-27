# Tasks: native fault application

**Feature**: `007-native-fault-kernel` | **Plan**: [plan.md](plan.md)

<!--
Scope: stages S1 to S4 of the plan (stories 1 to 5), T001 to T033.
Code base: branch `feat/007-native-fault-kernel`, worktree
`/home/lgrose/dev/LoopStructural-wt/feat-007-native-fault-kernel` (master 6874364c plus
907acf42, the 94-fault benchmark `benchmarks/fault_evaluation/many_faults.py`).
Python: `../feat-003-fault-evaluation/.venv` with
`PYTHONPATH=$PWD:$PWD/packages/loop_common/src:$PWD/packages/loop_interpolation/src`
(see HANDOFF.md), or a new venv. Use `OMP_NUM_THREADS=1` with `pytest -n`.
Measure on the T030 computer (Intel i7-12700H, Linux) with no other work running.
Python test command: `pytest <file> -q` from the repository root.
Rust commands run in `packages/loop_interpolation_core`; one cargo command at a time.
The crate is built with no features, with `parallel`, with `python` and for wasm
(`.github/workflows/native.yml`): fault code must compile without `parallel` (sequential
loop) and must not use `pyo3` outside `src/python.rs`.
Python tests that need the kernel: build it first in the active virtual environment with
`maturin develop --release -m packages/loop_interpolation_core/Cargo.toml`. These tests
skip with `pytest.importorskip("loop_interpolation_core")` and must not skip in the task
check.
Common rules for every Python task: `from __future__ import annotations` first;
numpy-style docstrings; type hints on public names; Python 3.9 syntax (no `match`);
line length 100; keep string quotes as written (Black skip-string-normalization);
`ruff check LoopStructural packages` finds nothing new.
Common rules for every Rust task: `cargo fmt --check` and
`cargo clippy --all-targets -- -D warnings` pass; no new crate dependency.
-->

## Orchestrator notes

- Run this file with the spec-orchestrator skill. It has the rules for the
  worktree, the sub agents, the checks and the commits.
- Waves: 11. Tiers: [S] = simple model, [M] = standard model.
- Do not start a task with "Blocked by".

## Phase 1: Setup

- [x] T001 [Setup] [S] Record the decision to add the native fault path in `ROADMAP.md` (FR-015).
  - Files: `ROADMAP.md`
  - Needs: none
  - Wave: 1
  - Done when: `ROADMAP.md` has the section `### Native fault application (spec 007)` under `## Decisions`, and `git diff --stat` shows only `ROADMAP.md`.
  - Context: Add a section after the last `###` section of `## Decisions`
    (before `## Stage sequence`). Content, 5 to 10 lines: the native path is
    an optional Rust kernel in `loop_interpolation_core`; the Python path
    stays the reference and the fallback; on by default when the kernel is
    installed and `FAULT_KERNEL_VERSION` matches; off with
    `LoopStructuralConfig.native_faults = False` or the environment variable
    `LOOPSTRUCTURAL_NATIVE_FAULTS=0` (the setting wins); tolerance 1e-8 x the
    model box (spec 003 FR-004); stages S1 (restore order once), S2 (native
    restore), S3 (native values and gradients), S4 (wheels). Link
    `specs/007-native-fault-kernel/spec.md`. ASD-STE100 style.

## Phase 2: Stage S1, restore order once (story 1)

- [x] T002 [US1] [M] Measure the cut faults of the 94-fault model before stage S1 (SC-001 baseline).
  - Files: `specs/007-native-fault-kernel/profile-s1.md`
  - Needs: none
  - Wave: 1
  - Done when: `profile-s1.md` exists with the baseline numbers at 500 points and a `cProfile` top-20 table (by cumulative time) of one `evaluate_model` call at 500 points.
  - Context: Baseline (Python path, master): copy from
    `/home/lgrose/dev/LoopStructural-wt/_runs007/many_faults_python.json`
    (`cases["500"]`: `cutting_mean_ms` 40.0, `cutting_max_ms` 117.4,
    `isolated_mean_ms` 2.65, `evaluate_model_seconds` 1.78). The bit-identity
    reference is `_runs007/many_faults_python_ref.npz`. Profile: in a script
    in `$CLAUDE_JOB_DIR/tmp` or `/tmp`, call `many_faults.build()` and
    `many_faults.points(500)` (`benchmarks/fault_evaluation/many_faults.py`;
    the build takes about 7 minutes), then `cProfile` one
    `model.evaluate_model(pts)`. Report the share of time in `with_cutters`,
    `cutters_first`, `_cutters`, `_evaluate_frame_legacy` and
    `restore_through` (`LoopStructural/utils/_fault_block.py`,
    `LoopStructural/modelling/features/fault/_fault_segment.py`). Record the
    computer and the commit. Change no code.

- [x] T003 [US1] [M] Test: the restore order is computed once for each state of the faults (FR-001, FR-002).
  - Files: `tests/fault_evaluation/test_restore_order.py`
  - Needs: none
  - Wave: 1
  - Done when: `pytest tests/fault_evaluation/test_restore_order.py -q` fails (no `restore_order`), then passes after T004.
  - Context: New function (T004) in `LoopStructural/utils/_fault_block.py`:
    `restore_order(faults, *, reverse=False, owner=None) -> list`. With
    `owner`, it keeps `owner._restore_order = (key, chain, ordered)`.
    Tests: (1) on `build_model("cutting_faults")` (`tests/fault_evaluation/models.py`),
    evaluate a faulted feature twice; with `unittest.mock.patch` on
    `LoopStructural.utils._fault_block.cutters_first` (wraps=original), the
    second call does not call it. (2) Story 1 scenario 2: add a new cutter to
    a fault builder (`fault.builder.builders[i].faults.append(other)`) after
    the first call; the next `restore_order` result differs and puts `other`
    before the fault. (3) Bit-identical: `capture(model, points)`
    (`tests/fault_evaluation/golden.py`) for `cutting_faults` and `f4_u3`
    with the cache equals a run with `owner=None`
    (`np.testing.assert_array_equal`, `equal_nan=True`). (4) `reverse=True`
    gives the reversed order rule of `restore_through`.

- [ ] T004 [US1] [M] Add `restore_order` with the order cache and pass `owner` from the two callers (D1a).
  - Files: `LoopStructural/utils/_fault_block.py`, `LoopStructural/modelling/features/_base_geological_feature.py`, `LoopStructural/modelling/features/_evaluation.py`
  - Needs: T003
  - Wave: 2
  - Done when: `pytest tests/fault_evaluation -q` passes.
  - Context: `restore_order(faults, *, reverse=False, owner=None) -> list`.
    Without owner: `cutters_first(with_cutters(faults))`; if reverse,
    `cutters_first(with_cutters(faults)[::-1])[::-1]`. With owner: compute
    `chain = with_cutters(faults)` and
    `key = (bool(reverse), tuple((id(f), tuple(id(c) for c in _cutters(f))) for f in chain))`;
    if `getattr(owner, "_restore_order", None)` has an equal key, return its
    `ordered`; else compute and set `owner._restore_order = (key, chain, ordered)`
    in one assignment (thread-safe). `restore_through(points, faults, *,
    reverse=False, block=None, restored_through=(), owner=None)` uses
    `restore_order` in place of its two order lines; nothing else changes.
    `BaseFeature._apply_faults` passes `owner=self`. `EvaluationContext.restore`
    passes `owner=feature` in both `restore_through` calls. Do not change the
    builder calls.

- [ ] T005 [US1] [M] Measure the cut faults again after D1a and compare with SC-001.
  - Files: `specs/007-native-fault-kernel/profile-s1.md`
  - Needs: T002, T004
  - Wave: 3
  - Done when: `profile-s1.md` has an "After D1a" section with the numbers of T002 and one line "SC-001: pass" or "SC-001: fail (<mean> ms)".
  - Context: Run
    `python -m benchmarks.fault_evaluation.many_faults --output /tmp/many_faults_d1a.json --points 500 --compare /home/lgrose/dev/LoopStructural-wt/_runs007/many_faults_python_ref.npz`
    (about 7 minutes; exit 1 means a model id or group value differs bit for
    bit). SC-001: `cases["500"]["cutting_mean_ms"]` at most 10 ms (baseline
    40.0) and exit code 0. Also run the T002 profile again and give the new
    shares. Change no code.

- [ ] T006 [US1] [M] Share one restoration among the members of a cut frame, only if T005 reports that SC-001 fails (D1b, OP3).
  - Files: `LoopStructural/modelling/features/fault/_fault_segment.py`, `tests/fault_evaluation/test_restore_order.py`, `specs/007-native-fault-kernel/profile-s1.md`
  - Needs: T005
  - Wave: 4
  - Done when: if T005 says "SC-001: pass", `profile-s1.md` has the line "D1b: not needed" and no code changes; else `pytest tests/fault_evaluation -q` passes, the `many_faults --compare` run of T005 exits with 0, and `profile-s1.md` has an "After D1b" section with "SC-001: pass" or "fail".
  - Context: In `FaultSegment._evaluate_frame_legacy`: if the 3 members have
    the same non-empty fault list (each item compared with `is`), compute
    `restored = restore_through(points, m0.faults, block=block,
    restored_through=restored_through, owner=m0)` once, and call each
    member (`evaluate_value`, and `evaluate_gradient(ignore_regions=True)`
    of member 1) with `restored` and
    `restored_through=tuple(restored_through) + tuple(restore_order(m0.faults, owner=m0))`.
    The masks must still use the original points; if this cannot be done
    without a change of a member method, stop and report. Add a test to
    `test_restore_order.py`: `capture` of `cutting_faults` is bit-identical
    with and without D1b. If a difference shows, revert D1b and record
    "D1b: dropped (not bit-identical)".

- [ ] T007 [US1] [S] Record stage S1 in `ROADMAP.md` (FR-015).
  - Files: `ROADMAP.md`
  - Needs: T001, T006
  - Wave: 5
  - Done when: `ROADMAP.md` `## Status log` has a dated line for spec 007 stage S1 with the SC-001 result from `profile-s1.md`.
  - Context: Add one line at the end of `## Status log`, format as the lines
    above it: date 2026-MM-DD of the day, "spec 007 S1: restore order
    computed once for each state of the faults (`restore_order`)", the
    SC-001 numbers (before and after) and "D1b: used / not needed / dropped"
    from `specs/007-native-fault-kernel/profile-s1.md`.

## Phase 3: Foundation, Rust module skeleton

- [x] T008 [Found] [M] Add the Rust module `fault` with all types and function stubs (D2a).
  - Files: `packages/loop_interpolation_core/src/fault/`, `packages/loop_interpolation_core/src/lib.rs`
  - Needs: none
  - Wave: 1
  - Done when: `cargo check --all-targets` and `cargo check --all-targets --features parallel` pass in `packages/loop_interpolation_core`.
  - Context: `lib.rs`: add `pub mod fault;`. Files in `src/fault/`:
    `mod.rs` (`pub struct FaultChain { pub faults: Vec<Fault>, pub block: BlockSettings }`,
    `pub struct Fault { pub frame: [GridField; 3], pub cutters: Vec<usize>, pub regions: Vec<Region>, pub profile: Profile, pub displacement: f64, pub steps: usize }`,
    `pub struct BlockSettings { pub side_epsilon: f64, pub extension: Extension, pub exact: f64, pub limit: f64, pub step: f64 }`,
    `pub enum Extension { Smooth, Mirror }`, `pub enum Side { Free, Absent, Hanging, Foot }`),
    `grid_field.rs` (`GridField` with fields of plan D2a, `value`, `gradient`),
    `profile.rs` (`Profile1`, `Profile`), `region.rs` (`Region` enum with
    `Sign` and `Unconformity`), `restore.rs` (`apply`, `restore`,
    `frame_value` on `FaultChain`), `nearest.rs` (`KdTree::new(&[[f64;3]])`,
    `nearest(&self, p) -> usize`), `gradient.rs` (`feature_gradient`),
    `pool.rs` (`#[cfg(feature = "parallel")] pub fn pool(threads: usize) -> Arc<rayon::ThreadPool>`).
    Exact signatures: plan.md D2a. Stub bodies: `unimplemented!("T0xx")`
    with the id of the task that writes the body. Derive `Debug, Clone` on
    data types. Doc comment on each public item.

## Phase 4: Stage S2, native restoration (stories 2, 4 and 5)

- [x] T009 [US2] [M] Find why the proof of concept differs from the Python path (C-007, OP9).
  - Files: `specs/007-native-fault-kernel/verification.md`
  - Needs: none
  - Wave: 1
  - Done when: `verification.md` has a section "C-007: proof-of-concept difference" with a table: item, Python behaviour (file:line), proof-of-concept behaviour (file:line), effect on the result.
  - Context: Read only, do not change: the proof of concept in
    `/home/lgrose/dev/LoopStructural-wt/_runs003/rust_poc/` (`README.md`,
    `loop-snapshot/`, `pyext/`, `export_snapshot.py`) and the Python reference:
    `LoopStructural/modelling/features/fault/_fault_segment.py`
    (`apply_to_points`, `_evaluate_frame`, `_displacement_fraction`),
    `.../fault/_fault_function.py` (`CubicFunction.__call__`,
    `_limit_value`, `continue_profile`, `block_displacement`, `_saturate`,
    `_taylor`), `LoopStructural/utils/_fault_block.py`. Check at least: the
    step direction of #425, `BLOCK_EXTENSION`, `_limit_value`, the local
    frame (`_project_points`, rotation), f32 against f64, the order of the
    steps and of the faults. If the folder is not there, write
    "not available" and list only the risks.

- [ ] T010 [US2] [M] Port the structured grid lookup to `grid_field.rs`, with Rust tests.
  - Files: `packages/loop_interpolation_core/src/fault/grid_field.rs`, `packages/loop_interpolation_core/tests/fault_grid_field.rs`
  - Needs: T008
  - Wave: 2
  - Done when: `cargo test --release --test fault_grid_field` fails against the stub, then passes, also with `--features parallel`.
  - Context: Port the grid methods `_lookup_points`, `trilinear`,
    `_values_from_coefficients`, `_evaluate_gradient_from_lookup`
    (`packages/loop_common/src/loop_common/supports/_3d_structured_grid.py`,
    `packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py`)
    and the interpolator methods `_project_points` and `_reproject_vectors`
    (world to local frame with `local_origin` and `local_rotation`; `None`
    means identity). `value(p)` NaN if a coordinate is NaN or the point is
    outside the grid; `gradient(p)` the same, returned in world axes.
    Tests: a field linear in x, y, z (exact value and gradient), a rotated
    local frame, a NaN point, an outside point, and 5 points with expected
    values copied from a Python run of the same grid (write the Python
    snippet in a comment of the test).

- [ ] T011 [US2] [M] Port the fault displacement profiles to `profile.rs`, with Rust tests.
  - Files: `packages/loop_interpolation_core/src/fault/profile.rs`, `packages/loop_interpolation_core/tests/fault_profile.rs`
  - Needs: T008
  - Wave: 2
  - Done when: `cargo test --release --test fault_profile` fails against the stub, then passes.
  - Context: `Profile1 { Ones, Zeros, SmoothPeak, Cubic { w: [f64;4], lim: [f64;2], min: f64, max: f64, value_points: Vec<(f64,f64)> }, Composite { positive: Box<Profile1>, negative: Box<Profile1> } }`,
    `Profile { scale, gx, gy, gz }`. Port from
    `LoopStructural/modelling/features/fault/_fault_function.py`:
    `smooth_peak`, `CubicFunction.__call__` (limits through `_limit_value`,
    then floor and ceiling), `Composite.__call__` (0 at v == 0),
    `FaultDisplacement.__call__` (`scale * gx * gy * gz`), `_saturate`,
    `_taylor`, `continue_profile` (including the `allclose` check at
    `side * [BLOCK_SIDE_EPSILON, BLOCK_EXTENSION_STEP]` with atol 1e-12) and
    `block_displacement` (Smooth and Mirror). Tests: values of
    `BaseFault.fault_displacement` and `BaseFault3D.fault_displacement` at
    gx, gy, gz in {-1.5, -1, -0.5, 0, 0.3, 1, 2}, and `block_displacement`
    for side +1 and -1, compared with values from Python to 1e-14.

- [ ] T012 [US2] [M] Implement the native restoration, regions and thread pool (`restore.rs`, `region.rs`, `pool.rs`).
  - Files: `packages/loop_interpolation_core/src/fault/restore.rs`, `packages/loop_interpolation_core/src/fault/region.rs`, `packages/loop_interpolation_core/src/fault/pool.rs`, `packages/loop_interpolation_core/src/fault/mod.rs`, `packages/loop_interpolation_core/tests/fault_restore.rs`
  - Needs: T010, T011
  - Wave: 3
  - Done when: `cargo test --release --test fault_restore` passes, also with `--features parallel`.
  - Context: Port `FaultSegment.apply_to_points` for one point (`force_side`,
    `_displacement_fraction`, `steps` Euler steps along the normalised
    gradient of frame member 1, `reverse`, block sides, NaN rules of #413,
    and #425: a NaN gradient of a cut frame uses the plain grid gradient of
    member 1) and the loop of `restore_through` (skip faults in
    `restored`, mark each fault done after it). Frame values of a fault with
    `cutters` not all in `restored`: restore the point through them first
    (`_evaluate_frame` semantics). `Region::Sign`: port
    `BaseSignRegion._calculate_value_and_distance` and
    `PositiveRegion`/`NegativeRegion.__call__` (`LoopStructural/utils/regions.py`).
    `Region::Unconformity`: port `UnconformityFeature.evaluate`
    (`.../features/_unconformity_feature.py`): `<= value` if `sign`, else
    `>= value`, NaN gives false. `pool(0)` = all cores; cache by thread
    count. Tests: one fault with constant displacement (exact shift along
    the slip direction), reverse undoes restore to 1e-9, NaN in gives NaN
    out, a block side 0 moves nothing, empty input.

- [ ] T013 [US2] [M] Add the PyO3 classes `FaultChain` and `GridField` and the constant `FAULT_KERNEL_VERSION = 1` (D2b).
  - Files: `packages/loop_interpolation_core/src/python.rs`
  - Needs: T012
  - Wave: 4
  - Done when: `cargo clippy --all-targets --features python -- -D warnings` passes, and after `maturin develop --release -m packages/loop_interpolation_core/Cargo.toml`, `python -c "import loop_interpolation_core as c; assert c.FAULT_KERNEL_VERSION == 1; c.FaultChain; c.GridField"` passes.
  - Context: Register in `#[pymodule] fn loop_interpolation_core` next to
    `PyProblem`. `FaultChain(spec: dict)`: parse the dict of plan D2b
    (keys `version`, `block`, `faults[*].displacement/steps/frame/cutters/regions/profile`;
    region `kind` = positive / negative / unconformity; profile `kind` =
    ones / zeros / smooth_peak / cubic / composite); `version != 1` or a
    bad key raises `ValueError`. `restore(points (N,3) f64, order i64,
    reverse bool, sides i8, restored bool, threads int) -> (N,3)`,
    `frame_value(fault, coordinate, points, sides, restored, threads) -> (N,)`.
    `GridField(spec)` with `values(points, threads) -> (N,)`,
    `gradients(points, threads) -> (N,3)`. Copy the points to a `Vec`, then
    run the work in `py.detach(...)` (as `PyProblem.solve` does). Internal
    error: `RuntimeError` with a message, no partial result. Empty input:
    empty output.

- [x] T014 [US4] [S] Test: the native setting, the variable, the version check and the thread count (FR-006, FR-006a, FR-006b, FR-012).
  - Files: `tests/fault_evaluation/test_native_settings.py`
  - Needs: none
  - Wave: 1
  - Done when: `pytest tests/fault_evaluation/test_native_settings.py -q` fails (no module `_native_faults`), then passes after T015.
  - Context: Module `LoopStructural.modelling.features._native_faults` with
    `REQUIRED_FAULT_KERNEL_VERSION = 1`, `native_enabled() -> bool`,
    `native_threads() -> int`. Settings `LoopStructuralConfig.native_faults`
    (None default) and `.native_threads` (None default). Use a fake module
    in `sys.modules["loop_interpolation_core"]` (`types.SimpleNamespace` with
    `FAULT_KERNEL_VERSION`) and `monkeypatch`; reset the module cache of
    `_native_faults` between tests (`importlib.reload`). Cases: variable
    `LOOPSTRUCTURAL_NATIVE_FAULTS` in {"0","false","OFF","no"} gives False;
    unset gives True; setting True with variable "0" gives True (setting
    wins); setting False gives False; import fails (`sys.modules[...] = None`)
    gives False and no warning; version 2 gives False and exactly one
    warning that contains "2" and "1" over two calls (`caplog`);
    `native_threads()` is 0 for None and 1 for 1.

- [ ] T015 [US4] [M] Add the settings and `native_enabled`, `native_threads` (D2c, D2e).
  - Files: `LoopStructural/__init__.py`, `LoopStructural/modelling/features/_native_faults.py`
  - Needs: T014
  - Wave: 2
  - Done when: `pytest tests/fault_evaluation/test_native_settings.py tests/unit/test_public_api_contract.py -q` passes.
  - Context: `LoopStructuralConfig` (dataclass in `LoopStructural/__init__.py`):
    add `native_faults: bool | None = None` and `native_threads: int | None = None`
    with numpy docstring entries; the file must keep Python 3.9 syntax
    (use `Optional[...]` if the file has no `from __future__ import annotations`).
    New module `_native_faults.py`: `REQUIRED_FAULT_KERNEL_VERSION: int = 1`.
    `native_enabled()`: setting if not None; else off if
    `os.environ.get("LOOPSTRUCTURAL_NATIVE_FAULTS", "").lower()` in
    {"0","false","off","no"}; else on. If on: import
    `loop_interpolation_core` once (keep the result in a module variable);
    import error: off, no warning; `getattr(core, "FAULT_KERNEL_VERSION", None)`
    not equal to the required version: off and one warning for each session
    (module flag) through `LoopStructural.utils.getLogger(__name__)`:
    `"native fault kernel version {found} found, LoopStructural needs version {required}; using the Python path"`.
    `native_threads()`: the setting, or 0 if None.

- [ ] T016 [US2] [M] Add `chain_for`, `restore` and the state token, and call them from the Python path (D2c, D2d).
  - Files: `LoopStructural/modelling/features/_native_faults.py`, `LoopStructural/modelling/features/_base_geological_feature.py`, `LoopStructural/modelling/features/_evaluation.py`
  - Needs: T004, T017, T018
  - Wave: 6
  - Done when: `pytest tests/fault_evaluation -q` passes with the kernel built (`test_native_faults.py` and `test_native_safety.py` do not skip).
  - Context: Follow plan.md D2c and D2d exactly (read those two sections).
    `chain_for(feature) -> tuple[FaultChain, list] | None` builds the spec
    dict of plan D2b from the faults of
    `restore_order(feature.faults, owner=feature)`; support checks as in
    D2c (FaultSegment; 3 members on `StructuredGrid` with
    `_has_uniform_lookup()`; equal member fault lists; regions Sign or
    Unconformity only; `FaultDisplacement` with `Ones`, `Zeros`,
    `CubicFunction` (call `solve()` first), `Composite`, or `smooth_peak` by
    `is`). Cache `feature._native_chain = (token, refs, chain, faults)`;
    token and refs as D2c. `restore(feature, points, *, reverse=False,
    block=None, restored_through=())` returns None for the Python path.
    `_apply_faults`: try native first. `_evaluation.py`: new private
    `_restore_direct(feature, points, reverse)` used by both calls.

- [ ] T017 [US2] [M] Test: native restoration agrees with the Python path on the golden models (story 2, FR-003, FR-005, FR-005a, FR-005c).
  - Files: `tests/fault_evaluation/test_native_faults.py`
  - Needs: T013, T015
  - Wave: 5
  - Done when: `pytest tests/fault_evaluation/test_native_faults.py -q` fails (no `_native_faults.restore`), then passes after T016.
  - Context: `pytest.importorskip("loop_interpolation_core")`. Models:
    `f1_u0`, `f4_u0`, `f4_u3`, `cutting_faults`, `plugin`, `stopped_fault`,
    `faulted_fold` (`build_model`, `sample_points(name, 2000, seed=1)` in
    `tests/fault_evaluation/models.py`). For each faulted feature: restored
    points with `_native_faults.restore(feature, pts)` against
    `restore_through(pts, feature.faults, owner=None)` with
    `LoopStructuralConfig.native_faults = False`. Tolerance
    `1e-8 * box_size` (as `assert_matches` in `golden.py`); NaN rows equal.
    Also `block={fault.name: +1}`, `{name: -1}`, `{name: 0}`, `reverse=True`
    and `restored_through=(feature.faults[0],)`. Record the largest
    difference: `record_property("max_native_diff", value)`. Story 4
    scenarios 1 to 3 with the real kernel: native used (the function returns
    an array), off by setting, off by variable.

- [ ] T018 [US5] [M] Test: the native path is safe to use (story 5 scenarios 2 to 7, edge cases).
  - Files: `tests/fault_evaluation/test_native_safety.py`
  - Needs: T013, T015
  - Wave: 5
  - Done when: `pytest tests/fault_evaluation/test_native_safety.py -q` fails before T016 (the native path is not used), then passes after T016.
  - Context: `pytest.importorskip("loop_interpolation_core")`. Model
    `build_model("f4_u0")`. (2) A fault whose `faultfunction` is a plain
    Python function: `chain_for` returns None and the values equal the
    native-off values. (3) Second call: no file opened (patch
    `builtins.open` to fail) and `feature._native_chain[2]` is the same
    object. (4) Change `fault.displacement`; the next values equal the
    Python path for the new state. (5) A thread that counts in a loop makes
    progress during a native `restore` on 10 million points. (6) Two threads
    evaluate one model; each result equals the result alone. (7)
    `LoopStructuralConfig.native_threads = 1` gives the same result as 0.
    Edge: empty points; a NaN coordinate; `FaultChain({"version": 99})`
    raises `ValueError` and the model state does not change. Reset the
    settings in a fixture.

- [ ] T019 [US4] [S] Add the option `--native-faults {on,off}` to the benchmark runner (FR-013).
  - Files: `benchmarks/fault_evaluation/run.py`
  - Needs: T015
  - Wave: 3
  - Done when: `python -m benchmarks.fault_evaluation.run --help` lists `--native-faults`, and `pytest tests/fault_evaluation/test_benchmark_run.py tests/fault_evaluation/test_benchmark_suites.py -q` passes.
  - Context: In `main(argv)`, after the `--suite` argument:
    `parser.add_argument("--native-faults", choices=("on", "off"), default="off", help="Use the native fault kernel (on) or the Python path (off, default).")`.
    After `parse_args`: `LoopStructuralConfig.native_faults = options.native_faults == "on"`
    (import from `LoopStructural`). Record the value in the result: where
    `result["suite"] = suite` is set, the `suite` dict gets
    `"native_faults": options.native_faults`; pass it through the existing
    call chain to that function. Schema version stays 1.

- [ ] T020 [US4] [S] Run the fault evaluation tests with the native kernel in CI (FR-013, SC-006).
  - Files: `.github/workflows/native.yml`
  - Needs: T013
  - Wave: 5
  - Done when: `python -c "import yaml,sys; yaml.safe_load(open('.github/workflows/native.yml'))"` passes and the `wheel` job has the new lines.
  - Context: Job `wheel`, step "Install, build and test": after
    `python -m pytest packages/loop_interpolation/tests -q` add
    `python -m pip install -e .` and
    `python -m pytest tests/fault_evaluation -q`. Add to both `paths` lists:
    `'LoopStructural/modelling/features/_native_faults.py'` and
    `'LoopStructural/utils/_fault_block.py'`. Change nothing else.

- [ ] T021 [US2] [M] Measure SC-002, SC-003 and SC-007 and record the native differences (FR-005c).
  - Files: `specs/007-native-fault-kernel/verification.md`
  - Needs: T009, T016, T019, T020
  - Wave: 7
  - Done when: `verification.md` has a section "Stage S2" with a table for SC-002, SC-003 and SC-007 (target, measured, pass or fail) and the largest native difference for each golden model.
  - Context: Computer of spec 003 T030 (Intel i7-12700H, 20 threads,
    Linux); if another computer is used, say so. SC-002: for `f1_u0`,
    `f4_u0`, `f4_u3`, `cutting_faults`, `plugin` at 10 000 points, native
    restore with `native_threads = 1` takes at most 20 % of the Python
    restore time. SC-003: with all 20 threads, at most 3 %. SC-007: for 10
    million points, the extra peak memory of a native restore is within
    spec 003 FR-010 (`acceptance.MEMORY_ALLOWANCE_BYTES` in
    `benchmarks/fault_evaluation/acceptance.py`). Differences: from
    `pytest tests/fault_evaluation/test_native_faults.py -q --junitxml=...`
    (`max_native_diff`). Change no code. If a target fails, record it; do
    not change a limit.

- [ ] T022 [US2] [S] Record stage S2 in `ROADMAP.md` (FR-015).
  - Files: `ROADMAP.md`
  - Needs: T007, T021
  - Wave: 8
  - Done when: `ROADMAP.md` `## Status log` has a dated line for spec 007 stage S2 with the SC-002, SC-003 and SC-007 results.
  - Context: One line at the end of `## Status log`, as the lines above it:
    "spec 007 S2: native fault restoration in `loop_interpolation_core`
    (`FaultChain`, `FAULT_KERNEL_VERSION = 1`)", with the results from the
    "Stage S2" section of `specs/007-native-fault-kernel/verification.md`.

## Phase 5: Stage S3, native values and gradients (story 3)

- [ ] T023 [US3] [M] Implement the exact nearest-point search `KdTree` (OP5).
  - Files: `packages/loop_interpolation_core/src/fault/nearest.rs`, `packages/loop_interpolation_core/tests/fault_nearest.rs`
  - Needs: T008
  - Wave: 2
  - Done when: `cargo test --release --test fault_nearest` fails against the stub, then passes.
  - Context: `pub struct KdTree`, `KdTree::new(points: &[[f64; 3]]) -> KdTree`,
    `fn nearest(&self, p: [f64; 3]) -> usize` (index into `points`).
    Euclidean distance, exact. For equal distances, the lowest index. Points
    with a NaN coordinate are not in the tree (the caller gives only known
    points). Panic-free for 1 point; `new(&[])` is allowed and `nearest`
    then returns `usize::MAX`. Tests: 2000 random points (fixed seed, a
    small LCG in the test, no new crate) against a brute-force search for
    500 queries; equal distances give the lowest index; 1 point.

- [ ] T024 [US3] [M] Implement the native tetrahedron gradient of a faulted feature (`gradient.rs`, OP5, FR-004).
  - Files: `packages/loop_interpolation_core/src/fault/gradient.rs`, `packages/loop_interpolation_core/tests/fault_gradient.rs`
  - Needs: T012, T023
  - Wave: 4
  - Done when: `cargo test --release --test fault_gradient` passes, also with `--features parallel`.
  - Context: `feature_gradient(chain, field: &GridField, points, order, sides, restored, jumping: &[bool], scale: f64) -> (Vec<[f64;3]>, usize)`.
    Port the faulted branch of `GeologicalFeature.evaluate_gradient`
    (`LoopStructural/modelling/features/_geological_feature.py`):
    `regular_tetraherdron_for_points` and `gradient_from_tetrahedron`
    (`packages/loop_common/src/loop_common/math/_maths.py`); side check:
    for each fault with `jumping`, `frame_value(fault, 0, corner, sides,
    no faults restored)` of all corners of the batch; NaN values take the
    value of the nearest known corner of the same batch (`KdTree` on the
    restored corner points, as `fillnan='nearest'`); straddle = not all 4
    > 0 and not all 4 < 0. At most 10 rounds: halve the scale, rebuild only
    the straddling tetrahedra, check again. Restore the corners with
    `restore`, read `field`, NaN for points that still straddle. The count is
    the number of points refined at least once. The batch of the nearest
    fill is the whole set of straddling points of the round, as in Python.
    Tests: a feature with no fault crossing equals the grid gradient to
    1e-10; a point on a fault gives NaN; count > 0 near a fault.

- [ ] T025 [US3] [M] Add the binding `FaultChain.feature_gradient` (D2b).
  - Files: `packages/loop_interpolation_core/src/python.rs`
  - Needs: T013, T024
  - Wave: 5
  - Done when: `cargo clippy --all-targets --features python -- -D warnings` passes and, after `maturin develop --release -m packages/loop_interpolation_core/Cargo.toml`, `python -c "import loop_interpolation_core as c; c.FaultChain.feature_gradient"` passes.
  - Context: `feature_gradient(self, field: GridField, points: (N,3) f64,
    order: i64, sides: i8, restored: bool, jumping: bool, scale: float,
    threads: int) -> tuple[np.ndarray (N,3), int]`. Same rules as `restore`
    of T013: copy the points, `py.detach`, `ValueError` for bad shapes
    (`jumping`, `sides` and `restored` must have one item for each fault),
    `RuntimeError` for an internal error, empty input gives `(empty, 0)`.

- [ ] T026 [US3] [M] Test: native values, gradients and model ids agree with the Python path (story 3, FR-004, FR-005, FR-005a, FR-005b).
  - Files: `tests/fault_evaluation/test_native_field.py`
  - Needs: T016
  - Wave: 7
  - Done when: `pytest tests/fault_evaluation/test_native_field.py -q` fails (no `field_for`), then passes after T027.
  - Context: `pytest.importorskip("loop_interpolation_core")`. For the
    golden models of T017 and `onlap_domain`: `capture(model, points)`
    (`tests/fault_evaluation/golden.py`) with native on against native off;
    compare with `assert_matches(name, actual, expected, box_size)`.
    Model ids (FR-005b): equal, except at points where the Python value is
    within `1e-8 * box_size` of a unit boundary. Story 3 scenario 2: a point
    on a fault of `f1_u0` gives a NaN gradient with native on. The warning
    "points are too close to a fault" is logged once for the call
    (`caplog`). A fault frame with an unconformity region (`onlap_domain`
    or a model with `add_unconformity` before a fault) uses the native path:
    `_native_faults.chain_for(feature)` is not None. Record
    `record_property("max_native_diff", value)`.

- [ ] T027 [US3] [M] Add `field_for` and use the native field and gradient in `GeologicalFeature` (D3).
  - Files: `LoopStructural/modelling/features/_native_faults.py`, `LoopStructural/modelling/features/_geological_feature.py`
  - Needs: T025, T026
  - Wave: 8
  - Done when: `pytest tests/fault_evaluation -q` passes with the kernel built.
  - Context: Follow plan.md D3 exactly (read it). `field_for(feature) -> GridField | None`,
    cache `feature._native_field = (token, refs, field)`, checks as
    `_shared_frame_members` for one interpolator. `evaluate_value`:
    `field.values(evaluation_points[mask, :], native_threads())` after
    `interpolator._check_frame()` and `interpolator.update()`.
    `evaluate_gradient` without faults: `field.gradients(...)`. With faults:
    if `field_for` and `chain_for` are not None, call
    `chain.feature_gradient(field, pos, order, sides, restored, jumping, element_scale_parameter, native_threads())`;
    rows outside the mask are NaN; log the Python warning once if the count
    > 0. Else keep the Python loop. Do not change the signatures (FR-007).

- [ ] T028 [US3] [M] Measure SC-004 and SC-005 with the native path.
  - Files: `specs/007-native-fault-kernel/verification.md`
  - Needs: T021, T027
  - Wave: 9
  - Done when: `verification.md` has a section "Stage S3" with SC-004 (measured time, target 0.5 s) and SC-005 (the golden test result with the native path on).
  - Context: SC-004: run
    `python -m benchmarks.fault_evaluation.many_faults --output /tmp/many_faults_native.json --points 5000`
    with the kernel built and the native path on (no
    `LOOPSTRUCTURAL_NATIVE_FAULTS`), on the T030 computer.
    `cases["5000"]["evaluate_model_seconds"]` at most 0.5 s (baseline 2.09 s
    in `/home/lgrose/dev/LoopStructural-wt/_runs007/many_faults_python.json`).
    Do not use `--compare` (native results are not bit-identical, C-001). SC-005: run
    `pytest tests/fault_evaluation/test_golden.py -q` with the native path on
    (no `LOOPSTRUCTURAL_NATIVE_FAULTS`), and record the result and the
    largest difference. Change no code; if a target fails, record it.

- [ ] T029 [US3] [S] Record stage S3 in `ROADMAP.md` (FR-015).
  - Files: `ROADMAP.md`
  - Needs: T022, T028
  - Wave: 10
  - Done when: `ROADMAP.md` `## Status log` has a dated line for spec 007 stage S3 with the SC-004 and SC-005 results.
  - Context: One line at the end of `## Status log`: "spec 007 S3: native
    field values and gradients (`GridField`, `FaultChain.feature_gradient`)",
    with the results from the "Stage S3" section of
    `specs/007-native-fault-kernel/verification.md`.

## Phase 6: Stage S4, prebuilt packages (story 4)

- [ ] T030 [US4] [S] Add a fault kernel smoke test to the wheel workflow (FR-014).
  - Files: `.github/workflows/core-wheels.yml`
  - Needs: T013
  - Wave: 5
  - Done when: `python -c "import yaml; yaml.safe_load(open('.github/workflows/core-wheels.yml'))"` passes and the smoke script has the fault check.
  - Context: In the step "Install the wheel and solve a small problem",
    add to the Python heredoc after the solve check: make
    `core.FaultChain(spec)` with one fault: `displacement` 10.0, `steps` 10,
    three frame grids of 21x21x21 nodes, origin 0, step 1, no local frame,
    node values x - 10 (member 0), y (member 1), z (member 2); `cutters` [],
    `regions` [], profile `{"scale": 1.0, "gx": {"kind": "ones"}, "gy": {"kind": "ones"}, "gz": {"kind": "ones"}}`,
    block settings of plan D2b. Restore the points (5, 5, 5) and
    (12, 4, 8) with `order=[0]`, `sides=[2]`, `restored=[False]`,
    `threads=0`; expected (5, 15, 5) and (12, 14, 8) (the step direction is
    the gradient of member 1, +y), to 1e-12. Change nothing else.

- [ ] T031 [US4] [S] Describe the fault kernel in the package README.
  - Files: `packages/loop_interpolation_core/README.md`
  - Needs: T025
  - Wave: 6
  - Done when: the README has a section "Fault kernel" and `git diff --stat` shows only that file.
  - Context: 10 to 20 lines, ASD-STE100 style: what the kernel does
    (restores points through the faults of a LoopStructural feature, and
    evaluates values and gradients); it is used by LoopStructural
    automatically; `FAULT_KERNEL_VERSION = 1`; the classes `FaultChain`
    (`restore`, `frame_value`, `feature_gradient`) and `GridField`
    (`values`, `gradients`); the switch `LoopStructuralConfig.native_faults`
    and `LOOPSTRUCTURAL_NATIVE_FAULTS=0`; the thread count
    `LoopStructuralConfig.native_threads` (None = all cores). Do not
    document the spec dict in detail; point to
    `specs/007-native-fault-kernel/plan.md` D2b.

- [ ] T032 [US4] [M] Check the prebuilt packages on each platform (SC-008).
  - Files: `specs/007-native-fault-kernel/verification.md`
  - Needs: T028, T030
  - Wave: 10
  - Done when: `verification.md` has a section "Stage S4" with one row for each platform of FR-014 (Linux x86_64, Linux aarch64, macOS arm64, macOS x86_64, Windows x64): wheel built, installed on a clean environment, golden tests with the native path on: pass or fail or "not checked".
  - Context: Use the artifacts of the `Native wheels` workflow
    (`.github/workflows/core-wheels.yml`; run it with `gh workflow run` only
    if the maintainer agrees, else use the last run on the branch). CI
    installs and smoke tests 4 platforms; macOS x86_64 is cross-compiled
    and not run in CI: mark it "maintainer check needed" unless a result
    exists. For the golden tests: `pip install loopstructural` from the
    branch plus the wheel, then `pytest tests/fault_evaluation/test_golden.py -q`.
    Change no code.

## Phase 7: Finish

- [ ] T033 [Finish] [M] Run the full test suite with the native path on, off and not installed, and record stage S4 (SC-006, FR-015).
  - Files: `specs/007-native-fault-kernel/verification.md`, `ROADMAP.md`
  - Needs: T006, T029, T031, T032
  - Wave: 11
  - Done when: `verification.md` has a section "SC-006" with the result of the three `pytest` runs, and `ROADMAP.md` `## Status log` has a dated line for spec 007 stage S4.
  - Context: Three runs from the repository root: (1) kernel built
    (`maturin develop --release -m packages/loop_interpolation_core/Cargo.toml`),
    `pytest -q`; (2) the same with `LOOPSTRUCTURAL_NATIVE_FAULTS=0 pytest -q`;
    (3) `pip uninstall -y loop-interpolation-core`, then `pytest -q`.
    Also run `ruff check LoopStructural packages`. Record pass counts and
    any failure with its test id. ROADMAP line: "spec 007 S4: prebuilt
    native packages with the fault kernel" with the SC-008 rows of
    `verification.md`. Change no code.

## Wave table

| Wave | Tasks | Tier [S] | Tier [M] |
|---|---|---|---|
| 1 | T001, T002, T003, T008, T009, T014 | 2 | 4 |
| 2 | T004, T010, T011, T015, T023 | 0 | 5 |
| 3 | T005, T012, T019 | 1 | 2 |
| 4 | T006, T013, T024 | 0 | 3 |
| 5 | T007, T017, T018, T020, T025, T030 | 3 | 3 |
| 6 | T016, T031 | 1 | 1 |
| 7 | T021, T026 | 0 | 2 |
| 8 | T022, T027 | 1 | 1 |
| 9 | T028 | 0 | 1 |
| 10 | T029, T032 | 1 | 1 |
| 11 | T033 | 0 | 1 |

## Critical path

T008 -> T010 -> T012 -> T013 -> T017 -> T016 -> T026 -> T027 -> T028 -> T029 -> T033

## Dependencies

- S1 (T002 to T007): T003 test, then T004 code. T002 and T005 measure before
  and after (the 94-fault benchmark `many_faults.py`); T006 runs D1b
  only if T005 reports that SC-001 fails.
- Rust (T008 to T013): T008 makes the skeleton. The grid (T010), the
  profiles (T011) and the nearest search (T023) run in parallel. T012 needs
  T010 and T011. T013 (PyO3) needs T012.
- Python S2: T014 test, then T015 settings (no kernel needed). T017 and
  T018 tests need T013 and T015. T016 makes them pass and needs T004.
  T019 and T020 select the path in benchmarks and CI. T021 measures.
- S3: T024 needs T012 and T023; T025 binds it; T026 test needs T016; T027
  makes it pass. T028 measures.
- S4: T030 smoke test needs T013; T031 README needs T025; T032 checks the
  platforms. T033 runs the full suite three ways and closes the spec.
- `ROADMAP.md` is a chain: T001, T007, T022, T029, T033.
  `verification.md` is a chain: T009, T021, T028, T032, T033.
