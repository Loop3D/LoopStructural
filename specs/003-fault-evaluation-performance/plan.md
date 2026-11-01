# Implementation plan: Faster evaluation of faulted models

**Feature**: `003-fault-evaluation-performance` | **Date**: 2026-11-01
**Spec**: [spec.md](spec.md) | **State**: Planned; decision gates remain open

## Open points

| Id | Question and evidence | Resolution rule | Tasks that wait |
|---|---|---|---|
| OP1 | What point count and call count apply to gradient and nested-fault work in FR-008 and FR-010? A gradient restores 4N corners, while its region masks evaluate N original points. Different features can refine to different tetrahedron sizes. A fault frame can itself apply faults at each integration step. Thus F total calls for all these different arrays is not established by sharing restored N-point arrays. | Ask the user before gradient changes. Proposed interpretation: one restoration per distinct ordered fault set AND distinct point batch, counting recursive frame batches separately; memory is 3 floats per actual restored point, excluding unchanged baseline work arrays. This is a proposal, not an accepted change. If the user retains the literal bound, use spec-reviewer/spec-planner to design a compatible algorithm before this task runs. Record the answer in Answers; use spec-cascade for changed requirements. | T021, T030, T031, T032 |
| OP2 | Can the complete benchmark matrix satisfy SC-008? There are at least 240 synthetic cases (16 models x 5 functions x 3 sizes), before special cases and repetitions. Research reports 8.65 s for all groups at only 5,000 points for one model. | T011 measures the full legacy run with one timed round per case after warmup. T030 measures the final suite. Never remove cases or reduce point counts to claim success. If either final runtime or run-to-run noise fails, report the failure and ask for an explicit spec decision. Timing success is not assumed. | T031, T032, if measurement fails |

## Summary

First freeze the current geology with stored solved fields, output references,
and a reproducible benchmark baseline. Next share restored points within an
evaluation and between calls. Last, reduce repeated structured-grid lookup
work inside one fault application. Do not change the solver, displacement law,
fault order, integration step count, or surface extraction algorithm.

## Technical context

- Python library, minimum Python 3.9. NumPy, SciPy, pandas and scikit-image are
  existing runtime dependencies. Add pytest-benchmark to the development extra
  only. No Rust work or Cargo command is required for this spec.
- `BaseFeature._apply_faults` in
  `LoopStructural/modelling/features/_base_geological_feature.py` loops over
  faults in order. `_calculate_mask` calls regions on original coordinates.
- `GeologicalFeature.evaluate_value` and `evaluate_gradient` in
  `LoopStructural/modelling/features/_geological_feature.py` build first.
  Gradients restore tetrahedra and refine them near faults.
- `UnconformityFeature` shares its parent's faults and builder, but has no
  regions. Its callable evaluates values on the original points.
- `GeologicalModel.evaluate_model` and `evaluate_model_gradient` in
  `LoopStructural/modelling/core/geological_model.py` copy input coordinates
  and loop over stratigraphic groups. Keep the unit priority and world/local
  coordinate handling unchanged.
- `FaultSegment.apply_to_points` already reuses initial coordinates in the
  first step and accepts `precomputed_gx`. Do not count this existing saving
  as new work. `use_threads` stays off.
- Builders expose `_up_to_date`, `update`, `update_build_arguments`,
  `set_not_up_to_date` and `add_data_from_data_frame`. Solved arrays are `c`
  on the discrete interpolator. Regions and fault lists are mutable.
- `StructuredGrid.evaluate_value`, `evaluate_gradient`,
  `position_to_local_coordinates` and `get_element_gradient_for_location`
  repeat cell lookup work in
  `packages/loop_common/src/loop_common/supports/_3d_structured_grid.py`.
- Read `research.md`, `measure_baseline.py`, and
  `data/plugin_model/plugin_model.py`. The plugin loader sets stored grids,
  solved fields, builder state and fault lists. It preserves world coordinates.
- Existing regression suites include `tests/unit/modelling/test_faults_segment.py`,
  `test_unconformity_feature.py`, `test_geological_feature.py`,
  `test_fault_cycle_detection.py`, `test_public_api_contract.py`, and the
  loop_common and loop_interpolation package tests.

## Planning inputs and worktree

The requested templates, `specs/README.md`, and
`docs/development/coding-tools.md` do not exist. Use the 001 and 004 plan/task
structure as the local precedent. Do not create an invented spec index or
provider policy. Tier S and M describe task size; the implementer must state
which runtime models it can use before workers start.

Use the existing worktree `C:/Users/lachl/dev/LoopStructural-native`, branch
`spec/003-fault-evaluation-performance`. It was clean at planning time. Do not
switch the main checkout or replace this existing branch with another branch.
The spec text is identical in the main checkout and this branch.

The approved constitution and toolkit update currently belong to the main
checkout `C:/Users/lachl/dev/LoopStructural`. Before implementation, carry those
reviewed setup changes to this worktree and run the strict doctor there. Do not
write production code until this prerequisite passes. The source CLI can be
run with `PYTHONPATH=C:/Users/lachl/dev/dev-tools/src` and
`python -m dev_tools doctor --require-constitution`; this is a local workaround
for the broken installed uv launcher, not a new committed project command.

No dependency on specs 001, 002 or 004 is required. Solvers stay unchanged.
Spec 004 also touches folded-feature tests; coordinate file changes if it runs
at the same time. This plan uses its own golden fixture directory.

## Constitution check

| Rule | Result | How this plan agrees |
|---|---|---|
| Python geological modelling and workspace packages | Pass | Changes stay in the modelling library, support package, tests and developer tools. |
| Supported test, lint and format commands | Pass | Use pytest, Ruff and Black. Do not run unrelated Cargo checks. |
| pytest and Ruff before merge | Pass | T032 runs both; do not waive existing failures. Report baseline failures separately. |
| NumPy docstrings, type hints, Python 3.9 lint target | Pass | New interfaces use compatible typing and local style; no new runtime dependency. |
| ROADMAP.md is the architecture source | Pass | T031 records the completed evaluation stage and its limits. |
| Independently testable, reversible stages | Pass | Golden tests, benchmark tools, call sharing, persistent cache and grid changes have separate checkpoints. |
| GeologicalModel public API boundary | Pass | Preserve all existing signatures and results; cache controls are additive. Run API contract tests. |
| QGIS plugin stays separate | Pass | Use the stored plugin model without QGIS; do not edit the plugin. |
| Spec, plan and task documents | Pass | This plan and tasks accompany the clarified spec. OP1 changes need a spec decision. |
| Scoped package releases | Pass | No new package; commit root and loop_common changes with their correct scopes. |
| Amendments need maintainer approval | Pass | Constitution accepted in chat on 2026-11-01; no further policy is introduced here. |

## Design

### D1. Golden models and immutable references

Create `tests/fault_evaluation/models.py` with
`model_names() -> list[str]`, `build_model(name: str) -> GeologicalModel`, and
`sample_points(name: str, count: int, seed: int = 1) -> np.ndarray`.
The 16 synthetic names are `f{F}_u{U}` for F in 0,1,2,4 and U in 0,1,2,3.
Use the data recipe of measure_baseline.py without importing it: that script
patches FaultSegment globally and runs measurements at import time. Identify
the oldest group from creation order and assert its identity in the model test.
Add names `stopped_fault`, `cutting_faults`, `onlap_domain`, `faulted_fold`,
and `plugin`. Each special recipe must construct the named relationships and
assert them. Include a stratigraphic column for every applicable model.

Put fixed recipes in `tests/fault_evaluation/data/recipes.json`, and frozen
fields/grids in `fields.npz`. Keep the supplied plugin files in place and load
them explicitly by path. Do not change those input fields during optimization.
Each synthetic stored field includes support type, grid geometry, local
transform, coefficient array and builder state. The folded case records its
fold frame and folded feature fields. Never solve in a normal golden test.

`golden.py` defines `capture(model, points) -> dict[str, np.ndarray]` and
`assert_matches(name, actual, expected, box_size) -> None`.
Store each model's fixed test coordinates, model ids, feature values and
gradients, restored points for every fault set, and surface summaries in
`tests/fault_evaluation/data/<name>.npz`. Use 128 fixed interior points plus
explicit points near fault surfaces, unconformities, unit boundaries and box
faces. Add empty, NaN and outside points in edge tests. Use a fixed 20x20x20
surface sampling box and fixed levels per recipe. Record vertex count, bounds
and triangle area for each surface. An unexpected empty surface is a failure.

Use `1e-8 * max(box maximum - box origin)` as absolute value/gradient/point
tolerance; compare NaN masks exactly. Unit ids must agree except points within
that value tolerance of a declared boundary; check the scalar boundary distance,
not the integer id distance. Surface vertex counts agree exactly; bounds and
area use rtol=1e-6 with an explicit zero handling rule (zero must stay zero).
Failure messages name model, output, maximum difference and failing point index.

`python -m tests.fault_evaluation.regenerate --write` is the only write command.
Without --write it validates and exits without changing a file. Normal tests
must never update references. `manifest.json` records input and output SHA256,
baseline source revision, field-generation revision, dependency versions,
coordinate transforms, seeds and levels. Record the actual legacy worktree
revision, not the historical research revision unless they are the same.
Binary output hashes and readable summary diffs make regeneration reviewable.
Capture before changing evaluation code. Test sensitivity by temporarily
perturbing fault movement by 1% in a scoped monkeypatch; leave no patch active.

### D2. Benchmark tools and baseline gate

Create a developer-only package `benchmarks/fault_evaluation` with `run.py`,
`compare.py`, `instrumentation.py` and `acceptance.py`. Add package __init__.py
files. pytest's existing testpaths=[tests] keeps benchmark tests out of CI.
The suite file is `benchmarks/fault_evaluation/bench_cases.py`; it is selected
explicitly by the runner. Never add an autouse runtime counter or CI job.

`python -m benchmarks.fault_evaluation.run --output PATH --rounds 1` calls
pytest-benchmark with each of the 21 models, five operations (oldest value,
all group values, model ids, oldest gradient, model gradient), and 5,000,
50,000 and 500,000 points. Add 20,000-point plugin and single-fault cases,
5,000-point SC-002/003 cases, and warm second calls. Reuse stored fields;
exclude model construction from timed evaluation. Time cold calls after
clearing the cache, then warm calls separately. Warmup must not turn a cold
case into a cache hit. Include the cost of hashing and invalidation in timing.

Run time, counters and peak RSS in separate fresh-process passes so wrappers
and memory tracing do not distort timings. `instrumentation.py` patches
FaultSegment.apply_to_points and the selected interpolator field evaluation
methods only inside its context manager, then restores the original methods
even on failure. Record total and top-level counts and point counts separately;
total count is the spec metric until OP1 is resolved. Peak RSS is process high
water minus the built-model idle RSS, sampled consistently with platform units
recorded; do not use tracemalloc alone to claim NumPy peak memory coverage.
Use native process APIs via stdlib ctypes/resource, not a runtime dependency.

Result schema v1 has case_id, model, operation, points, cold/warm state, rounds,
time samples and median seconds, fault calls, field calls, processed point
counts, peak extra RSS bytes, code revision and dirty state, input hash,
Python/NumPy/SciPy/LoopStructural versions, OS, CPU and RAM. Reject duplicate
ids and incompatible or missing cases when comparing. Preserve pytest-benchmark
data in the JSON. `compare BASE NEW --output PATH` writes ratios, noise marks
above 1.10, fault/unconformity scaling tables and provenance differences.
Report undefined ratios for zero baselines explicitly. Do not silently skip
cases or compare warm new results with cold baseline results.

T011 writes baseline.json and baseline_repeat.json plus a readable baseline
report. It records complete runtime and noise. T011 is a hard gate: all golden
checks must pass and all baseline cases must exist before any evaluation edit.
If the run fails, capture the exact failing case; do not substitute research
numbers for a machine-readable measured baseline. SC-007/008 are measured
requirements, not promised outcomes of planning.

### D3. Evaluation scope and restoration sharing

Add private `_evaluation.py` with `EvaluationContext` and
`evaluation_scope(model=None)` as a context manager backed by ContextVar.
Existing public signatures stay unchanged. A scope owns original point tokens,
ordered fault-set keys and a dictionary of immutable restored arrays. Nested
calls reuse the active scope; a top-level call resets it in finally. A scope
for a different model is separate. Standalone features use a temporary scope.
Wrap feature value/gradient and model evaluation entry points, so regions,
fault frames and fold frames can participate without passing new arguments.

`restore(feature, points, reverse=False) -> np.ndarray` is the sole private
cache entry. Its key includes exact point content/shape/dtype, ordered fault
identities, dependency state, reverse, and faults_enabled. Map repeated original
arrays to a scope token; do not depend on object identity between public calls.
Treat reordered rows, sliced arrays, mutated arrays and tetrahedron arrays as
different point batches unless their exact content is equal. Read original
points for masks; never feed restored points into region predicates. Preserve
fault order including reverse behavior as it exists today; do not fix geology
bugs as part of the performance change. Do not cache partial fault prefixes.

Returned restored arrays are read-only internally and never exposed for caller
mutation. Audit each consumer; copy only where a consumer writes. An exception
must not publish a partial restoration. Build all required reachable builders
before the first cache lookup, preserving first-use lazy building. If a build
or solve occurs during a scope, discard its state and restart lookup with the
new revision. Keep cycle checks in the faults setter.

### D4. Persistent cache and invalidation

`_evaluation_cache.py` defines `RestorationCache` with `enabled: bool`,
`max_bytes: int`, `clear() -> None`, `invalidate() -> None`,
`get(key) -> np.ndarray | None` and `put(key, restored) -> None`.
Expose it as additive `GeologicalModel.evaluation_cache`, with enabled=True
and max_bytes=64*1024*1024. max_bytes=0 or enabled=False stops retention;
negative limits raise ValueError; shrinking the limit evicts immediately.
Use OrderedDict LRU and account retained ndarray.nbytes exactly. Entries larger
than the limit are not retained. Do not retain a second copy in the call scope.
Document metadata and bounded hashing workspace separately from stored arrays.

Use a model RLock for lazy builds, cache publication, invalidation and full
top-level evaluation. This first version serializes evaluations of one model;
different models can run concurrently. Nested evaluations reuse the lock.
Do not share a ContextVar scope across threads. Keep fault worker threads off.

Use a monotonic model revision for known mutation paths: BaseBuilder dirty,
build/update and argument changes; GeologicalFeatureBuilder data/geometry/build;
DiscreteInterpolator solve/update; feature faults, regions, faults_enabled;
fault displacement, faultfunction and steps; model data, column and feature
relationships. Direct assignment must retain its previous behavior.
Add property setters or invalidation calls at these existing boundaries.

Revision alone is insufficient: public lists, dictionaries and arrays are
mutable. Before a persistent lookup compute a streaming SHA256 dependency
fingerprint, including ordered reachable faults/frame features/parents,
regions and region thresholds/sign, builder data/build arguments/dirty state,
field coefficients and support transforms/geometry, displacement/function/steps,
and faults_enabled. Use cycle-safe traversal and bounded chunks, without a full
array copy. Detect input point mutation by content digest on every outer call.
Use hashes only for keys; do not retain input arrays. Full fingerprint and
hashing cost is included in warm-call measurements.

If any dependency is opaque (a user callable with mutable global/closure state,
an unsupported support or custom interpolator), bypass persistent reuse for
that reachable fault set. Do not guess that an object id proves it unchanged.
Call sharing may remain only where the dependency is known stable within that
call. A custom impure region must retain its legacy evaluation behavior.
Tests include direct/in-place edits of points, c, lists, build arguments and
region state, in addition to ordinary setters. Changed-then-restored state is
safe because output depends on the full current state. Clear old entries on
revision change to avoid retained obsolete data. Never cache masks or geology
outputs: only restored coordinates.

### D5. Gradient decision and thread/memory checks

Preserve tetrahedron construction, near-fault refinement, element scale and
gradient reconstruction. Share equal point batches across groups; never reuse
an N-point restoration for a 4N-point tetrahedron. OP1 blocks implementation
of a claimed total F-call gradient bound. The gradient tests must show counts
for original points, corners, refinement and nested fault-frame calls, so the
decision has concrete evidence. Do not mark FR-008 satisfied from top-level
value counts alone. If OP1 changes the spec, revise this plan before coding.

Use small deterministic unit tests for retained-byte accounting and separate
local subprocess acceptance cases for 10 million points and peak RSS. Exclude
persistent cache retention when measuring call-only extra memory; measure both
modes. Record original input size, restored batch sizes, distinct sets, retained
bytes and process peak. Do not claim FR-010 from dictionary length alone.
Concurrent tests compare two same-model calls with serial references, first
build happens once, independent models remain independent, and an exception
cleans scope state and releases locks.

### D6. Cheaper structured-grid evaluation

After the baseline gate, add private StructuredGrid
`_lookup_points(pos) -> (corners, inside, local_coords)` and
`_gradient_weights(local_coords) -> np.ndarray`. Compute the cell index once,
derive corners and local coordinates from it, and preserve corner order,
boundary clipping, outside values and NaN behavior. Route evaluate_value and
evaluate_gradient through the common lookup. Preserve summation order and the
existing scale factors. Test the current implementation as the reference before
changing these functions. Do not change other support types.

Add a private FaultSegment `_evaluate_frame(points, need_gradient=False)` fast
path only when frame coordinates use StructuredGrid supports with equal geometry
and transforms, no child faults, and equivalent region masks. Reuse the lookup
for all three values and coordinate-1 gradient. Keep public method signatures,
precomputed_gx, initial mask, all 10 steps and their arithmetic order. Fall back
to the original feature calls for folded frames, different supports/transforms,
child faults, fillnan behavior or incompatible regions. No cached step lookup
survives a point movement. Test fallback as well as the fast path. Measure the
SC-009 target; if it fails, report the measured shortfall and replan, rather than
changing steps, tolerances or displacements to obtain a passing result.

## Test strategy and checkpoints

1. Golden harness tests fail first because fixture APIs are absent, then pass.
   The geology golden checks themselves must pass on unchanged evaluation code.
   A 1% fault perturbation must fail; test NaN masks and failure messages too.
2. All legacy golden cases run under 120 seconds; record actual elapsed time.
3. Benchmark tools have small stub models for fast runner/schema/counter tests.
   Their test files live under tests; they do not run the full benchmarks in CI.
4. Freeze complete baseline and repeat results before T012-T027 production work.
5. Add failing tests before each context/cache/invalidation/gradient/grid change.
   The orchestrator observes failure; workers do not run competing build commands.
6. At each implementation wave run relevant pytest filters once, plus Ruff on
   changed production files. Existing unrelated failures are reported, not hidden.
7. Final: Black check, Ruff, root pytest, loop_common pytest and
   loop_interpolation pytest. Run API contract tests and all golden cases.
8. Local acceptance runs compare complete before/after results for SC-002 through
   SC-010, and measure 10-million-point memory. These runs are not CI tests.
   A missed performance target keeps acceptance open.

## Requirement coverage

| Requirement | Tasks |
|---|---|
| FR-001, FR-001a | T001, T002, T005, T006 |
| FR-002, FR-003, FR-004 | T003, T004, T005, T006 |
| FR-005, SC-001 | T006, T030, T032 |
| FR-006 | T006, T015, T021, T027, T032 |
| FR-007 | T014, T015, T018, T019 |
| FR-008 | T014, T015, T020, T021; OP1 unresolved for gradients/nested calls |
| FR-009 | T016, T017, T018, T019, T022, T023 |
| FR-010 | T012, T013, T018, T019, T022, T023, T030; OP1 for gradient budget |
| FR-010a | T022, T023 |
| FR-011 | T014, T015, T019, T032 |
| FR-012, FR-013, FR-014, FR-015 | T007, T008, T009, T010, T011 |
| FR-016 | T011 |
| SC-002, SC-003, SC-004, SC-005 | T028, T029, T030 |
| SC-006 | T022, T023, T028, T029, T030 |
| SC-007, SC-008 | T011, T028, T029, T030; OP2 measurement gate |
| SC-009 | T024, T025, T026, T027, T030 |
| SC-010 | T016, T017, T018, T019, T028, T029, T030 |

## Handover and release

No protocol, schema, generated solver fixtures or Rust package changes.
The loop_common change is a separate package release concern. Run its tests
and describe it in that package's changelog workflow. Keep all existing plugin
imports and model evaluation signatures. Baseline and golden regeneration are
explicit developer commands. Do not edit stored references after optimization
to hide a regression. Do not push, commit or open a PR during planning.

## Answers

- 2026-11-01: The user approved the constitution draft in the main checkout.
  The draft marker was removed and the stricter pytest/Ruff merge rule activated.
- OP1, 2026-11-01: The user accepted a limit per distinct ordered fault set and distinct point batch. Count recursive frame batches separately. Measure memory per actual restored point, including tetrahedron corners, excluding unchanged baseline work arrays. This resolves the interpretation for T021 and dependent tasks; the measured targets still apply.
- OP2: Pending measurements. The 15-minute target and 10% noise limit remain.
