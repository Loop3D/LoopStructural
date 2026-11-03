# Implementation plan: Faster evaluation of faulted models

**Feature**: `003-fault-evaluation-performance` | **Date**: 2026-11-01
**Spec**: [spec.md](spec.md) | **State**: Replanned 2026-11-03 for the fault-set cache key (D4a, T040-T043) after the failed T030; T001-T029 and T033-T039 are done; T030 is blocked

## Open points

| Id | Question and evidence | Resolution rule | Tasks that wait |
|---|---|---|---|
| OP1 (resolved) | Gradient and nested-fault count/memory interpretation. | The user accepted one restoration per distinct ordered fault set and distinct point batch on 2026-11-01. Count recursive frame batches separately and memory per actual restored point, including corners. See Answers. | none |
| OP2 | Can the standard suite meet its time budget (about 1 hour before this feature, SC-008 15 minutes after it) and the SC-007 noise limit? The first run of the 860-case full matrix did not finish in 4 hours and was stopped on 2026-11-02. 500,000-point cases were 90 % of the time measured. | The spec now defines a standard suite (C-006 to C-010). T011 measures the standard suite twice before any performance edit. T030 measures it again after. If the baseline run of the standard suite takes more than 90 minutes, or if SC-007 or SC-008 fail, report the measured result and ask for a spec decision. Never remove cases from the standard suite to claim success. | T031, T032, if measurement fails |
| OP3 | SC-008 (15 minutes for the standard suite). T030 measured 3,017 s (baseline 3,271 s). The sum of the candidate medians is only 334 s. About 2,000 s is work that is not evaluation: counters pass, memory pass, warm priming and process starts. Faster faults cannot reach 900 s. | Needs a spec decision. Options: (a) change the SC-008 target to a measured value; (b) cut the fixed cost of the passes (for example, fewer process starts); (c) make the standard suite smaller. (b) and (c) change FR-012a. T043 splits the wall time by pass, so the user decides with numbers. Recommendation: (b), then (a) for the part that stays. Never remove cases to claim success. | T030, T031 |
| OP4 | SC-006 (peak memory at most 120 % of the baseline). T030 shows many cases above 1.2. Several have a very small baseline memory, so a few MB of cache or hash workspace gives a large ratio. | T043 lists each failing case with the baseline bytes, the candidate bytes, and the cause (cache retention, hash workspace, or other). The user then decides whether the rule needs an absolute allowance for small baselines. This is a spec change. Do not change the 1.20 limit in code without that answer. | T030 |
| OP5 | Can SC-010 (warm `evaluate_model` at most 10 % of the cold call) be met? T030 diagnosis: a warm call with no fault calls and no tokens takes about 0.017 s. The cold call is about 0.22 s after the fix of D4a. The limit is then about 0.022 s, so the key work may use about 5 ms. | T042 measures the key time and the warm and cold times. If the ratio is above 0.10 after D4a, report the numbers and ask the user. Do not cache the key per revision: that hides in-place edits (D4). | T030 |

## Summary

First test and repair the existing gradient output-array defect in a separate
correctness commit (T033, T034). Then freeze solved fields, golden outputs and
the benchmark baseline before any performance change. Next share restored
points within an evaluation and between calls. Last, reduce repeated
structured-grid lookup work inside one fault application. Do not change the solver, displacement law,
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
provider policy. The existing `.agents/model-policy.json` maps Codex simple
workers to `gpt-6-luna` at low effort and standard workers to `gpt-6-luna` at
medium effort. Its worker limit is 5; this runtime has 3 available worker slots
besides the main agent. Resolve runtime availability again before starting work.

Use the existing worktree `C:/Users/lachl/dev/LoopStructural-native`, branch
`spec/003-fault-evaluation-performance`. It was clean at planning time. Do not
switch the main checkout or replace this existing branch with another branch.
The spec text is identical in the main checkout and this branch.

The approved constitution and toolkit setup were carried into this worktree
in commit `d22db4f1`; strict doctor passes. T001 was completed in `c5775418`,
with a correction for actual fold-frame coordinate dependencies in `667ca332`.
T002 has uncommitted fixture helper drafts but no generated fields archive.
Preserve those files; do not treat T002 as complete or repeat T001's red check.
The new test/fix tasks run before T002 resumes. The source CLI can be
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

### D0. Repair the gradient output before baseline capture

The first T002 field-export run failed in
`GeologicalFeature.evaluate_gradient` at the assignment to `v[mask, :]`.
The method first allocates `v` as N by 3. Its tetrahedron refinement loop then
assigns N by 4 fault-coordinate values to the same variable. Reconstruction
produces N by 3 gradients and raises a shape mismatch. A two-point linear
field with an identity fault reproduces this without a solver: shape (2,3)
cannot be assigned into shape (2,4). The failed exporter produced no fields.npz.

T033 writes the standalone regression in
`tests/unit/modelling/test_faulted_gradient_output.py`. It must not import
the unfinished golden fixture package or require saved fields. Use a local
BaseFeature test double that implements evaluate_value, evaluate_gradient,
copy and get_data, supplies coordinate 0 through __getitem__, and returns
unchanged points from apply_to_points. Use a field double with
`evaluate_value(points) = points @ [2.0, -3.0, 0.5] + 7.0`. At two finite
points away from a fault boundary the expected result has shape (2,3) and
each row [2.0, -3.0, 0.5], with rtol=0 and atol=1e-12. Test one and two
identity faults, one masked point remaining NaN, and one controlled refinement
iteration (mixed coordinate-0 corner signs on the first call, positive on
the next). Assert the loop refines once and returns the same gradient.
These cases must fail first with the existing output shape error.

T034 changes only the temporary name inside the refinement loop:
`v = f[0].evaluate_value(...).reshape(...)` becomes `fault_values = ...`,
and the two sign tests in `flag` use fault_values. Keep `v` as the original
gradient output. Do not change the tetrahedra, refinement condition, scale,
fault order, restored coordinates, solver, integration steps, masks or
public signature. No cache or performance work belongs in this task.
Check the diff and commit this correctness fix separately before any capture.

The pre-performance baseline is the code after this correctness fix and
before caching/grid optimization. Record three immutable revisions in the
golden manifest and both benchmark result files:
`prerequisite_parent_revision`, `prerequisite_fix_revision`, and
`baseline_revision`. Record the original exception and the regression results
in `gradient-prerequisite.md`; the fix revision is the commit that contains
only the reviewed fix and its report, after the failing test commit.
The raw parent has no valid affected-gradient result. Do not claim numerical
parity for that absent output, replace it with zeros, omit the case, or use
an optimized run as the baseline. Where supplied pre-fix results do exist
(the plugin values and restored coordinates), they must still match.
Golden and timing results must say that this prerequisite was applied.
This is a correctness prerequisite to FR-002/FR-006 and US4, not a claimed
speed improvement. Keep the faulted-gradient regression in all later checks.

### D1. Golden models and immutable references

Create `tests/fault_evaluation/models.py` with
`model_names() -> list[str]`, `build_model(name: str) -> GeologicalModel`, and
`sample_points(name: str, count: int, seed: int = 1) -> np.ndarray`.
The 16 synthetic names are `f{F}_u{U}` for F in 0,1,2,4 and U in 0,1,2,3.
Coordinate contract and reference epochs (spec 005, 2026-11-03): public
model inputs, fixed query points, restored points and public evaluation outputs
are in world coordinates. Arrays at a support or solver boundary are in that
support's computational coordinates; convert once at the boundary. The mapping
is the support transform. Spec 005 owns repairs to frame construction and
world/local conversion. This spec fingerprints the current restoration
transform to prevent stale cache use. Spec 005's coordinate owner rebuilds
dependent local state after an accepted transform change or rejects the change
before mutation. Existing golden and benchmark manifests remain valid
historical records. If a future reviewed correction changes a reference, append
its correctness revision and affected-reference record; do not rewrite prior
provenance. Performance comparisons require matching cases, run state,
correctness reference epoch and equivalent machine/environment provenance.
Report correctness changes separately from measured performance. T004's
existing capture remains valid for its recorded reference epoch.

Use the data recipe of measure_baseline.py without importing it: that script
patches FaultSegment globally and runs measurements at import time. Identify
the oldest group from creation order and assert its identity in the model test.
Add names `stopped_fault`, `cutting_faults`, `onlap_domain`, `faulted_fold`,
and `plugin`. Each special recipe must construct the named relationships and
assert them. Include a stratigraphic column for every applicable model.

Put fixed recipes in `tests/fault_evaluation/data/recipes.json`, and frozen
fields/grids in `fields.npz`. Keep the supplied plugin files in place and load
them explicitly by path. Do not change those input fields during optimization.
Each synthetic stored field includes support type, grid geometry, the explicit
support-local-to-world transform and its inverse world-to-local mapping,
coefficient array and builder state. Do not infer or apply an unspecified
transform twice. If a future field input has no recoverable transform, require a
supplied transform or explicit identity declaration; otherwise fail before
publishing a partial field. Preserve existing captured fields and manifests as
valid historical records. Spec 005 owns shared frame and coordinate-conversion
repairs. The folded case records its
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
must never update references. `manifest.json` records input and output SHA256, baseline source revision, the
three D0 revision fields, field-generation revision, dependency versions,
coordinate transforms, seeds and levels. Preserve this as a valid historical
record. If a future reviewed correction changes references, append a distinct
correctness revision and affected-file record; keep performance baseline
revision and provenance separate. Record the actual corrected pre-performance
worktree revision, not the historical research revision unless they are the
same. Do not compare across correctness epochs as if only performance changed.
Binary output hashes and readable summary diffs make regeneration reviewable.
Capture after T034 and before any performance change. Use a scoped monkeypatch
to perturb fault movement by 1% for the sensitivity test; leave no patch active.

### D2. Benchmark tools and baseline gate

Create a developer-only package `benchmarks/fault_evaluation` with `run.py`,
`compare.py`, `instrumentation.py` and `acceptance.py`. Add package __init__.py
files. pytest's existing testpaths=[tests] keeps benchmark tests out of CI.
The suite file is `benchmarks/fault_evaluation/bench_cases.py`; it is selected
explicitly by the runner. Never add an autouse runtime counter or CI job.

`python -m benchmarks.fault_evaluation.run --output PATH [--rounds 1]
[--suite standard|full] [--work-dir PATH] [--resume]` calls pytest-benchmark.
The runner keeps `--suite standard` as the default (FR-012b). The runner
already shows progress, saves each finished case to a checkpoint file and can
resume (commits 9610f6a1 and 290edf3d). Put `--output` and `--work-dir`
outside the repository for a real run, because the dirty flag counts untracked
files, and use `--resume` after a stop.

Case lists (all in `run.py`):

- `enumerate_cases()` keeps its name and returns the FULL suite: 21 models, five
  operations (oldest value, all group values, model ids, oldest gradient, model
  gradient), cold and warm states, and 5,000, 50,000, 100,000 and 500,000
  points, plus the 20,000-point cases of `plugin` and `f1_u0`. That is 860
  cases. The order is stable. The tests of T007 depend on this function.
- `standard_cases()` returns the STANDARD suite, a subset of the full list in the
  same order: (a) all 210 cases at 5,000 points; (b) the 20 cases at 20,000
  points; (c) 24 cases at 50,000 and 100,000 points, 12 for each size, taken
  from a fixed list of case ids `STANDARD_LARGE_CASES` in `run.py`. T036 draws the
  list once with `numpy.random.default_rng(2003)` from the pool of cases of
  these models: `f0_u0..f0_u3`, `f1_u0..f1_u3`, `stopped_fault`, `onlap_domain`
  and `faulted_fold`. The code stores the ids as literal text, because the
  random stream of NumPy can change between versions. That is 254 cases.
- `select_cases(suite: str) -> list[dict]` returns `standard_cases()` for
  `"standard"` and `enumerate_cases()` for `"full"`, and raises ValueError for
  any other name. `main` uses it.

Every case that SC-002 to SC-010 use is in the standard suite: `f4_u3` and
`f4_u0` at 5,000 points (SC-002 to SC-004, SC-010), `plugin` at 20,000 points
(SC-005) and `f1_u0` at 20,000 points (SC-009).

Reuse stored fields; exclude model construction from timed evaluation. Time cold
calls after clearing the cache, then warm calls separately. A warm case primes
the cache with one untimed call before the timed call, but only if the model
has an `evaluation_cache`. On the code before this feature the model has none,
so the priming call is useless and the runner skips it: this makes a warm
legacy case cost one call, not two, and changes no measured value. Warmup must
not turn a cold case into a cache hit. Include the cost of hashing and
invalidation in timing.

Short cases (SC-007): the timing pass runs every case with one round. When it
finishes, the runner repeats each case with a median below 1.0 s for five
rounds in a second step, and the new samples replace the old ones for that
case. The counters and memory passes keep one round. The comparison and the
acceptance check use the rule in D2b.

Measured costs on 2026-11-02, on the development computer, for one cold call at
5,000 points: `f4_u3` 42 s for the five operations, `cutting_faults` 164 s,
`plugin` 11 s, `f4_u0` 9 s, `f1_u0` 2 s, other models under 3 s. The cold sum
over all 21 models is about 320 s. A case costs one call in the cold state and
two in the warm state on the code before this feature (one call with the
priming skip), in each of three passes. The standard suite should take about
50 minutes for the whole run.

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
counts, peak extra RSS bytes, code revision and dirty state, the three D0
revision fields, input hash, Python/NumPy/SciPy/LoopStructural versions, OS, CPU and RAM. Reject duplicate
ids and incompatible or missing cases when comparing. Preserve pytest-benchmark
data in the JSON. `compare BASE NEW --output PATH` writes ratios, noise marks
above 1.10, fault/unconformity scaling tables and provenance differences.
Report undefined ratios for zero baselines explicitly. Do not silently skip
cases or compare warm new results with cold baseline results.

T011 writes baseline.json and baseline_repeat.json (the standard suite, run
twice) plus a readable baseline report. It records complete runtime and noise.
T011 is a hard gate: all golden checks must pass and all baseline cases must
exist before any performance edit. The only earlier evaluation edit is the
isolated T034 correctness fix. If the run fails, capture the exact failing
case; do not substitute research numbers for a machine-readable measured
baseline. SC-007/008 are measured requirements, not promised outcomes of
planning. Run T011 from a clean committed tree, with `--output` and
`--work-dir` outside the repository, then copy the result files into the spec
folder. Start no other task and edit no tracked file while a run is active.

### D2b. How the acceptance check treats noise and memory

SC-006 compares the peak memory of each case of the standard suite. SC-007
uses the second run (`new["repeat"]`): for each case whose baseline median is
1.0 s or more, the larger ratio of the two runs must be below 1.10. A case with a
shorter baseline median is reported with its noise ratio, in a separate list
`short_cases`, and never fails SC-007. A missing repeat file or a zero median
still fails. SC-008 compares `suite_wall_seconds` of the standard suite with 900
s. The old acceptance test of T028 treats every case the same way for SC-007;
T039 changes it.

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

### D4a. Cache key from the fault set (replan of 2026-11-03)

This section corrects D4. D4 puts the fingerprint of the whole feature and the
state of the whole model in the key. T030 shows two defects. The 4 groups of
`f4_u3` have the same ordered faults but 4 different feature fingerprints, so
separate top-level calls never share an entry. The model-wide token costs 74 %
of a warm `evaluate_model`. All other rules of D4 stay (bounded LRU, content
hash of the points, opaque bypass, in-place edit detection, no cached masks).
The transform owner in spec 005 defines whether a transform mutation is accepted
and rebuilt or rejected before mutation. D4a only fingerprints the current
restoration transform and prevents a stale cache hit; it does not implement
coordinate-state mutation or rebuilding.

What a restoration depends on. `FaultSegment.apply_to_points` and
`_evaluate_frame` (`_fault_segment.py`) read only the fault itself: its
frame members (fields, supports, builder data and build arguments, the member
regions that make the mask, the faults of the frame), `displacement`, `steps`,
`faultfunction`, and the bounding-box transform of the interpolator. They never
read the feature that owns the fault list. So the key uses the ordered fault set
and nothing of the root feature (its data, build arguments, regions, parent,
`value` and `sign`).

Names and signatures (all private, no public API change):

- `fault_set_fingerprint(faults: Sequence[BaseFeature], build: bool = True) ->
  bytes | None` in `_evaluation_cache.py`. With `build`, call
  `_build_reachable(fault)` for each fault in order. Then make one
  `_StateHasher`, write `f"S{len(faults)};"`, and call
  `state.node(fault, True)` for each fault in order. A value that is not a
  feature, or `_Opaque`, gives None. Do not change `dependency_fingerprint`
  or its tests; the new function reuses `_StateHasher`.
- `_update_interpolator` also writes the current bounding-box transform:
  `box = getattr(interpolator, "bounding_box", None)`; with a box, write
  `np.asarray(box.world_to_local_matrix, dtype=float)` with `_encode`; without
  a box, write `b"Nbox;"`. `_project_points` uses this transform, and the old
  fingerprint did not read it. This also strengthens `dependency_fingerprint`.
- `_update_frame` is rewritten to loop over `frame.items()`, with no `iloc`
  call. It keeps the column names, dtypes, shape and object-column hashing.
  The digest bytes may change; only equality inside one process matters.
- `RestorationCache.sync_revision(revision: int) -> None`: under the lock, if
  `revision` differs from the stored revision (start value None), call
  `_clear_locked()` and store it. It is for internal use; document that.
- `EvaluationContext.restore` changes as follows. `_fingerprint` takes the
  tuple of faults and memoises by `tuple(id(f) for f in faults)` in
  `self._fingerprints` (the scope holds the faults, so ids stay valid). It calls
  `fault_set_fingerprint(faults)` under `cache.lock`. Delete
  `_state_of_model` and `self._model_state`. Read
  `revision = int(getattr(self.model, "evaluation_revision", 0))` after the
  fingerprint, because the fingerprint builds dirty builders. Call
  `cache.sync_revision(revision)`. The persistent key is
  `("persistent", digest, points.shape, points.dtype.str, bool(reverse),
  fingerprint, revision)`. After the faults are applied, read the revision
  again. If it changed, keep the result in the scope only and do not call
  `cache.put` (D3: a build or solve during a scope voids the entry).
- The scope key (`_fault_state`) stays as it is.

Result. Groups with equal ordered faults share one entry, also across top-level
calls. The cost of the key is one fingerprint for each distinct ordered fault
tuple for each scope, not one for each feature plus the whole model. The
revision term is O(1). It also clears reuse when anything in the model changes
through a setter, for example `GeologicalModel.update()` (it always bumps the
revision). This is conservative and the documentation (T031) must say so.
Edits in place to the root feature (its data, regions) do not change a
restoration, so they no longer invalidate; edits in place to anything that a
fault reaches still change the fingerprint.

Tests that change (T040): `Mutation` in `test_invalidation.py` gets a field
`touches_faults: bool = True`, and the four root-region mutations
(`region_add`, `region_remove_inplace`, `region_threshold`, `region_sign`) get
False. `test_mutated_model_gives_a_fresh_changed_result` then asserts equal
results always, equal fault-call counts when `touches_faults`, and cached
calls not above the cache-off calls otherwise.
`test_persistent.py::test_model_change_restores_again_and_matches_the_changed_reference`
asserts `calls > 0` when `touches_faults or mutation.revision`, and result
equality always. All other T016/T017/T018 tests stay.

Measured gates for T042 (report them; a miss is a failed task):

- `f4_u3`, 5,000 points, cold, cache on: `all_group_values` and `model_ids`
  each make 4 fault calls (SC-003).
- `f4_u3` cold `model_ids` with the cache on is at most 1.10 times the cold
  call with the cache off (the cache must not slow a first call).
- Warm `model_ids`: 0 fault calls; the time of `fault_set_fingerprint` for
  `f4_u3` is reported; the target is at most 5 ms (OP5).
- `onlap_domain`, 5,000 points, `model_gradient` cold: report the time against
  the baseline (0.92 s). It was 2.52 s on 2026-11-03. No success criterion
  uses it; the number goes into measurements.md.

### D5. Gradient decision and thread/memory checks

Preserve tetrahedron construction, near-fault refinement, element scale and
gradient reconstruction. Share equal point batches across groups; never reuse
an N-point restoration for a 4N-point tetrahedron. Apply the accepted OP1
rule: one restoration per distinct ordered fault set and distinct point batch.
Tests must show counts for original points, corners, refinement and nested
fault-frame calls. Do not claim a total F-call bound across different batches
or mark FR-008 satisfied from top-level value counts alone. Retain the D0
output-array regression while changing gradient reuse.

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

1. T001's helper-API red check is complete. Next T033 must fail on the
   demonstrated gradient shape error and T034 must make it pass. Then resume
   T002 and the golden harness. The geology golden checks must pass on the
   corrected pre-performance baseline; the raw failure stays recorded in D0.
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

The implementation preflight passed 19 existing unit/API tests. Ruff reported
474 existing findings before any production change. This remains a final merge
gate failure until addressed; do not suppress rules or call it a pass. The D0
fix must add no lint findings. Do not mix unrelated lint cleanup into D0.

## Requirement coverage

| Requirement | Tasks |
|---|---|
| FR-001, FR-001a | T001, T002, T005, T006 |
| FR-002, FR-003, FR-004 | T033, T034, T003, T004, T005, T006; D0 separates the correctness fix from capture; C-006 permits only reviewed affected-reference changes and requires separate correctness provenance |
| FR-005, SC-001 | T033, T034, T003-T006, T030, T032; pre-performance baseline as defined in D0; C-006 records correctness and performance epochs separately and requires comparable provenance |
| FR-006 | T033, T034, T006, T015, T021, T027, T032; C-006 keeps correctness and performance epochs distinct and preserves unaffected results and numeric limits |
| FR-007 | T014, T015, T018, T019, T040, T041 |
| FR-008 | T014, T015, T020, T021, T040, T041; apply the accepted OP1 batch rule |
| FR-009 | T016, T017, T018, T019, T022, T023, T040, T041 (the coordinate owner controls accepted/rejected transform mutation; D4a fingerprints the current transform and prevents stale cache reuse) |
| FR-010 | T012, T013, T018, T019, T022, T023, T030; accepted OP1 actual-point budget |
| FR-010a | T022, T023 |
| FR-011 | T014, T015, T019, T032 |
| FR-012, FR-013, FR-014, FR-015 | T007, T008, T009, T010, T011 |
| FR-012a, FR-012b | T035, T036, T037, T038, T011 |
| FR-016 | T011 |
| SC-002, SC-003, SC-004, SC-005 | T028, T029, T040, T041, T042, T030 |
| SC-006 | T022, T023, T028, T039, T029, T043, T030; OP4 |
| SC-007, SC-008 | T011, T028, T037, T038, T039, T029, T043, T030; OP2 measurement gate; OP3 gap (no plan can reach 900 s without a spec decision) |
| SC-009 | T024, T025, T026, T027, T030 |
| SC-010 | T016, T017, T018, T019, T028, T029, T040, T041, T042, T030; OP5 |

## Handover and release

Spec 005 owns frame construction defects, shared world/local conversion
repairs, and mutation policy for coordinate-owner state, including persistence
metadata and copy behavior. This spec owns world-coordinate evaluation outputs
and fingerprinting restoration transforms so the cache cannot reuse a stale
result. Record separate correctness and performance epochs only when a future
reviewed correction changes references; preserve current manifest and reference
records. Keep benchmark comparisons within a matching correctness epoch. Preserve existing timing, memory, count and numeric limits;
report actual failures without claiming a gate passed.

Private names added by D4a: `fault_set_fingerprint`,
`RestorationCache.sync_revision`, the test field `Mutation.touches_faults`. No
public name changes. `dependency_fingerprint` stays and keeps its tests.

No protocol, schema, generated solver fixtures or Rust package changes.
The loop_common change is a separate package release concern. Run its tests
and describe it in that package's changelog workflow. Keep all existing plugin
imports and model evaluation signatures. Baseline and golden regeneration are
explicit developer commands. Do not edit stored references after optimization
to hide a regression. Do not push, commit or open a PR during planning.

## Resumption order

Keep all existing task IDs and T001's completed checkbox. The newly allocated
T033 and T034 are appended to tasks.md so no existing ID is reassigned.
They execute before T002 through dependencies, not by their numeric position.
T032 remains the final scheduled task and depends on all other tasks.
Next: T033 regression test -> T034 isolated fix -> finish T002's draft fixture
loader/exporter -> T003-T011 capture and baseline -> performance work.
No build, solve or test is required to make this planning edit.

## Replan of 2026-11-02 (standard suite)

New task ids T035-T039 are appended; no id is reused and the done tasks T001-T010
and T028 keep their `[x]`. Order: T035 test suite selection -> T036 implement
suite selection and the priming skip -> T037 test short-case repeats -> T038
implement the repeats -> T011 baseline of the standard suite. T039 (acceptance
test change) runs beside T036 and T037, and T029 needs it. Waves from T011
onward move by four; the dependencies of T012-T032 are unchanged except where
tasks.md says so.

## Replan of 2026-11-03 (fault-set cache key)

Cause: T030 failed SC-003, SC-004, SC-006, SC-007, SC-008 and SC-010 (see the
Answers section for the numbers and the diagnosis). The user chose to fix the
cache design. D4a gives the design. Changes of the first replan stay.

New task ids T040-T043 are appended; no id is reused. All done tasks keep
their `[x]`. T030 stays `[!]` until T041, T042 and the decisions for OP3 and
OP4 are recorded; then it is set to `[ ]` and run again. Order:

1. T040 (tests, red): fault-set key, revision term, publish guard, bounding-box
   term, `_update_frame` equality, and the two changed existing tests.
   Files: `tests/fault_evaluation/test_fault_set_key.py` (new),
   `tests/fault_evaluation/test_invalidation.py`,
   `tests/fault_evaluation/test_persistent.py`.
2. T041 (code): D4a. Files: `LoopStructural/modelling/features/_evaluation_cache.py`,
   `LoopStructural/modelling/features/_evaluation.py`. Needs T040. Done when
   `pytest tests/fault_evaluation -q` passes, the golden tests stay under 120 s,
   and Ruff finds nothing new in the two files.
3. T042 (measure, read-only for code): run the D4a gates on the clean tree with
   the runner of D2 for the cases `f4_u3`, `f4_u0` and `onlap_domain` at 5,000
   points, cold and warm, all five operations; write
   `specs/003-fault-evaluation-performance/profile-d4a.md` with fault calls,
   key time, cold and warm times, and the `onlap_domain` numbers. Needs T041.
   If SC-003 count (4) fails, stop and replan.
4. T043 (diagnosis, read-only for code): from the raw results in
   `C:/Users/lachl/dev/ls-bench-003/cand/` write
   `specs/003-fault-evaluation-performance/t030-diagnosis.md`: (a) SC-006,
   each failing case with baseline bytes, candidate bytes and cause; (b) SC-008,
   the wall time split by pass (timing, repeat of short cases, counters,
   memory, warm priming, process starts) and the two `onlap_domain` 100,000-point
   gradient cases. It makes the facts for OP3 and OP4. Needs T029 only, so it
   can run beside T040. It changes no tracked file other than the report.
5. T030 runs again after T041, T042 and the answers to OP3 and OP4. `acceptance.py`
   is not changed: SC-010 divides the warm call by the cold call of the same
   candidate run, and that is what the spec says ("the first call").
6. T031 must also document D4a, including the conservative revision term.

The orchestrator must not edit `specs/003-fault-evaluation-performance/spec.md`.
OP3 and OP4 need the spec-interview or spec-reviewer skill and then the
spec-cascade skill.

## Answers

- 2026-11-01: The user approved the constitution draft in the main checkout.
  The draft marker was removed and the stricter pytest/Ruff merge rule activated.
- OP1, 2026-11-01: The user accepted a limit per distinct ordered fault set and distinct point batch. Count recursive frame batches separately. Measure memory per actual restored point, including tetrahedron corners, excluding unchanged baseline work arrays. This resolves the interpretation for T021 and dependent tasks; the measured targets still apply.
- OP2: Pending measurements. The 15-minute target and 10% noise limit remain.
  Both now apply to the standard suite (see the 2026-11-02 answers below).

- 2026-11-01: After the exporter failed on the existing gradient shape defect,
  the user requested spec-planner to add the correctness prerequisite. T033/T034
  define the regression and narrow fix. No production change was made in this
  planning run; golden capture and timing wait for that fix to pass.

- 2026-11-02: The user instructed the implementer to keep working through fixture repair retries, overriding the skill's two-failure stop limit for those repairs. Keep tests and review requirements unchanged.

- 2026-11-02: The first run of the full 650-case suite did not finish in 4
  hours and the user stopped it. The user asked for a suite of about 1 hour.
  Through the spec-reviewer skill the user chose: a standard suite of about 1
  hour on the code before this feature, with 12 fixed random cases at 50,000
  and 100,000 points each (cheaper models only), a 10% noise limit for cases
  of 1 s or more, and a 15-minute limit after the change (spec C-006 to C-010).
  The user also chose that every time is measured and none is inferred.
- OP2, answered 2026-11-03: the baseline of the standard suite failed SC-007 (35 of
  124 cases of 1.0 s or more differ by 10 % or more between the two runs; see
  baseline.md). The user chose to continue with the cache and grid tasks and to
  decide at T030, after a new measurement on a quiet computer. SC-007 is not
  changed. A failed SC-007 at T030 still keeps acceptance open.
- 2026-11-02, planner decisions that follow from this and need no user answer:
  the runner skips the useless priming call when the model has no cache;
  short cases repeat for five rounds in a second step; the list of large cases
  is stored as literal ids. The runner correction (progress output, checkpoint
  file, `--resume`, `--work-dir`, `suite_wall_seconds`) is already committed.
- T030, 2026-11-03: the first acceptance measurement failed. Code at revision
  c8e3f9f1 (T001 to T029 done). Two runs of the standard suite (254 cases),
  outside the repository, clean tree. The raw results are in
  `C:/Users/lachl/dev/ls-bench-003/cand/` (`candidate.json`,
  `candidate_repeat.json`, `acceptance_cand.json`, logs). Passed: SC-002,
  SC-005, SC-009. Failed: SC-003 (16 fault calls, not 4; time ratio 0.214),
  SC-004 (6.16 against 1.5), SC-006 (many cases above 1.2; several have a very
  small baseline memory, check the rule), SC-007 (cases above 1.10; the same
  noise problem as the baseline), SC-008 (3,017 s against 900 s; baseline
  3,271 s), SC-010 (0.48 against 0.10).
- T030 diagnosis, 2026-11-03 (read-only profile of `f4_u3`, 5,000 points):
  (a) The persistent key uses `dependency_fingerprint(feature)` of the whole
  feature, with its own data and build arguments. The 4 groups have the same
  ordered faults but 4 different fingerprints, so separate top-level calls
  never share an entry (16 fault calls with the cache on and off). (b) The
  model-wide state token `_state_of_model` is rebuilt in every scope. It costs
  0.137 s of a 0.185 s warm `evaluate_model` (74 %), mostly pandas `iloc`
  calls in `_update_frame`. The cache makes cold `evaluate_model` 0.135 s
  slower (0.355 s against 0.220 s). (c) SC-010 in `acceptance.py` divides the
  warm case by the candidate COLD case (0.213 s / 0.445 s), not by the
  baseline. Against the baseline (6.4 s) the warm ratio is already 0.033.
  Warm time with zero fault calls and no tokens is about 0.017 s. (d) SC-008:
  the sum of the candidate medians is 334 s (baseline 943 s), but the suite
  takes 3,017 s. About 2,000 s is not evaluation (counters pass, memory pass,
  warm priming, process starts). The two `onlap_domain` 100,000-point gradient
  cases take about 128 s per pass and are almost unchanged. Faster faults alone
  cannot reach 900 s. (e) `onlap_domain` is slower than the baseline at 5,000
  points (`model_gradient` cold 2.52 s against 0.92 s); the token and
  fingerprint cost is the probable cause (not profiled).
- T030 decision, 2026-11-03: the user chose to fix the cache design. Proposed
  by the diagnosis, for the spec-planner skill to turn into tasks: (A) the
  persistent key comes from the ordered fault set and the state of those faults
  only (not the root data, build arguments or regions of the feature, which
  cannot change a fault restoration; check `FaultSegment` for this); this needs
  a correction of D4 and a change of the T016/T017 invalidation tests that
  expect root data edits to invalidate; (B) the model-state token must be O(1)
  per call, for example cached per `(model, evaluation_revision)`, or removed
  if A covers it; every in-place edit that no setter sees must then be decided
  in D4; (C) optional: cache `dependency_fingerprint` per revision; hash
  `_update_frame` through `frame.items()`. SC-003 and SC-004 should then be
  reachable; SC-010 needs a warm call of about 0.02 s. The SC-008 choice (change
  the target or shorten the suite) is not yet decided by the user. T030 stays
  `[!]` until the new tasks are done. The implementer stops here: a plan change
  needs the spec-planner skill.
- Pickle defect, 2026-11-03: the orchestrator found it in the check before the
  T040 wave. `RestorationCache` (T013, T023) keeps a `threading.RLock` and has no
  `__getstate__`. After a `dill` round trip (`GeologicalModel.to_file`, then
  `from_file`) the lock cannot be acquired, so the first evaluation of a loaded
  model waits forever in `evaluation_scope`
  (`tests/integration/test_p2fold_model.py::test_p2fold_save_load` hangs). The
  user chose to add the fix to T040 and T041: T040 adds a red regression test,
  and T041 adds `RestorationCache.__getstate__` and `__setstate__`. The state
  keeps `enabled` and `max_bytes`, drops the entries, and makes a new `RLock`
  on load. A loaded cache starts empty, with no stored revision.
- T040 transform test, 2026-11-04 (orchestrator, from the cascade text): spec
  005 has no plan and selects no transform policy (its FR-020 allows rebuild or
  reject before mutation). Today a transform change is a direct in-place edit.
  T040 group (5) therefore tests only the rule that both policies need: after a
  change of the bounding-box transform, both fingerprints change, the cache
  gives no stale hit, and the cache-on result equals the cache-off result. The
  test of the selected owner policy belongs to spec 005 (its SC-016).
- D4 tests after D4a, 2026-11-04: T041 made the D4a code. Four existing tests
  still asserted the old D4 key. D4a keys each nested fault set on its own
  faults, so after an in-place edit of a fault list (`faults_inplace_reverse`,
  `faults_inplace_pop`) or of one frame (`frame_c_inplace_child`) the nested
  sets that did not change are reused correctly: the results equal cache-off,
  with fewer fault calls (4 vs 8, 3 vs 7, 20 vs 22). The user chose:
  `test_mutated_model_gives_a_fresh_changed_result` always asserts equal results
  and cached calls <= cache-off calls, and when `touches_faults` also cached
  calls > 0. `test_opaque_dependency_bypasses_persistent_entries` puts the
  opaque region on a fault of `cutting_faults`, not on the root group (D4a does
  not read root regions; they change the mask, never the restoration). It keeps
  the check of a hidden change of the closure state. T041 code stays as it is.
  Detail (orchestrator, from D4a): `FaultSegment.regions` is read only in
  `evaluate_gradient`; the restoration reads the regions of the frame members
  `fault[0..2]`. The opaque closure therefore goes on
  `_fault(model, 0)[0].regions` of `cutting_faults`.
