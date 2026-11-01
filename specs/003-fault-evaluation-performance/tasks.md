# Tasks: Faster evaluation of faulted models

**Input**: spec.md and plan.md in `specs/003-fault-evaluation-performance/`
**State**: Replanned; T001 complete. Run T033/T034 before resuming T002.

## Orchestrator notes

- **Worktree**: use the existing `C:/Users/lachl/dev/LoopStructural-native`,
  branch `spec/003-fault-evaluation-performance`. Do not switch the main branch.
  Approved setup is present and strict doctor passes. Use the existing
  fixture drafts and keep T001 complete. See plan.md, Planning inputs and worktree.
- **Rules for each sub agent**: read the approved constitution and plan.md;
  change only Files; do not commit, change branches or run competing checks.
  The orchestrator runs tests and records the expected failure before the next
  implementation task. Use ASD-STE100 comments and docs.
- **Tiers**: 3 [S], 31 [M]. Resolve the actual worker models from
  `.agents/model-policy.json` and available runtime tools. Codex workers use
  gpt-6-luna (low effort for S; medium for M). Limit is 5 workers, capped by
  runtime slots (currently 3 besides the orchestrator). Report unavailable models.
- **Waves**: 26. Use parallel workers only for ready tasks with disjoint Files.
- **After each wave**: one relevant pytest run plus Ruff for changed production
  files; new test tasks must fail for the stated reason. No Cargo command is
  appropriate. Staged tests may fail at named missing private APIs.
- **Checkpoints**: T033 fails on the N-by-4/N-by-3 output error; T034 fixes
  only that defect before T002 resumes. T006 golden suite passes in <120s and
  rejects a 1% perturbation; T011 freezes the corrected pre-performance baseline;
  T019 verifies value sharing/invalidation; T021 applies the accepted OP1 rule;
  T030 verifies measured targets; T032 runs the full project/package checks.
- **Capture gate**: D0/T034 is the only evaluation edit allowed before T005/T011.
  Record prerequisite parent/fix/baseline revisions; never capture after a
  performance change or remove a failing gradient case. T001 stays complete.
- **Task IDs**: T033/T034 are appended to preserve existing IDs. Schedule from
  Needs and Wave. Next task is T033, not T002. T032 is the final scheduled task.
- **Python environment**: the shared environment is
  C:/Users/lachl/dev/LoopStructural/.venv/Scripts/python.exe. Set PYTHONPATH to
  this worktree root plus packages/loop_common/src and packages/loop_interpolation/src
  (absolute paths), so editable imports do not resolve to the main checkout.
  Pytest needs write access to its worktree cache. Black and Ruff are installed.
- **Known lint baseline**: 474 existing Ruff findings. Keep the final gate strict;
  do not suppress findings or mix unrelated lint cleanup into the D0 fix.
- **Blocked tasks**: do not start a task with unresolved Blocked by. A proposed
  interpretation in plan.md is not a user answer. Record answers before work.
- **Before done**: Black check, Ruff, root pytest, both affected package suites,
  complete golden coverage and the local acceptance report must pass.
- **Commits**: Conventional Commits, one passed implementation wave per commit.
  Do not commit, push or open a PR during planning. Do not regenerate optimized
  golden outputs to hide a difference.

## Phase 1: Golden tests

- [x] T001 [US1] [M] Test the golden model catalogue and stored-field loader.
  - Do: Test D1 model_names, build_model and sample_points for all 21 names. Assert fault order, group age, unconformity relationships, all special-case topology, stratigraphic columns, deterministic points and no solver calls. Validate plugin coordinates against qgis_check.npz. Name each missing API in the failure.
  - Files: `tests/fault_evaluation/test_models.py`
  - Needs: none
  - Wave: 1
  - Done when: python -m pytest tests/fault_evaluation/test_models.py -q fails for absent model APIs, then passes after T002.

- [ ] T002 [US1] [M] Create deterministic golden models and saved fields.
  - Do: Implement the exact D1 catalogue/signatures. Extract the synthetic recipe without importing measure_baseline.py. After T034, resume the existing drafts and store solved grids/coefficients and fold state on the corrected pre-performance code. No field archive exists yet. Preserve actual frame-coordinate ownership and skip analytical coordinates; do not repeat the completed T001 test task. Load the supplied plugin data by path. Do not edit input plugin files or tests.
  - Files: `tests/fault_evaluation/__init__.py`, `tests/fault_evaluation/models.py`, `tests/fault_evaluation/data/recipes.json`, `tests/fault_evaluation/data/fields.npz`
  - Needs: T034
  - Wave: 4
  - Done when: python -m pytest tests/fault_evaluation/test_models.py -q passes.

- [ ] T003 [US1] [M] Test golden capture, comparisons and reference write protection.
  - Do: Test D1 capture keys and shapes, all outputs, exact NaN masks, box-scaled tolerances, boundary-id exception, surface vertex/bounds/area checks including zero, and model/output/max-error messages. Test missing reference and non-write mode. A mismatching fault result must be rejected.
  - Files: `tests/fault_evaluation/test_golden_tools.py`
  - Needs: T002
  - Wave: 5
  - Done when: python -m pytest tests/fault_evaluation/test_golden_tools.py -q fails for absent capture/comparison/regeneration APIs, then passes after T004.

- [ ] T004 [US1] [M] Implement golden capture and explicit regeneration.
  - Do: Implement D1 capture and assert_matches. The regenerate module exports main(argv=None), loads frozen fields, writes only with --write, and records revision, hashes, versions, seeds, transforms, surface levels, prerequisite_parent_revision, prerequisite_fix_revision and baseline_revision from D0. No evaluation or solver change. Do not modify the tests.
  - Files: `tests/fault_evaluation/golden.py`, `tests/fault_evaluation/regenerate.py`
  - Needs: T003
  - Wave: 6
  - Done when: python -m pytest tests/fault_evaluation/test_golden_tools.py -q passes.

- [ ] T005 [US1] [M] Freeze all legacy golden outputs and provenance.
  - Do: Run the sole explicit D1 write command once. Record the actual pre-performance source revision, dirty-state details and all three D0 revision fields; identify the known failing raw parent separately. Confirm plugin outputs match qgis_check.npz. Inspect readable summary/hashes, surface counts and NaN masks. The isolated T034 correctness fix is the sole earlier evaluation change; no performance code may change before this task and T011 pass.
  - Files: `tests/fault_evaluation/data/f0_u0.npz`, `tests/fault_evaluation/data/f0_u1.npz`, `tests/fault_evaluation/data/f0_u2.npz`, `tests/fault_evaluation/data/f0_u3.npz`, `tests/fault_evaluation/data/f1_u0.npz`, `tests/fault_evaluation/data/f1_u1.npz`, `tests/fault_evaluation/data/f1_u2.npz`, `tests/fault_evaluation/data/f1_u3.npz`, `tests/fault_evaluation/data/f2_u0.npz`, `tests/fault_evaluation/data/f2_u1.npz`, `tests/fault_evaluation/data/f2_u2.npz`, `tests/fault_evaluation/data/f2_u3.npz`, `tests/fault_evaluation/data/f4_u0.npz`, `tests/fault_evaluation/data/f4_u1.npz`, `tests/fault_evaluation/data/f4_u2.npz`, `tests/fault_evaluation/data/f4_u3.npz`, `tests/fault_evaluation/data/stopped_fault.npz`, `tests/fault_evaluation/data/cutting_faults.npz`, `tests/fault_evaluation/data/onlap_domain.npz`, `tests/fault_evaluation/data/faulted_fold.npz`, `tests/fault_evaluation/data/plugin.npz`, `tests/fault_evaluation/data/manifest.json`
  - Needs: T004
  - Wave: 7
  - Done when: python -m tests.fault_evaluation.regenerate --write completes after T034 on the corrected pre-performance baseline and records all 21 models with D0 provenance.

- [ ] T006 [US1] [M] Add the passing legacy golden suite and prove sensitivity.
  - Do: Load every D1 output and reference; never regenerate in tests. Cover stopped/cutting faults, onlap/domain, faulted fold, reverse, ignore_regions, faults_enabled, zero displacement, empty/outside/NaN points and first-use build. The sensitivity test expects assert_matches to raise under a temporary 1% fault-movement monkeypatch. Document that geology tests pass before optimization; unlike new feature tests these must pass on the corrected pre-performance baseline defined in D0.
  - Files: `tests/fault_evaluation/test_golden.py`, `tests/fault_evaluation/test_edges.py`, `tests/fault_evaluation/legacy-check.md`
  - Needs: T005
  - Wave: 8
  - Done when: python -m pytest tests/fault_evaluation -q passes in less than 120 seconds; the scoped 1% movement perturbation is rejected.

## Phase 2: Benchmark baseline

- [ ] T007 [US3] [M] Test benchmark runner, schema, counters and isolation.
  - Do: Use tiny fake models; do not run real benchmarks in CI. Test D2 full case enumeration, cold/warm isolation, separate process passes, schema fields, exceptional unpatch, normal runtime imports without wrappers, explicit suite selection and valid memory units. Assert all 16 F/U combinations, special models and five functions at all three sizes.
  - Files: `tests/fault_evaluation/test_benchmark_run.py`
  - Needs: T006
  - Wave: 9
  - Done when: python -m pytest tests/fault_evaluation/test_benchmark_run.py -q fails for absent benchmark modules, then passes after T008.

- [ ] T008 [P] [US3] [M] Implement benchmark runner and opt-in instrumentation.
  - Do: Implement D2 run/main and case generation using pytest-benchmark. Add it to the dev extra only and update the lock with the repository tool. No production counters. Use cold cache clears outside timed regions and include lookup costs. Record real field/fault calls, provenance, RSS and raw timing samples in schema v1. Keep benchmark suite out of normal pytest collection.
  - Files: `pyproject.toml`, `uv.lock`, `benchmarks/__init__.py`, `benchmarks/fault_evaluation/__init__.py`, `benchmarks/fault_evaluation/run.py`, `benchmarks/fault_evaluation/bench_cases.py`, `benchmarks/fault_evaluation/instrumentation.py`
  - Needs: T007
  - Wave: 10
  - Done when: python -m pytest tests/fault_evaluation/test_benchmark_run.py -q passes.

- [ ] T009 [P] [US3] [S] Test benchmark comparison ratios and missing-case errors.
  - Do: Test compare_results(base: dict, new: dict) -> dict with fixed v1 sample cases: time 2 to 3 gives ratio 1.5 and slower=True; 2 to 2.1 is within 10%; zero denominator gives an explicit undefined ratio. Missing/duplicate cases, different point counts or cold/warm states raise ValueError. Check F/U scaling table and provenance warnings.
  - Files: `tests/fault_evaluation/test_benchmark_compare.py`
  - Needs: T007
  - Wave: 10
  - Done when: python -m pytest tests/fault_evaluation/test_benchmark_compare.py -q fails for absent compare_results, then passes after T010.

- [ ] T010 [US3] [S] Implement benchmark result comparison.
  - Do: Implement D2 compare_results(base: dict, new: dict) -> dict and main(argv=None). Parse base/new paths and --output. Emit time and count ratios, >1.10 noise marks, F/U scaling and provenance differences. Use only stdlib. Do not change tests.
  - Files: `benchmarks/fault_evaluation/compare.py`
  - Needs: T009, T008
  - Wave: 11
  - Done when: python -m pytest tests/fault_evaluation/test_benchmark_compare.py -q passes.

- [ ] T011 [P] [US3] [M] Record the complete pre-change benchmark baseline twice.
  - Do: Run D2 runner twice after T034 but before any performance edit. Store machine identity, actual revision, input hashes, full times/counts/peak memory, total duration and the three D0 revision fields. Raw failing gradients have no numeric baseline; use the corrected pre-performance revision for complete comparisons. Run the compare CLI. Record SC-007/008 failures honestly and activate OP2 if they fail. Missing/failed cases block this task; do not replace with research numbers or a reduced suite.
  - Files: `specs/003-fault-evaluation-performance/baseline.json`, `specs/003-fault-evaluation-performance/baseline_repeat.json`, `specs/003-fault-evaluation-performance/baseline.md`
  - Needs: T006, T010
  - Wave: 12
  - Done when: Both complete legacy benchmark result files exist, have identical case ids and are accepted by the compare command; baseline.md records runtime and noise.

## Phase 3: Evaluation sharing and cache

- [ ] T012 [P] [Found] [M] Test restoration cache storage and exact byte limits.
  - Do: Test D4 constructor/defaults, get/put, LRU access and eviction, overwrite accounting, max_bytes shrink, negative limit, zero/disabled modes, oversize entries, clear/invalidate, immutable arrays and exception-safe publication. Count ndarray.nbytes and reject a second retained copy.
  - Files: `tests/fault_evaluation/test_cache_storage.py`
  - Needs: T011
  - Wave: 13
  - Done when: python -m pytest tests/fault_evaluation/test_cache_storage.py -q fails for absent RestorationCache, then passes after T013.

- [ ] T013 [P] [Found] [M] Implement bounded restoration cache storage.
  - Do: Implement D4 RestorationCache using OrderedDict and RLock with exact array-byte accounting. Store restored arrays read-only. No input arrays, masks or final outputs are retained. Keep metadata accounting explicit. Do not yet add public evaluation reuse.
  - Files: `LoopStructural/modelling/features/_evaluation_cache.py`
  - Needs: T012
  - Wave: 14
  - Done when: python -m pytest tests/fault_evaluation/test_cache_storage.py -q passes.

- [ ] T014 [US2] [M] Test one-call restoration sharing and public behavior.
  - Do: Test D3 nested value/region calls and evaluate_model for F=4,U=3: four top-level calls with identical fault sets, correct original mask coordinates and golden values. Test different/reordered sets, reverse, disabled faults, nested scopes, standalone features, child-fault batches, opaque callables and exception cleanup. Check public signatures. Keep total nested counts visible and apply the accepted OP1 rule from Answers.
  - Files: `tests/fault_evaluation/test_context.py`
  - Needs: T013
  - Wave: 15
  - Done when: python -m pytest tests/fault_evaluation/test_context.py -q fails on repeated legacy restoration counts, then passes after T015.

- [ ] T015 [US2] [M] Share value restorations inside an evaluation scope.
  - Do: Implement D3 EvaluationContext, evaluation_scope and restore. Wrap value/model entry points, build before lookup, use original points for masks, preserve ordering/priority/world coordinates and reset ContextVar in finally. Model initializes an additive cache property but persistent reuse stays off until T019. Do not alter the gradient algorithm or tests.
  - Files: `LoopStructural/modelling/features/_evaluation.py`, `LoopStructural/modelling/features/_base_geological_feature.py`, `LoopStructural/modelling/features/_geological_feature.py`, `LoopStructural/modelling/core/geological_model.py`
  - Needs: T014
  - Wave: 16
  - Done when: python -m pytest tests/fault_evaluation/test_context.py tests/fault_evaluation/test_golden.py -q passes.

- [ ] T016 [US2] [M] Test complete cache invalidation and direct mutations.
  - Do: Test D4 setters plus direct/in-place mutations of points, coefficients, data, support transforms, build arguments, fault order, regions/thresholds/sign, faults_enabled, displacement, steps, faultfunction, build and solve. Test dependency frames and opaque callables. Compare every second result to cache-disabled changed-model evaluation and assert no stale cache lookup is used. Test dirty lazy-build occurs exactly once.
  - Files: `tests/fault_evaluation/test_invalidation.py`
  - Needs: T015
  - Wave: 17
  - Done when: python -m pytest tests/fault_evaluation/test_invalidation.py -q fails for absent revision/fingerprint invalidation, then passes after T017.

- [ ] T017 [US2] [M] Implement revision hooks and conservative state fingerprints.
  - Do: Implement D4 monotonic revision hooks and bounded-memory dependency fingerprints. Keep property assignment semantics and detect bypass mutations. Invalidate on build/solve/data/geometry/relationship changes. Traverse faults, frames and region parents cycle-safely; bypass persistent entries for opaque dependencies. Do not monkeypatch global NumPy or change solver math. Do not modify tests.
  - Files: `LoopStructural/modelling/features/_evaluation_cache.py`, `LoopStructural/modelling/features/_base_geological_feature.py`, `LoopStructural/modelling/features/fault/_fault_segment.py`, `LoopStructural/modelling/features/builders/_base_builder.py`, `LoopStructural/modelling/features/builders/_geological_feature_builder.py`, `packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py`, `LoopStructural/modelling/core/geological_model.py`
  - Needs: T016
  - Wave: 18
  - Done when: python -m pytest tests/fault_evaluation/test_invalidation.py -q passes.

- [ ] T018 [US2] [M] Test reuse between calls and additive controls.
  - Do: Test D4 model.evaluation_cache enabled/max_bytes/clear, identical copies, strided rows, row reorder, input mutation, different models/sets/reverse, oversized results, LRU and scope/cache shared storage. Second same-grid evaluate_model must make zero fault calls; changes must restore again and match changed references. Assert full lookup hashing occurs outside no timed shortcut.
  - Files: `tests/fault_evaluation/test_persistent.py`
  - Needs: T017
  - Wave: 19
  - Done when: python -m pytest tests/fault_evaluation/test_persistent.py -q fails while cross-call reuse is off, then passes after T019.

- [ ] T019 [US2] [M] Enable persistent restoration reuse with safe invalidation.
  - Do: Connect D3 scopes to D4 model cache, default enabled 64 MiB, content/state keys and conservative bypass. Compute fingerprints after builds and include lookup work in real call path. Reuse one stored array between scope and persistent cache. Clear on revision changes and do not publish partial work. Preserve existing public call signatures.
  - Files: `LoopStructural/modelling/features/_evaluation.py`, `LoopStructural/modelling/features/_evaluation_cache.py`, `LoopStructural/modelling/core/geological_model.py`
  - Needs: T018
  - Wave: 20
  - Done when: python -m pytest tests/fault_evaluation/test_persistent.py tests/fault_evaluation/test_invalidation.py tests/fault_evaluation/test_golden.py -q passes.

- [ ] T020 [P] [US4] [M] Test gradient batch counts and legacy numerical behavior.
  - Do: Test original mask batches, 4N corners, different element scales, tetrahedron refinement, repeated groups, different fault sets, child-fault frames, reverse/disabled faults and faulted fold frames. Assert golden gradients and model gradient group priority. Produce explicit batch/count/memory diagnostics under the accepted OP1 distinct-batch rule and retain the D0 output-array regression.
  - Files: `tests/fault_evaluation/test_gradient_context.py`
  - Needs: T019
  - Wave: 21
  - Done when: python -m pytest tests/fault_evaluation/test_gradient_context.py -q exposes repeated equal-batch work, then passes after T021 under the accepted OP1 rule.

- [ ] T021 [P] [US4] [M] Share equal gradient batches after the count decision.
  - Do: Implement D5 with the accepted OP1 count and memory rule recorded in Answers. Reuse identical original/corner batches; retain tetrahedron refinement and gradient arithmetic. If the decision requires a new algorithm, stop for a revised plan. Do not weaken tests or claim F total calls from a top-level-only counter.
  - Files: `LoopStructural/modelling/features/_evaluation.py`, `LoopStructural/modelling/features/_geological_feature.py`, `LoopStructural/modelling/core/geological_model.py`
  - Needs: T020
  - Wave: 22
  - Done when: python -m pytest tests/fault_evaluation/test_gradient_context.py tests/fault_evaluation/test_golden.py -q passes under the recorded OP1 decision.

- [ ] T022 [P] [US2] [M] Test simultaneous calls, cleanup and memory accounting.
  - Do: Test D4/D5 two threads evaluating one model against serial references, first-build once, separate models and context isolation, exception release and invalidation under lock. Unit tests use small batches and exact retained bytes. The local large RSS case is added by T029, not run in normal pytest. Test no duplicate scope/cache arrays and bounded hashing workspace.
  - Files: `tests/fault_evaluation/test_concurrency.py`, `tests/fault_evaluation/test_memory.py`
  - Needs: T019
  - Wave: 21
  - Done when: python -m pytest tests/fault_evaluation/test_concurrency.py tests/fault_evaluation/test_memory.py -q fails on missing synchronization/accounting behavior, then passes after T023.

- [ ] T023 [US2] [M] Complete locking and bounded-memory integration.
  - Do: Use per-model reentrant lock for full outer evaluation, builds, revisions and cache publication; nested calls reuse it. ContextVar tokens reset on exceptions. Enforce byte limits and bounded fingerprint chunks; keep faults use_threads off. Different models must not share a global lock. Do not change tests.
  - Files: `LoopStructural/modelling/features/_evaluation.py`, `LoopStructural/modelling/features/_evaluation_cache.py`, `LoopStructural/modelling/core/geological_model.py`
  - Needs: T022, T021
  - Wave: 23
  - Done when: python -m pytest tests/fault_evaluation/test_concurrency.py tests/fault_evaluation/test_memory.py -q passes.

## Phase 4: Fault-step cost

- [ ] T024 [P] [US5] [M] Test structured-grid lookup equivalence and boundaries.
  - Do: Test D6 exact private helper tuple and weights against frozen legacy calculations for random/interior/corner/face/outside/NaN/empty points and non-unit grids. Compare evaluate_value/evaluate_gradient with old arithmetic. Assert one cell-index lookup per call and unchanged scale/corner ordering. No timing assertions in unit tests.
  - Files: `packages/loop_common/tests/test_evaluation_lookup.py`
  - Needs: T011
  - Wave: 13
  - Done when: python -m pytest packages/loop_common/tests/test_evaluation_lookup.py -q fails for absent _lookup_points/_gradient_weights, then passes after T025.

- [ ] T025 [P] [US5] [M] Remove repeated structured-grid cell searches.
  - Do: Implement D6 _lookup_points and _gradient_weights and route value/gradient evaluation through them. Preserve boundary and NaN behavior, summation order and scale factors. Leave other support types and solvers unchanged. Do not alter tests or tolerances.
  - Files: `packages/loop_common/src/loop_common/supports/_3d_structured_grid.py`
  - Needs: T024
  - Wave: 14
  - Done when: python -m pytest packages/loop_common/tests/test_evaluation_lookup.py -q passes.

- [ ] T026 [P] [US5] [M] Test fault-frame shared lookup and fallback.
  - Do: Test D6 common-grid frame values and coordinate-1 gradients, initial/precomputed_gx path, all steps, masks and restored-point equivalence. Test differing geometry/transforms, child faults, opaque regions and alternate supports use the legacy fallback. Assert shared lookup counts and unchanged integration step count.
  - Files: `tests/fault_evaluation/test_fault_frame_lookup.py`
  - Needs: T025, T019
  - Wave: 21
  - Done when: python -m pytest tests/fault_evaluation/test_fault_frame_lookup.py -q fails for absent _evaluate_frame fast path, then passes after T027.

- [ ] T027 [P] [US5] [M] Reuse frame grid lookup within each fault step.
  - Do: Implement D6 private _evaluate_frame with a guarded common-geometry fast path and legacy fallback. Use the existing transform/mask semantics; preserve every step and floating-point operation order required by the golden tolerances. Do not change displacement functions, default steps or thread setting.
  - Files: `LoopStructural/modelling/features/fault/_fault_segment.py`, `packages/loop_common/src/loop_common/supports/_3d_structured_grid.py`
  - Needs: T026
  - Wave: 22
  - Done when: python -m pytest tests/fault_evaluation/test_fault_frame_lookup.py tests/fault_evaluation/test_golden.py -q passes.

## Phase 5: Acceptance tools and final checks

- [ ] T028 [P] [Found] [S] Test numerical performance acceptance reporting.
  - Do: Test check_acceptance(base: dict, new: dict) -> dict with fixed schema-v1 files. Check SC-002 <=0.40, SC-003 <=0.15, SC-004 <=1.5, SC-005 <=0.50, SC-006 <=1.20, SC-007 <10%, SC-008 <900s, SC-009 <=0.50, SC-010 zero calls and <=0.10. Missing cases, incompatible provenance and missing gradient decision are explicit failures. Never turn a failed target into a warning.
  - Files: `tests/fault_evaluation/test_acceptance.py`
  - Needs: T010
  - Wave: 12
  - Done when: python -m pytest tests/fault_evaluation/test_acceptance.py -q fails for absent check_acceptance, then passes after T029.

- [ ] T029 [P] [Found] [M] Implement local acceptance and large-memory measurements.
  - Do: Implement D2/D5 check_acceptance and CLI: --baseline PATH --candidate PATH --repeat PATH --output PATH. Add explicit --memory-points 10000000 local subprocess mode using frozen models and cold/warm/off cache modes. Emit actual RSS, restored batch sizes, retention and FR-010 accounting plus every SC result. Never run this large case via normal pytest.
  - Files: `benchmarks/fault_evaluation/acceptance.py`
  - Needs: T028, T008
  - Wave: 13
  - Done when: python -m pytest tests/fault_evaluation/test_acceptance.py -q passes.

- [ ] T030 [US3] [M] Measure full before/after results and all success criteria.
  - Do: Run the full D2 suite twice with identical baseline case ids and machine/environment. Run D5 large-memory cases and acceptance CLI. Include gradients, total recursive counts and warm invalidation. Record raw provenance, peak memory and complete elapsed times. A target miss is a failed task, not permission to remove a case or alter golden references. OP2 requires a user decision only if runtime/noise targets fail; it does not block starting measurement.
  - Files: `specs/003-fault-evaluation-performance/candidate.json`, `specs/003-fault-evaluation-performance/candidate_repeat.json`, `specs/003-fault-evaluation-performance/acceptance.json`, `specs/003-fault-evaluation-performance/measurements.md`
  - Needs: T021, T023, T027, T029
  - Wave: 24
  - Done when: The acceptance command exits 0 for every SC, all golden tests pass under 120 seconds, and measurements.md records the 10-million-point memory result.

- [ ] T031 [Finish] [M] Document cache controls, benchmark workflow and achieved limits.
  - Do: Document D1 regeneration/review, D2 run/compare/acceptance, model.evaluation_cache enabled/max_bytes/clear, invalidation/opaque fallback, serialized same-model calls, memory accounting and measured limits. Update ROADMAP only for stages actually complete. Map every FR and SC to exact test node or measured case and acceptance result. Do not create specs/README.md because no index exists.
  - Files: `doc/advanced/fault_evaluation.rst`, `ROADMAP.md`, `specs/003-fault-evaluation-performance/coverage.md`
  - Needs: T030
  - Wave: 25
  - Done when: Documentation contains runnable commands, all cache controls, accepted decision text, measured results and a requirement-to-test/measurement table with no unproved pass.

- [ ] T032 [Finish] [M] Run final project and package checks and record the result.
  - Do: Run python -m black --check LoopStructural tests/fault_evaluation benchmarks/fault_evaluation packages/loop_common/src/loop_common/supports/_3d_structured_grid.py; python -m ruff check LoopStructural packages; python -m pytest; python -m pytest packages/loop_common/tests; python -m pytest packages/loop_interpolation/tests. Use the prepared development environment. Record failures exactly. No Cargo command is needed. Do not claim done with any failed task or unresolved acceptance gate.
  - Files: `specs/003-fault-evaluation-performance/validation.md`
  - Needs: T001, T002, T003, T004, T005, T006, T007, T008, T009, T010, T011, T012, T013, T014, T015, T016, T017, T018, T019, T020, T021, T022, T023, T024, T025, T026, T027, T028, T029, T030, T031, T033, T034
  - Wave: 26
  - Done when: Black, Ruff, root pytest, loop_common pytest and loop_interpolation pytest all exit 0; validation.md records commands, revisions, golden runtime and acceptance status.

## Added foundation prerequisite (execute before T002)

New IDs follow the previous allocation. Needs and Wave put these tasks before
field export; existing IDs and completed tasks are unchanged.

- [x] T033 [Found] [M] Test the faulted-gradient output shape without a solver.
  - Do: Follow D0. In one new test module, define the BaseFeature identity-fault double with evaluate_value, evaluate_gradient, __getitem__, apply_to_points, copy and get_data; define the linear field f(xyz)=xyz @ [2,-3,0.5]+7. Use two finite points, builder=None, explicit element scale 1 and one or two identity faults. Assert shape (2,3), gradient [2,-3,0.5] at rtol=0/atol=1e-12, a region-masked row remains NaN, and a controlled mixed-sign first coordinate evaluation causes one refinement then returns the same gradient. No solves, saved fields, golden helper import or production change. The current bug must raise the demonstrated shape mismatch; do not assert that exception as the expected final behavior.
  - Files: `tests/unit/modelling/test_faulted_gradient_output.py`
  - Needs: T001
  - Wave: 2
  - Done when: `python -m pytest tests/unit/modelling/test_faulted_gradient_output.py -q` fails first with the N-by-4/N-by-3 shape mismatch, then passes after T034.

- [x] T034 [Found] [M] Preserve the N-by-3 gradient output during refinement.
  - Do: Follow D0. In GeologicalFeature.evaluate_gradient, rename only the temporary v assigned by f[0].evaluate_value(...).reshape(N,4) to fault_values, and use fault_values in both sign comparisons that form flag. Keep the original v output allocation and final masked assignment intact. Do not change geometry, masks, fault order, refinement, steps, interpolation arithmetic, solver or signatures. Do not edit tests. Record the raw exception, prerequisite parent revision, exact test commands/results and narrow diff in gradient-prerequisite.md. Commit this correctness change separately before field/golden/timing capture; T004/T005/T011 record its commit as prerequisite_fix_revision, their actual capture revision as baseline_revision and its parent as prerequisite_parent_revision. Do not call this a performance improvement or claim successful raw gradient parity.
  - Files: `LoopStructural/modelling/features/_geological_feature.py`, `specs/003-fault-evaluation-performance/gradient-prerequisite.md`
  - Needs: T033
  - Wave: 3
  - Done when: `python -m pytest tests/unit/modelling/test_faulted_gradient_output.py tests/unit/modelling/test_geological_feature.py tests/unit/modelling/test_unconformity_feature.py tests/unit/modelling/test_faults_segment.py tests/unit/test_public_api_contract.py -q` passes; the report records the raw failure and the correctness-fix provenance.

## Dependencies

T001 is complete. T033 depends on it, proves the existing gradient output
defect without a solver, and must fail before T034 repairs it. T002 waits for
T034 and resumes its saved drafts. T003-T006 then capture and check all golden
outputs on the corrected pre-performance revision; T011 freezes both full
benchmark baseline runs before any cache or grid optimization. D0 is the sole
earlier evaluation edit. The accepted OP1 batch rule applies to gradients.
OP2 is resolved by measurement; a failed runtime/noise target triggers a user
decision and never counts as passing acceptance. Cache and grid tasks keep
their prior dependency branches. T030 joins them, T031 documents passing
results, and T032 is the final scheduled task and depends on every other task.

## Wave table

| Wave | Tasks | Tier [S] | Tier [M] |
|---|---|---|---|
| 1 | T001 | 0 | 1 |
| 2 | T033 | 0 | 1 |
| 3 | T034 | 0 | 1 |
| 4 | T002 | 0 | 1 |
| 5 | T003 | 0 | 1 |
| 6 | T004 | 0 | 1 |
| 7 | T005 | 0 | 1 |
| 8 | T006 | 0 | 1 |
| 9 | T007 | 0 | 1 |
| 10 | T008, T009 | 1 | 1 |
| 11 | T010 | 1 | 0 |
| 12 | T011, T028 | 1 | 1 |
| 13 | T012, T024, T029 | 0 | 3 |
| 14 | T013, T025 | 0 | 2 |
| 15 | T014 | 0 | 1 |
| 16 | T015 | 0 | 1 |
| 17 | T016 | 0 | 1 |
| 18 | T017 | 0 | 1 |
| 19 | T018 | 0 | 1 |
| 20 | T019 | 0 | 1 |
| 21 | T020, T022, T026 | 0 | 3 |
| 22 | T021, T027 | 0 | 2 |
| 23 | T023 | 0 | 1 |
| 24 | T030 | 0 | 1 |
| 25 | T031 | 0 | 1 |
| 26 | T032 | 0 | 1 |

## Critical path

T001 -> T033 -> T034 -> T002 -> T003 -> T004 -> T005 -> T006 -> T007 -> T009 -> T010 -> T011 -> T012 -> T013 -> T014 -> T015 -> T016 -> T017 -> T018 -> T019 -> T020 -> T021 -> T023 -> T030 -> T031 -> T032

26 sequential waves; T001 is already complete. Performance acceptance can extend this path.
