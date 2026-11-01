# Tasks: Faster evaluation of faulted models

**Input**: spec.md and plan.md in `specs/003-fault-evaluation-performance/`
**State**: Ready for staged implementation after worktree setup; OP1 blocks gradients.

## Orchestrator notes

- **Worktree**: use the existing `C:/Users/lachl/dev/LoopStructural-native`,
  branch `spec/003-fault-evaluation-performance`. Do not switch the main branch.
  Bring the approved main-checkout setup into this worktree and pass strict
  doctor before any production work. See plan.md, Planning inputs and worktree.
- **Rules for each sub agent**: read the approved constitution and plan.md;
  change only Files; do not commit, change branches or run competing checks.
  The orchestrator runs tests and records the expected failure before the next
  implementation task. Use ASD-STE100 comments and docs.
- **Tiers**: 3 [S], 29 [M]. Resolve the actual worker models from
  available runtime tools. No coding-tools.md provider policy exists. Report
  unavailable capabilities; do not invent a model mapping.
- **Waves**: 24. Use parallel workers only for ready tasks with disjoint Files.
- **After each wave**: one relevant pytest run plus Ruff for changed production
  files; new test tasks must fail for the stated reason. No Cargo command is
  appropriate. Staged tests may fail at named missing private APIs.
- **Checkpoints**: T006 legacy golden suite passes in <120s and rejects a 1%
  perturbation; T011 full legacy baseline is frozen before evaluation edits;
  T019 value sharing/invalidation passes; T021 waits for OP1; T030 verifies all
  measured targets; T032 runs the full project/package checks.
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
  - Do: Implement the exact D1 catalogue/signatures. Extract the synthetic recipe without importing measure_baseline.py. Store solved grids/coefficients and fold state on unchanged evaluation code. Load the supplied plugin data by path. Do not edit input plugin files or tests.
  - Files: `tests/fault_evaluation/__init__.py`, `tests/fault_evaluation/models.py`, `tests/fault_evaluation/data/recipes.json`, `tests/fault_evaluation/data/fields.npz`
  - Needs: T001
  - Wave: 2
  - Done when: python -m pytest tests/fault_evaluation/test_models.py -q passes.

- [ ] T003 [US1] [M] Test golden capture, comparisons and reference write protection.
  - Do: Test D1 capture keys and shapes, all outputs, exact NaN masks, box-scaled tolerances, boundary-id exception, surface vertex/bounds/area checks including zero, and model/output/max-error messages. Test missing reference and non-write mode. A mismatching fault result must be rejected.
  - Files: `tests/fault_evaluation/test_golden_tools.py`
  - Needs: T002
  - Wave: 3
  - Done when: python -m pytest tests/fault_evaluation/test_golden_tools.py -q fails for absent capture/comparison/regeneration APIs, then passes after T004.

- [ ] T004 [US1] [M] Implement golden capture and explicit regeneration.
  - Do: Implement D1 capture and assert_matches. The regenerate module exports main(argv=None), loads frozen fields, writes only with --write, and records revision, hashes, versions, seeds, transforms and surface levels. No evaluation or solver change. Do not modify the tests.
  - Files: `tests/fault_evaluation/golden.py`, `tests/fault_evaluation/regenerate.py`
  - Needs: T003
  - Wave: 4
  - Done when: python -m pytest tests/fault_evaluation/test_golden_tools.py -q passes.

- [ ] T005 [US1] [M] Freeze all legacy golden outputs and provenance.
  - Do: Run the sole explicit D1 write command once. Record actual source revision and dirty-state details. Confirm plugin outputs match qgis_check.npz. Inspect readable summary/hashes, surface counts and NaN masks. No numerical code may change before this task and T011 pass.
  - Files: `tests/fault_evaluation/data/f0_u0.npz`, `tests/fault_evaluation/data/f0_u1.npz`, `tests/fault_evaluation/data/f0_u2.npz`, `tests/fault_evaluation/data/f0_u3.npz`, `tests/fault_evaluation/data/f1_u0.npz`, `tests/fault_evaluation/data/f1_u1.npz`, `tests/fault_evaluation/data/f1_u2.npz`, `tests/fault_evaluation/data/f1_u3.npz`, `tests/fault_evaluation/data/f2_u0.npz`, `tests/fault_evaluation/data/f2_u1.npz`, `tests/fault_evaluation/data/f2_u2.npz`, `tests/fault_evaluation/data/f2_u3.npz`, `tests/fault_evaluation/data/f4_u0.npz`, `tests/fault_evaluation/data/f4_u1.npz`, `tests/fault_evaluation/data/f4_u2.npz`, `tests/fault_evaluation/data/f4_u3.npz`, `tests/fault_evaluation/data/stopped_fault.npz`, `tests/fault_evaluation/data/cutting_faults.npz`, `tests/fault_evaluation/data/onlap_domain.npz`, `tests/fault_evaluation/data/faulted_fold.npz`, `tests/fault_evaluation/data/plugin.npz`, `tests/fault_evaluation/data/manifest.json`
  - Needs: T004
  - Wave: 5
  - Done when: python -m tests.fault_evaluation.regenerate --write completes on unchanged evaluation code and records all 21 models.

- [ ] T006 [US1] [M] Add the passing legacy golden suite and prove sensitivity.
  - Do: Load every D1 output and reference; never regenerate in tests. Cover stopped/cutting faults, onlap/domain, faulted fold, reverse, ignore_regions, faults_enabled, zero displacement, empty/outside/NaN points and first-use build. The sensitivity test expects assert_matches to raise under a temporary 1% fault-movement monkeypatch. Document that geology tests pass before optimization; unlike new feature tests these must not fail on legacy code.
  - Files: `tests/fault_evaluation/test_golden.py`, `tests/fault_evaluation/test_edges.py`, `tests/fault_evaluation/legacy-check.md`
  - Needs: T005
  - Wave: 6
  - Done when: python -m pytest tests/fault_evaluation -q passes in less than 120 seconds; the scoped 1% movement perturbation is rejected.

## Phase 2: Benchmark baseline

- [ ] T007 [US3] [M] Test benchmark runner, schema, counters and isolation.
  - Do: Use tiny fake models; do not run real benchmarks in CI. Test D2 full case enumeration, cold/warm isolation, separate process passes, schema fields, exceptional unpatch, normal runtime imports without wrappers, explicit suite selection and valid memory units. Assert all 16 F/U combinations, special models and five functions at all three sizes.
  - Files: `tests/fault_evaluation/test_benchmark_run.py`
  - Needs: T006
  - Wave: 7
  - Done when: python -m pytest tests/fault_evaluation/test_benchmark_run.py -q fails for absent benchmark modules, then passes after T008.

- [ ] T008 [P] [US3] [M] Implement benchmark runner and opt-in instrumentation.
  - Do: Implement D2 run/main and case generation using pytest-benchmark. Add it to the dev extra only and update the lock with the repository tool. No production counters. Use cold cache clears outside timed regions and include lookup costs. Record real field/fault calls, provenance, RSS and raw timing samples in schema v1. Keep benchmark suite out of normal pytest collection.
  - Files: `pyproject.toml`, `uv.lock`, `benchmarks/__init__.py`, `benchmarks/fault_evaluation/__init__.py`, `benchmarks/fault_evaluation/run.py`, `benchmarks/fault_evaluation/bench_cases.py`, `benchmarks/fault_evaluation/instrumentation.py`
  - Needs: T007
  - Wave: 8
  - Done when: python -m pytest tests/fault_evaluation/test_benchmark_run.py -q passes.

- [ ] T009 [P] [US3] [S] Test benchmark comparison ratios and missing-case errors.
  - Do: Test compare_results(base: dict, new: dict) -> dict with fixed v1 sample cases: time 2 to 3 gives ratio 1.5 and slower=True; 2 to 2.1 is within 10%; zero denominator gives an explicit undefined ratio. Missing/duplicate cases, different point counts or cold/warm states raise ValueError. Check F/U scaling table and provenance warnings.
  - Files: `tests/fault_evaluation/test_benchmark_compare.py`
  - Needs: T007
  - Wave: 8
  - Done when: python -m pytest tests/fault_evaluation/test_benchmark_compare.py -q fails for absent compare_results, then passes after T010.

- [ ] T010 [US3] [S] Implement benchmark result comparison.
  - Do: Implement D2 compare_results(base: dict, new: dict) -> dict and main(argv=None). Parse base/new paths and --output. Emit time and count ratios, >1.10 noise marks, F/U scaling and provenance differences. Use only stdlib. Do not change tests.
  - Files: `benchmarks/fault_evaluation/compare.py`
  - Needs: T009, T008
  - Wave: 9
  - Done when: python -m pytest tests/fault_evaluation/test_benchmark_compare.py -q passes.

- [ ] T011 [P] [US3] [M] Record the complete pre-change benchmark baseline twice.
  - Do: Run D2 runner twice before any evaluation edit. Store machine identity, actual revision, input hashes, full times/counts/peak memory and total duration. Run the compare CLI. Record SC-007/008 failures honestly and activate OP2 if they fail. Missing/failed cases block this task; do not replace with research numbers or a reduced suite.
  - Files: `specs/003-fault-evaluation-performance/baseline.json`, `specs/003-fault-evaluation-performance/baseline_repeat.json`, `specs/003-fault-evaluation-performance/baseline.md`
  - Needs: T006, T010
  - Wave: 10
  - Done when: Both complete legacy benchmark result files exist, have identical case ids and are accepted by the compare command; baseline.md records runtime and noise.

## Phase 3: Evaluation sharing and cache

- [ ] T012 [P] [Found] [M] Test restoration cache storage and exact byte limits.
  - Do: Test D4 constructor/defaults, get/put, LRU access and eviction, overwrite accounting, max_bytes shrink, negative limit, zero/disabled modes, oversize entries, clear/invalidate, immutable arrays and exception-safe publication. Count ndarray.nbytes and reject a second retained copy.
  - Files: `tests/fault_evaluation/test_cache_storage.py`
  - Needs: T011
  - Wave: 11
  - Done when: python -m pytest tests/fault_evaluation/test_cache_storage.py -q fails for absent RestorationCache, then passes after T013.

- [ ] T013 [P] [Found] [M] Implement bounded restoration cache storage.
  - Do: Implement D4 RestorationCache using OrderedDict and RLock with exact array-byte accounting. Store restored arrays read-only. No input arrays, masks or final outputs are retained. Keep metadata accounting explicit. Do not yet add public evaluation reuse.
  - Files: `LoopStructural/modelling/features/_evaluation_cache.py`
  - Needs: T012
  - Wave: 12
  - Done when: python -m pytest tests/fault_evaluation/test_cache_storage.py -q passes.

- [ ] T014 [US2] [M] Test one-call restoration sharing and public behavior.
  - Do: Test D3 nested value/region calls and evaluate_model for F=4,U=3: four top-level calls with identical fault sets, correct original mask coordinates and golden values. Test different/reordered sets, reverse, disabled faults, nested scopes, standalone features, child-fault batches, opaque callables and exception cleanup. Check public signatures. Keep total nested counts visible; do not claim OP1 solved.
  - Files: `tests/fault_evaluation/test_context.py`
  - Needs: T013
  - Wave: 13
  - Done when: python -m pytest tests/fault_evaluation/test_context.py -q fails on repeated legacy restoration counts, then passes after T015.

- [ ] T015 [US2] [M] Share value restorations inside an evaluation scope.
  - Do: Implement D3 EvaluationContext, evaluation_scope and restore. Wrap value/model entry points, build before lookup, use original points for masks, preserve ordering/priority/world coordinates and reset ContextVar in finally. Model initializes an additive cache property but persistent reuse stays off until T019. Do not alter the gradient algorithm or tests.
  - Files: `LoopStructural/modelling/features/_evaluation.py`, `LoopStructural/modelling/features/_base_geological_feature.py`, `LoopStructural/modelling/features/_geological_feature.py`, `LoopStructural/modelling/core/geological_model.py`
  - Needs: T014
  - Wave: 14
  - Done when: python -m pytest tests/fault_evaluation/test_context.py tests/fault_evaluation/test_golden.py -q passes.

- [ ] T016 [US2] [M] Test complete cache invalidation and direct mutations.
  - Do: Test D4 setters plus direct/in-place mutations of points, coefficients, data, support transforms, build arguments, fault order, regions/thresholds/sign, faults_enabled, displacement, steps, faultfunction, build and solve. Test dependency frames and opaque callables. Compare every second result to cache-disabled changed-model evaluation and assert no stale cache lookup is used. Test dirty lazy-build occurs exactly once.
  - Files: `tests/fault_evaluation/test_invalidation.py`
  - Needs: T015
  - Wave: 15
  - Done when: python -m pytest tests/fault_evaluation/test_invalidation.py -q fails for absent revision/fingerprint invalidation, then passes after T017.

- [ ] T017 [US2] [M] Implement revision hooks and conservative state fingerprints.
  - Do: Implement D4 monotonic revision hooks and bounded-memory dependency fingerprints. Keep property assignment semantics and detect bypass mutations. Invalidate on build/solve/data/geometry/relationship changes. Traverse faults, frames and region parents cycle-safely; bypass persistent entries for opaque dependencies. Do not monkeypatch global NumPy or change solver math. Do not modify tests.
  - Files: `LoopStructural/modelling/features/_evaluation_cache.py`, `LoopStructural/modelling/features/_base_geological_feature.py`, `LoopStructural/modelling/features/fault/_fault_segment.py`, `LoopStructural/modelling/features/builders/_base_builder.py`, `LoopStructural/modelling/features/builders/_geological_feature_builder.py`, `packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py`, `LoopStructural/modelling/core/geological_model.py`
  - Needs: T016
  - Wave: 16
  - Done when: python -m pytest tests/fault_evaluation/test_invalidation.py -q passes.

- [ ] T018 [US2] [M] Test reuse between calls and additive controls.
  - Do: Test D4 model.evaluation_cache enabled/max_bytes/clear, identical copies, strided rows, row reorder, input mutation, different models/sets/reverse, oversized results, LRU and scope/cache shared storage. Second same-grid evaluate_model must make zero fault calls; changes must restore again and match changed references. Assert full lookup hashing occurs outside no timed shortcut.
  - Files: `tests/fault_evaluation/test_persistent.py`
  - Needs: T017
  - Wave: 17
  - Done when: python -m pytest tests/fault_evaluation/test_persistent.py -q fails while cross-call reuse is off, then passes after T019.

- [ ] T019 [US2] [M] Enable persistent restoration reuse with safe invalidation.
  - Do: Connect D3 scopes to D4 model cache, default enabled 64 MiB, content/state keys and conservative bypass. Compute fingerprints after builds and include lookup work in real call path. Reuse one stored array between scope and persistent cache. Clear on revision changes and do not publish partial work. Preserve existing public call signatures.
  - Files: `LoopStructural/modelling/features/_evaluation.py`, `LoopStructural/modelling/features/_evaluation_cache.py`, `LoopStructural/modelling/core/geological_model.py`
  - Needs: T018
  - Wave: 18
  - Done when: python -m pytest tests/fault_evaluation/test_persistent.py tests/fault_evaluation/test_invalidation.py tests/fault_evaluation/test_golden.py -q passes.

- [ ] T020 [P] [US4] [M] Test gradient batch counts and legacy numerical behavior.
  - Do: Test original mask batches, 4N corners, different element scales, tetrahedron refinement, repeated groups, different fault sets, child-fault frames, reverse/disabled faults and faulted fold frames. Assert golden gradients and model gradient group priority. Produce explicit batch/count/memory diagnostics for OP1; do not impose an unapproved weaker bound.
  - Files: `tests/fault_evaluation/test_gradient_context.py`
  - Needs: T019
  - Wave: 19
  - Done when: python -m pytest tests/fault_evaluation/test_gradient_context.py -q exposes repeated equal-batch work, then passes after T021 under the accepted OP1 rule.

- [ ] T021 [P] [US4] [M] Share equal gradient batches after the count decision.
  - Do: Implement D5 with the exact accepted count and memory rule recorded in Answers. Reuse identical original/corner batches; retain tetrahedron refinement and gradient arithmetic. If the decision requires a new algorithm, stop for a revised plan. Do not weaken tests or claim F total calls from a top-level-only counter.
  - Files: `LoopStructural/modelling/features/_evaluation.py`, `LoopStructural/modelling/features/_geological_feature.py`, `LoopStructural/modelling/core/geological_model.py`
  - Needs: T020
  - Wave: 20
  - Done when: python -m pytest tests/fault_evaluation/test_gradient_context.py tests/fault_evaluation/test_golden.py -q passes under the recorded OP1 decision.
  - Blocked by: OP1

- [ ] T022 [P] [US2] [M] Test simultaneous calls, cleanup and memory accounting.
  - Do: Test D4/D5 two threads evaluating one model against serial references, first-build once, separate models and context isolation, exception release and invalidation under lock. Unit tests use small batches and exact retained bytes. The local large RSS case is added by T029, not run in normal pytest. Test no duplicate scope/cache arrays and bounded hashing workspace.
  - Files: `tests/fault_evaluation/test_concurrency.py`, `tests/fault_evaluation/test_memory.py`
  - Needs: T019
  - Wave: 19
  - Done when: python -m pytest tests/fault_evaluation/test_concurrency.py tests/fault_evaluation/test_memory.py -q fails on missing synchronization/accounting behavior, then passes after T023.

- [ ] T023 [US2] [M] Complete locking and bounded-memory integration.
  - Do: Use per-model reentrant lock for full outer evaluation, builds, revisions and cache publication; nested calls reuse it. ContextVar tokens reset on exceptions. Enforce byte limits and bounded fingerprint chunks; keep faults use_threads off. Different models must not share a global lock. Do not change tests.
  - Files: `LoopStructural/modelling/features/_evaluation.py`, `LoopStructural/modelling/features/_evaluation_cache.py`, `LoopStructural/modelling/core/geological_model.py`
  - Needs: T022, T021
  - Wave: 21
  - Done when: python -m pytest tests/fault_evaluation/test_concurrency.py tests/fault_evaluation/test_memory.py -q passes.

## Phase 4: Fault-step cost

- [ ] T024 [P] [US5] [M] Test structured-grid lookup equivalence and boundaries.
  - Do: Test D6 exact private helper tuple and weights against frozen legacy calculations for random/interior/corner/face/outside/NaN/empty points and non-unit grids. Compare evaluate_value/evaluate_gradient with old arithmetic. Assert one cell-index lookup per call and unchanged scale/corner ordering. No timing assertions in unit tests.
  - Files: `packages/loop_common/tests/test_evaluation_lookup.py`
  - Needs: T011
  - Wave: 11
  - Done when: python -m pytest packages/loop_common/tests/test_evaluation_lookup.py -q fails for absent _lookup_points/_gradient_weights, then passes after T025.

- [ ] T025 [P] [US5] [M] Remove repeated structured-grid cell searches.
  - Do: Implement D6 _lookup_points and _gradient_weights and route value/gradient evaluation through them. Preserve boundary and NaN behavior, summation order and scale factors. Leave other support types and solvers unchanged. Do not alter tests or tolerances.
  - Files: `packages/loop_common/src/loop_common/supports/_3d_structured_grid.py`
  - Needs: T024
  - Wave: 12
  - Done when: python -m pytest packages/loop_common/tests/test_evaluation_lookup.py -q passes.

- [ ] T026 [P] [US5] [M] Test fault-frame shared lookup and fallback.
  - Do: Test D6 common-grid frame values and coordinate-1 gradients, initial/precomputed_gx path, all steps, masks and restored-point equivalence. Test differing geometry/transforms, child faults, opaque regions and alternate supports use the legacy fallback. Assert shared lookup counts and unchanged integration step count.
  - Files: `tests/fault_evaluation/test_fault_frame_lookup.py`
  - Needs: T025, T019
  - Wave: 19
  - Done when: python -m pytest tests/fault_evaluation/test_fault_frame_lookup.py -q fails for absent _evaluate_frame fast path, then passes after T027.

- [ ] T027 [P] [US5] [M] Reuse frame grid lookup within each fault step.
  - Do: Implement D6 private _evaluate_frame with a guarded common-geometry fast path and legacy fallback. Use the existing transform/mask semantics; preserve every step and floating-point operation order required by the golden tolerances. Do not change displacement functions, default steps or thread setting.
  - Files: `LoopStructural/modelling/features/fault/_fault_segment.py`, `packages/loop_common/src/loop_common/supports/_3d_structured_grid.py`
  - Needs: T026
  - Wave: 20
  - Done when: python -m pytest tests/fault_evaluation/test_fault_frame_lookup.py tests/fault_evaluation/test_golden.py -q passes.

## Phase 5: Acceptance tools and final checks

- [ ] T028 [P] [Found] [S] Test numerical performance acceptance reporting.
  - Do: Test check_acceptance(base: dict, new: dict) -> dict with fixed schema-v1 files. Check SC-002 <=0.40, SC-003 <=0.15, SC-004 <=1.5, SC-005 <=0.50, SC-006 <=1.20, SC-007 <10%, SC-008 <900s, SC-009 <=0.50, SC-010 zero calls and <=0.10. Missing cases, incompatible provenance and missing gradient decision are explicit failures. Never turn a failed target into a warning.
  - Files: `tests/fault_evaluation/test_acceptance.py`
  - Needs: T010
  - Wave: 10
  - Done when: python -m pytest tests/fault_evaluation/test_acceptance.py -q fails for absent check_acceptance, then passes after T029.

- [ ] T029 [P] [Found] [M] Implement local acceptance and large-memory measurements.
  - Do: Implement D2/D5 check_acceptance and CLI: --baseline PATH --candidate PATH --repeat PATH --output PATH. Add explicit --memory-points 10000000 local subprocess mode using frozen models and cold/warm/off cache modes. Emit actual RSS, restored batch sizes, retention and FR-010 accounting plus every SC result. Never run this large case via normal pytest.
  - Files: `benchmarks/fault_evaluation/acceptance.py`
  - Needs: T028, T008
  - Wave: 11
  - Done when: python -m pytest tests/fault_evaluation/test_acceptance.py -q passes.

- [ ] T030 [US3] [M] Measure full before/after results and all success criteria.
  - Do: Run the full D2 suite twice with identical baseline case ids and machine/environment. Run D5 large-memory cases and acceptance CLI. Include gradients, total recursive counts and warm invalidation. Record raw provenance, peak memory and complete elapsed times. A target miss is a failed task, not permission to remove a case or alter golden references. OP2 requires a user decision if runtime/noise targets fail.
  - Files: `specs/003-fault-evaluation-performance/candidate.json`, `specs/003-fault-evaluation-performance/candidate_repeat.json`, `specs/003-fault-evaluation-performance/acceptance.json`, `specs/003-fault-evaluation-performance/measurements.md`
  - Needs: T021, T023, T027, T029
  - Wave: 22
  - Done when: The acceptance command exits 0 for every SC, all golden tests pass under 120 seconds, and measurements.md records the 10-million-point memory result.
  - Blocked by: OP1; OP2 if measured runtime/noise fails

- [ ] T031 [Finish] [M] Document cache controls, benchmark workflow and achieved limits.
  - Do: Document D1 regeneration/review, D2 run/compare/acceptance, model.evaluation_cache enabled/max_bytes/clear, invalidation/opaque fallback, serialized same-model calls, memory accounting and measured limits. Update ROADMAP only for stages actually complete. Map every FR and SC to exact test node or measured case and acceptance result. Do not create specs/README.md because no index exists.
  - Files: `doc/advanced/fault_evaluation.rst`, `ROADMAP.md`, `specs/003-fault-evaluation-performance/coverage.md`
  - Needs: T030
  - Wave: 23
  - Done when: Documentation contains runnable commands, all cache controls, accepted decision text, measured results and a requirement-to-test/measurement table with no unproved pass.
  - Blocked by: OP1; OP2 if unresolved

- [ ] T032 [Finish] [M] Run final project and package checks and record the result.
  - Do: Run python -m black --check LoopStructural tests/fault_evaluation benchmarks/fault_evaluation packages/loop_common/src/loop_common/supports/_3d_structured_grid.py; python -m ruff check LoopStructural packages; python -m pytest; python -m pytest packages/loop_common/tests; python -m pytest packages/loop_interpolation/tests. Use the prepared development environment. Record failures exactly. No Cargo command is needed. Do not claim done with any failed task or unresolved acceptance gate.
  - Files: `specs/003-fault-evaluation-performance/validation.md`
  - Needs: T001, T002, T003, T004, T005, T006, T007, T008, T009, T010, T011, T012, T013, T014, T015, T016, T017, T018, T019, T020, T021, T022, T023, T024, T025, T026, T027, T028, T029, T030, T031
  - Wave: 24
  - Done when: Black, Ruff, root pytest, loop_common pytest and loop_interpolation pytest all exit 0; validation.md records commands, revisions, golden runtime and acceptance status.
  - Blocked by: OP1; OP2 if unresolved

## Dependencies

Golden model and tool tests precede the bodies. Golden outputs must be made
and verified on unchanged evaluation code. Benchmark tools then make the full
legacy baseline. T011 is the common gate for all performance code changes.
Cache storage precedes one-call scopes, invalidation and cross-call reuse.
Gradient tests and the OP1 decision precede gradient changes. Concurrency and
grid lookup work use separate file lists where they can run together. Frame
lookup needs both the grid helper and the evaluation scope. Acceptance tools
are foundation tasks and can run before the optimized code. Final measurement
joins all implementation branches. Documentation follows passing measurement;
the final validation task depends on every other task.

## Wave table

| Wave | Tasks | Tier [S] | Tier [M] |
|---|---|---|---|
| 1 | T001 | 0 | 1 |
| 2 | T002 | 0 | 1 |
| 3 | T003 | 0 | 1 |
| 4 | T004 | 0 | 1 |
| 5 | T005 | 0 | 1 |
| 6 | T006 | 0 | 1 |
| 7 | T007 | 0 | 1 |
| 8 | T008, T009 | 1 | 1 |
| 9 | T010 | 1 | 0 |
| 10 | T011, T028 | 1 | 1 |
| 11 | T012, T024, T029 | 0 | 3 |
| 12 | T013, T025 | 0 | 2 |
| 13 | T014 | 0 | 1 |
| 14 | T015 | 0 | 1 |
| 15 | T016 | 0 | 1 |
| 16 | T017 | 0 | 1 |
| 17 | T018 | 0 | 1 |
| 18 | T019 | 0 | 1 |
| 19 | T020, T022, T026 | 0 | 3 |
| 20 | T021, T027 | 0 | 2 |
| 21 | T023 | 0 | 1 |
| 22 | T030 | 0 | 1 |
| 23 | T031 | 0 | 1 |
| 24 | T032 | 0 | 1 |

## Critical path

T001 -> T002 -> T003 -> T004 -> T005 -> T006 -> T007 -> T009 -> T010 -> T011 -> T012 -> T013 -> T014 -> T015 -> T016 -> T017 -> T018 -> T019 -> T020 -> T021 -> T023 -> T030 -> T031 -> T032

24 sequential waves. OP1 and measured acceptance failures can extend this path.
