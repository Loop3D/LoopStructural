# Tasks: Faster evaluation of faulted models

**Input**: spec.md and plan.md in `specs/003-fault-evaluation-performance/`
**State**: Replanned 2026-11-18 for the reduced scope (spec C-012 to C-016, plan D7 to D11). T001-T029 and T033-T045 are done. Next: T046 (merge master), then T047 and T049, then T048, T050, T051. T030 then runs again, with no block.

## Orchestrator notes

- Run this file with the spec-orchestrator skill. It has the rules for the
  worktree, the sub agents, the checks and the commits.
- Waves: 40. Tiers: [S] = simple model, [M] = standard model.
- Do not start a task with "Blocked by".
- **Merge (T046)**: before you send T046, run `git merge --no-commit master`
  in the worktree. Do not rebase. The sub agent resolves the conflict only in
  the Files of T046. You make the merge commit after the check passes. The
  run-state keeps its old waves; the next wave starts at T046.
- **Worktree**: use the existing `C:/Users/lachl/dev/LoopStructural-native`,
  branch `spec/003-fault-evaluation-performance`. Do not switch the main branch.
  Approved setup is present and strict doctor passes. Use the existing
  fixture drafts and keep T001 complete. See plan.md, Planning inputs and worktree.
- **Rules for each sub agent**: read the approved constitution and plan.md;
  change only Files; do not commit, change branches or run competing checks.
  The orchestrator runs tests and records the expected failure before the next
  implementation task. Use ASD-STE100 comments and docs.
- **Tiers**: 4 [S], 49 [M]. Resolve the actual worker models from
  `.agents/model-policy.json` and available runtime tools. Codex workers use
  gpt-6-luna (low effort for S; medium for M). Limit is 5 workers, capped by
  runtime slots (currently 3 besides the orchestrator). Report unavailable models.
- **Waves**: 40. Use parallel workers only for ready tasks with disjoint Files.
- **Correctness and performance epochs**: if a future analytic correction
  changes references, append provenance for the reviewed affected results.
  Keep current references and manifests as valid historical records. Preserve
  all numeric acceptance and performance limits; report measured failures.
  Compare performance only across matching correctness references and
  equivalent case/run provenance.
- **Coordinate handover**: spec 005 owns frame construction, shared world/local
  repairs, and coordinate-state mutation/rebuild policy. After T048 this spec
  keeps no restored points between calls, so no stale transform result can
  occur. Do not add a
  whole-spec dependency on 005; its planning gate waits for these document
  amendments only.
- **Long benchmark runs (T011, T030)**: run them from a clean committed tree, with
  `--output` and `--work-dir` as absolute paths OUTSIDE the repository (for
  example `C:/Users/lachl/dev/ls-bench-003/`). Use `--suite standard` (the default)
  and `--resume` after a stop. Do not start another task or edit a tracked file
  while a run is active; an edit would make the baseline invalid. Copy the result
  files into the spec folder when the runs finish. Stop and report if one run of
  the standard suite takes more than 90 minutes. Never reduce the case list.
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
- **Task IDs**: T033/T034, T035-T039, T040-T043, T044/T045 and T046-T051 are appended to preserve existing IDs.
  Schedule from Needs and Wave, not from the number. T032 is the final scheduled task.
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

- [x] T002 [US1] [M] Create deterministic golden models and saved fields.
  - Do: Implement the exact D1 catalogue/signatures. Extract the synthetic recipe without importing measure_baseline.py. After T034, resume the existing drafts and store solved grids/coefficients and fold state on the corrected pre-performance code. No field archive exists yet. Preserve actual frame-coordinate ownership and skip analytical coordinates; do not repeat the completed T001 test task. Load the supplied plugin data by path. Do not edit input plugin files or tests.
  - Files: `tests/fault_evaluation/__init__.py`, `tests/fault_evaluation/models.py`, `tests/fault_evaluation/data/recipes.json`, `tests/fault_evaluation/data/fields.npz`
  - Needs: T034
  - Wave: 4
  - Done when: python -m pytest tests/fault_evaluation/test_models.py -q passes.

- [x] T003 [US1] [M] Test golden capture, comparisons and reference write protection.
  - Do: Test D1 capture keys and shapes, all outputs, exact NaN masks, box-scaled tolerances, boundary-id exception, surface vertex/bounds/area checks including zero, and model/output/max-error messages. Test missing reference and non-write mode. A mismatching fault result must be rejected.
  - Files: `tests/fault_evaluation/test_golden_tools.py`
  - Needs: T002
  - Wave: 5
  - Done when: python -m pytest tests/fault_evaluation/test_golden_tools.py -q fails for absent capture/comparison/regeneration APIs, then passes after T004.

- [x] T004 [US1] [M] Implement golden capture and explicit regeneration.
  - Do: Implement D1 capture and assert_matches. The regenerate module exports main(argv=None), loads frozen fields, writes only with --write, and records revision, hashes, versions, seeds, transforms, surface levels, prerequisite_parent_revision, prerequisite_fix_revision and baseline_revision from D0. No evaluation or solver change. Do not modify the tests.
  - Files: `tests/fault_evaluation/golden.py`, `tests/fault_evaluation/regenerate.py`
  - Needs: T003
  - Wave: 6
  - Done when: python -m pytest tests/fault_evaluation/test_golden_tools.py -q passes.

- [x] T005 [US1] [M] Freeze all legacy golden outputs and provenance.
  - Do: Run the sole explicit D1 write command once. Record the actual pre-performance source revision, dirty-state details and all three D0 revision fields; identify the known failing raw parent separately. Confirm plugin outputs match qgis_check.npz. Inspect readable summary/hashes, surface counts and NaN masks. The isolated T034 correctness fix is the sole earlier evaluation change; no performance code may change before this task and T011 pass.
  - Files: `tests/fault_evaluation/data/f0_u0.npz`, `tests/fault_evaluation/data/f0_u1.npz`, `tests/fault_evaluation/data/f0_u2.npz`, `tests/fault_evaluation/data/f0_u3.npz`, `tests/fault_evaluation/data/f1_u0.npz`, `tests/fault_evaluation/data/f1_u1.npz`, `tests/fault_evaluation/data/f1_u2.npz`, `tests/fault_evaluation/data/f1_u3.npz`, `tests/fault_evaluation/data/f2_u0.npz`, `tests/fault_evaluation/data/f2_u1.npz`, `tests/fault_evaluation/data/f2_u2.npz`, `tests/fault_evaluation/data/f2_u3.npz`, `tests/fault_evaluation/data/f4_u0.npz`, `tests/fault_evaluation/data/f4_u1.npz`, `tests/fault_evaluation/data/f4_u2.npz`, `tests/fault_evaluation/data/f4_u3.npz`, `tests/fault_evaluation/data/stopped_fault.npz`, `tests/fault_evaluation/data/cutting_faults.npz`, `tests/fault_evaluation/data/onlap_domain.npz`, `tests/fault_evaluation/data/faulted_fold.npz`, `tests/fault_evaluation/data/plugin.npz`, `tests/fault_evaluation/data/manifest.json`
  - Needs: T004
  - Wave: 7
  - Done when: python -m tests.fault_evaluation.regenerate --write completes after T034 on the corrected pre-performance baseline and records all 21 models with D0 provenance.

- [x] T006 [US1] [M] Add the passing legacy golden suite and prove sensitivity.
  - Do: Load every D1 output and reference; never regenerate in tests. Cover stopped/cutting faults, onlap/domain, faulted fold, reverse, ignore_regions, faults_enabled, zero displacement, empty/outside/NaN points and first-use build. The sensitivity test expects assert_matches to raise under a temporary 1% fault-movement monkeypatch. Document that geology tests pass before optimization; unlike new feature tests these must pass on the corrected pre-performance baseline defined in D0.
  - Files: `tests/fault_evaluation/test_golden.py`, `tests/fault_evaluation/test_edges.py`, `tests/fault_evaluation/legacy-check.md`
  - Needs: T005
  - Wave: 8
  - Done when: python -m pytest tests/fault_evaluation -q passes in less than 120 seconds; the scoped 1% movement perturbation is rejected.

## Phase 2: Benchmark baseline

- [x] T007 [US3] [M] Test benchmark runner, schema, counters and isolation.
  - Do: Use tiny fake models; do not run real benchmarks in CI. Test D2 full case enumeration, cold/warm isolation, separate process passes, schema fields, exceptional unpatch, normal runtime imports without wrappers, explicit suite selection and valid memory units. Assert all 16 F/U combinations, special models and five functions at all three sizes.
  - Files: `tests/fault_evaluation/test_benchmark_run.py`
  - Needs: T006
  - Wave: 9
  - Done when: python -m pytest tests/fault_evaluation/test_benchmark_run.py -q fails for absent benchmark modules, then passes after T008.

- [x] T008 [P] [US3] [M] Implement benchmark runner and opt-in instrumentation.
  - Do: Implement D2 run/main and case generation using pytest-benchmark. Add it to the dev extra only and update the lock with the repository tool. No production counters. Use cold cache clears outside timed regions and include lookup costs. Record real field/fault calls, provenance, RSS and raw timing samples in schema v1. Keep benchmark suite out of normal pytest collection.
  - Files: `pyproject.toml`, `uv.lock`, `benchmarks/__init__.py`, `benchmarks/fault_evaluation/__init__.py`, `benchmarks/fault_evaluation/run.py`, `benchmarks/fault_evaluation/bench_cases.py`, `benchmarks/fault_evaluation/instrumentation.py`
  - Needs: T007
  - Wave: 10
  - Done when: python -m pytest tests/fault_evaluation/test_benchmark_run.py -q passes.

- [x] T009 [P] [US3] [S] Test benchmark comparison ratios and missing-case errors.
  - Do: Test compare_results(base: dict, new: dict) -> dict with fixed v1 sample cases: time 2 to 3 gives ratio 1.5 and slower=True; 2 to 2.1 is within 10%; zero denominator gives an explicit undefined ratio. Missing/duplicate cases, different point counts or cold/warm states raise ValueError. Check F/U scaling table and provenance warnings.
  - Files: `tests/fault_evaluation/test_benchmark_compare.py`
  - Needs: T007
  - Wave: 10
  - Done when: python -m pytest tests/fault_evaluation/test_benchmark_compare.py -q fails for absent compare_results, then passes after T010.

- [x] T010 [US3] [S] Implement benchmark result comparison.
  - Do: Implement D2 compare_results(base: dict, new: dict) -> dict and main(argv=None). Parse base/new paths and --output. Emit time and count ratios, >1.10 noise marks, F/U scaling and provenance differences. Use only stdlib. Do not change tests.
  - Files: `benchmarks/fault_evaluation/compare.py`
  - Needs: T009, T008
  - Wave: 11
  - Done when: python -m pytest tests/fault_evaluation/test_benchmark_compare.py -q passes.

- [x] T011 [P] [US3] [M] Record the complete pre-change baseline of the standard suite twice.
  - Do: CHANGED by the replan of 2026-11-02: the earlier text asked for the full 650-case matrix, which did not finish in 4 hours (stopped, no result file). Now run the STANDARD suite (254 cases, D2) twice, after T038 and before any performance edit, as the notes on long runs say: `python -m benchmarks.fault_evaluation.run --output <dir>/baseline.json --work-dir <dir>/run1.work` and the same with `baseline_repeat.json` and `run2.work`; `<dir>` is an absolute path outside the repository. Use `--resume` after a stop. The runs must start from a clean committed tree. Store machine identity, actual revision, input hashes, full times/counts/peak memory, `suite_wall_seconds` and the three D0 revision fields. Raw failing gradients have no numeric baseline; use the corrected pre-performance revision. Copy the two files into the spec folder. Run the compare CLI. Write baseline.md with the runtime of each run, the cases below 1 s, and the SC-007 noise (D2b). Record SC-007/008 failures honestly and activate OP2 if they fail. Missing/failed cases block this task; do not replace with research numbers or a reduced suite. If one run takes more than 90 minutes, stop and report.
  - Files: `specs/003-fault-evaluation-performance/baseline.json`, `specs/003-fault-evaluation-performance/baseline_repeat.json`, `specs/003-fault-evaluation-performance/baseline.md`
  - Needs: T006, T010, T038
  - Wave: 16
  - Done when: Both complete standard-suite result files exist, have identical case ids (254) and are accepted by the compare command; baseline.md records runtime and noise.

## Phase 3: Evaluation sharing and cache

- [x] T012 [P] [Found] [M] Test restoration cache storage and exact byte limits.
  - Do: Test D4 constructor/defaults, get/put, LRU access and eviction, overwrite accounting, max_bytes shrink, negative limit, zero/disabled modes, oversize entries, clear/invalidate, immutable arrays and exception-safe publication. Count ndarray.nbytes and reject a second retained copy.
  - Files: `tests/fault_evaluation/test_cache_storage.py`
  - Needs: T011
  - Wave: 17
  - Done when: python -m pytest tests/fault_evaluation/test_cache_storage.py -q fails for absent RestorationCache, then passes after T013.

- [x] T013 [P] [Found] [M] Implement bounded restoration cache storage.
  - Do: Implement D4 RestorationCache using OrderedDict and RLock with exact array-byte accounting. Store restored arrays read-only. No input arrays, masks or final outputs are retained. Keep metadata accounting explicit. Do not yet add public evaluation reuse.
  - Files: `LoopStructural/modelling/features/_evaluation_cache.py`
  - Needs: T012
  - Wave: 18
  - Done when: python -m pytest tests/fault_evaluation/test_cache_storage.py -q passes.

- [x] T014 [US2] [M] Test one-call restoration sharing and public behavior.
  - Do: Test D3 nested value/region calls and evaluate_model for F=4,U=3: four top-level calls with identical fault sets, correct original mask coordinates and golden values. Test different/reordered sets, reverse, disabled faults, nested scopes, standalone features, child-fault batches, opaque callables and exception cleanup. Check public signatures. Keep total nested counts visible and apply the accepted OP1 rule from Answers.
  - Files: `tests/fault_evaluation/test_context.py`
  - Needs: T013
  - Wave: 19
  - Done when: python -m pytest tests/fault_evaluation/test_context.py -q fails on repeated legacy restoration counts, then passes after T015.

- [x] T015 [US2] [M] Share value restorations inside an evaluation scope.
  - Do: Implement D3 EvaluationContext, evaluation_scope and restore. Wrap value/model entry points, build before lookup, use original points for masks, preserve ordering/priority/world coordinates and reset ContextVar in finally. Model initializes an additive cache property but persistent reuse stays off until T019. Do not alter the gradient algorithm or tests.
  - Files: `LoopStructural/modelling/features/_evaluation.py`, `LoopStructural/modelling/features/_base_geological_feature.py`, `LoopStructural/modelling/features/_geological_feature.py`, `LoopStructural/modelling/core/geological_model.py`
  - Needs: T014
  - Wave: 20
  - Done when: python -m pytest tests/fault_evaluation/test_context.py tests/fault_evaluation/test_golden.py -q passes.

- [x] T016 [US2] [M] Test complete cache invalidation and direct mutations.
  - Do: Test D4 setters plus direct/in-place mutations of points, coefficients, data, support transforms, build arguments, fault order, regions/thresholds/sign, faults_enabled, displacement, steps, faultfunction, build and solve. Test dependency frames and opaque callables. Compare every second result to cache-disabled changed-model evaluation and assert no stale cache lookup is used. Test dirty lazy-build occurs exactly once.
  - Files: `tests/fault_evaluation/test_invalidation.py`
  - Needs: T015
  - Wave: 21
  - Done when: python -m pytest tests/fault_evaluation/test_invalidation.py -q fails for absent revision/fingerprint invalidation, then passes after T017.

- [x] T017 [US2] [M] Implement revision hooks and conservative state fingerprints.
  - Do: Implement D4 monotonic revision hooks and bounded-memory dependency fingerprints. Keep property assignment semantics and detect bypass mutations. Invalidate on build/solve/data/geometry/relationship changes. Traverse faults, frames and region parents cycle-safely; bypass persistent entries for opaque dependencies. Do not monkeypatch global NumPy or change solver math. Do not modify tests.
  - Files: `LoopStructural/modelling/features/_evaluation_cache.py`, `LoopStructural/modelling/features/_base_geological_feature.py`, `LoopStructural/modelling/features/fault/_fault_segment.py`, `LoopStructural/modelling/features/builders/_base_builder.py`, `LoopStructural/modelling/features/builders/_geological_feature_builder.py`, `packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py`, `LoopStructural/modelling/core/geological_model.py`
  - Needs: T016
  - Wave: 22
  - Done when: python -m pytest tests/fault_evaluation/test_invalidation.py -q passes.

- [x] T018 [US2] [M] Test reuse between calls and additive controls.
  - Do: Test D4 model.evaluation_cache enabled/max_bytes/clear, identical copies, strided rows, row reorder, input mutation, different models/sets/reverse, oversized results, LRU and scope/cache shared storage. Second same-grid evaluate_model must make zero fault calls; changes must restore again and match changed references. Assert full lookup hashing occurs outside no timed shortcut.
  - Files: `tests/fault_evaluation/test_persistent.py`
  - Needs: T017
  - Wave: 23
  - Done when: python -m pytest tests/fault_evaluation/test_persistent.py -q fails while cross-call reuse is off, then passes after T019.

- [x] T019 [US2] [M] Enable persistent restoration reuse with safe invalidation.
  - Do: Connect D3 scopes to D4 model cache, default enabled 64 MiB, content/state keys and conservative bypass. Compute fingerprints after builds and include lookup work in real call path. Reuse one stored array between scope and persistent cache. Clear on revision changes and do not publish partial work. Preserve existing public call signatures.
  - Files: `LoopStructural/modelling/features/_evaluation.py`, `LoopStructural/modelling/features/_evaluation_cache.py`, `LoopStructural/modelling/core/geological_model.py`
  - Needs: T018
  - Wave: 24
  - Done when: python -m pytest tests/fault_evaluation/test_persistent.py tests/fault_evaluation/test_invalidation.py tests/fault_evaluation/test_golden.py -q passes.

- [x] T020 [P] [US4] [M] Test gradient batch counts and legacy numerical behavior.
  - Do: Test original mask batches, 4N corners, different element scales, tetrahedron refinement, repeated groups, different fault sets, child-fault frames, reverse/disabled faults and faulted fold frames. Assert golden gradients and model gradient group priority. Produce explicit batch/count/memory diagnostics under the accepted OP1 distinct-batch rule and retain the D0 output-array regression.
  - Files: `tests/fault_evaluation/test_gradient_context.py`
  - Needs: T019
  - Wave: 25
  - Done when: python -m pytest tests/fault_evaluation/test_gradient_context.py -q exposes repeated equal-batch work, then passes after T021 under the accepted OP1 rule.

- [x] T021 [P] [US4] [M] Share equal gradient batches after the count decision.
  - Do: Implement D5 with the accepted OP1 count and memory rule recorded in Answers. Reuse identical original/corner batches; retain tetrahedron refinement and gradient arithmetic. If the decision requires a new algorithm, stop for a revised plan. Do not weaken tests or claim F total calls from a top-level-only counter.
  - Files: `LoopStructural/modelling/features/_evaluation.py`, `LoopStructural/modelling/features/_geological_feature.py`, `LoopStructural/modelling/core/geological_model.py`
  - Needs: T020
  - Wave: 26
  - Done when: python -m pytest tests/fault_evaluation/test_gradient_context.py tests/fault_evaluation/test_golden.py -q passes under the recorded OP1 decision.

- [x] T022 [P] [US2] [M] Test simultaneous calls, cleanup and memory accounting.
  - Do: Test D4/D5 two threads evaluating one model against serial references, first-build once, separate models and context isolation, exception release and invalidation under lock. Unit tests use small batches and exact retained bytes. The local large RSS case is added by T029, not run in normal pytest. Test no duplicate scope/cache arrays and bounded hashing workspace.
  - Files: `tests/fault_evaluation/test_concurrency.py`, `tests/fault_evaluation/test_memory.py`
  - Needs: T019
  - Wave: 25
  - Done when: python -m pytest tests/fault_evaluation/test_concurrency.py tests/fault_evaluation/test_memory.py -q fails on missing synchronization/accounting behavior, then passes after T023.

- [x] T023 [US2] [M] Complete locking and bounded-memory integration.
  - Do: Use per-model reentrant lock for full outer evaluation, builds, revisions and cache publication; nested calls reuse it. ContextVar tokens reset on exceptions. Enforce byte limits and bounded fingerprint chunks; keep faults use_threads off. Different models must not share a global lock. Do not change tests.
  - Files: `LoopStructural/modelling/features/_evaluation.py`, `LoopStructural/modelling/features/_evaluation_cache.py`, `LoopStructural/modelling/core/geological_model.py`
  - Needs: T022, T021
  - Wave: 27
  - Done when: python -m pytest tests/fault_evaluation/test_concurrency.py tests/fault_evaluation/test_memory.py -q passes.

## Phase 4: Fault-step cost

- [x] T024 [P] [US5] [M] Test structured-grid lookup equivalence and boundaries.
  - Do: Test D6 exact private helper tuple and weights against frozen legacy calculations for random/interior/corner/face/outside/NaN/empty points and non-unit grids. Compare evaluate_value/evaluate_gradient with old arithmetic. Assert one cell-index lookup per call and unchanged scale/corner ordering. No timing assertions in unit tests.
  - Files: `packages/loop_common/tests/test_evaluation_lookup.py`
  - Needs: T011
  - Wave: 17
  - Done when: python -m pytest packages/loop_common/tests/test_evaluation_lookup.py -q fails for absent _lookup_points/_gradient_weights, then passes after T025.

- [x] T025 [P] [US5] [M] Remove repeated structured-grid cell searches.
  - Do: Implement D6 _lookup_points and _gradient_weights and route value/gradient evaluation through them. Preserve boundary and NaN behavior, summation order and scale factors. Leave other support types and solvers unchanged. Do not alter tests or tolerances.
  - Files: `packages/loop_common/src/loop_common/supports/_3d_structured_grid.py`
  - Needs: T024
  - Wave: 18
  - Done when: python -m pytest packages/loop_common/tests/test_evaluation_lookup.py -q passes.

- [x] T026 [P] [US5] [M] Test fault-frame shared lookup and fallback.
  - Do: Test D6 common-grid frame values and coordinate-1 gradients, initial/precomputed_gx path, all steps, masks and restored-point equivalence. Test differing geometry/transforms, child faults, opaque regions and alternate supports use the legacy fallback. Assert shared lookup counts and unchanged integration step count.
  - Files: `tests/fault_evaluation/test_fault_frame_lookup.py`
  - Needs: T025, T019
  - Wave: 25
  - Done when: python -m pytest tests/fault_evaluation/test_fault_frame_lookup.py -q fails for absent _evaluate_frame fast path, then passes after T027.

- [x] T027 [P] [US5] [M] Reuse frame grid lookup within each fault step.
  - Do: Implement D6 private _evaluate_frame with a guarded common-geometry fast path and legacy fallback. Use the existing transform/mask semantics; preserve every step and floating-point operation order required by the golden tolerances. Do not change displacement functions, default steps or thread setting.
  - Files: `LoopStructural/modelling/features/fault/_fault_segment.py`, `packages/loop_common/src/loop_common/supports/_3d_structured_grid.py`
  - Needs: T026
  - Wave: 26
  - Done when: python -m pytest tests/fault_evaluation/test_fault_frame_lookup.py tests/fault_evaluation/test_golden.py -q passes.

## Phase 5: Acceptance tools and final checks

- [x] T028 [P] [Found] [S] Test numerical performance acceptance reporting.
  - Do: Test check_acceptance(base: dict, new: dict) -> dict with fixed schema-v1 files. Check SC-002 <=0.40, SC-003 <=0.15, SC-004 <=1.5, SC-005 <=0.50, SC-006 <=1.20, SC-007 <10%, SC-008 <900s, SC-009 <=0.50, SC-010 zero calls and <=0.10. Missing cases, incompatible provenance and missing gradient decision are explicit failures. Never turn a failed target into a warning.
  - Files: `tests/fault_evaluation/test_acceptance.py`
  - Needs: T010
  - Wave: 12
  - Done when: python -m pytest tests/fault_evaluation/test_acceptance.py -q fails for absent check_acceptance, then passes after T029.

- [x] T029 [P] [Found] [M] Implement local acceptance and large-memory measurements.
  - Do: Implement D2/D5 check_acceptance and CLI: --baseline PATH --candidate PATH --repeat PATH --output PATH. Follow D2b for SC-006/007/008: the memory limit covers the standard suite; the 10 % noise limit covers only cases with a baseline median of 1.0 s or more; shorter cases go in `short_cases` and never fail SC-007 (T039 changes the tests to say this). Add explicit --memory-points 10000000 local subprocess mode using frozen models and cold/warm/off cache modes. Emit actual RSS, restored batch sizes, retention and FR-010 accounting plus every SC result. Never run this large case via normal pytest.
  - Files: `benchmarks/fault_evaluation/acceptance.py`
  - Needs: T028, T008, T039
  - Wave: 14
  - Done when: python -m pytest tests/fault_evaluation/test_acceptance.py -q passes.

- [ ] T030 [US3] [M] Measure full before/after results and all success criteria.
  - Do: SECOND RUN, after the replan of 2026-11-18 (plan D7 to D11). The first run (revision c8e3f9f1, 2026-11-03) failed SC-003, SC-004, SC-006, SC-007, SC-008 and SC-010; its raw files stay in `C:/Users/lachl/dev/ls-bench-003/cand/`. Use a new output folder. Do not change `acceptance.py`: T050 made it apply the D10 rules (`cold_view` of the baseline, SC-006 passes at 1.20 times the baseline or at most `MEMORY_ALLOWANCE_BYTES` = 16 MiB more, no SC-010 row, time ratios on `baseline.json`, C-015). The standard suite now has 131 cold cases (D9, C-014, C-016); the counter pass runs only the 5,000 and 20,000-point cases, the memory pass only the 5,000-point cases. SC-010 is removed (C-012); write that in measurements.md, and put the T051 numbers there. Run the STANDARD D2 suite twice with identical baseline case ids and machine/environment, in the way the notes on long runs say (outside the repository, clean tree, `suite_wall_seconds` recorded for SC-008). The full suite is not required. Run D5 large-memory cases (cold mode only) and acceptance CLI. Include gradients and total recursive counts. Record raw provenance, peak memory and complete elapsed times. A target miss is a failed task, not permission to remove a case or alter golden references. OP2 requires a user decision only if runtime/noise targets fail; it does not block starting measurement.
  - Files: `specs/003-fault-evaluation-performance/candidate.json`, `specs/003-fault-evaluation-performance/candidate_repeat.json`, `specs/003-fault-evaluation-performance/acceptance.json`, `specs/003-fault-evaluation-performance/measurements.md`
  - Needs: T027, T029, T042, T043, T051, T053
  - Wave: 38
  - Done when: The acceptance command exits 0 for every SC, all golden tests pass under 120 seconds, and measurements.md records the 10-million-point memory result.

- [ ] T031 [Finish] [M] Document cache controls, benchmark workflow and achieved limits.
  - Do: Document D1 regeneration/review, D2 run/compare/acceptance (including `--suite standard|full`, `--work-dir`, `--resume`, and why the work directory and outputs belong outside the repository), and the state after the replan of 2026-11-18: no cache between calls and no lock (C-012; two calls on one model run at the same time); the per-call scope of D3 (restored points are shared inside one call and released at its end); the standard suite of D9 with its pass point sets (`PASS_POINTS`: timing all cases, counters 5,000 and 20,000 points, memory 5,000 points) and the 16 large ids (C-016); the SC-006 rule of D10 (120 % or +16 MiB, C-013); time ratios on `baseline.json` (C-015); call-only memory accounting and measured limits. Do not document `evaluation_cache`, `evaluation_revision` or the fault-set key: T048 removed them. Update ROADMAP only for stages actually complete. Map every FR and SC to exact test node or measured case and acceptance result. Do not create specs/README.md because no index exists.
  - Files: `doc/advanced/fault_evaluation.rst`, `ROADMAP.md`, `specs/003-fault-evaluation-performance/coverage.md`
  - Needs: T030
  - Wave: 39
  - Done when: Documentation contains runnable commands, all cache controls, accepted decision text, measured results and a requirement-to-test/measurement table with no unproved pass.

- [ ] T032 [Finish] [M] Run final project and package checks and record the result.
  - Do: Run python -m black --check LoopStructural tests/fault_evaluation benchmarks/fault_evaluation packages/loop_common/src/loop_common/supports/_3d_structured_grid.py; python -m ruff check LoopStructural packages; python -m pytest; python -m pytest packages/loop_common/tests; python -m pytest packages/loop_interpolation/tests. Use the prepared development environment. Record failures exactly. No Cargo command is needed. Do not claim done with any failed task or unresolved acceptance gate.
  - Files: `specs/003-fault-evaluation-performance/validation.md`
  - Needs: T001, T002, T003, T004, T005, T006, T007, T008, T009, T010, T011, T012, T013, T014, T015, T016, T017, T018, T019, T020, T021, T022, T023, T024, T025, T026, T027, T028, T029, T030, T031, T033, T034, T035, T036, T037, T038, T039, T040, T041, T042, T043, T044, T045, T046, T047, T048, T049, T050, T051, T052, T053
  - Wave: 40
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

## Added standard-suite tasks (replan of 2026-11-02)

The first run of the full matrix did not finish in 4 hours. The spec now defines a
standard suite of about 1 hour (FR-012a, FR-012b, C-006 to C-010; plan D2, D2b).
These tasks add it. The done tasks T007, T008 and T010 are not reopened: T035 and
T036 change the case lists they created, and T039 changes the tests of T028.
The runner correction of 9610f6a1 and 290edf3d (progress, checkpoint file,
`--resume`, `--work-dir`, `suite_wall_seconds`) is already in the code.

- [x] T035 [US3] [M] Test the standard and full suite selection and the priming skip.
  - Do: CHANGES to done work: T007 tests the case list with `enumerate_cases()`; that function stays and now returns the full suite of 860 cases (sizes 5,000, 50,000, 100,000, 500,000, plus the 20,000-point cases of `plugin` and `f1_u0`). In a new file test D2: `standard_cases()` returns 254 cases (210 at 5,000 points, 20 at 20,000, 24 in `STANDARD_LARGE_CASES`: 12 at 50,000 and 12 at 100,000 points); every standard id is in the full list, in the same order; the 24 large ids use only the models `f0_u0..f0_u3`, `f1_u0..f1_u3`, `stopped_fault`, `onlap_domain`, `faulted_fold`; the ids are literal text; the standard list contains every case that SC-002 to SC-010 use (`f4_u3` and `f4_u0` at 5,000 points, cold and warm, all five operations; `plugin` and `f1_u0` at 20,000 points); `select_cases("standard")`, `select_cases("full")` and ValueError for another name; `main` has `--suite` with default `standard`, and `--suite full` selects the full list. Test the priming rule with a fake model: a warm case on a model with no `evaluation_cache` makes one call in each round (no priming call); with a cache object it primes once and clears the cache before each round; a cold case clears the cache and never primes. Also edit the two existing tests that replace `enumerate_cases` with a fake list and call `main` with no suite option (`test_benchmark_run.py` near the fake `enumerate_cases`, `test_benchmark_resume.py` near the fake `enumerate_cases`): they must pass `--suite full`. Change nothing else in those two files.
  - Files: `tests/fault_evaluation/test_benchmark_suites.py`, `tests/fault_evaluation/test_benchmark_run.py`, `tests/fault_evaluation/test_benchmark_resume.py`
  - Needs: T010
  - Wave: 12
  - Done when: `python -m pytest tests/fault_evaluation/test_benchmark_suites.py tests/fault_evaluation/test_benchmark_run.py tests/fault_evaluation/test_benchmark_resume.py -q` fails for absent `standard_cases`, `select_cases` and the `--suite` option, then passes after T036.

- [x] T036 [US3] [M] Implement the standard suite, the suite option and the priming skip.
  - Do: Implement D2. Add 100,000 to the sizes of the full suite. Add `STANDARD_LARGE_CASES` (a tuple of 24 literal case ids): draw it once, with `numpy.random.default_rng(2003)`, 12 cases for each of the sizes 50,000 and 100,000, from the pool of cases of the models named in D2; write the ids as text and keep a comment that says how they were drawn. Add `standard_cases()`, `select_cases(suite)` and `--suite {standard,full}` (default `standard`) to `main`. Keep `enumerate_cases()` as the full suite. In `bench_cases.run_case`, skip the priming call of a warm case when the model has no `evaluation_cache`. Do not change tests. Do not add counters or a production change.
  - Files: `benchmarks/fault_evaluation/run.py`, `benchmarks/fault_evaluation/bench_cases.py`
  - Needs: T035
  - Wave: 13
  - Done when: `python -m pytest tests/fault_evaluation/test_benchmark_suites.py tests/fault_evaluation/test_benchmark_run.py tests/fault_evaluation/test_benchmark_resume.py tests/fault_evaluation/test_benchmark_compare.py -q` passes.

- [x] T039 [P] [Found] [M] Change the acceptance tests for the SC-007 rule of short cases.
  - Do: CHANGES to done work: T028 treats every case in the same way for SC-007 (the larger ratio between the two runs must be below 1.10). Follow D2b: only cases with a baseline median of 1.0 s or more are checked. Add tests: a case of 1.0 s or more with a noise ratio of 1.10 or more fails SC-007 and names the case; a case with a baseline median below 1.0 s and a noise ratio of 1.5 does not fail SC-007 and appears in `result["short_cases"]` with its noise ratio; a missing repeat document or a zero median still fails. Change a fixture or an existing assertion only where the new rule makes it wrong, and say which in the final report. Keep all other tests of T028 as they are.
  - Files: `tests/fault_evaluation/test_acceptance.py`
  - Needs: T028
  - Wave: 13
  - Done when: `python -m pytest tests/fault_evaluation/test_acceptance.py -q` fails for absent `check_acceptance` and the new tests name the `short_cases` rule, then passes after T029.

- [x] T037 [US3] [M] Test the repeat of short cases in the timing pass.
  - Do: Test D2 (short cases) with a fake launcher and no real benchmark. Constants in `run.py`: `SHORT_CASE_SECONDS = 1.0` and `SHORT_CASE_ROUNDS = 5`. After the timing pass, the runner repeats each case with a median below 1.0 s for 5 rounds, and the new samples replace the old ones in the checkpoint rows and the result; a case of 1.0 s or more is not repeated; the counters and memory passes are not repeated; the recorded `rounds` and `time_samples_seconds` of a repeated case are those of the 5 rounds; a stop between the two steps and a `--resume` do not repeat a finished case a second time and do not lose one; the repeat adds its time to `suite_wall_seconds`.
  - Files: `tests/fault_evaluation/test_benchmark_repeat.py`
  - Needs: T036
  - Wave: 14
  - Done when: `python -m pytest tests/fault_evaluation/test_benchmark_repeat.py -q` fails for the absent constants and repeat step, then passes after T038.

- [x] T038 [US3] [M] Repeat short cases in the timing pass.
  - Do: Implement the repeat step of D2 in `run.py`, with the constants named in T037. Use the checkpoint file and the resume rules that exist. Keep schema v1. Do not change tests, `bench_cases.py` or the other passes.
  - Files: `benchmarks/fault_evaluation/run.py`
  - Needs: T037
  - Wave: 15
  - Done when: `python -m pytest tests/fault_evaluation/test_benchmark_repeat.py tests/fault_evaluation/test_benchmark_suites.py tests/fault_evaluation/test_benchmark_run.py tests/fault_evaluation/test_benchmark_resume.py -q` passes.

## Added fault-set cache key tasks (replan of 2026-11-03)

The first T030 measurement failed SC-003, SC-004, SC-006, SC-007, SC-008 and
SC-010 (plan.md, Answers, "T030, 2026-11-03" and "T030 diagnosis"). The user
chose to fix the cache design. Plan D4a corrects D4: the persistent key comes
from the ordered fault set, not from the root feature and the whole model.
The done tasks T016 to T019 are not reopened: T040 changes two of their tests,
and T041 changes the code of T019 and T023. T043 gives the facts for the spec
decisions OP3 and OP4. T030 waits for T042, T043 and those decisions.

- [x] T040 [US2] [M] Test the fault-set cache key, the revision term and the publish guard.
  - Do: Write the D4a tests in a new file and change two existing tests. For the world-to-local transform case: if the coordinate owner accepts and rebuilds the change, assert the current-transform fingerprints change and cache-on restoration agrees with cache-off; if the owner rejects before mutation, assert model state, fingerprints, and results stay unchanged. Do not require both policies. This cache task only prevents stale reuse of the previous frame result. New tests: (1) `fault_set_fingerprint` gives equal bytes for the 4 groups of `f4_u3` (same ordered faults, different root data); a different fault order or subset gives different bytes; an opaque region on a fault gives None. (2) Cross-call sharing on `f4_u3` with the cache on: after `group_0.evaluate_value(points)`, a separate top-level `group_1.evaluate_value(points)` makes 0 fault calls; a cold `model.evaluate_model(points)` (cache cleared) makes exactly 4 fault calls (SC-003). (3) `RestorationCache.sync_revision`: the first call stores the revision and keeps entries put after it; the same revision keeps entries; a different revision clears all entries (`retained_bytes == 0`). (4) Publish guard: monkeypatch one fault's `apply_to_points` so that it bumps `model.evaluation_revision` (through a setter) and then calls the original; after the call the cache holds no entry for that key, and the result equals the cache-off result. (5) Transform policy: for an accepted and rebuilt transform change, the bounding-box transform changes both fingerprints and cache-on restoration agrees with cache-off. For a rejected-before-mutation change, the transform, fingerprints, and results remain unchanged. Test the policy selected by the coordinate owner; do not require both. (6) `_update_frame`: equal DataFrames give equal digests; a changed value, column name, dtype or object-column value gives a different digest. Tests (5) and (6) can pass before T041 only where the current code already reads that state; say which in the report. (7) Pickle defect (plan.md, Answers, "Pickle defect"): `dill.loads(dill.dumps(RestorationCache(enabled=False, max_bytes=1234)))` keeps `enabled` and `max_bytes`, has `retained_bytes == 0`, and its `lock.acquire(timeout=2)` returns True; a dill round trip of a built `f1_u0` model (`dill.dumps`/`dill.loads`, or `to_file`/`from_file` in `tmp_path`) then `evaluate_model(points)` returns the same ids as before the save. Run this test in a daemon thread with `join(timeout=60)` and fail if it is still alive, so that a hang fails and does not block pytest. CHANGES to done tests: see Context.
  - Files: `tests/fault_evaluation/test_fault_set_key.py`, `tests/fault_evaluation/test_invalidation.py`, `tests/fault_evaluation/test_persistent.py`
  - Needs: T023
  - Wave: 28
  - Done when: `python -m pytest tests/fault_evaluation/test_fault_set_key.py tests/fault_evaluation/test_invalidation.py tests/fault_evaluation/test_persistent.py -q` fails for the absent `fault_set_fingerprint` and `sync_revision`, for the cross-call fault count and for the lock after a dill round trip, then passes after T041.
  - Context: Code under test (do not change it): `LoopStructural/modelling/features/_evaluation_cache.py` (`RestorationCache`, `dependency_fingerprint(feature, build=True) -> bytes | None`, `_update_frame(hasher, frame)`, `_StateHasher`) and `_evaluation.py` (`EvaluationContext.restore`). New private names from T041: `fault_set_fingerprint(faults: Sequence[BaseFeature], build: bool = True) -> bytes | None` in `_evaluation_cache.py`; `RestorationCache.sync_revision(revision: int) -> None`. Persistent key after T041: `("persistent", digest, points.shape, points.dtype.str, bool(reverse), fingerprint, revision)`. Helpers: `tests.fault_evaluation.models.build_model`, `sample_points`; `benchmarks.fault_evaluation.instrumentation.EvaluationInstrumentation` (`snapshot()["fault_calls_total"]`); `model.evaluation_cache`. In `test_invalidation.py`: add the field `touches_faults: bool = True` to `Mutation` (after `effect`), and set it False for `region_add`, `region_remove_inplace`, `region_threshold`, `region_sign`. In `test_mutated_model_gives_a_fresh_changed_result`: always assert equal results and cached calls <= cache-off calls; when `touches_faults`, also assert cached calls > 0 (plan.md, Answers, "D4 tests after D4a"). In `test_persistent.py::test_opaque_dependency_bypasses_persistent_entries`: put the opaque region on a fault of `cutting_faults`, not on the root group; keep the hidden closure-change check. In `test_persistent.py::test_model_change_restores_again_and_matches_the_changed_reference`: assert `calls > 0` only when `mutation.touches_faults or mutation.revision`; always assert result equality. Change nothing else in those two files. Rules: pytest; NumPy docstrings and type hints (Ruff `D`, `ANN`); Python 3.9 target. Environment: `C:/Users/lachl/dev/LoopStructural/.venv/Scripts/python.exe`, PYTHONPATH = worktree root, `packages/loop_common/src`, `packages/loop_interpolation/src` (absolute).

- [x] T041 [US2] [M] Make the persistent key from the ordered fault set (D4a).
  - Do: Implement D4a. Fingerprint the current restoration transform and prevent reuse of a stale restoration result after its value changes. The coordinate owner controls whether a transform mutation is accepted and rebuilt or rejected before mutation; do not implement that owner policy here. In `_evaluation_cache.py`: add `fault_set_fingerprint` (with `build`, call `_build_reachable(fault)` for each fault in order; then one `_StateHasher`, write `f"S{len(faults)};"`, call `state.node(fault, True)` for each fault; a value that is not a feature or `_Opaque` gives None); make `_update_interpolator` also write the bounding-box transform; rewrite `_update_frame` to loop over `frame.items()` with no `iloc`; add `RestorationCache.sync_revision`. In `_evaluation.py`: key the persistent entry on the fault tuple and the revision, delete `_state_of_model` and `self._model_state`, and do not publish when the revision changes during the fault application. Do not change `dependency_fingerprint` or the scope key `_fault_state`. Do not change tests. No public API change.
  - Files: `LoopStructural/modelling/features/_evaluation_cache.py`, `LoopStructural/modelling/features/_evaluation.py`
  - Needs: T040
  - Wave: 29
  - Done when: `python -m pytest tests/fault_evaluation -q` passes, the golden tests take less than 120 s, and `ruff check LoopStructural/modelling/features/_evaluation_cache.py LoopStructural/modelling/features/_evaluation.py` shows no finding that was not there before.
  - Context: `_update_interpolator`: `box = getattr(interpolator, "bounding_box", None)`; with a box, `_encode(hasher, np.asarray(box.world_to_local_matrix, dtype=float))`; without a box, `hasher.update(b"Nbox;")`. `_update_frame` keeps the shape record, column names, dtypes and the `pd.util.hash_pandas_object(column, index=False)` hash for object columns. `sync_revision(revision)`: under `self._lock`, if `revision` differs from the stored value (start None), call `self._clear_locked()` and store it; docstring says it is for internal use. `EvaluationContext._fingerprint(faults: tuple, cache)` memoises in `self._fingerprints` by `tuple(id(f) for f in faults)` (store the faults with the digest, so the ids stay valid) and calls `fault_set_fingerprint(faults)` under `cache.lock`. In `restore`: read `revision = int(getattr(self.model, "evaluation_revision", 0))` AFTER the fingerprint (the fingerprint builds dirty builders), call `cache.sync_revision(revision)`, then use the key `("persistent", digest, points.shape, points.dtype.str, bool(reverse), fingerprint, revision)`. After the faults are applied, read the revision again; if it changed, keep the result in the scope only and do not call `cache.put`. Pickle defect (plan.md, Answers): add `RestorationCache.__getstate__` (keep `enabled` and `max_bytes`; drop the lock, the entries, the byte count and the stored revision) and `__setstate__` (make a new `threading.RLock`, empty entries, revision None). Rules: Ruff `D` and `ANN`, NumPy docstrings, Python 3.9 target; ASD-STE100 comments. Environment: as in T040.

- [x] T042 [US3] [M] Measure the D4a gates and write the profile report.
  - Do: SECOND RUN, after T045 (plan D4b). The first run (598ef21e, 2026-11-04) failed gate G4: `fault_set_fingerprint` took about 8 ms. Keep the first run in `profile-d4a.md` and add a new section for this run. Report again the OP5 ratio warm / cold `evaluate_model` for `f4_u3`, `f4_u0` and `onlap_domain`; the `onlap_domain` ratio (0.284 in the first run) has no gate here: if it stays above 0.10, report it, and T030 decides it (SC-010). Read-only for tracked code. On the clean committed tree, measure `f4_u3`, `f4_u0` and `onlap_domain` at 5,000 points, cold and warm, all five operations, with the D2 tools. Put any script in the scratchpad, outside the repository, and delete it at the end. Write `profile-d4a.md` with: fault calls (total) for each case; the cold `f4_u3` `model_ids` time with the cache on and with `model.evaluation_cache.enabled = False`, and their ratio; the warm `model_ids` fault calls; the median time of `fault_set_fingerprint` for the 4 faults of `f4_u3` (with `build=False`, 20 repeats); warm and cold `evaluate_model` times and their ratio (OP5); the `onlap_domain` `model_gradient` cold time against the baseline 0.92 s. Give the revision, the Python path and the commands. Mark each gate pass or fail. A miss is a failed task. If the SC-003 count is not 4, stop and report: the plan must change.
  - Files: `specs/003-fault-evaluation-performance/profile-d4a.md`
  - Needs: T045
  - Wave: 32
  - Done when: `profile-d4a.md` shows, for `f4_u3` at 5,000 points: cold `all_group_values` and `model_ids` each 4 fault calls; cold `model_ids` cache-on time at most 1.10 times cache-off time; warm `model_ids` 0 fault calls; the `fault_set_fingerprint` time (target at most 5 ms) and the `onlap_domain` number reported.
  - Context: Case ids are `"{model}:{operation}:{points}:{state}"`; operations `oldest_value`, `all_group_values`, `model_ids`, `oldest_gradient`, `model_gradient`. `benchmarks.fault_evaluation.run.enumerate_cases()` gives the full case list; filter it by model and `points == 5000`. `benchmarks.fault_evaluation.bench_cases.run_case(benchmark, model, points, case, rounds=...)` clears the cache, primes a warm case and times one call; the points come from `bench_cases._points(model, count)`. Count calls with `benchmarks.fault_evaluation.instrumentation.EvaluationInstrumentation` (`snapshot()["fault_calls_total"]`). `fault_set_fingerprint(faults, build=False)` is in `LoopStructural/modelling/features/_evaluation_cache.py`; the faults of a group are `feature.faults`. The `--suite` runner has no case filter: do not change it. Use `time.perf_counter`, the median of 5 rounds for each time. Environment: as in T040. Rule: every time is measured, none is inferred (plan.md, Answers 2026-11-02).

- [x] T043 [P] [US3] [M] Split the failed T030 results by cause for the OP3 and OP4 decisions.
  - Do: Read-only for tracked code. From the raw results in `C:/Users/lachl/dev/ls-bench-003/cand/` (`candidate.json`, `candidate_repeat.json`, `acceptance_cand.json`, `run1.log`, `run2.log`, the `run1.work` and `run2.work` checkpoints) and the spec baseline (`specs/003-fault-evaluation-performance/baseline.json`), write `t030-diagnosis.md`. Part (a), SC-006: a table of each case above 1.20, with the baseline peak extra RSS bytes, the candidate bytes, the ratio and the probable cause (cache retention, hash workspace, or other), and the count of failing cases with a baseline below 10 MB. Part (b), SC-008: the wall time of each run split by pass (timing, repeat of short cases, counters, memory), with the time for warm priming and for process starts where the data shows it, and the time of the two `onlap_domain` 100,000-point gradient cases in each pass. Say which numbers are measured and which are derived. Do not run the suite again, change the raw files or edit `plan.md` or `spec.md`.
  - Files: `specs/003-fault-evaluation-performance/t030-diagnosis.md`
  - Needs: T029
  - Wave: 15
  - Done when: `t030-diagnosis.md` has the SC-006 table with one row for each case of `acceptance_cand.json` above 1.20, and the SC-008 table whose pass times add up to the measured `suite_wall_seconds` of each run within 1 %.
  - Context: Result schema v1 (`benchmarks/fault_evaluation/run.py`): rows have `case_id`, `model`, `operation`, `points`, `state`, `rounds`, `time_samples_seconds`, median seconds, fault and field calls, peak extra RSS bytes, and the provenance fields; the document has `suite_wall_seconds`. The acceptance output (`benchmarks/fault_evaluation/acceptance.py`, `_sc006`, `_sc008`) has one row for each SC. The SC-006 limit is 1.20 and the SC-008 limit is 900 s. Known values: candidate 3,017 s, baseline 3,271 s, sum of candidate medians 334 s (plan.md, Answers, "T030 diagnosis"). These facts serve OP3 (SC-008 target) and OP4 (SC-006 for small baselines); do not propose a resolution as fact.

## Added faster frame digest tasks (replan of 2026-11-04)

T042 failed its key-time gate: `fault_set_fingerprint` for the 4 faults of
`f4_u3` took about 8 ms (target 5 ms; plan.md, Answers, "T042, 2026-11-04").
The cause is the pandas cost for each column of the 12 small builder frames.
Plan D4b changes only the body of the private `_update_frame`: a small frame is
read with one call. T040, T041 and T043 stay done. T042 runs again after T045.

- [x] T044 [US2] [M] Test the D4b frame digest for small and large frames and in-place edits.
  - Do: Change section (6) of the test file only. Add two helpers after `_frame()`: `_small_frame_cells() -> int` reads `_SMALL_FRAME_CELLS` from the module `_evaluation_cache` with `getattr(module, "_SMALL_FRAME_CELLS", None)` and asserts that it is not None, with the message "D4b needs _SMALL_FRAME_CELLS"; `_sized_frame(rows: int) -> pd.DataFrame` gives the columns of `_frame()` with `rows` rows (see Context). Add a size parameter `"small"` / `"large"` (a string, so that collection does not need the constant; compute the rows inside the test). Parametrize `test_equal_frames_have_equal_digests` and `test_changed_frame_has_another_digest` also over the size, and use `_sized_frame`. Keep the 4 changes (`value`, `column_name`, `dtype`, `object_value`); the value change and the object-value change edit row 0 of a copy. New `test_in_place_edit_changes_the_digest`, over the two sizes and two edits (`frame.iloc[1, 0] = 99.25`, `frame.loc[1, "name"] = "edited"`): take the digest, edit the same object, assert that the digest changed. New `test_object_values_of_other_types_have_other_digests`: two small frames with one column of dtype `object`, values `[1, "a"]` and `[1.0, "a"]`; the digests differ. New `test_small_frame_does_not_hash_each_column`: with `monkeypatch`, replace `pd.util.hash_pandas_object` and `pd.DataFrame.items` with functions that raise `AssertionError("per-column path used for a small frame")`; then `_frame_digest(_frame())` must end with no error. No time test. Change no other section and no code.
  - Files: `tests/fault_evaluation/test_fault_set_key.py`
  - Needs: T041
  - Wave: 30
  - Done when: `pytest tests/fault_evaluation/test_fault_set_key.py -q` fails only for `test_small_frame_does_not_hash_each_column` and for each test that calls `_small_frame_cells()` (the constant is absent), and `pytest tests/fault_evaluation -q` has no other failure; then all pass after T045.
  - Context: Existing helpers in section (6): `_frame_digest(frame: pd.DataFrame) -> bytes` (a `hashlib.sha256()` given to `_update_frame`, then `digest()`) and `_frame()` (`X` float64 `[1.0, 2.0, 3.0]`, `n` int64 `[1, 2, 3]`, `name` text `["a", "b", "c"]`). The file imports `_update_frame` from `LoopStructural.modelling.features._evaluation_cache`; import the module itself for `getattr` (for example `from LoopStructural.modelling.features import _evaluation_cache`). `_sized_frame(rows)`: `X = np.arange(rows) * 0.5`, `n = np.arange(rows, dtype=np.int64)`, `name = [f"s{i}" for i in range(rows)]`. Sizes: `small` = 3 rows; `large` = `_small_frame_cells() // 3 + 1` rows (above the limit). The large path keeps the column loop and `pd.util.hash_pandas_object` for object columns; the small path calls `repr(frame.to_numpy().tolist())` once. The current code (before T045) calls both patched functions for the text column, so the new no-column test is red. Rules: pytest; NumPy docstrings and type hints (Ruff `D`, `ANN`); Python 3.9 target; ASD-STE100 docstrings. Environment: as in T040.

- [x] T045 [US2] [S] Read a small DataFrame with one call in the frame digest (D4b).
  - Do: In `_evaluation_cache.py`, add the module constant `_SMALL_FRAME_CELLS = 2048` after `_CHUNK_BYTES`, with the comment: a frame with at most this number of cells is read with one call; below this size, the cost of each column is larger than the cost of each cell. In `_update_frame(hasher: _Updatable, frame: object) -> None`, keep the signature, the docstring purpose, `import pandas as pd` and the `raise _Opaque` check for a value that is not a `pd.DataFrame`. Replace the rest of the body with the code in Context, exactly. Do not cache a digest for each frame or each revision (OP5). Change no other function and no test. No public API change.
  - Files: `LoopStructural/modelling/features/_evaluation_cache.py`
  - Needs: T044
  - Wave: 31
  - Done when: `pytest tests/fault_evaluation -q` passes, the golden tests take less than 120 s, and `ruff check LoopStructural/modelling/features/_evaluation_cache.py` shows no finding that was not there before.
  - Context: New body after the type check: `dtypes = frame.dtypes.tolist()`; `texts = {dtype: str(dtype) for dtype in set(dtypes)}`; `hasher.update(f"D{frame.shape};".encode())`; `_encode(hasher, [repr(name) for name in frame.columns.tolist()])`; `_encode(hasher, [texts[dtype] for dtype in dtypes])`; `if frame.size <= _SMALL_FRAME_CELLS:` then the comment "One call reads all cells. ``repr`` of a float is exact, and it keeps the type of each value in an object column (1 is not 1.0 or "1").", then `hasher.update(b"V;")`, `_encode(hasher, repr(frame.to_numpy().tolist()))`, `return`. After the `if`: `for item in frame.items():` `column = item[1]` (user answer, plan.md, Answers, "T045, 2026-11-04": the form `for _, column` gives Ruff PERF102); `values = column.to_numpy()`; `if values.dtype.hasobject: values = pd.util.hash_pandas_object(column, index=False).to_numpy()`; `_update_array(hasher, values)`. The index stays out of the digest. Existing helpers in the file: `_encode(hasher, value)`, `_update_array(hasher, array)`, `_Opaque`. The digest bytes change; only equality in one process matters, so no stored file changes. Rules: Ruff `D` and `ANN`; Python 3.9 target; ASD-STE100 comments. Environment: as in T040.

## Added reduced-scope tasks (replan of 2026-11-18)

A review of the branch and of T030 changed the spec (C-012 to C-016): no
cache between calls and no lock, SC-006 passes at 120 % or +16 MiB, a shorter
standard suite, time ratios on `baseline.json`, 16 cold large cases. Plan D7
to D11 give the design. All done tasks keep `[x]`; T048 removes their cache
code. OP3 and OP4 are resolved, OP5 is closed, so T030 has no block now. In
these tasks, "Environment" means the Python environment of the orchestrator
notes, and "Ruff finds nothing new" means `ruff check LoopStructural packages`
shows no finding that the baseline of 474 findings does not have.

- [x] T046 [Found] [M] Resolve the master merge and move the value lookup to the geometry class (D8).
  - Do: The orchestrator has started `git merge --no-commit master`. The only conflict is in `packages/loop_common/src/loop_common/supports/_3d_structured_grid.py`. Master (#367, #369) moved `evaluate_value`, `trilinear`, `position_to_local_coordinates`, `position_to_dof_coefs` and `position_to_cell_corners` into `StructuredGrid3DGeometry` in `packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py`; master's `evaluate_value` also accepts a property array of shape `(n_nodes, k)`. Resolve it as in Context. Add two tests to `test_evaluation_lookup.py`: `test_lookup_value_matches_master_path_for_multi_column_property` (a `(n_nodes, 3)` property array; the fast path gives the same values as the master path) and `test_rectilinear_grid_uses_the_master_path` (a `RectilinearGrid` returns False from `_has_uniform_lookup`). Do not change other files of the merge, and do not commit, abort or reset the merge. If a golden test fails, stop and report: only C-006 can change a reference.
  - Files: `packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py`, `packages/loop_common/src/loop_common/supports/_3d_structured_grid.py`, `packages/loop_common/tests/test_evaluation_lookup.py`
  - Needs: T045
  - Wave: 32
  - Done when: `git diff --name-only --diff-filter=U` is empty; `pytest packages/loop_common/tests` and `pytest tests/fault_evaluation` pass (the golden tests in less than 120 s); the two new lookup tests pass; `git diff HEAD -- tests/fault_evaluation/data` is empty; Ruff finds nothing new.
  - Context: Add to `StructuredGrid3DGeometry` (private): `_has_uniform_lookup(self) -> bool` (True when `type(self)` does not override `position_to_cell_index`, `position_to_local_coordinates` or `position_to_cell_corners` of `StructuredGrid3DGeometry`); `_lookup_points(self, pos) -> tuple[np.ndarray, np.ndarray, np.ndarray]` (the branch body: one `position_to_cell_index` call; corners with -1 for outside rows; local coordinates clipped to 0..1); staticmethod `_values_from_coefficients(idc, inside, coefs, property_array)` (the master weighted sum, for 1-D and `(n_nodes, k)` arrays); `_evaluate_value_from_lookup(self, lookup, property_array) -> np.ndarray`. `evaluate_value` uses `_lookup_points` and `trilinear` when `_has_uniform_lookup()` is True, else the master path; keep the master checks of the property array shape. `StructuredGrid` (support file) keeps the master text plus the branch gradient methods `_gradient_weights`, `_evaluate_gradient_from_lookup`, the fast path in `evaluate_gradient`, and `get_element_gradient_for_location` through `_gradient_weights`; its `_has_uniform_lookup` overrides the geometry method and also requires that `get_element_gradient_for_location` is the method of `StructuredGrid`. Move the branch value methods out of the support file. `_fault_segment._shared_frame_members` stays as it is (`isinstance(support, StructuredGrid)` and `support._has_uniform_lookup()`). The master tests `test_structured_support_reference.py`, `test_property_grid.py` and `test_bounding_box_structured_grid.py` must pass. Rules: NumPy docstrings and type hints (Ruff `D`, `ANN`); Python 3.9 target; ASD-STE100 comments. Environment: as in the orchestrator notes.

- [x] T047 [US2] [M] Change the evaluation tests for no cache between calls (D7, red first).
  - Do: Delete `test_cache_storage.py`, `test_persistent.py`, `test_invalidation.py`, `test_fault_set_key.py` and `test_fingerprint_cycle.py` (each tests only removed code). In `test_concurrency.py`, remove the tests of the model lock and of the cache (for example `test_held_lock_blocks_evaluation_also_with_the_cache_disabled`, `test_nested_and_repeated_calls_reuse_the_lock_in_one_thread`, `test_thread_holding_the_lock_of_one_model_does_not_block_another_model`, `test_failed_call_publishes_no_entry`, and the lock and entry parts of the other tests); keep the tests that compare two threads with serial references; add `test_two_calls_on_one_model_run_at_the_same_time` (Context). In `test_memory.py`, remove the retained-bytes, cache-mode and fingerprint tests; keep the call-only FR-010 accounting and the bounded digest workspace tests of `_digest`; add `test_no_restoration_stays_after_the_call` (Context). In `test_context.py`, `test_gradient_context.py` and `test_fault_frame_lookup.py`, remove only the lines that set or read `evaluation_cache`; do not change the asserted counts. Add the new file `test_no_retention.py` (Context). Do not change production code or other tests.
  - Files: `tests/fault_evaluation/test_cache_storage.py`, `tests/fault_evaluation/test_persistent.py`, `tests/fault_evaluation/test_invalidation.py`, `tests/fault_evaluation/test_fault_set_key.py`, `tests/fault_evaluation/test_fingerprint_cycle.py`, `tests/fault_evaluation/test_concurrency.py`, `tests/fault_evaluation/test_memory.py`, `tests/fault_evaluation/test_context.py`, `tests/fault_evaluation/test_gradient_context.py`, `tests/fault_evaluation/test_fault_frame_lookup.py`, `tests/fault_evaluation/test_no_retention.py`
  - Needs: T046
  - Wave: 33
  - Done when: `pytest tests/fault_evaluation/test_concurrency.py tests/fault_evaluation/test_memory.py tests/fault_evaluation/test_context.py tests/fault_evaluation/test_gradient_context.py tests/fault_evaluation/test_fault_frame_lookup.py tests/fault_evaluation/test_no_retention.py` fails for exactly 3 tests, for the reasons in Context, then passes after T048.
  - Context: `test_concurrency.py` and `test_memory.py` import helpers from `tests.fault_evaluation.test_invalidation` (for example `_fault`, `_group`, `_points`), and `test_memory.py` imports from `_evaluation_cache`. Copy each helper that a kept test still uses into the file that uses it, and remove these imports. `_digest` is in `_evaluation.py` and stays. `test_two_calls_on_one_model_run_at_the_same_time`: patch `FaultSegment.apply_to_points` with a wrapper that calls `threading.Barrier(2, timeout=10).wait()` once for each thread and then calls the original; two threads call `model.evaluate_model(points)` on the golden model `f1_u0`; the test fails if the barrier breaks (`threading.BrokenBarrierError`); both results equal the serial result. `test_no_restoration_stays_after_the_call`: after `evaluate_model`, `_evaluation._ACTIVE.get()` is None, and after `gc.collect()` no object in `gc.get_objects()` is an `EvaluationContext`. `test_no_retention.py` (FR-009, user story 2 scenario 6): for each change, evaluate `evaluate_model(points)`, make the change, evaluate again, build a second model with the same change from the start, and assert equal results within the FR-004 tolerance of `tests.fault_evaluation.golden`; the changes are `fault.displacement *= 1.5` on the first fault of `f1_u0`, `group.faults_enabled = False` on the oldest group of `f1_u1`, and `fault.steps = 5` on the first fault of `f1_u0`. No transform test (spec 005 owns it, its SC-016). `test_repeat_call_applies_the_faults_again`: two equal `evaluate_model(points)` calls on `f1_u0` make the same `snapshot()["fault_calls_total"]` count (`benchmarks.fault_evaluation.instrumentation.EvaluationInstrumentation`) and equal results. `test_model_has_no_cache_attributes`: `GeologicalModel` has no attribute `evaluation_cache` or `evaluation_revision`. Helpers: `tests.fault_evaluation.models.build_model(name)`, `sample_points(name, count, seed=1)`. Expected red before T048: the barrier test (the lock serialises the calls), the repeat test (the cache gives 0 fault calls in the second call) and the attribute test; the three change tests pass already. Rules: pytest; NumPy docstrings and type hints (Ruff `D`, `ANN`); Python 3.9 target. Environment: as in the orchestrator notes.

- [x] T048 [US2] [M] Remove the cache between calls and the model lock (D7).
  - Do: Implement D7 (Context). Delete `_evaluation_cache.py`. Keep the scope for each call (D3), the gradient batch sharing (D5) and the fault-frame lookup (D6). Do not change tests or `benchmarks/`: T049 and T050 change the benchmark files and their tests. No change of results. The removed names were never released, so the released API does not change.
  - Files: `LoopStructural/modelling/features/_evaluation_cache.py`, `LoopStructural/modelling/features/_evaluation.py`, `LoopStructural/modelling/core/geological_model.py`, `LoopStructural/modelling/features/_base_geological_feature.py`, `LoopStructural/modelling/features/builders/_base_builder.py`, `LoopStructural/modelling/features/builders/_geological_feature_builder.py`, `LoopStructural/modelling/features/fault/_fault_segment.py`, `packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py`
  - Needs: T047
  - Wave: 34
  - Done when: `pytest tests/fault_evaluation --ignore=tests/fault_evaluation/test_benchmark_run.py --ignore=tests/fault_evaluation/test_benchmark_suites.py --ignore=tests/fault_evaluation/test_acceptance.py --ignore=tests/fault_evaluation/test_benchmark_compare.py` passes (the golden tests in less than 120 s); `pytest packages/loop_interpolation/tests` passes; `git diff master -- LoopStructural/modelling/features/builders/_base_builder.py LoopStructural/modelling/features/builders/_geological_feature_builder.py packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py` is empty; `git grep -n -E "RestorationCache|fault_set_fingerprint|dependency_fingerprint|bump_revision|values_differ|_SMALL_FRAME_CELLS|evaluation_cache|evaluation_revision|_revision_hook" -- LoopStructural packages` finds nothing; Ruff finds nothing new.
  - Context: `_evaluation.py`: remove the import of `fault_set_fingerprint` and of the `RestorationCache` type; `EvaluationContext.__init__(self, model=None)` keeps `self.model` and `self._restorations: dict`, and loses `_fingerprints`, `_published`, `_opaque`; `restore(feature, points, reverse=False)` keeps the `faults_enabled` check, the read-only view for no faults, the `builder.up_to_date()` loop, `_digest`, the scope key `(digest, points.shape, points.dtype.str, tuple(_fault_state(f) for f in faults), bool(reverse))` (remove the last literal `True`), the scope lookup, the fault loop, the copy when the result shares memory with `points`, `flags.writeable = False` and the store after all faults; remove `_persistent_cache`, `_fingerprint`, `_withdraw` and the persistent path; `evaluation_scope(model=None)` keeps the nested-scope rule and becomes `context = EvaluationContext(model); token = _ACTIVE.set(context); try: yield context; finally: _ACTIVE.reset(token)`; the module docstring says the scope keeps entries for one call only. `geological_model.py`: remove the `RestorationCache` import, `_evaluation_revision`, `_evaluation_cache`, the properties `evaluation_revision` and `evaluation_cache`, `_bump_evaluation_revision` and each call; keep each `with evaluation_scope(self):`. `_base_geological_feature.py`: remove `bump_revision`; restore the plain attribute `self.faults_enabled = True` in `__init__` (no property); keep `evaluation_scope` in `evaluate_normalised_value`. The two builder files: remove `bump_revision` and the `_revision_hook` assignment. `_fault_segment.py`: remove `bump_revision` and `values_differ`, the `displacement` and `steps` properties (plain attributes again), the `bump_revision` calls in `set_displacement` and the `try`/`finally` of the `faultfunction` setter (master body); keep `_shared_frame_members`, `_evaluate_frame`, `_evaluate_frame_legacy` and the changed `apply_to_points`. `_discrete_interpolator.py`: remove the class attribute `_revision_hook` and its call in `solve_system`. Use `git show master:<path>` to see the master text. Rules: NumPy docstrings and type hints (Ruff `D`, `ANN`); Python 3.9 target; `GeologicalModel` public signatures stay; ASD-STE100 comments. Environment: as in the orchestrator notes.

- [x] T049 [US3] [M] Test the shorter suite and the new acceptance rules (D9, D10, red first).
  - Do: Change the four benchmark test files for D9 and D10 (Context). Remove the tests of warm cases, of the priming call, of `evaluation_cache`, of `_cache_fields`, of SC-010 and of the `"warm"` and `"off"` memory modes. Change the case counts and the fixtures that the new rules make wrong, and say which in the report. Add the new tests of Context. Do not change production code or other tests.
  - Files: `tests/fault_evaluation/test_benchmark_run.py`, `tests/fault_evaluation/test_benchmark_suites.py`, `tests/fault_evaluation/test_acceptance.py`, `tests/fault_evaluation/test_benchmark_compare.py`
  - Needs: T046
  - Wave: 33
  - Done when: `pytest tests/fault_evaluation/test_benchmark_run.py tests/fault_evaluation/test_benchmark_suites.py tests/fault_evaluation/test_acceptance.py tests/fault_evaluation/test_benchmark_compare.py` fails only in new or changed tests, for the absent D9 and D10 names and rules; the tests that this task did not touch pass; then all pass after T050.
  - Context: Module `benchmarks.fault_evaluation.run`: `STATES == ("cold",)`; `enumerate_cases()` gives 430 cases; `standard_cases()` gives 131 (105 at 5,000 points, 10 at 20,000, 9 at 50,000, 7 at 100,000); `STANDARD_LARGE_CASES` is exactly the 16 ids of plan D9 in this order: `f0_u0:all_group_values:100000:cold`, `f0_u1:oldest_value:100000:cold`, `f0_u1:model_ids:50000:cold`, `f0_u2:oldest_value:50000:cold`, `f0_u3:model_ids:100000:cold`, `f1_u0:all_group_values:100000:cold`, `f1_u1:oldest_value:50000:cold`, `f1_u2:oldest_value:100000:cold`, `f1_u2:oldest_gradient:50000:cold`, `f1_u3:oldest_value:50000:cold`, `f1_u3:all_group_values:100000:cold`, `stopped_fault:model_ids:50000:cold`, `stopped_fault:oldest_gradient:50000:cold`, `onlap_domain:model_ids:50000:cold`, `faulted_fold:oldest_value:50000:cold`, `faulted_fold:all_group_values:100000:cold`; no `onlap_domain` gradient case at 100,000 points in the standard list. `PASS_POINTS == {"timing": None, "counters": (5_000, 20_000), "rss": (5_000,)}`; `pass_cases(cases, pass_kind, suite)`: `"full"` gives `cases`, `"standard"` filters by `points` (None = all), ValueError for an unknown pass or suite; `main` gives each `_run_fresh_process_pass` call the `pass_cases` list of its pass. `_merge`: a case absent from the counter pass gets None for each of `COUNTER_FIELDS`; a case absent from the memory pass gets `peak_extra_rss_bytes` None and `rss_unit` `"bytes"`; top-level `"suite"`. `validate_result` accepts None only for `COUNTER_FIELDS` and `peak_extra_rss_bytes`, and an optional `"suite"`; `SCHEMA_VERSION` stays 1. `cold_view(document)`: a copy where a `:warm` row with no `:cold` row of the same case becomes `:cold` (`case_id` and `state`), and other `:warm` rows go. `bench_cases.run_case` never primes and never reads `evaluation_cache`; `warmup_rounds` stays 0. Module `acceptance`: `check_acceptance(base, new)` applies `cold_view` to `base`; no `"SC-010"` row; `MEMORY_MODES == ("cold",)`; `MEMORY_ALLOWANCE_BYTES == 16 * 1024 * 1024`; SC-006 passes a case when `new <= 1.20 * base` or `new - base <= MEMORY_ALLOWANCE_BYTES` (1.5 x a 1 MB baseline passes; 1.5 x a 100 MB baseline fails; a 0 baseline uses only the allowance); a candidate None that is not in the memory pass is not a failure; a criterion that needs a None value fails with the reason "case X was not measured in the Y pass". Module `compare`: a base with `:warm` ids goes through `cold_view`; a None value gives the ratio text "not measured", no failure and no noise mark. Rules: pytest; NumPy docstrings and type hints (Ruff `D`, `ANN`); Python 3.9 target. Environment: as in the orchestrator notes.

- [x] T050 [US3] [M] Implement the shorter suite and the new acceptance rules (D9, D10).
  - Do: Implement D9 and D10 so that the T049 tests pass (Context). Keep the comment on the draw of 2026-11-02 above `STANDARD_LARGE_CASES` and add the rule of C-016 (merge the 24 old ids without the state, remove the 5 gradient cases at 100,000 points). The `cold_view` docstring says: use it only on a file of the code before this feature (T011); that code had no cache and the runner skipped the priming call, so its warm case was one cold call (C-016). Keep `_time_ratio_row` and `_sc009` on the baseline file (C-015) and change no limit. Do not change tests or files outside `benchmarks/fault_evaluation/`.
  - Files: `benchmarks/fault_evaluation/run.py`, `benchmarks/fault_evaluation/bench_cases.py`, `benchmarks/fault_evaluation/acceptance.py`, `benchmarks/fault_evaluation/compare.py`
  - Needs: T049, T048
  - Wave: 35
  - Done when: `pytest tests/fault_evaluation` passes; `git grep -n -E "RestorationCache|fault_set_fingerprint|dependency_fingerprint|bump_revision|values_differ|_SMALL_FRAME_CELLS|evaluation_cache|evaluation_revision|_revision_hook" -- LoopStructural packages benchmarks tests` finds nothing, except the test lines that assert that `evaluation_cache` is absent or never read (`test_benchmark_run.py`, `test_benchmark_suites.py`, `test_no_retention.py`; plan.md, Answers, 2026-11-18, T050); Ruff finds nothing new and `ruff check benchmarks/fault_evaluation` shows no finding that was not there before.
  - Context: `run.py`: `STATES = ("cold",)` (ids keep `:cold`; `state` stays in schema v1); the 16 `STANDARD_LARGE_CASES` ids of plan D9 (also in the Context of T049), in the order of the full list; `PASS_POINTS: dict[str, tuple[int, ...] | None] = {"timing": None, "counters": (5_000, 20_000), "rss": (5_000,)}`; `pass_cases(cases: list[dict[str, Any]], pass_kind: str, suite: str) -> list[dict[str, Any]]`; `main` calls `_run_fresh_process_pass(pass_kind, pass_cases(cases, pass_kind, options.suite), ...)`; `_merge` gives None for each of `COUNTER_FIELDS` of a case that the counter pass did not run, None for `peak_extra_rss_bytes` (with `rss_unit` `"bytes"`) of a case that the memory pass did not run, and the top-level key `"suite"`; `validate_result` accepts None only for those fields, and an optional `"suite"`; `SCHEMA_VERSION` stays 1; `cold_view(document: dict[str, Any]) -> dict[str, Any]`. `bench_cases.py`: remove the warm path of `run_case` (cache clear, priming call, `evaluation_cache`); keep `warmup_rounds=0`. `acceptance.py`: `check_acceptance` applies `run.cold_view` to `base`; remove `_sc010`, `IDS_COLD_F4U3`, `IDS_WARM_F4U3`, `WARM_F4U3`, the `"SC-010"` row, `_cache_fields`, the cache lines of `_memory_worker`; `MEMORY_MODES = ("cold",)`; `MEMORY_ALLOWANCE_BYTES = 16 * 1024 * 1024`; `_sc006(base, new)`: compare cases where both peaks are not None; pass at `new <= 1.20 * base` or `new - base <= MEMORY_ALLOWANCE_BYTES`; the row value is the worst ratio, the reason lists each failing case with both values and the row gives the number of compared cases as the integer key `"compared_cases"` (cases with a memory value in both documents; plan.md, Answers, 2026-11-18, review of wave 33); `_check_cases` accepts None counter and memory fields when the criterion does not read them; a needed None value fails with "case X was not measured in the Y pass". `compare.py`: apply `run.cold_view` to a base with `:warm` ids; None on either side gives "not measured", never a failure or a noise mark. Rules: stdlib only in `compare.py`; NumPy docstrings and type hints (Ruff `D`, `ANN`); Python 3.9 target. Environment: as in the orchestrator notes.

- [x] T051 [US3] [M] Measure the gates on the merged tree with no cache (D11 step 4).
  - Do: Read-only for tracked code. On the clean committed tree after T048 and T050, measure `f4_u3`, `f4_u0` and `onlap_domain` at 5,000 points, cold, all five operations, with the D2 tools. Put any script in the scratchpad, outside the repository, and delete it at the end. Add a new top section to `profile-d4a.md` ("Third run (<revision>, no cache, after the master merge)"), above the second run, and keep the old sections. Report: G1, fault calls of `f4_u3` cold `model_ids` (target 4, SC-003 after spec C-017), and the `all_group_values` count for information (16 expected: 4 separate calls, no reuse between calls); G2 (new), the cold `f4_u3` `model_ids` median time (target at most 1.10 times 0.2208 s, the cache-on time of the second run of T042); G5, the `onlap_domain` `model_gradient` cold time against the baseline 0.92 s (report only); the fault calls and median times of every measured case; and that `pytest tests/fault_evaluation` (the golden tests) passes on this tree, with its time. G3 and G4 are gone (no cache, no fingerprint). Give the revision, the Python path and the commands. Mark each gate pass or fail. A miss is a failed task. If the `model_ids` count is not 4, stop and report: the plan must change. The first try (2026-11-18, revision 3865c80b) stopped on the `all_group_values` count of 16; the user then chose spec C-017 (plan.md, Answers).
  - Files: `specs/003-fault-evaluation-performance/profile-d4a.md`
  - Needs: T048, T050
  - Wave: 36
  - Done when: `profile-d4a.md` has the third-run section with the revision, G1 (4 for `model_ids`), G2 at most 0.2429 s, the G5 time, and the golden test result, each marked pass or fail.
  - Context: Case ids are `"{model}:{operation}:{points}:cold"`; operations `oldest_value`, `all_group_values`, `model_ids`, `oldest_gradient`, `model_gradient`. `benchmarks.fault_evaluation.run.enumerate_cases()` gives the full case list; filter it by model and `points == 5000`. `benchmarks.fault_evaluation.bench_cases.run_case(benchmark, model, points, case, rounds=...)` times one call; the points come from `bench_cases._points(model, count)`. Count calls with `benchmarks.fault_evaluation.instrumentation.EvaluationInstrumentation` (`snapshot()["fault_calls_total"]`). The runner has no case filter: do not change it. Use `time.perf_counter` and the median of 5 rounds for each time. Rule: every time is measured, none is inferred (plan.md, Answers 2026-11-02). Environment: as in the orchestrator notes.

- [ ] T052 [US3] [M] Test SC-003 and SC-004 on one evaluate_model call (spec C-017, red first).
  - Do: In `tests/fault_evaluation/test_acceptance.py`, change the fixtures and tests of SC-003 and SC-004 for spec C-017 (Context). Add the rows of `f4_u3:model_ids:5000:cold` and `f4_u0:model_ids:5000:cold` to the base and candidate fixtures if they are absent. Keep the SC-003 time test on `all_group_values`. Change the SC-003 count test so that it sets the count of the `model_ids` case. Add: an `all_group_values` count of 16 does not fail SC-003 when the `model_ids` count is 4; a `model_ids` count of 5 fails SC-003 and the reason names the `model_ids` id; SC-004 uses the `model_ids` times (a slow `f4_u3:all_group_values` does not fail SC-004; a `model_ids` ratio above 1.5 fails it); a `model_ids` case absent from the candidate fails SC-003 and SC-004 with the id in the reason. Do not change other tests or production code.
  - Files: `tests/fault_evaluation/test_acceptance.py`
  - Needs: T050
  - Wave: 36
  - Done when: `pytest tests/fault_evaluation/test_acceptance.py -q` fails only in the new or changed SC-003 and SC-004 tests, because `acceptance.py` still reads `all_group_values` for the count and for SC-004; then all pass after T053.
  - Context: `benchmarks/fault_evaluation/acceptance.py` has `ALL_F4U3 = "f4_u3:all_group_values:5000:cold"`, `ALL_F4U0 = "f4_u0:all_group_values:5000:cold"`, `_time_ratio_row(base, new, case_id, limit, calls)` (time ratio and exact count on one case) and `_sc004(new)` (ratio of the two `all_group_values` medians, limit 1.5). After T053: new constants `IDS_F4U3 = "f4_u3:model_ids:5000:cold"` and `IDS_F4U0 = "f4_u0:model_ids:5000:cold"`; SC-003 = time ratio of `ALL_F4U3` at most 0.15 of the baseline AND `fault_calls_total` of `IDS_F4U3` equal to 4; SC-004 = median of `IDS_F4U3` / median of `IDS_F4U0` at most 1.5. Measured values for realistic fixtures (T051 first try, 2026-11-18): `model_ids` f4_u3 0.192 s and 4 calls, f4_u0 0.186 s and 4 calls; `all_group_values` f4_u3 0.733 s and 16 calls. Baseline (`baseline.json`): `model_ids` f4_u3 6.34 s and 52 calls, f4_u0 0.48 s. Use the helpers of the file (`_set`, `_row`, the base and candidate fixtures). Rules: pytest; NumPy docstrings and type hints (Ruff `D`, `ANN`); Python 3.9 target. Environment: as in the orchestrator notes.

- [ ] T053 [US3] [M] Measure the SC-003 count and SC-004 on one evaluate_model call (spec C-017).
  - Do: In `acceptance.py`, add `IDS_F4U3` and `IDS_F4U0` (Context of T052). Change the SC-003 row so that the time ratio stays on `ALL_F4U3` and the count check reads `IDS_F4U3`; a missing or unmeasured `IDS_F4U3` fails with its id in the reason (the existing "was not measured in the counters pass" form). Change `_sc004` to the two `model_ids` cases. Change the docstrings. Do not change tests or other files.
  - Files: `benchmarks/fault_evaluation/acceptance.py`
  - Needs: T052
  - Wave: 37
  - Done when: `pytest tests/fault_evaluation -q` passes; `ruff check benchmarks/fault_evaluation` is clean.
  - Context: see the Context of T052. Keep the limits (0.15, 4 calls, 1.5) and the row form `_row(passed, value, limit, reason)`. The value of the SC-003 row stays the time ratio. Rules: NumPy docstrings and type hints (Ruff `D`, `ANN`); Python 3.9 target; ASD-STE100 comments. Environment: as in the orchestrator notes.

## Dependencies

T001 is complete. T033 depends on it, proves the existing gradient output
defect without a solver, and must fail before T034 repairs it. T002 waits for
T034 and resumes its saved drafts. T003-T006 then capture and check all golden
outputs on the corrected pre-performance revision. The benchmark tools T007-T010
exist. T035 -> T036 add the standard suite and T037 -> T038 add the repeat of
short cases; T039 (acceptance tests) runs beside them and T029 needs it. T011
freezes both baseline runs of the standard suite before any cache or grid
optimization, so it needs T038. D0 is the sole earlier evaluation edit. The
accepted OP1 batch rule applies to gradients. OP2 is resolved by measurement; a
failed runtime/noise target triggers a user decision and never counts as
passing acceptance. Cache and grid tasks keep their prior dependency branches
and start after T011. The first T030 run failed (2026-11-03). T040 -> T041 fix the D4a fault-set key;
T044 -> T045 make the frame digest faster (D4b), and T042 then measures D4a; T043 splits the failed results by
cause and can run at any time after T029. The replan of 2026-11-18 adds
T046 (merge master, D8), then T047 -> T048 (remove the cache, D7) beside
T049 -> T050 (shorter suite, acceptance rules, D9, D10; T050 also needs T048),
then T051 (gates on the merged tree). T052 -> T053 move the SC-003 count and
SC-004 to `model_ids` (spec C-017). T030 joins T042, T043, T051 and T053; OP3 and
OP4 are resolved, so it has no block. T031 documents passing
results, and T032 is the final scheduled task and depends on every other task.
Spec 005 owns frame/world-local repairs; its planning gate waits for these documentation amendments only, so no implementation dependency on 005 is
added here. If a future reviewed correction changes references, append its
provenance and keep current historical records. Keep performance limits and
unaffected baselines fixed; compare equivalent correctness and run provenance.

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
| 12 | T028, T035 | 1 | 1 |
| 13 | T036, T039 | 0 | 2 |
| 14 | T029, T037 | 0 | 2 |
| 15 | T038, T043 | 0 | 2 |
| 16 | T011 | 0 | 1 |
| 17 | T012, T024 | 0 | 2 |
| 18 | T013, T025 | 0 | 2 |
| 19 | T014 | 0 | 1 |
| 20 | T015 | 0 | 1 |
| 21 | T016 | 0 | 1 |
| 22 | T017 | 0 | 1 |
| 23 | T018 | 0 | 1 |
| 24 | T019 | 0 | 1 |
| 25 | T020, T022, T026 | 0 | 3 |
| 26 | T021, T027 | 0 | 2 |
| 27 | T023 | 0 | 1 |
| 28 | T040 | 0 | 1 |
| 29 | T041 | 0 | 1 |
| 30 | T044 | 0 | 1 |
| 31 | T045 | 1 | 0 |
| 32 | T042, T046 | 0 | 2 |
| 33 | T047, T049 | 0 | 2 |
| 34 | T048 | 0 | 1 |
| 35 | T050 | 0 | 1 |
| 36 | T051, T052 | 0 | 2 |
| 37 | T053 | 0 | 1 |
| 38 | T030 | 0 | 1 |
| 39 | T031 | 0 | 1 |
| 40 | T032 | 0 | 1 |

## Critical path

T001 -> T033 -> T034 -> T002 -> T003 -> T004 -> T005 -> T006 -> T007 -> T008 -> T010 -> T035 -> T036 -> T037 -> T038 -> T011 -> T012 -> T013 -> T014 -> T015 -> T016 -> T017 -> T018 -> T019 -> T020 -> T021 -> T023 -> T040 -> T041 -> T044 -> T045 -> T046 -> T047 -> T048 -> T050 -> T052 -> T053 -> T030 -> T031 -> T032

40 sequential waves. T001-T029 and T033-T050 are complete. T051 failed its first try on the SC-003 count (spec C-017 answers it); T052 and T053 move the SC-003 count and SC-004 to `model_ids`. No task is blocked.
