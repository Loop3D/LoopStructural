# Profile report: D4a gates (T042)

Result: FAILED. One gate misses its target (G4, `fault_set_fingerprint` time).

## Setup

- Revision: 598ef21e98861875717d34fbe98e304d3ee46374 (tree clean).
- Python: `C:/Users/lachl/dev/LoopStructural/.venv/Scripts/python.exe` (3.13.12).
- Machine: Windows 11 (10.0.26200), Intel64 Family 6 Model 154, 20 logical CPUs. No other heavy work ran.
- PYTHONPATH: worktree root; `packages/loop_common/src`; `packages/loop_interpolation/src`.
- Method: a script in the scratchpad (outside the repository, deleted at the end). It called
  `benchmarks.fault_evaluation.bench_cases.run_case` with a small fake `benchmark` object.
  The fake object has `extra_info = {}` and `pedantic(call, setup, rounds, warmup_rounds)`.
  `pedantic` runs `setup()`, then times `call()` with `time.perf_counter`, for each round.
  So each time is exactly the time that `run_case` times.
- Cases: `enumerate_cases()` filtered by model (`f4_u3`, `f4_u0`, `onlap_domain`) and `points == 5000`. 30 cases.
- Time: 5 rounds for each case. The table gives the median.
- Fault calls: one extra run of the case with `EvaluationInstrumentation` and `rounds=1`.
  `after_setup=counters.reset`. The value is `snapshot()["fault_calls_total"]`.
  The timing runs did not use the instrumentation.
- Each model was built once with `tests.fault_evaluation.models.build_model`, outside the timed part.
- Command: `python m.py` (the scratchpad script), run from the scratchpad folder.
- Every time is measured. None is inferred.

## Fault calls and median time, 5,000 points

| Case | Fault calls | Median (s) |
|---|---|---|
| f4_u3 oldest_value cold | 4 | 0.2389 |
| f4_u3 oldest_value warm | 0 | 0.0147 |
| f4_u3 all_group_values cold | 4 | 0.2565 |
| f4_u3 all_group_values warm | 0 | 0.0521 |
| f4_u3 model_ids cold | 4 | 0.2378 |
| f4_u3 model_ids warm | 0 | 0.0259 |
| f4_u3 oldest_gradient cold | 8 | 1.7519 |
| f4_u3 oldest_gradient warm | 0 | 0.3873 |
| f4_u3 model_gradient cold | 8 | 2.8562 |
| f4_u3 model_gradient warm | 0 | 1.5138 |
| f4_u0 oldest_value cold | 4 | 0.2324 |
| f4_u0 oldest_value warm | 0 | 0.0098 |
| f4_u0 all_group_values cold | 4 | 0.2667 |
| f4_u0 all_group_values warm | 0 | 0.0108 |
| f4_u0 model_ids cold | 4 | 0.2618 |
| f4_u0 model_ids warm | 0 | 0.0104 |
| f4_u0 oldest_gradient cold | 4 | 1.7437 |
| f4_u0 oldest_gradient warm | 0 | 0.4822 |
| f4_u0 model_gradient cold | 4 | 1.5863 |
| f4_u0 model_gradient warm | 0 | 0.4105 |
| onlap_domain oldest_value cold | 1 | 0.0334 |
| onlap_domain oldest_value warm | 0 | 0.0084 |
| onlap_domain all_group_values cold | 1 | 0.0425 |
| onlap_domain all_group_values warm | 0 | 0.0144 |
| onlap_domain model_ids cold | 1 | 0.0356 |
| onlap_domain model_ids warm | 0 | 0.0101 |
| onlap_domain oldest_gradient cold | 2 | 0.3343 |
| onlap_domain oldest_gradient warm | 0 | 0.2292 |
| onlap_domain model_gradient cold | 2 | 0.5555 |
| onlap_domain model_gradient warm | 0 | 0.4454 |

Samples (s) of the key cases:

- f4_u3 model_ids cold, cache on: 0.2670, 0.2728, 0.2378, 0.2290, 0.2316.
- f4_u3 model_ids cold, cache off: 0.2223, 0.2196, 0.3192, 0.2194, 0.2355 (fault calls 4).
- f4_u3 model_ids warm: 0.0303, 0.0241, 0.0259, 0.0257, 0.0267.
- onlap_domain model_gradient cold: 0.5555, 0.5486, 0.5582, 0.5474, 0.5564.

An earlier run of the same script gave a cache-off median of 0.2547 s (samples 0.2547, 0.3171, 0.2372, 0.3172, 0.2147).
The cache-on and cache-off times are close, and the noise between rounds is about 10 to 40 percent.
Both runs give a ratio below 1.10 (0.93 and 1.07).

## Gates

| Gate | Measured | Target | Result |
|---|---|---|---|
| G1 SC-003 count: f4_u3 cold all_group_values and model_ids | 4 and 4 | 4 each | pass |
| G2 cold model_ids, cache on / cache off (f4_u3) | 0.2378 / 0.2223 = 1.07 | at most 1.10 | pass |
| G3 warm model_ids fault calls (f4_u3) | 0 | 0 | pass |
| G4 `fault_set_fingerprint` median, 4 faults of f4_u3, `build=False`, 20 repeats | 7.57 to 8.07 ms (group_0 8.07, group_1 8.02, group_2 7.96, group_3 7.57) | at most 5 ms | FAIL |
| G5 onlap_domain model_gradient cold | 0.5555 s | below the baseline 0.92 s | pass |

Fingerprint minimum per call was 6.6 to 6.9 ms. It did not reach 5 ms in any group.

## OP5: warm evaluate_model / cold evaluate_model

| Model | Cold (s) | Warm (s) | Ratio warm / cold |
|---|---|---|---|
| f4_u3 | 0.2378 | 0.0259 | 0.109 |
| f4_u0 | 0.2618 | 0.0104 | 0.040 |
| onlap_domain | 0.0356 | 0.0101 | 0.284 |

## Verdict

Gate G4 fails: the fingerprint costs about 8 ms, not 5 ms or less. The task is failed.
The warm `f4_u3` `model_ids` time (25.9 ms) holds about 8 ms of fingerprint time.
