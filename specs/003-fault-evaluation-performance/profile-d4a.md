# Profile report: D4a gates (T042)

Overall result: PASSED in the second run (45980adc, after D4b). All five gates pass.
The first run (598ef21e) failed gate G4. G4 now measures 3.50 ms at most (target 5 ms).
The `onlap_domain` OP5 ratio is 0.254. It is above 0.10 and has no gate here. T030 decides it (SC-010).

## Second run (45980adc, after D4b)

### Setup

- Revision: 45980adc34cbe13acb6a00ab4889860ed572a942 (`git status --short` printed nothing before and after).
- Python: `C:/Users/lachl/dev/LoopStructural/.venv/Scripts/python.exe` (3.13.12).
- PYTHONPATH: worktree root; `packages/loop_common/src`; `packages/loop_interpolation/src`.
- Machine: same as the first run (Windows 11, 20 logical CPUs). No other heavy work ran.
- Method: same as the first run. A scratchpad script (outside the repository, deleted at the end) gave a fake
  `benchmark` object (`extra_info`, `pedantic`) to `benchmarks.fault_evaluation.bench_cases.run_case`.
  Cases: `enumerate_cases()` filtered by model (`f4_u3`, `f4_u0`, `onlap_domain`) and `points == 5000` (30 cases).
  Time: median of 5 rounds, `time.perf_counter`. Fault calls: one extra run with `EvaluationInstrumentation`,
  `rounds=1`, `after_setup=reset`, value `snapshot()["fault_calls_total"]`.
  Cache off: `model.evaluation_cache.enabled = False`, then the cold `model_ids` case.
  Fingerprint: `fault_set_fingerprint(group.faults, build=False)`, 20 repeats for each group, median and minimum.
- Command: `python m.py` (the scratchpad script), run from the scratchpad folder.
- Every time is measured. None is inferred.

### Fault calls and median time, 5,000 points

| Case | Fault calls | Median (s) |
|---|---|---|
| f4_u3 oldest_value cold | 4 | 0.2114 |
| f4_u3 oldest_value warm | 0 | 0.0088 |
| f4_u3 all_group_values cold | 4 | 0.2276 |
| f4_u3 all_group_values warm | 0 | 0.0309 |
| f4_u3 model_ids cold | 4 | 0.2208 |
| f4_u3 model_ids warm | 0 | 0.0205 |
| f4_u3 oldest_gradient cold | 8 | 1.8037 |
| f4_u3 oldest_gradient warm | 0 | 0.3785 |
| f4_u3 model_gradient cold | 8 | 2.8261 |
| f4_u3 model_gradient warm | 0 | 1.5080 |
| f4_u0 oldest_value cold | 4 | 0.2051 |
| f4_u0 oldest_value warm | 0 | 0.0048 |
| f4_u0 all_group_values cold | 4 | 0.1921 |
| f4_u0 all_group_values warm | 0 | 0.0050 |
| f4_u0 model_ids cold | 4 | 0.2116 |
| f4_u0 model_ids warm | 0 | 0.0051 |
| f4_u0 oldest_gradient cold | 4 | 1.4983 |
| f4_u0 oldest_gradient warm | 0 | 0.3902 |
| f4_u0 model_gradient cold | 4 | 1.4992 |
| f4_u0 model_gradient warm | 0 | 0.3509 |
| onlap_domain oldest_value cold | 1 | 0.0283 |
| onlap_domain oldest_value warm | 0 | 0.0065 |
| onlap_domain all_group_values cold | 1 | 0.0364 |
| onlap_domain all_group_values warm | 0 | 0.0117 |
| onlap_domain model_ids cold | 1 | 0.0338 |
| onlap_domain model_ids warm | 0 | 0.0086 |
| onlap_domain oldest_gradient cold | 2 | 0.3303 |
| onlap_domain oldest_gradient warm | 0 | 0.2235 |
| onlap_domain model_gradient cold | 2 | 0.5580 |
| onlap_domain model_gradient warm | 0 | 0.4208 |

Samples (s) of the key cases:

- f4_u3 model_ids cold, cache on: 0.2171, 0.2295, 0.2262, 0.2160, 0.2208.
- f4_u3 model_ids cold, cache off: 0.2117, 0.2019, 0.2068, 0.2226, 0.2142 (median 0.2117, fault calls 4).
- f4_u3 model_ids warm: 0.0187, 0.0319, 0.0184, 0.0205, 0.0221.
- onlap_domain model_gradient cold: 0.6537, 0.5901, 0.5553, 0.5580, 0.5449.

### Fingerprint time, `fault_set_fingerprint(faults, build=False)`, 20 repeats

| Model and group | Faults | Median (ms) | Minimum (ms) |
|---|---|---|---|
| f4_u3 group_0 | 4 | 3.502 | 2.675 |
| f4_u3 group_1 | 4 | 2.759 | 2.585 |
| f4_u3 group_2 | 4 | 2.872 | 2.615 |
| f4_u3 group_3 | 4 | 2.643 | 2.538 |
| f4_u0 group_0 | 4 | 2.551 | 2.395 |
| onlap_domain group_1 | 1 | 2.275 | 1.956 |
| onlap_domain group_2 | 1 | 2.201 | 1.933 |

Groups with 0 faults take 0.001 ms and are left out of the table (`onlap_domain` group_0). The group_0 median of `f4_u3` is
the first call in the process, and it is higher than the other groups. Its minimum is in line with the others.
The worst median is 3.50 ms, so the target of 5 ms holds in every group.

### Gates

| Gate | Measured | Target | Result |
|---|---|---|---|
| G1 SC-003 count: f4_u3 cold all_group_values and model_ids | 4 and 4 | 4 each | pass |
| G2 cold model_ids, cache on / cache off (f4_u3) | 0.2208 / 0.2117 = 1.043 | at most 1.10 | pass |
| G3 warm model_ids fault calls (f4_u3) | 0 | 0 | pass |
| G4 `fault_set_fingerprint`, 4 faults of f4_u3, `build=False`, 20 repeats | median 2.64 to 3.50 ms (worst 3.50); minimum 2.54 to 2.68 ms | at most 5 ms | pass |
| G5 onlap_domain model_gradient cold | 0.5580 s | below the baseline 0.92 s | pass |

### OP5: warm evaluate_model / cold evaluate_model

| Model | Cold (s) | Warm (s) | Ratio warm / cold |
|---|---|---|---|
| f4_u3 | 0.2208 | 0.0205 | 0.093 |
| f4_u0 | 0.2116 | 0.0051 | 0.024 |
| onlap_domain | 0.0338 | 0.0086 | 0.254 |

The `onlap_domain` ratio stays above 0.10. It has no gate here. It is reported for T030 (SC-010).

### Comparison with the first run

| Item | First run (598ef21e) | Second run (45980adc) |
|---|---|---|
| G4 fingerprint median, f4_u3 | 7.57 to 8.07 ms | 2.64 to 3.50 ms |
| G4 fingerprint minimum, f4_u3 | 6.6 to 6.9 ms | 2.54 to 2.68 ms |
| OP5 f4_u3 | 0.109 | 0.093 |
| OP5 f4_u0 | 0.040 | 0.024 |
| OP5 onlap_domain | 0.284 | 0.254 |
| f4_u3 model_ids warm | 25.9 ms | 20.5 ms |

D4b (T045) cut the fingerprint time by about a factor of 2.3 to 2.9. The warm times fell by 5 ms or more for `f4_u3`.
Fault counts did not change (G1 and G3 still pass).

## First run (598ef21e)

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
