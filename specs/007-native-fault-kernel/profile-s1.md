# Profile of stage S1 (spec 007)

This file gives the cut-fault numbers of the 94-fault model before stage S1 (SC-001 baseline), and a profile of one `evaluate_model` call at 500 points.

- Commit: `d979cff40551c7b55851f20effff81cc3b70662b` (worktree `feat-007-native-fault-kernel`). No code was changed for the measurement.
- Computer: 12th Gen Intel(R) Core(TM) i7-12700H, 20 threads (default thread settings), Ubuntu 26.04.1 LTS, kernel 7.0.0-44-generic, Python 3.13.13.
- Limit: other agents can run work on the same computer at the same time. The times are not exact. Use them to compare shares, not to set tight limits.
- Scripts: `many_faults.build()` and `many_faults.points(500)`, run from the scratch folder `tmp/t002` (`prof.py`).

## Baseline (Python path, master)

Source: `_runs007/many_faults_python.json`, case `"500"`. The bit-identity reference is `_runs007/many_faults_python_ref.npz`.

| Value (500 points) | Result |
|---|---|
| `evaluate_model_seconds` | 1.78 s |
| `fault_calls` | 94 |
| `cutting_mean_ms` | 40.0 ms |
| `cutting_max_ms` | 117.4 ms |
| `isolated_mean_ms` | 2.65 ms |

## Profile before S1

Calls: the first `evaluate_model` call was a warm-up. The second call, without the profiler, took 1.73 s. The third call, with `cProfile`, took 8.17 s. The profiled call is the third call. The profiler makes the run about 4.7 times slower, because the code makes 5.3 million small `is_in` calls. The shares below are therefore only a guide. Use the times without the profiler for SC-001.

The cumulative time of a function includes the time of the functions it calls. Some functions call themselves (recursion), so the shares can add up to more than 100 percent and must not be added together.

| Function | ncalls (primitive) | cumulative time (s) | share of 8.17 s |
|---|---|---|---|
| `_calculate_mask` | 2868 (2) | 8.098 | 99.1 % |
| `restore_through` | 2557 (1) | 8.022 | 98.1 % |
| `apply_to_points` | 94 | 8.001 | 97.9 % |
| `_evaluate_frame` | 1034 | 7.908 | 96.7 % |
| `_evaluate_frame_legacy` | 698 | 7.808 | 95.5 % |
| `cutters_first` | 2557 | 6.114 | 74.8 % |
| `_cutters` | 296025 | 5.587 | 68.4 % |
| `is_in` | 5308363 | 4.903 | 60.0 % |
| `up_to_date` (`_base_builder.py`) | 716108 (4696) | 0.667 | 8.2 % |
| `up_to_date` (`_fault_builder.py`) | 235948 (4912) | 0.647 | 7.9 % |
| `with_cutters` | 2557 | 0.596 | 7.3 % |
| `_shared_frame_members` | 754 | 0.026 | 0.3 % |

Top 20 by cumulative time:

| Function | ncalls | tottime (s) | cumtime (s) |
|---|---|---|---|
| `_api_registry.py:25(wrapper)` | 1 | 0.000 | 8.174 |
| `geological_model.py:1235(evaluate_model)` | 1 | 0.000 | 8.174 |
| `_geological_feature.py:38(wrapper)` | 2560 (2) | 0.006 | 8.174 |
| `_geological_feature.py:153(evaluate_value)` | 1900 (2) | 0.012 | 8.174 |
| `_geological_feature.py:143(_restore_points)` | 2560 (4) | 0.001 | 8.168 |
| `_evaluation.py:222(restore)` | 4 | 0.000 | 8.168 |
| `_evaluation.py:102(restore)` | 4 | 0.000 | 8.168 |
| `_base_geological_feature.py:236(_calculate_mask)` | 2868 (2) | 0.004 | 8.098 |
| `_unconformity_feature.py:89(__call__)` | 2 | 0.000 | 8.098 |
| `_unconformity_feature.py:65(evaluate)` | 2 | 0.000 | 8.098 |
| `_fault_block.py:210(restore_through)` | 2557 (1) | 0.007 | 8.022 |
| `_fault_segment.py:648(apply_to_points)` | 94 | 0.028 | 8.001 |
| `_fault_segment.py:355(_evaluate_frame)` | 1034 | 0.007 | 7.908 |
| `_fault_segment.py:446(_evaluate_frame_legacy)` | 698 | 0.004 | 7.808 |
| `_base_geological_feature.py:272(_apply_faults)` | 2556 | 0.002 | 6.932 |
| `_fault_block.py:189(cutters_first)` | 2557 | 0.173 | 6.114 |
| `_fault_block.py:151(_cutters)` | 296025 | 1.264 | 5.587 |
| `_fault_block.py:146(is_in)` | 5308363 | 1.996 | 4.903 |
| `~:0(<built-in method builtins.any>)` | 5261591 (5578846) | 1.670 | 2.802 |
| `_fault_block.py:148(<genexpr>)` | 20193440 | 1.299 | 1.299 |

The call to `is_in` and its generator, `_cutters`, `cutters_first` and `any` are in `LoopStructural/utils/_fault_block.py`. `_calculate_mask` is in `_base_geological_feature.py`. The next row after the top 20 is `evaluate_gradient` (1.978 s).

## Where the time goes

The restore order work is the largest part: `cutters_first` takes 74.8 % of the profiled time, and nearly all of it is in `_cutters` and `is_in` (5.3 million calls). Together with their `any` calls and generators, the self time of these functions is about 6.4 s of 8.17 s. `with_cutters` takes only 7.3 % and `_shared_frame_members` less than 1 %. The frame evaluation work (`_evaluate_frame_legacy`) is high in the cumulative column (95.5 %), but this is because it calls `restore_through` again for each cutting fault; its own time is very small. The `up_to_date` checks add about 16 % in two functions that overlap with the rest. Plan OP3: S1 must remove the `is_in` and `_cutters` loops first, because they are the main cost; the cost of the frame evaluation itself stays small in this profile.

## After D1a

Run: the benchmark and the profile ran one after the other, with the code of D1a (T004). No code was changed for the measurement. Limit: other agents can build Rust code on the same computer at the same time, so the times are not exact.

Command: `python -m benchmarks.fault_evaluation.many_faults --output many_faults_d1a.json --points 500 2000 5000 --compare _runs007/many_faults_python_ref.npz`. The reference file holds arrays for 500, 2000 and 5000 points, and `compare()` marks an array that is missing from the run as a difference. Therefore the run used all three point counts. SC-001 uses the 500 case.

| Value (500 points) | Baseline | After D1a |
|---|---|---|
| `evaluate_model_seconds` | 1.78 s | 0.730 s |
| `fault_calls` | 94 | 94 |
| `cutting_mean_ms` | 40.0 ms | 13.6 ms |
| `cutting_max_ms` | 117.4 ms | 25.4 ms |
| `isolated_mean_ms` | 2.65 ms | 2.66 ms |

Other point counts: 2000 points: cutting mean 15.8 ms, max 27.7 ms, isolated mean 2.75 ms, `evaluate_model` 0.808 s. 5000 points: cutting mean 20.1 ms, max 32.9 ms, isolated mean 3.42 ms, `evaluate_model` 1.019 s. The model build took 195 s.

`--compare` exit code: 1. Message: `differ: ['group_0_2000', 'group_0_500', 'group_0_5000', 'group_1_2000', 'group_1_500', 'group_1_5000']`. The arrays `model_ids_*` are bit-identical. The values of the two groups differ from the reference at all three point counts. This needs a check by the owner of D1a (the difference can come from a change in the order of the fault restore).

SC-001: fail (13.6 ms)

### Profile after D1a

Calls: one warm-up `evaluate_model` call, then one call without the profiler (0.678 s), then one call with `cProfile` (2.456 s; before D1a: 1.73 s and 8.17 s). The shares are of 2.456 s. As before, cumulative times include called functions, so the shares must not be added together.

| Function | ncalls | cumulative time (s) | share of 2.456 s | share before D1a |
|---|---|---|---|---|
| `_calculate_mask` | 2868 | 2.381 | 96.9 % | 99.1 % |
| `restore_through` (now `_fault_block.py:255`) | 2557 | 2.306 | 93.9 % | 98.1 % |
| `apply_to_points` | 94 | 2.302 | 93.7 % | 97.9 % |
| `_evaluate_frame` | 1034 | 2.213 | 90.1 % | 96.7 % |
| `_evaluate_frame_legacy` | 698 | 2.115 | 86.1 % | 95.5 % |
| `_apply_faults` | 2556 | 1.287 | 52.4 % | 84.8 % |
| `restore_order` (replaces `cutters_first`) | 2557 | 1.051 | 42.8 % | 74.8 % (`cutters_first`) |
| `is_in` | 790406 | 0.936 | 38.1 % | 60.0 % (5308363 calls) |
| `_cutters` | 52428 | 0.834 | 34.0 % | 68.4 % (296025 calls) |
| `up_to_date` (`_base_builder.py`) | 716108 | 0.649 | 26.4 % | 8.2 % |
| `up_to_date` (`_fault_builder.py`) | 235948 | 0.629 | 25.6 % | 7.9 % |
| `with_cutters` | 2557 | 0.589 | 24.0 % | 7.3 % |
| `any` (built-in) | 791078 | 0.580 | 23.6 % | n/a |

The `is_in` and `_cutters` work is still the largest single part (about 38 % of the profiled time), but it fell from 5.3 million to 0.79 million `is_in` calls. The `up_to_date` checks and `with_cutters` now have a larger share, because the total is smaller (their absolute times are about the same).
