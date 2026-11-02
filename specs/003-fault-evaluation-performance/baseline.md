# Baseline of the standard benchmark suite (T011)

This file records the benchmark baseline of the code before the performance
changes of spec 003. It was measured on 2026-11-02. The numbers come from
`baseline.json` and `baseline_repeat.json` in this folder.

## What was measured

- **Suite**: the standard suite (254 cases): all 210 cases at 5,000 points, the
  20 cases at 20,000 points, and the 24 fixed cases at 50,000 and 100,000 points
  (plan D2). Both result files have the same 254 case ids.
- **Code**: revision `617cfea703f25a8ca3c1fc092752c7bc9ec48b27`, a clean tree
  (dirty flag false). The corrected gradient code of T034 is included. No
  performance change is included.
- **Computer**: Intel64 Family 6 Model 154, 34.0 GB RAM, Windows 11
  (10.0.26200). Python 3.13.12, NumPy 2.4.6, SciPy 1.18.0, LoopStructural 1.11.0.
- **Commands** (output and work directories were outside the repository, and the
  result files were copied here):

  ```
  python -m benchmarks.fault_evaluation.run --suite standard --output <dir>/baseline.json --work-dir <dir>/run1.work
  python -m benchmarks.fault_evaluation.run --suite standard --output <dir>/baseline_repeat.json --work-dir <dir>/run2.work
  python -m benchmarks.fault_evaluation.compare baseline.json baseline_repeat.json --output <path>
  ```

- **Rounds**: every case runs one round first. A timing case with a median below
  1.0 s runs again for 5 rounds (130 of 254 cases). 124 cases have a median of 1.0
  s or more.

## Runtime of the suite

| Run | Timing pass | Counters pass | Memory pass | Whole run |
|---|---|---|---|---|
| 1 (`baseline.json`) | 1,219 s | 969 s | 1,083 s | 3,271 s (54.5 min) |
| 2 (`baseline_repeat.json`) | 1,280 s | 1,057 s | 1,088 s | 3,426 s (57.1 min) |

Both runs are inside the budget of about 1 hour (C-006). SC-008 (under 15 minutes)
is a target for the code after this feature. It is not tested here.

## Values that the success criteria use (run 1, cold, median seconds)

| Case | Time | Fault calls | Field calls | Fault points | Peak extra memory |
|---|---|---|---|---|---|
| `f4_u3:oldest_value:5000` | 1.963 s | 16 | 644 | 80,000 | 2.0 MB |
| `f4_u3:all_group_values:5000` | 6.478 s | 52 | 2,093 | 260,000 | 5.5 MB |
| `f4_u3:model_ids:5000` | 6.337 s | 52 | 2,093 | 260,000 | 6.3 MB |
| `f4_u0:all_group_values:5000` | 0.505 s | 4 | 161 | 20,000 | 7.9 MB |
| `plugin:oldest_value:20000` | 4.397 s | 12 | 502 | 240,000 | 12.9 MB |
| `f1_u0:oldest_value:20000` | 0.623 s | 1 | 41 | 20,000 | 14.1 MB |

- The fault-call counts agree with the spec: 16 for the oldest group and 52 for all
  groups of the model with 4 faults and 3 unconformities.
- The times on this computer are lower than the numbers of `research.md`
  (1.96 s against 2.55 s, 6.48 s against 8.65 s). The success criteria use ratios
  against this baseline, so this does not change a target. For SC-004 the
  baseline ratio of `f4_u3` to `f4_u0` (all groups, 5,000 points) is 12.8
  (spec: 12.9).
- A warm case on the code before this feature takes the same time as a cold case
  (`f4_u3:model_ids:5000`: 6.337 s cold, 6.419 s warm), as expected, because there
  is no cache. The runner does not make a priming call for a model with no cache.

## Noise between the two runs (SC-007)

The limit is a difference of less than 10 % (ratio below 1.10) for each case with a
baseline time of 1.0 s or more (spec SC-007, plan D2b).

| Group | Cases | Median noise ratio | Largest | Cases at 1.10 or more |
|---|---|---|---|---|
| Baseline time 1.0 s or more | 124 | 1.053 | 1.442 | **35** |
| Baseline time below 1.0 s (reported only) | 130 | 1.097 | 1.457 | 64 |

**SC-007 FAILS on this baseline: 35 of 124 checked cases differ by 10 % or more.**
The largest are `f1_u2:oldest_gradient:50000:warm` (11.8 s and 17.0 s, ratio 1.44),
`cutting_faults:model_gradient:5000:warm` (44.6 s and 61.0 s, 1.37) and
`f1_u2:model_gradient:5000:warm` (2.5 s and 3.4 s, 1.37).

Observations, not conclusions:

- Run 2 is slower in most cases. The sum of the medians is 943 s in run 1 and 1,044
  s in run 2 (11 % more), and the whole run took 4.7 % longer. The direction is the
  same in most groups, which points to a drift of the computer (heat, background
  load) more than to random noise in single cases. This record cannot prove the
  cause.
- A case of 1.0 s or more has one timed round in the first step, so it has no
  averaging. Only a case that the first step measured below 1.0 s has 5 rounds.
- The compare tool marks 98 cases as slower than the noise limit. The count ratios
  (fault and field calls) are exactly 1.0 for every case.

The scaling tables of the compare tool (median ratio of run 2 to run 1, by group):

| Unconformities (F = 0) | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| Ratio | 1.08 | 1.11 | 1.15 | 1.09 |

| Faults (U = 0) | 0 | 1 | 2 | 4 |
|---|---|---|---|---|
| Ratio | 1.08 | 1.10 | 1.00 | 1.00 |

## Decision needed (OP2)

SC-007 fails, so plan open point OP2 is active. The plan says to report the measured
result and to ask for an explicit spec decision. Nothing was removed from the suite
and no number was changed to make the check pass. Options for the decision:

1. Keep SC-007 as written and repeat the baseline on a quiet computer, with the
   same code, to see if the drift goes away.
2. Measure the cases of 1.0 s or more with more rounds (for example 3), which
   lengthens the suite.
3. Change SC-007 to a looser limit, or to a comparison of the sum of the medians.

Work on the next tasks does not depend on this decision. The final acceptance
check (T030) does.

## Files

- `baseline.json`: run 1. `baseline_repeat.json`: run 2.
- The work directories (checkpoint files and the raw benchmark data) stay outside
  the repository, in `C:/Users/lachl/dev/ls-bench-003/std/`.
