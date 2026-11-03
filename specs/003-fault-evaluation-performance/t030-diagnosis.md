# T030 diagnosis, split by cause (T043)

This report supports the decisions OP3 (SC-008 target) and OP4 (SC-006 for
small baselines). It proposes no resolution. It only splits the failed T030
results by cause.

Labels: **M** means measured (the number is in a raw file). **D** means
derived (computed from M numbers, or a rule that this report applies).
A cause that has no profile is "probable".

Sources (all read only):

- `C:/Users/lachl/dev/ls-bench-003/cand/candidate.json` (run 1) and
  `candidate_repeat.json` (run 2). Both have `suite_wall_seconds`.
- `acceptance_cand.json` (the SC-006 and SC-008 rows).
- `run1.work` and `run2.work`: `meta.json` (`wall_seconds` of each pass) and
  the `*.checkpoint.jsonl` files (`elapsed_seconds` and the samples of each case).
- `run1.log` and `run2.log` (the pytest line "N passed in X s" of each launch).
- `specs/003-fault-evaluation-performance/baseline.json`.

Note on provenance (M): the baseline ran on Python 3.13.12. The candidate ran
on Python 3.14.4. `acceptance_cand.json` reports this difference as a
provenance problem. Both candidate runs used code revision `c8e3f9f1`.

## Part (a): SC-006 (peak extra RSS, limit 1.20)

### Result

`SC-006` fails (M). The reason text of `acceptance_cand.json` names 45 cases:

- 43 cases have a ratio above 1.20. The highest ratio is 45.6.
- 2 cases have a baseline of 0 bytes and a candidate value above 0 bytes.
  SC-006 has no ratio for them, but it counts them as failures.

Failing cases with a baseline below 10 MB (10,000,000 bytes) (D, counted from
the table): **34 of 45**. This is 32 of the 43 ratio cases and both
zero-baseline cases. The 11 other failing cases have a baseline of 10 MB or
more. The highest ratio among them is 1.5564
(`plugin:oldest_gradient:5000:cold`). The count is the same with a limit of
10 MiB.

### Rule for the probable cause (D)

The suite does not profile memory by source. The cause column uses this rule.
Treat it as a hypothesis, not as a finding.

1. If the extra bytes (candidate minus baseline) are below 1 MiB
   (1,048,576), the cause is **other**. RSS values are multiples of 4,096
   bytes. A change of some pages can come from the allocator or from noise.
2. Otherwise, for a **cold** case, the cause is **cache retention**. The cache
   is empty at the start. The timed call fills it, and the sampler sees the
   stored entry.
3. Otherwise, for a **warm** case, the cause is **hash workspace**. The
   priming call fills the cache before the RSS sampler starts (`run_case`
   calls `after_setup` last). The stored entry is then not in the extra bytes.
   The temporary arrays of the key and token work in the timed call stay as
   the probable cause.

The rule uses only the state, the extra bytes and the model name. It does not
check the operation or the code path. For the `f0_*` models (no fault calls
in the baseline or the candidate), the label comes from the rule alone.

Count by cause (D): cache retention 19, other 16, hash workspace 10
(total 45). Among the 34 cases with a baseline below 10 MB: cache retention
13, other 16, hash workspace 5.

The column "Fault calls" shows the baseline and the candidate value of
`fault_calls_total` (M). A candidate value of 0 for a warm case means that the
fault work was a cache hit.

### Table (one row for each failing case of `acceptance_cand.json`)

| Case | Baseline peak extra RSS (bytes, M) | Candidate (bytes, M) | Ratio (D) | Extra bytes (D) | Fault calls base / cand (M) | Probable cause (D) |
|---|---:|---:|---:|---:|---:|---|
| f0_u0:oldest_value:5000:cold | 1,085,440 | 1,359,872 | 1.2528 | 274,432 | 0 / 0 | other |
| f0_u0:all_group_values:5000:warm | 4,096 | 147,456 | 36.0000 | 143,360 | 0 / 0 | other |
| f0_u0:model_ids:5000:cold | 2,383,872 | 2,977,792 | 1.2491 | 593,920 | 0 / 0 | other |
| f0_u1:oldest_value:100000:warm | 37,838,848 | 54,001,664 | 1.4271 | 16,162,816 | 0 / 0 | hash workspace |
| f0_u1:model_ids:5000:cold | 4,096 | 28,672 | 7.0000 | 24,576 | 0 / 0 | other |
| f0_u1:model_ids:50000:warm | 16,773,120 | 22,990,848 | 1.3707 | 6,217,728 | 0 / 0 | hash workspace |
| f0_u1:oldest_gradient:5000:cold | 1,134,592 | 1,368,064 | 1.2058 | 233,472 | 0 / 0 | other |
| f0_u2:oldest_value:5000:cold | 8,192 | 155,648 | 19.0000 | 147,456 | 0 / 0 | other |
| f0_u2:all_group_values:5000:cold | 2,613,248 | 3,338,240 | 1.2774 | 724,992 | 0 / 0 | other |
| f0_u2:model_ids:5000:warm | 4,096 | 28,672 | 7.0000 | 24,576 | 0 / 0 | other |
| f0_u2:model_gradient:5000:cold | 122,880 | 225,280 | 1.8333 | 102,400 | 0 / 0 | other |
| f0_u3:model_ids:5000:cold | 24,576 | 319,488 | 13.0000 | 294,912 | 0 / 0 | other |
| f0_u3:model_ids:100000:warm | 41,807,872 | 54,038,528 | 1.2925 | 12,230,656 | 0 / 0 | hash workspace |
| f0_u3:oldest_gradient:5000:cold | 163,840 | 2,891,776 | 17.6500 | 2,727,936 | 0 / 0 | cache retention |
| f0_u3:model_gradient:5000:cold | 20,480 | 933,888 | 45.6000 | 913,408 | 0 / 0 | other |
| f4_u0:model_gradient:5000:cold | 17,010,688 | 21,667,840 | 1.2738 | 4,657,152 | 4 / 4 | cache retention |
| f4_u3:oldest_value:5000:cold | 1,986,560 | 4,427,776 | 2.2289 | 2,441,216 | 16 / 4 | cache retention |
| stopped_fault:oldest_value:5000:cold | 61,440 | 2,555,904 | 41.6000 | 2,494,464 | 1 / 1 | cache retention |
| stopped_fault:model_ids:5000:cold | 28,672 | 49,152 | 1.7143 | 20,480 | 1 / 1 | other |
| stopped_fault:oldest_gradient:5000:cold | 5,922,816 | 10,801,152 | 1.8237 | 4,878,336 | 1 / 1 | cache retention |
| stopped_fault:oldest_gradient:5000:warm | 6,901,760 | 12,095,488 | 1.7525 | 5,193,728 | 1 / 0 | hash workspace |
| stopped_fault:model_gradient:5000:cold | 5,246,976 | 7,798,784 | 1.4863 | 2,551,808 | 1 / 1 | cache retention |
| cutting_faults:oldest_value:5000:cold | 16,277,504 | 20,570,112 | 1.2637 | 4,292,608 | 42 / 22 | cache retention |
| cutting_faults:all_group_values:5000:cold | 18,223,104 | 25,538,560 | 1.4014 | 7,315,456 | 42 / 22 | cache retention |
| onlap_domain:oldest_value:5000:cold | 348,160 | 1,458,176 | 4.1882 | 1,110,016 | 2 / 1 | cache retention |
| onlap_domain:all_group_values:5000:warm | 24,576 | 32,768 | 1.3333 | 8,192 | 4 / 0 | other |
| onlap_domain:model_ids:5000:cold | 20,480 | 32,768 | 1.6000 | 12,288 | 4 / 1 | other |
| onlap_domain:model_ids:50000:warm | 22,892,544 | 28,590,080 | 1.2489 | 5,697,536 | 4 / 0 | hash workspace |
| onlap_domain:oldest_gradient:5000:cold | 4,497,408 | 14,020,608 | 3.1175 | 9,523,200 | 2 / 2 | cache retention |
| onlap_domain:oldest_gradient:5000:warm | 3,186,688 | 14,815,232 | 4.6491 | 11,628,544 | 2 / 0 | hash workspace |
| onlap_domain:oldest_gradient:100000:cold | 213,581,824 | 286,228,480 | 1.3401 | 72,646,656 | 2 / 2 | cache retention |
| onlap_domain:model_gradient:5000:cold | 5,267,456 | 12,357,632 | 2.3460 | 7,090,176 | 4 / 2 | cache retention |
| onlap_domain:model_gradient:5000:warm | 6,467,584 | 14,635,008 | 2.2628 | 8,167,424 | 4 / 0 | hash workspace |
| faulted_fold:oldest_value:5000:warm | 7,946,240 | 10,334,208 | 1.3005 | 2,387,968 | 1 / 1 | hash workspace |
| faulted_fold:all_group_values:5000:warm | 3,874,816 | 15,978,496 | 4.1237 | 12,103,680 | 2 / 1 | hash workspace |
| faulted_fold:model_ids:5000:cold | 3,072,000 | 5,783,552 | 1.8827 | 2,711,552 | 1 / 1 | cache retention |
| faulted_fold:oldest_gradient:5000:cold | 42,700,800 | 52,666,368 | 1.2334 | 9,965,568 | 1 / 1 | cache retention |
| faulted_fold:oldest_gradient:5000:warm | 40,345,600 | 52,117,504 | 1.2918 | 11,771,904 | 1 / 1 | hash workspace |
| plugin:oldest_value:5000:cold | 2,351,104 | 6,877,184 | 2.9251 | 4,526,080 | 12 / 6 | cache retention |
| plugin:all_group_values:5000:cold | 4,075,520 | 10,326,016 | 2.5337 | 6,250,496 | 12 / 6 | cache retention |
| plugin:model_ids:5000:cold | 1,982,464 | 5,328,896 | 2.6880 | 3,346,432 | 12 / 6 | cache retention |
| plugin:oldest_gradient:5000:cold | 15,605,760 | 24,289,280 | 1.5564 | 8,683,520 | 12 / 12 | cache retention |
| plugin:model_gradient:5000:cold | 9,281,536 | 19,865,600 | 2.1403 | 10,584,064 | 12 / 12 | cache retention |
| f0_u1:all_group_values:5000:warm | 0 | 32,768 | n/a (baseline 0) | 32,768 | 0 / 0 | other |
| f0_u2:oldest_value:5000:warm | 0 | 28,672 | n/a (baseline 0) | 28,672 | 0 / 0 | other |
The candidate values come from `candidate.json` (run 1). The memory pass
runs each case once (1 round). Run 2 (`candidate_repeat.json`) is not used for
SC-006.

## Part (b): SC-008 (suite wall time, limit 900 s)

### Result

`SC-008` fails (M). `suite_wall_seconds` is 3,016.748 s for run 1 and
2,779.902 s for run 2. The baseline value is 3,270.599 s. The sum of the case
`median_seconds` (D) is 334.27 s for run 1 (the value that plan.md names),
365.70 s for run 2, and 943.15 s for the baseline.

### How the runner takes the time (M)

`run.py` runs three passes in this order: timing, counters, memory (`rss`).
Each launch runs in a new pytest process. The timing pass has two launches:
one for all 254 cases (1 round), and one "repeat" launch for the short cases
(5 rounds). The runner adds the wall time of each pass into `meta.json`, and
`suite_wall_seconds` is the sum of the three values.

For each case, `elapsed_seconds` in the checkpoint file covers the whole test
body: point sampling, model build (in the first case of a model in a launch),
cache clear, warm priming, and the timed calls.

### Table 1: pass times (M, from `meta.json`)

| Pass | Run 1 (s) | Run 2 (s) |
|---|---:|---:|
| timing (main launch, repeat launch, process starts) | 1,568.134 | 1,543.143 |
| counters | 813.964 | 622.054 |
| rss (memory) | 634.650 | 614.706 |
| **Sum of the passes (D)** | **3,016.748** | **2,779.902** |
| `suite_wall_seconds` (M) | 3,016.748 | 2,779.902 |
| Difference (D) | 0.000 (0.00 %) | 0.000 (0.00 %) |

### Table 2: split of each pass by launch and by activity

Columns:

- "Session" is the pytest time from the log (M).
- "Timed calls" is the sum of all time samples in the checkpoint file (M).
- "Warm non-timed" is the `elapsed_seconds` of the warm cases minus their
  timed samples (D). It is mainly the warm priming call, plus the cache clear.
- "Cold non-timed" is the same value for the cold cases (D). It holds point
  sampling and model build. A cold case has no priming call.
- "Session overhead" is the session time minus the sum of `elapsed_seconds`
  (D). It holds collection, the benchmark statistics and teardown.
- "Outside pytest" is the pass wall time minus the session times (D). It holds
  the process start and the imports.

Run 1 (all values in seconds):

| Pass and launch | Session (M) | Timed calls (M) | Warm non-timed (D) | Cold non-timed (D) | Session overhead (D) |
|---|---:|---:|---:|---:|---:|
| timing, main launch (254 cases, 1 round) | 568.58 | 322.6 | 238.2 | 1.5 | 6.2 |
| timing, repeat launch (217 short cases, 5 rounds) | 996.49 | 278.9 | 708.3 | 1.8 | 7.5 |
| timing, outside pytest (2 process starts) | 3.06 (D) | | | | |
| counters | 811.98 | 433.0 | 368.9 | 2.2 | 7.9 |
| counters, outside pytest | 1.98 (D) | | | | |
| rss | 631.97 | 348.8 | 271.8 | 1.9 | 9.4 |
| rss, outside pytest | 2.68 (D) | | | | |
| **Total** | **3,016.75** | | | | |

Run 2 (all values in seconds):

| Pass and launch | Session (M) | Timed calls (M) | Warm non-timed (D) | Cold non-timed (D) | Session overhead (D) |
|---|---:|---:|---:|---:|---:|
| timing, main launch (254 cases, 1 round) | 645.65 | 367.7 | 271.0 | 1.8 | 5.2 |
| timing, repeat launch (212 short cases, 5 rounds) | 894.12 | 228.3 | 658.0 | 1.5 | 6.2 |
| timing, outside pytest (2 process starts) | 3.37 (D) | | | | |
| counters | 620.41 | 351.2 | 261.9 | 1.6 | 5.6 |
| counters, outside pytest | 1.64 (D) | | | | |
| rss | 612.53 | 348.7 | 256.2 | 1.9 | 5.7 |
| rss, outside pytest | 2.18 (D) | | | | |
| **Total** | **2,779.90** | | | | |

Check against `suite_wall_seconds` (D): the session times plus the outside
times add up to 3,016.75 s (run 1) and 2,779.90 s (run 2). The difference
from the measured values is below 0.01 s, so it is inside the 1 % limit. In
each row, the four split columns add up to the session value within 0.1 s
(rounding).

### What the data shows and does not show

- Warm priming (D): the "warm non-timed" time is 238.2 + 708.3 + 368.9 +
  271.8 = 1,587.2 s in run 1 (52.6 % of the suite). It is 271.0 + 658.0 +
  261.9 + 256.2 = 1,447.1 s in run 2 (52.1 %). This is the largest part of the
  suite. The repeat launch of the timing pass costs 996.49 s in run 1 and
  894.12 s in run 2. Of this, 708.3 s and 658.0 s is warm non-timed time. The
  data does not split priming from the cache clear. The two are in one value.
- Process starts (D): 3.06 s, 1.98 s and 2.68 s outside pytest in run 1 (7.7 s
  in total). The timing pass has two launches. The time of one launch is not
  in the data.
- Model build: the data has no separate model build time. It is in the "cold
  non-timed" time (1.5 to 2.2 s for each launch, an upper bound) and in the
  first case of each model.
- Timed calls (M): in the four launches of run 1, the timed calls add up to
  322.6 + 278.9 + 433.0 + 348.8 = 1,383.3 s (45.9 % of the suite). The
  counters pass has instrumented timed calls, so those values are not
  comparable to the timing pass.
- The 334 s value in plan.md is the sum of `median_seconds`. It is not a pass
  time. The timed calls of the main timing launch (322.6 s) are close to it.
- plan.md says "about 2,000 s is not evaluation". The difference between
  `suite_wall_seconds` and the sum of the medians is 2,682.5 s (D, run 1). The
  tables above show the split of this time. For the warm cases, the data
  cannot split the non-timed time between the priming call and other work.

### Table 3: the two `onlap_domain` 100,000-point gradient cases

The two cases are `onlap_domain:oldest_gradient:100000:cold` and
`onlap_domain:model_gradient:100000:warm` (M). They are the only cases of this
model with 100,000 points and a gradient operation. The first has no priming.
The second has one priming call, so its non-timed time is mostly that call
(D). The two cases are not in the repeat launch.

Each cell shows "elapsed / timed" in seconds. "Elapsed" is `elapsed_seconds`
(M). "Timed" is the sum of the samples (M). The last two columns are sums (D).

| Run and pass | `oldest_gradient:100000:cold` | `model_gradient:100000:warm` | Both, elapsed | Both, timed |
|---|---:|---:|---:|---:|
| Run 1, timing | 43.20 / 43.19 | 167.52 / 85.32 | 210.72 | 128.51 |
| Run 1, counters | 52.61 / 52.60 | 222.61 / 101.02 | 275.22 | 153.62 |
| Run 1, rss | 50.00 / 50.00 | 192.14 / 87.13 | 242.14 | 137.13 |
| Run 1, total | | | 728.08 (24.1 % of 3,016.75) | 419.26 |
| Run 2, timing | 50.49 / 50.48 | 178.55 / 91.51 | 229.04 | 141.99 |
| Run 2, counters | 46.03 / 46.02 | 179.92 / 90.38 | 225.95 | 136.40 |
| Run 2, rss | 47.12 / 47.11 | 180.13 / 89.07 | 227.25 | 136.18 |
| Run 2, total | | | 682.24 (24.5 % of 2,779.90) | 414.57 |

Baseline medians for the same two cases (M, `baseline.json`): 45.51 s and
92.11 s (137.62 s in total). The candidate medians in run 1 are 43.19 s and
85.32 s (128.51 s in total). The change is small. The "about 128 s per pass"
figure in plan.md is the timed value of the run 1 timing pass.

The two cases take 24 % of the suite wall time in both runs (D). In the run 1
timing pass, the priming call of the warm case takes 82.2 s (D, elapsed minus
timed).

## Open points for OP3 and OP4 (not resolved here)

- OP3: the timed calls of run 1 add up to 1,383.3 s in four launches. The
  non-timed time is 1,633.4 s (D). A target of 900 s is below both the timed
  total and the non-timed total, so the data shows that faster fault work
  alone cannot meet it. This report does not choose a new target or a method.
- OP4: 34 of the 45 failing SC-006 cases have a baseline below 10 MB. In 16 of
  the 45 cases the extra bytes are below 1 MiB. Two cases have a baseline of
  0 bytes, so SC-006 has no ratio for them. Whether a small baseline needs a
  different limit is a decision for the user.
