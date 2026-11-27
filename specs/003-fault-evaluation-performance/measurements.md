# Measurements of T030 (final run, 2026-11-27)

T030 measures the code after spec 003 against the code before it, on one
computer, with the standard suite run twice for each side.

## What was measured

- **Computer**: Intel i7-12700H, 32 GB RAM, Linux 7.0, Python 3.13.13, NumPy
  2.4.6, SciPy 1.18.0. No other work ran during the runs.
- **Baseline (T062, spec C-019)**: branch `bench/003-linux-baseline-db5e65ac`
  at `23a37ee3`: master `db5e65ac` (the fixes #424, #425, #427, #428 and
  #429 merged, no spec 003 code) with the harness of this branch at
  `e88269b1`; no library change. Files `baseline_linux.json` and
  `baseline_linux_repeat.json` (906 s and 988 s).
- **Candidate**: `e88269b1` on `feat/003-fault-evaluation-performance`, clean
  tree. Files `candidate.json` and `candidate_repeat.json` (329 s and 327 s).
- **Harness**: standard suite; tracemalloc memory meter (C-020); a timing
  case below 1 s uses 5 rounds, a case from 1 s up to 5 s the median of 3
  rounds (C-022).
- **Commands** (outside the repository; see `baseline.md`):

  ```
  python -m benchmarks.fault_evaluation.run --suite standard --output <dir>/candidate.json --work-dir <dir>/cand1.work
  python -m benchmarks.fault_evaluation.acceptance --baseline baseline_linux.json --candidate candidate.json --repeat candidate_repeat.json --output acceptance.json
  python -m benchmarks.fault_evaluation.acceptance --memory-points 10000000 --output <dir>/memory.json
  ```

## Result: `acceptance.json` passes

| SC | Case | Baseline | Candidate | Value | Limit |
|---|---|---|---|---|---|
| SC-002 | `f4_u3:oldest_value:5000` | 1.418 s, 16 fault calls | 0.139 s, 4 | 0.098 | 0.40 |
| SC-003 | `f4_u3:all_group_values:5000` (time), `model_ids` (count) | 4.610 s; 52 | 0.551 s; 4 | 0.120 | 0.15 |
| SC-004 | `f4_u3` / `f4_u0` `model_ids` | | 0.147 s / 0.139 s | 1.06 | 1.5 |
| SC-005 | `plugin:oldest_value:20000` | 3.012 s | 0.652 s | 0.217 | 0.50 |
| SC-006 | worst: `stopped_fault:oldest_value:5000` | 2.83 MB | 3.44 MB | 1.215 (+0.6 MB) | 1.20 or +16 MiB |
| SC-007 | candidate run 1 against run 2 (worst `f4_u2:model_ids:5000`) | | | 1.071 | 1.10 |
| SC-008 | standard suite, candidate | | 329 s | 329 s | 900 s |
| SC-009 | `f1_u0:oldest_value:20000`, per fault point | 0.403 s | 0.140 s | 0.347 | 0.50 |

- SC-006 compares every case of the memory pass; the worst three pass by the
  16 MiB allowance.
- The baseline noise is at most 1.086 for its cases of 1 s or more.
- Other cases for information: `cutting_faults:model_gradient:5000` 1.08 s to
  0.56 s (9 fault calls on both sides: the master fixes already removed the
  nesting).
- SC-010 was removed (C-012).

## Earlier runs (superseded)

- Windows, 2026-11-03 and 2026-11-18: see `t030-diagnosis.md` and the plan
  answers.
- Linux, 2026-11-26 (candidate `6c2bbb03`, baseline master `a942177b` with
  #424): all criteria passed with the gate of that date; superseded by the
  master merges and the reference updates of C-021.
- Linux, 2026-11-27 (candidate `a558a1bd`): SC-007 failed twice on a
  different case of one round each (1.158, then 1.121), which led to C-022.

## Large memory check (FR-010, `memory_10m.json`)

`f1_u0`, `evaluate_model` on 10 000 000 points, cold: tracemalloc peak 7.11 GB.
The call restores one batch of 10 000 000 points for its one fault set (240 MB,
the FR-010 accounting). The rest of the peak is the work arrays of the
evaluation; this check has no baseline value.
