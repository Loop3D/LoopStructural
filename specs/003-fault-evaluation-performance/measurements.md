# Measurements of T030 (2026-11-26)

T030 measures the code after spec 003 against the code before it, on one
computer, with the standard suite run twice for each side.

## What was measured

- **Computer**: Intel i7-12700H, 32 GB RAM, Linux 7.0 (the same computer as the
  Windows `baseline.json`), Python 3.13.13, NumPy 2.4.6, SciPy 1.18.0. No other
  work ran during the runs.
- **Baseline (T062, spec C-019)**: branch `bench/003-linux-baseline` at
  `7b4d910c`: master `a942177b` with the fix of #424, plus the harness of this
  branch, no library change. Files `baseline_linux.json` and
  `baseline_linux_repeat.json` (846 s and 842 s).
- **Candidate**: `6c2bbb03` on `feat/003-fault-evaluation-performance`, clean
  tree. Files `candidate.json` and `candidate_repeat.json` (331 s and 332 s).
- **Memory meter**: tracemalloc peak (C-020), in all four files.
- **Commands** (outside the repository; see `baseline.md`):

  ```
  python -m benchmarks.fault_evaluation.run --suite standard --output <dir>/candidate.json --work-dir <dir>/cand1.work
  python -m benchmarks.fault_evaluation.acceptance --baseline baseline_linux.json --candidate candidate.json --repeat candidate_repeat.json --output acceptance.json
  python -m benchmarks.fault_evaluation.acceptance --memory-points 10000000 --output <dir>/memory.json
  ```

## Result: `acceptance.json` passes

| SC | Case | Baseline | Candidate | Value | Limit |
|---|---|---|---|---|---|
| SC-002 | `f4_u3:oldest_value:5000` | 1.327 s, 16 fault calls | 0.131 s, 4 | 0.098 | 0.40 |
| SC-003 | `f4_u3:all_group_values:5000` (time), `model_ids` (count) | 4.295 s; 52 | 0.524 s; 4 | 0.122 | 0.15 |
| SC-004 | `f4_u3` / `f4_u0` `model_ids` | | 0.140 s / 0.129 s | 1.09 | 1.5 |
| SC-005 | `plugin:oldest_value:20000` | 2.936 s | 0.633 s | 0.216 | 0.50 |
| SC-006 | worst: `stopped_fault:oldest_value:5000` | 2.83 MB | 3.44 MB | 1.215 (+0.6 MB) | 1.20 or +16 MiB |
| SC-007 | candidate run 1 against run 2 | | | 1.06 | 1.10 |
| SC-008 | standard suite, candidate | | 331 s | 331 s | 900 s |
| SC-009 | `f1_u0:oldest_value:20000`, per fault point | 0.393 s | 0.144 s | 0.365 | 0.50 |

- SC-006 compares 105 cases of the memory pass; every case passes, the worst
  three by the 16 MiB allowance (ratios 1.21 to 1.22, about 0.6 MB more).
- The baseline noise is at most 1.05 for the 53 cases of 1 s or more.
- SC-010 was removed (C-012). The first runs of T030 (Windows, 2026-11-03 and
  2026-11-18; Linux, 2026-11-26 before the review fixes) are superseded; see
  `t030-diagnosis.md` and the plan answers.
- Other cases for information: `cutting_faults:model_gradient:5000` 30.2 s to
  8.2 s (49 to 29 fault calls).

## Large memory check (FR-010, `memory_10m.json`)

`f1_u0`, `evaluate_model` on 10 000 000 points, cold: tracemalloc peak 7.10 GB
(about 710 bytes for each point). The call restores one batch of 10 000 000
points for its one fault set (240 MB, the FR-010 accounting). The rest of the
peak is the work arrays of the evaluation, which the code before this feature
also has; this check has no baseline value.
