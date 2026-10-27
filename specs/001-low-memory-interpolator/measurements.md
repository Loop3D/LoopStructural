# Measurements

## Baseline

bytes_per_dof=84.64394

Note: R1 is a 100x100x100 grid, 2000 value and 500 normal constraints, default regularisation, `solve(5000, 1e-8)`. The value is peak heap bytes divided by dof.

## After table cuts

bytes_per_dof=79.071076

## After CG state machine

bytes_per_dof=63.29426

Note: measured by the ignored test `session_r1` in `tests/memory.rs` (same R1 data, through `Problem::solve`). The test `baseline_r1` still measures the one-shot solver of the ported code, so it stays at the value under `After table cuts`. The value does not include the copy of the field that goes to JavaScript (8 bytes per dof). Same iteration count as the baseline (1054).
