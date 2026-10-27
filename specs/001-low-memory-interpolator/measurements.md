# Measurements

## Baseline

bytes_per_dof=84.64394

Note: R1 is a 100x100x100 grid, 2000 value and 500 normal constraints, default regularisation, `solve(5000, 1e-8)`. The value is peak heap bytes divided by dof.

## After table cuts

bytes_per_dof=79.071076

## After CG state machine

bytes_per_dof=63.29426

Note: measured by the ignored test `session_r1` in `tests/memory.rs` (same R1 data, through `Problem::solve`). The test `baseline_r1` still measures the one-shot solver of the ported code, so it stays at the value under `After table cuts`. The value does not include the copy of the field that goes to JavaScript (8 bytes per dof). Same iteration count as the baseline (1054).

## Inequality constraints (R2 and R3)

Measured with the ignored tests `r2` and `r3` in `tests/scale_inequality.rs` (release build, counting allocator, 100^3 grid, 2000 value and 500 gradient rows, active pair bounds).

| Problem | Pair bounds | Peak bytes | Ratio to R1 | Limit | Result |
|---|---|---|---|---|---|
| R1 (no bounds) | 0 | 63 195 988 | 1.000 | - | - |
| R2 | 200 | 63 258 388 | 1.001 | 1.2 (SC-002) | pass |
| R3 | 20 000 | 69 435 988 | 1.099 | 1.2 (SC-011) | pass |

The first run had a ratio of 1.280 for R3 (fail). Two changes fixed it: ADMM no longer keeps a copy of the data right-hand side (it rebuilds it each outer iteration), and the constraint rows are reserved before they are added (no `Vec` doubling). All runs converge in 1062 iterations.

## Reference problems

Measured on the development machine (Windows, release build, one run). Tests: `tests/reference_problems.rs`. Run the ignored tests with `cargo test --release --test reference_problems -- --ignored --nocapture --test-threads=1`.

| Item | Result | Limit | Result |
|---|---|---|---|
| SC-001: R1 (100^3, 2000 value and 500 gradient rows) peak | 63 195 988 bytes | 100 MB | pass |
| SC-007: R1 solve time (desktop, report only) | 53.5 s and 57.2 s | 60 s (Chromium, reference machine) | report only |
| SC-008: 250^3, single precision, 2 GB budget, peak | 456 131 092 bytes | 1.5 GB | pass |
| SC-008: time of that solve | 2907 s | none | report only |
| SC-012: peak against dof at 50^3, 100^3, 150^3 (R^2) | 0.99986 | 0.99 | pass |
| SC-009: estimate against peak, 10^3 to 200^3 (15 sizes) | 1.049 to 1.122 | 1.0 to 1.25 | pass |

SC-009 table (500 value and 125 gradient rows, 20 iterations):

| n | estimate | peak | ratio |
|---|---|---|---|
| 10 | 286 077 | 255 060 | 1.122 |
| 20 | 957 057 | 892 340 | 1.073 |
| 40 | 5 307 297 | 5 057 972 | 1.049 |
| 60 | 15 915 777 | 15 032 756 | 1.059 |
| 80 | 35 479 137 | 33 525 684 | 1.058 |
| 100 | 66 694 017 | 62 770 612 | 1.063 |
| 120 | 112 257 057 | 105 263 540 | 1.066 |
| 150 | 213 403 317 | 200 613 652 | 1.064 |
| 200 | 491 923 617 | 459 732 404 | 1.070 |

Notes:

- The memory model has three parts: about 53 bytes for each node, about 1020 bytes for each n^2 (shell and border scratch), and about 190 bytes for each stored constraint row. The first estimate had a low surface term and failed SC-009 from n = 120.
- The R1 time of 53 to 57 s is close to the 60 s limit of SC-007. The limit is for Chromium on the reference machine and is checked by T052.
- The field of R1 for the browser comparison is written to `packages/loop_interpolation_core/target/r1_field.f64` (8 000 000 bytes, 1062 iterations, converged).
