# Verification: native fault application

This file collects the evidence of spec 007.

## C-007: proof-of-concept difference

Sources: the proof of concept (PoC) in `/home/lgrose/dev/LoopStructural-wt/_runs003/rust_poc/` (read only), and the Python reference in this worktree. Paths of PoC files are relative to that folder. Line numbers are those of the files when this was written.

### Summary

The most likely cause is the linearisation of the gy and gz profiles in the exporter.

- The default gy and gz profile of the fault is `smooth_peak` (`_fault_function.py:12`). It is a plain function, not a `CubicFunction`.
- `export_snapshot.py:138` (`linearise`, line 60) does not export it as a formula. It samples it at 2001 points over the node range of the frame coordinate and the PoC interpolates between the samples.
- A piecewise linear curve has an error of about `h^2 * |f''| / 8` for sample step `h`. I tested this with `np.interp` on `smooth_peak`. The maximum error is 1.0e-6 for a span of 2, 4.0e-6 for a span of 4, and 8.9e-6 for a span of 6. The span is the node range of gy or gz.
- This is a relative error of 1e-6 to 1e-5 on the displacement. For a displacement of 40 to 120 m it gives about 4e-5 to 1e-3 m. The measured differences are 5.4e-5 m (f4_u0) and 1.9e-4 m (10 km box). The size and the scaling agree.
- The error exists with one fault: model f1_u0 has 1.4e-5 m (`poc_bench.json`). A step-order or fault-order cause would need two or more faults. A f32 cause is not possible (see the table).
- The cause was not proved by a run of the PoC. I tested the size of the error only.

Other items are far smaller (rounding, about 1e-16 relative) or change only special points (NaN, cutting faults, regions).

### Table

The items are in order of the likely size of their effect.

| # | Item | Python behaviour | PoC behaviour | Effect on the result |
|---|------|------------------|---------------|----------------------|
| 1 | Profile of gy and gz (`smooth_peak`, other plain functions) | Exact formula, `_fault_function.py:12-17`; used in `FaultDisplacement.__call__` (line 304) | Sampled at 2001 points and linearly interpolated: `export_snapshot.py:60-95`, `profile.rs:117` (`Linear::eval`) | Relative error 1e-6 to 1e-5 on the displacement, about 1e-4 m on 100 m. Probable main cause. Also, the curve has a kink at the sample points and at gy = +-1 |
| 2 | Slope and curvature of the continuation of a linearised profile | `_taylor`, `_fault_function.py:357-377`: finite differences of the exact function | Exact values at six points, saved as `stencil`: `export_snapshot.py:79-92`, `profile.rs:255-300` | None when the stencil is used. It affects only a block with `BLOCK_EXTENSION = "smooth"`; the restoration without a block does not use it |
| 3 | Fault with a frame cut by a younger fault (#421) | If the gradient is NaN, it uses the gradient of the frame without the cutting faults: `_fault_segment.py:716-724` | No fallback. A NaN gradient stays NaN: `fault.rs:198-207`. The comment in `README.md` says the frame is read at restored points | None on the NaN points that agree. A point on a cutting fault can differ (NaN in the PoC, finite in Python). The test models have no such point |
| 4 | Regions of the frame members | `_calculate_mask` sets values outside a region to NaN: `_fault_segment.py:422-427` | The exporter refuses a region that is not an erosional unconformity: `export_snapshot.py` (`check_regions`, `grid_to_dict`, line 195-210) | None on exported models. The PoC cannot read a model with an abutting fault |
| 5 | Point that is not finite | `v[nan_rows] = nan` for values, NaN gradient (#413): `_fault_segment.py:425-445` | A NaN local coordinate gives the value 0.0 and the gradient 0: `grid.rs:250`, `grid.rs:274-277` | The point does not move in the PoC. In Python the value is NaN, so the point also does not move (d = 0). The NaN points agreed in the tests |
| 6 | Step direction (#425) | The step moves along the gradient of frame coordinate 1 (`member1.evaluate_gradient`, `_fault_segment.py:460`), made unit length at each step, and the sign of `reverse` is applied to `d`: lines 732-741 | Same: `g1` is used (`fault.rs:148-158`), unit length, `d` multiplied by -1 for `reverse`: `fault.rs:196-212` | None. The direction agrees |
| 7 | Number of steps and the points that move | The set of moving points is fixed by `|d| > 0` before step 0: `_fault_segment.py:695`; all later steps keep it | `d0 == 0` returns the point before step 0, so the set is the same: `fault.rs:181-190` | None. A point that has `d = 0` at a later step does not move in either path |
| 8 | Order of the steps and of the faults | `restore_through`: `cutters_first(with_cutters(faults))`, `_fault_block.py:242-249` | The exporter writes the same order: `export_snapshot.py` (`cutters_first(with_cutters(...))`) | None. Each fault is applied once, in the same order |
| 9 | Local frame (`_project_points`, rotation) | `bounding_box.project` and `reproject_vectors`: `_geological_interpolator.py:253-275` | One affine map `local = A * world + b`, built by the exporter from `project` of a unit basis: `export_snapshot.py:213-218`, `grid.rs:152-158`. A rotated grid is refused | Rounding only, about 1e-16 relative. No effect at 1e-8 times the box size |
| 10 | Cell search and trilinear sum | `_lookup_points`, `_structured_grid_3d.py:515-547`: cell by `position_to_cell_index`, clip of local coordinates to [0, 1] | `lookup`, `grid.rs:164-215`: same cell rule, clamp to [0, 1], tolerance `32 * eps * max(extent, 1)` | Rounding only, about 1e-16 relative. The order of the 8 corner terms (`sum8`, `grid.rs:18`) can differ from `_values_from_coefficients` by 1 ulp |
| 11 | Cubic profile | `CubicFunction.__call__`, `_fault_function.py:163-171`, `v**3`; the limit value comes from `_limit_value` (line 152) | `Cubic::eval`, `profile.rs:79-98`, `powf(3.0)`; `limit_value` is exported from `_limit_value` | 1 ulp. The order (limit, then floor and ceiling) is the same |
| 12 | `BLOCK_EXTENSION`, `_saturate`, `continue_profile` | `_fault_function.py:345-431`; `_fault_block.py:39-61` | Constants are read from the manifest; `saturate` in `profile.rs:247` is the same formula | None, if the manifest holds the values 0.25, 0.5, 1e-3, 1e-9. This is an item to check at export |
| 13 | Single or double precision | All arrays are float64 | `<f8` in the files (`export_snapshot.py:177`), `f64` in all Rust code | None. No float32 exists in either path |
| 14 | Version of the fault code | Spec 007 code, with #421, #425, #427 and block extension | The crate was written for an older version (kiridashi `loop-snapshot`) and then changed by hand | Unknown. The risk is that a rule added after the copy is not in the crate. Items 3 and 12 are the known cases |

### Rules for the native kernel

1. Do not linearise a profile. Give the kernel the exact formula of each profile: `smooth_peak` as a polynomial with a zero value for |x| >= 1, `CubicFunction`, `Composite`, `Ones`, `Zeros`. If a profile has no known formula, refuse the model; do not sample it.
2. Keep the 8-term trilinear sum and the cell rule of Python (`_lookup_points`), including the clip of local coordinates to [0, 1].
3. Compute the step in float64 only, in the order of Python: the frame at the current point, `d` from the profile, unit gradient of coordinate 1, `newp += g * (1 / steps) * d`.
4. Fix the set of moving points before step 0 (`|d| > 0`). Do not change it in later steps.
5. Apply the faults in the order of `restore_through` (`cutters_first(with_cutters(...))`), and give each fault the list of faults already restored.
6. Use the fallback gradient of #421: if the gradient is NaN and the point is finite, use the gradient of the frame without the cutting faults. Do not copy the PoC behaviour here.
7. A point that is not finite gives NaN values and a NaN gradient, as in Python (#413). Do not return 0.0 for it.
8. Apply the region mask of each frame member (value NaN outside the region), or refuse the model. Do not ignore a region.
9. Read `BLOCK_SIDE_EPSILON`, `BLOCK_EXTENSION`, `BLOCK_EXTENSION_EXACT`, `BLOCK_EXTENSION_LIMIT` and `BLOCK_EXTENSION_STEP` from `_fault_block.py` at the time of use, and write them in the manifest. Do not copy the numbers into the Rust code.
10. Build the local frame from `bounding_box.project`, as the PoC does, but also accept a rotated grid, or refuse it with a clear error.
11. Check the kernel with a model that has a fault with `smooth_peak` in gy, at the golden tolerance of 1e-8 times the box size. A pass on a model with only `Ones` profiles proves nothing about item 1.
12. To confirm the cause (open point): export the model f1_u0 with 2001, 20001 and 200001 samples and compare with Python. If the maximum difference falls with the square of the step, item 1 is the cause.

## Stage S2

Measures SC-002, SC-003 and SC-007 and records the native differences (FR-005c). Task T021.

### Set-up

- Computer: the computer of spec 003 T030. Intel Core i7-12700H (12th Gen), 20 threads, Linux.
- Commit: `dcad71c78bc58501456311cd57a66b028e84f88d` (HEAD of the worktree). Other work was in the tree (an untracked test file of another task).
- Limit of the numbers: another agent could run tests at the same time (light load). The timings can be a little high. Each timing is the median of 7 repeats (31 repeats for runs below 10 ms), after two warm-up calls.
- Time method: for each model, 10 000 points (`sample_points(name, 10000, seed=1)`). For each faulted feature (`feature.faults` not empty, `faults_enabled`):
  - Python: `restore_through(pts, feature.faults, owner=feature)` with `native_faults = False`.
  - Native: `_native_faults.restore(feature, pts)` with `native_faults = True`, after one call that builds the FaultChain. `native_threads = 1` (SC-002) and `native_threads = 0` (all 20 threads, SC-003).
  - The time of a model is the sum over its faulted features. Ratio = native / Python.
  - Every timed feature gave the same result shape and no `None`: no feature fell back to Python and none was excluded.
- Memory method: `f4_u0`, first faulted feature (`group_0`), 10 000 000 points (`np.tile` of 10 000 sample points, float64 (N, 3), 240 MB), `native_threads = 0`. One restore in a fresh process for each path. The baseline is `VmRSS` after the points exist (and, for native, after the chain is built). The peak is `ru_maxrss` after the call. Extra = peak - baseline. The Python run took about 7.5 minutes, so the full 10 million points were used.
- Scratch scripts and results: `/home/lgrose/.claude/jobs/8f5b2ef3/tmp/t021/`.

### SC-002: native restore, 1 thread, at most 20 % of Python

| Model | Python ms | Native ms (1 thread) | Ratio | Target | Result |
|-------|-----------|----------------------|-------|--------|--------|
| f1_u0 | 76.4 | 22.9 | 0.299 | <= 0.20 | fail |
| f4_u0 | 314.0 | 92.8 | 0.295 | <= 0.20 | fail |
| f4_u3 | 2148.4 | 627.2 | 0.292 | <= 0.20 | fail |
| cutting_faults | 487.9 | 112.8 | 0.231 | <= 0.20 | fail |
| plugin | 725.6 | 191.6 | 0.264 | <= 0.20 | fail |

### SC-003: native restore, 20 threads, at most 3 % of Python

| Model | Python ms | Native ms (20 threads) | Ratio | Target | Result |
|-------|-----------|------------------------|-------|--------|--------|
| f1_u0 | 76.4 | 2.2 | 0.028 | <= 0.03 | pass |
| f4_u0 | 314.0 | 8.1 | 0.026 | <= 0.03 | pass |
| f4_u3 | 2148.4 | 57.4 | 0.027 | <= 0.03 | pass |
| cutting_faults | 487.9 | 11.2 | 0.023 | <= 0.03 | pass |
| plugin | 725.6 | 17.8 | 0.025 | <= 0.03 | pass |

The margin of SC-003 is small (0.023 to 0.028 against 0.03). A load from other work can change the result for `f1_u0`.

### SC-007: extra peak memory of a native restore, 10 million points

| Path | Baseline MB | Peak MB | Extra MB |
|------|-------------|---------|----------|
| Native | 479.2 | 1164.3 | 685.1 |
| Python | 476.0 | 6912.2 | 6430.4 |

- Rule (spec 003 FR-010, `acceptance.py` `_check_memory`, `limit = 1.20`): a new value fails only if it is above `1.20 * base` AND more than `MEMORY_ALLOWANCE_BYTES` (16 MiB) above `base`. Base = the Python extra memory.
- New (native) extra is 685.1 MB. The limit value is 1.20 * 6430.4 = 7716.5 MB. Native is below base.
- Result: pass. Native uses about 11 % of the extra memory of Python.
- Note: the native extra memory is about 2.9 times the size of the point array (240 MB). The reason was not examined.

### Largest native difference (`max_native_diff`)

Source: `OMP_NUM_THREADS=1 pytest tests/fault_evaluation/test_native_faults.py -q -p no:cacheprovider --junitxml=.../native.xml`. Result: 10 passed. The test uses 2000 points and compares many variants (default, reverse, block sides, `restored_through`). The box size is the largest box extent.

| Model | Max difference (m) | Box size (m) | Difference / box | Fallback cases |
|-------|--------------------|--------------|------------------|----------------|
| f1_u0 | 2.27e-13 | 1000 | 2.3e-16 | none |
| f4_u0 | 2.27e-13 | 1000 | 2.3e-16 | none |
| f4_u3 | 2.27e-13 | 1000 | 2.3e-16 | none |
| cutting_faults | 5.68e-13 | 1000 | 5.7e-16 | `group_0` with `reverse=True` uses Python |
| plugin | 9.31e-10 | 29 656.82 | 3.1e-14 | none |
| stopped_fault | 1.14e-13 | 1000 | 1.1e-16 | none |
| faulted_fold | 1.82e-12 | 10 000 | 1.8e-16 | none |

All differences are far below the golden tolerance of 1e-8 times the box size. The PoC difference of C-007 (5.4e-5 m) is not present.

### Result

- SC-002: FAIL for all 5 models. The ratio is 0.23 to 0.30 against a target of 0.20. The limit was not changed.
- SC-003: PASS for all 5 models. The ratio is 0.023 to 0.028 against a target of 0.03.
- SC-007: PASS. Native extra memory 685 MB, Python 6430 MB.
- FR-005c: the largest native difference is 9.3e-10 m (`plugin`), or 3.1e-14 of the box. For the 1000 m models it is 5.7e-13 m or less.
- One fallback case: `cutting_faults`, `group_0`, `reverse=True` (the native path returns `None`, Python is used). No fallback in the timed (non-reverse) calls.

## Stage S2 after T035

Task T035. It makes the 1-thread native restore faster, with the same results.

### What changed

- `FaultChain::apply` reads the frame point (the point restored through the cutters) one time for each step. Before, it read it up to four times for each step. It also keeps the displacement fraction of step 0 instead of computing it again.
- `FaultChain::restore` and `Region::Sign::contains` use a stack array (up to 64 faults) in place of `to_vec()` and `vec![false; n]`. There is no allocation for each call.
- The floating-point operations and their order are the same. No public signature changed.

### Profile

- `perf record` is not allowed on this computer (`perf_event_paranoid`), so there is no function profile.
- The costs came from the code: the two allocations for each point and the repeated `frame_point` call (a full restore through the cutters). Both are removed.
- Remaining costs (not measured by a profile): the grid look-ups of `GridField` (`lookup` is done again for `value` and for `gradient` of the same point), and the cutter restores of `Region::Unconformity` and `Region::Sign`. `grid_field.rs` is not in the Files of T035. If SC-002 must pass, add it to a new task.

### Set-up

- Same method as "Stage S2" (10 000 points, median of 7 or 31 repeats, `native_threads = 1` and `0`). Script copy: `/home/lgrose/.claude/jobs/8f5b2ef3/tmp/t035/time.py`.
- Limit of the numbers: the computer had a load average of about 11 (other agents ran tests). The Python times are about 1.8 times the T021 times, and the 20-thread run is hurt most by the load. Measure again on a quiet computer before you trust a pass or fail near the target.

### SC-002: native restore, 1 thread, at most 20 % of Python

| Model | Python ms | Native ms (1 thread) | Ratio | Target | Result |
|-------|-----------|----------------------|-------|--------|--------|
| f1_u0 | 137.6 | 37.0 | 0.269 | <= 0.20 | fail |
| f4_u0 | 752.2 | 144.8 | 0.192 | <= 0.20 | pass |
| f4_u3 | 4654.8 | 1050.0 | 0.226 | <= 0.20 | fail |
| cutting_faults | 907.7 | 184.7 | 0.204 | <= 0.20 | fail |
| plugin | 1228.6 | 288.0 | 0.234 | <= 0.20 | fail |

Native time of T021 (quiet computer) against now (loaded computer): the ratios went from 0.23-0.30 to 0.19-0.27.

### SC-003: native restore, 20 threads, at most 3 % of Python

| Model | Python ms | Native ms (20 threads) | Ratio | Target | Result |
|-------|-----------|------------------------|-------|--------|--------|
| f1_u0 | 137.6 | 6.2 | 0.045 | <= 0.03 | fail (load) |
| f4_u0 | 752.2 | 19.6 | 0.026 | <= 0.03 | pass |
| f4_u3 | 4654.8 | 136.6 | 0.029 | <= 0.03 | pass |
| cutting_faults | 907.7 | 21.9 | 0.024 | <= 0.03 | pass |
| plugin | 1228.6 | 32.1 | 0.026 | <= 0.03 | pass |

### Largest native difference (`max_native_diff`)

Source: `OMP_NUM_THREADS=1 pytest tests/fault_evaluation/test_native_faults.py --junitxml=...`. 10 passed.

| Model | Max difference (m) | Box size (m) | Difference / box |
|-------|--------------------|--------------|------------------|
| f1_u0 | 2.27e-13 | 1000 | 2.3e-16 |
| f4_u0 | 2.27e-13 | 1000 | 2.3e-16 |
| f4_u3 | 2.27e-13 | 1000 | 2.3e-16 |
| cutting_faults | 1.25e-12 | 1000 | 1.3e-15 |
| plugin | 9.31e-10 | 29 656.82 | 3.1e-14 |
| stopped_fault | 1.14e-13 | 1000 | 1.1e-16 |
| faulted_fold | 1.82e-12 | 10 000 | 1.8e-16 |

The `cutting_faults` value is higher than in T021 (5.68e-13). The Rust code computes the same operations; the Python side (sample points or solver) is the likely reason. This was not examined. All values are far below 1e-8 times the box.

### Result

- See the quiet re-measure below: SC-002 fails for all 5 models; SC-003 passes for all 5 models.
- Rust tests pass (`cargo test --release`, also `--features parallel`). Clippy passes for no features, `parallel`, `python`, `python,parallel`.
- Python tests pass: `test_native_faults.py`, `test_native_safety.py`, `test_golden.py` (43 passed).
- SC-002: 1 of 5 models passes (`f4_u0` 0.192), 4 fail (0.204 to 0.269), under load. Not met for all models. (measured under load; replaced by the quiet re-measure below)
- SC-003: 4 of 5 models pass; `f1_u0` (0.045) fails under load. T021 had 0.028 for this model on a quiet computer. (measured under load; replaced by the quiet re-measure below)
- FR-005c: the largest native difference is 9.3e-10 m (`plugin`), 3.1e-14 of the box.

### Quiet re-measure (orchestrator)

Measured on a quiet computer (load average below 3, nothing else running, 2026-11-27 15:53). Same script and method as T021: medians of 7 or more repeats, 10 000 points.

| Model | Python ms | Native 1 thread ms | Ratio 1 thread (SC-002, target 0.20) | Native 20 threads ms | Ratio 20 threads (SC-003, target 0.03) |
|---|---|---|---|---|---|
| f1_u0 | 77.7 | 23.1 | 0.298 fail | 2.1 | 0.0275 pass |
| f4_u0 | 305.1 | 88.2 | 0.289 fail | 8.3 | 0.0273 pass |
| f4_u3 | 2157.1 | 634.4 | 0.294 fail | 58.0 | 0.0269 pass |
| cutting_faults | 493.4 | 113.8 | 0.231 fail | 11.1 | 0.0225 pass |
| plugin | 726.3 | 194.8 | 0.268 fail | 17.8 | 0.0244 pass |

The ratios are the same as in T021 (Stage S2) within the noise. The changes of T035 did not make the 1-thread restoration measurably faster. The better ratios above this sub-section came from the load on the computer, which made the Python times longer. SC-002 fails for all 5 models. SC-003 passes for all 5 models.

## Stage S3

Measures SC-004 (many-faults model) and SC-005 (golden test with the native path on). Task T028.

### Set-up

- Computer: the computer of spec 003 T030. Intel Core i7-12700H, 20 threads, Linux. Load average was 1.8 to 3 during the runs (another agent edited ROADMAP.md only).
- Commit: `d12101b223acf7c758a410ea0d6cdc978450b3eb` (HEAD of the worktree, kernel built, no code changed by this task).
- Native path on: `LOOPSTRUCTURAL_NATIVE_FAULTS` not set. A separate check gave `native_enabled()` = True in this venv. In the 94-fault model, `chain_for(...)` is not `None` for 3 of 3 faulted features (the features in `model.features` that have `faults`).
- Method: the model of `benchmarks/fault_evaluation/many_faults.py` (94 faults, built once), points from `many_faults.points(n)`, `evaluate_model` timed on the first call for each point count. 500, 2000 and 5000 points. Native build time: 153 s. Python build time: 176 s.
- Problem found: `python -m benchmarks.fault_evaluation.many_faults` stops with the native path on. The script times the calls of `FaultSegment.apply_to_points`. The native path does not call it, so the lists of per-fault times are empty and `np.max` raises `ValueError`. No JSON is written. The benchmark script was not changed (this task changes no code). The native times come from a scratch script (`/home/lgrose/.claude/jobs/8f5b2ef3/tmp/t028/native_run.py`) that imports `build` and `points` from the benchmark and times `evaluate_model` in the same way. The JSON is `/home/lgrose/.claude/jobs/8f5b2ef3/tmp/t028/many_faults_native.json`. The first call of each point count is the value in the table; the median of 5 repeated calls is 0.053, 0.056 and 0.060 s.
- Python run: the benchmark script with `LOOPSTRUCTURAL_NATIVE_FAULTS=0`, same venv, `many_faults_python.json` in the same folder. This is the Python path after stage S1 of spec 007.
- Baselines: (a) same venv, spec 007 S1, quiet computer: 0.615 s at 500 points, 0.894 s at 5000 points (2.05 s at 5000 points before S1). (b) `_runs007/many_faults_python.json`: 2.09 s at 5000 points. This is the master Python path before spec 007.

### SC-004: `evaluate_model` at 5000 points, at most 0.5 s

| Points | Native s | Python s (same run conditions) | Ratio native / Python | Python s (S1 baseline) | Master s (_runs007) |
|---|---|---|---|---|---|
| 500 | 0.059 | 0.676 | 0.087 | 0.615 | - |
| 2000 | 0.056 | 0.742 | 0.076 | - | - |
| 5000 | 0.059 | 0.971 | 0.061 | 0.894 | 2.09 |

- Target: 0.5 s at 5000 points. Measured: 0.059 s. Pass. (The ratio to master, 2.09 s, is 0.028.)
- The native time is almost the same for 500, 2000 and 5000 points. The work is not limited by the point count in this range.
- Mean time of one fault call (Python run, 5000 points): cutting faults 18.3 ms (maximum 28.5 ms), isolated faults 3.83 ms. The native path does not call `apply_to_points`, so these means are not available for it.
- `fault_calls` is 94 in the Python run. It is not available for the native run (same reason).
- Native results are not bit-identical to Python (C-001), so `--compare` was not used.

### SC-005: golden test with the native path on

- Command: `OMP_NUM_THREADS=1 .venv/bin/python -m pytest tests/fault_evaluation/test_golden.py -q -n 8 -p no:cacheprovider`, no environment variable. Result: 23 passed, 0 failed (warnings only). The test prints no difference values.
- Largest field difference: `tests/fault_evaluation/test_native_field.py` gives 11 passed. The largest `max_native_diff` in the junit file is 1.16e-10 (the other values are 2.3e-13 to 1.8e-12).

### Result

- SC-004: pass. 0.059 s at 5000 points (target 0.5 s). Python in the same conditions: 0.971 s.
- SC-005: pass. 23 of 23 golden tests pass with the native path on.
- Largest native field difference: 1.16e-10.
- Open point: `benchmarks/fault_evaluation/many_faults.py` fails with the native path on (empty per-fault timer lists). It must write the JSON without per-fault times for the native path. A separate task is needed.
