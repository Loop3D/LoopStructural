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
