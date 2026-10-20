# Review: interpolation code hardening and matrix-free Rust port

Status: review only. No code was changed. No test was run. All findings come from reading the code.

Scope: shared path (`_discrete_interpolator.py`, `_solver_strategy.py`, `_solver_pipeline.py`), P1, P2, FD, DFI, fold setup.
Package: `packages/loop_interpolation/src/loop_interpolation/`
Related: `HANDOFF_p2_fold.md` (the P2 fold interpolator plan).
Update: step 1 of section 4 (golden tests) is done, see `packages/loop_interpolation/tests/golden/`. It confirmed bugs 3, 4 and 11 by running code, and found three new issues (F1 to F3), listed in its README. The "No test was run" statement above describes the original review only.

## 1. Corrections to the first review of the P2 fold handoff

1. The volume factor does nothing. `add_constraints_to_least_squares` divides each row of A and B by the row norm (`_discrete_interpolator.py:361-366`). A volume factor on both A and B cancels exactly. The prototype and DFI both multiply by volume. It has no effect.
2. The rule "weight `w/n_samples` per sample" is wrong. The weight multiplies the residual (`:697-698`). The energy is therefore proportional to w². With n samples of weight w/n, the energy is w²/n. Gauss4 has 1/4 of the centroid strength. The correct split is `w*sqrt(q_i)`, where the quadrature weights q_i sum to 1. Re-run the handoff table with this fix. The balance between fold constraints and regularisation changes.

## 2. Bugs found

### High (wrong results, no error)

| # | Where | Problem |
|---|---|---|
| 1 | `_discrete_interpolator.py:835` | `add_equality_block` is never called. Equality constraints have no effect. Users: `_fault_builder.py:693`, `_geological_feature_builder.py:336`. |
| 2 | `_solver_strategy.py` | Inequality constraints work only with ADMM. CG and LSMR build `Q` and ignore it. No warning. |
| 3 | `_solver_strategy.py:66-71` | CG returns `True` after an iteration-limit warning. The model is marked up to date. |
| 4 | `_p2interpolator.py:122,151,166` | P2 multiplies weights into A and B before row normalisation. `gpw`, `npw`, `tpw` have no effect in P2. |
| 5 | `_p2interpolator.py:210` | `np.all(neighbours > 3)` is wrong. Border neighbours are negative. It should be `>= 0`. Tets next to elements 0 to 3 are dropped by mistake. |
| 6 | `_discrete_interpolator.py:357-359` | One NaN row removes the whole constraint family (warning only). DFI does not filter NaN rows. The prototype does. Check S0 and S2 logs for "Constraints contain nan". |
| 7 | `_discrete_interpolator.py:349` | `np.tile(w, n)` gives a component-major layout. Rows are point-major. Per-point weights go to the wrong rows. No effect today because the weights are uniform. |

### Medium

| # | Where | Problem |
|---|---|---|
| 8 | `_discrete_fold_interpolator.py:192,200` | `fold_normalisation` scales B and is also the row weight. Any value other than 1 changes the target. |
| 9 | P1, P2, FD | Weight conventions differ. FD value constraints use the per-point weight column; P1 and P2 do not. FD splits `npw` by dimension; P1 and P2 do not. Default `npw`: FD 10, P1 1. The weight is a residual multiplier, not an energy weight, and no document says so. |
| 10 | P1 `:287`, FD `:379,513` | The caller's `vectors` and `points[:, 3:6]` are changed in place. FD hard-codes columns 3:6. |
| 11 | `_discrete_interpolator.py:1337` | `update()` calls `setup_interpolator()` with no arguments. `fold_weights`, `regularisation`, `directional_regularisation` are lost. A re-solve can differ from the first solve. |
| 12 | `_pre_solve` | `c` is set to NaN. A failed solve leaves it NaN with `up_to_date=False`. Every evaluation retries the solve. |
| 13 | `_post_solve` | Commented out. All constraint matrices stay in memory after the solve. This matters at 500k+ tets. |
| 14 | `:376-380` | The duplicate-name code calls `int()` on the second `_` token. A name such as `gradient_orthogonal` crashes on its second use. |

### Low

15. P2 `wtfunc` replaces `w` and receives element indices, not points (`_p2interpolator.py:274-275`).
16. P2 steepness: off-diagonal Hessian terms count once. A Frobenius norm counts them twice.
17. `evaluate_shape_derivatives` returns 3 values in P1 and 2 in P2.
18. `ridge_factor=1e-8` becomes 1e-16 in the normal equations. This is below float precision. The term has no effect.

## 3. Problems for a matrix-free Rust port

1. Row normalisation defines the weights. Each regularisation row is divided by its own norm. The result is not a mesh integral. It does not change with element size. The value-constraint weight changes by up to 2x across one element (the shape-function norm is between 0.5 and 1). Decide the intended energy first, for example the integral over each element. Implement that in Rust. Do not copy the current behaviour.
2. FD already has a partial matrix-free path (`_assemble_operator`, `_build_fused_cg_regularisation_operator`). It covers regularisation only. Data constraints stay as a sparse matrix (`AtA_data`). P1 and P2 have none. The fused FD kernel is not exact near the boundary, so the code patches the 2-cell shell. Reuse the design (stencil kernels, Jacobi preconditioner). Do not reuse the patch.
3. Both paths solve the normal equations `AtA x = Atb`. This squares the condition number. Use one of:
   - CG on the normal equations with a Jacobi preconditioner. Column scaling does this today. The diagonal of AtA is easy to compute matrix-free.
   - LSQR or LSMR on the operator (no squaring).
4. Tet meshes are the hard part. The operator is a sum over elements and faces. Edge jumps need face-neighbour tables. Define each operator as `apply(x)` and `apply_transpose(y)`. Test every pair with the dot-product check: <Ax, y> = <x, ATy>.
5. Point location and sampling are in Python. `evaluate_shape_derivatives` raises `LinAlgError` for points outside the mesh. Rust must return NaN or a mask. For inside points, it must give the same answer as Python.
6. DFI uses `rng.shuffle` and `[::step]`. Runs are not repeatable. Rust needs a deterministic rule.

## 4. Proposed order of work

1. Pin the current behaviour. Write golden tests (small meshes, fixed seeds) for P1, P2, FD and DFI. Save the matrix for each constraint family. These tests become the Rust parity tests.
2. Fix bugs 1 to 7 and 14. Fix each in its own commit. Write a failing test first. Keep the golden files to show exactly what changes.
3. Write the operator specification: one page per family. State the energy, the weight convention (w or w squared), and the row-normalisation rule. Settle items 9 and 10 here.
4. Add a matrix-free `LinearOperator` for P1 and P2 in Python, with the transpose test. Check it against the assembled matrix to about 1e-10.
5. Port to Rust behind the same operator interface. Use the assembled matrix as the oracle.
6. Build the P2 fold interpolator after step 4. It then needs only a quadrature loop. This avoids writing it twice.

## 5. Not checked

- Support code was not read in detail: shape functions, `element_size`, `get_quadrature_points`, P2 node ordering.
- `constraints.py` and the `get_*_constraints` code were not read.
- No test was run. Items 5, 6 and 7 come from reading the code. Confirm item 6 from the logs.
