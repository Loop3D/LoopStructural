# Phase 1: `P2FoldInterpolator`

## Goal
A working `P2FoldInterpolator(P2Interpolator)` with a quadrature sample rule. It is not yet registered in the factory (Phase 2).

New file: `packages/loop_interpolation/src/loop_interpolation/_p2_fold_interpolator.py`.
Templates: `_discrete_fold_interpolator.py`, `_fold_setup.py`, `_fd_fold_interpolator.py`.

## Corrections to the handoff (from the review, section 1)
- Do not multiply by tet volume. Row normalisation cancels it.
- Do not use weight `w / n_samples`. The weight multiplies the residual, so energy goes with w squared. Use `w * sqrt(q_i)`, where the quadrature weights `q_i` sum to 1. Then the total energy per tet does not depend on the rule.

## Design

1. **Quadrature.** Add two rules on the reference tet:
   - `centroid`: one point (1/4, 1/4, 1/4, 1/4), weight 1.
   - `gauss4`: four points with barycentric coordinates (a, a, a, b) and permutations. a = 0.1381966011250105, b = 0.5854101966249685. Each weight is 1/4.
   Check `get_quadrature_points` in the P2 support first. Reuse it if it matches. The review did not read it.
2. **Sample points.** Map the rule to every tet (or to every tet in the active region). Keep the tet index of each point.
3. **Shape derivatives.** `support.evaluate_shape_derivatives(points, elements=idx)`. Shape `(n, 3, 10)`. The P2 gradient is linear in a tet. One point does not control it. This is why `gauss4` beats `centroid`.
4. **Three constraint families**, as in `DiscreteFoldInterpolator`:
   - gradient . fold_direction = 0 (default weight 10)
   - gradient . fold_axis = 0 (default weight 10)
   - gradient . dgz = fold_norm (default weight 1)
   Get the shared vectors from `_fold_setup.py`. Get the norm target from `resolve_fold_norm_target`, so the sign matches DFI.
5. **Weights.** Row = `einsum(dN, vector)`. Target = `fold_norm` or 0. Weight = `w * sqrt(q_i)`. Pass the weight through `add_constraints_to_least_squares` so it is applied after normalisation.
6. **Regularisation.** Use the P2 default. Add no fold-specific terms (the handoff found no gain).
7. **Options.** `fold_weights["sample_rule"]` is `"gauss4"` (default) or `"centroid"`.
8. **NaN rows.** Filter with the Phase 0 code from bug 6.
9. **`update()`.** Must keep `fold_weights` (Phase 0, bug 11).
10. **Outside points.** `evaluate_shape_derivatives` returns NaN for points outside the mesh, not `LinAlgError`. Make the P1 and P2 return signatures match (review item 17).

## Risks
- The handoff numbers used the wrong weight rule. Results will differ. Phase 4 re-runs them.
- Per-point weights meet bug 7 (`np.tile` layout). If a per-sample weight vector goes through `add_constraints_to_least_squares`, check the row order.
- Memory: 4 samples per tet at 526k tets is 2.1M rows per family. Check memory at `nelements=1e5`.

## Done when
- The class builds, solves and evaluates on the `fold_sparse` golden case when you construct it by hand.
- The quadrature unit test passes.
- Both sample rules work, and the weight option changes the field.
