# Phase 1 tasks

Needs Phase 0 done.

- [ ] 1.1 Read the P2 support: `get_quadrature_points`, `evaluate_shape_derivatives`, node ordering, `element_size`
  - [ ] Note what already exists. Write findings in a comment at the top of the new file
- [ ] 1.2 Quadrature
  - [ ] Add the `centroid` and `gauss4` rules (reuse the support function if it matches)
  - [ ] Test: the rule integrates all quadratics exactly on the reference tet
  - [ ] Test: weights sum to 1
- [ ] 1.3 `evaluate_shape_derivatives` for outside points
  - [ ] Write a failing test: a point outside the mesh gives NaN
  - [ ] Fix. Remove the `LinAlgError`
  - [ ] Make P1 and P2 return the same number of values (review item 17). Update callers
- [ ] 1.4 Create `_p2_fold_interpolator.py` with class `P2FoldInterpolator(P2Interpolator)`
  - [ ] Constructor takes the support and keeps `fold`
  - [ ] Add `fold_weights` defaults: 10, 10, 1, `sample_rule="gauss4"`
- [ ] 1.5 `setup_interpolator` for the three constraint families
  - [ ] Sample points per tet, with the tet index
  - [ ] Row = `einsum(dN, vector)`, no volume factor
  - [ ] Weight `w * sqrt(q_i)`, applied after normalisation
  - [ ] Norm target from `resolve_fold_norm_target`
  - [ ] Filter NaN rows
- [ ] 1.6 Keep the settings on `update()` (check against the Phase 0 fix)
- [ ] 1.7 Check order of rows when per-sample weights are used (bug 7). Add a test with unequal weights
- [ ] 1.8 Smoke test by hand on the `fold_sparse` golden case. Compare `centroid` and `gauss4`
- [ ] 1.9 Memory check at about 5e5 tets. Record the peak memory and solve time
