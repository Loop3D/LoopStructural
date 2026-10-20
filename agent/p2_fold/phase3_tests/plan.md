# Phase 3: tests

## Goal
Pin the behaviour of `P2Fold` so that later work (the matrix-free operator and the Rust port) has an oracle.

Tests live in `packages/loop_interpolation/tests/` and in the golden suite `tests/golden/`.

## Tests

### Unit
1. Quadrature: exact for all quadratics on the reference tet. Weights sum to 1. (Written in Phase 1.)
2. Outside points give NaN. (Written in Phase 1.)
3. Weight rule: for the same fold data, the total energy of the fold family is the same for `centroid` and `gauss4` (within a small tolerance), because of `w * sqrt(q_i)`.
4. Weights have an effect: change `fold_weights` from 10 to 50. The field moves by more than a fixed limit. This catches bug 4 returning.
5. Fold norm sign: the sign matches `DiscreteFoldInterpolator` for the same fold event. Use `resolve_fold_norm_target`.
6. `update()` keeps `fold_weights` and `sample_rule`.
7. NaN rows: a fold event that gives NaN at some points still builds. The count of dropped rows is logged.
8. Dot-product test, only if you build the constraints as a `LinearOperator`: `<Ax, y> = <x, ATy>`.

### Golden
1. Add `P2Fold` to `tests/golden/engines.py`.
2. Run it on `fold_sparse` and `refolded_sparse`. Save the golden files (field and constraint fingerprints).
3. Set bounds in `test_golden.py::BOUNDS` about 1.5 times above what the code gives (see golden README).
4. Add `P2Fold` to the figures in `tests.golden.run`.

### Accuracy
1. On `fold_sparse`, the median gradient direction error of `P2Fold` is not worse than `DFI` at equal unknowns. Use a fixed seed and an explicit buffer for both.
2. Constant fold angle, no orientation data: the median angle to the expected direction is below a fixed limit.
3. `gauss4` is better than `centroid`. Check this on the golden cases. If it is not, record it and do not make it a hard test.

### Integration
1. `create_and_add_folded_fold_frame(interpolatortype="P2Fold")` does not fall back to DFI (no warning, `type(interpolator)` is `P2FoldInterpolator`).
2. Save and load keeps the type.

## Notes
- Mesh sizes in tests must be small. The P2 golden mesh is 8^3 cells (4.4k unknowns). Keep the whole suite near 20 s plus the new cases.
- Do not compare numbers between separate runs without a fixed seed.
- The refolded case has a 2 deg floor from the fold event model (see golden README). Do not set a bound below it.

## Done when
- All new tests pass on a clean checkout.
- Golden diff is reviewed.
- The golden README lists P2 fold as covered.
