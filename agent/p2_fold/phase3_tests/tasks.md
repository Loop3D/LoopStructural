# Phase 3 tasks

Needs Phase 2 done. Run from `packages/loop_interpolation`.

## Unit
- [ ] 3.1 Weight rule: fold-family energy is the same for `centroid` and `gauss4`
- [ ] 3.2 Weights have an effect: `fold_weights` 10 to 50 moves the field
- [ ] 3.3 Fold norm sign matches `DiscreteFoldInterpolator`
- [ ] 3.4 `update()` keeps `fold_weights` and `sample_rule`
- [ ] 3.5 NaN fold rows do not stop the build. Dropped rows are logged
- [ ] 3.6 (Optional) Dot-product test if constraints become a `LinearOperator`

## Golden
- [ ] 3.7 Add `P2Fold` to `tests/golden/engines.py`
- [ ] 3.8 Run `python -m tests.golden.run --case fold_sparse --engine P2Fold`. Read the figure
- [ ] 3.9 Run the same for `refolded_sparse`. Read the figure
- [ ] 3.10 Set `BOUNDS` in `test_golden.py` (about 1.5 times the observed values)
- [ ] 3.11 `python -m pytest tests/golden --update-golden`. Review the new files only
- [ ] 3.12 Update `tests/golden/README.md`: engine list, "Not covered yet", findings

## Accuracy
- [ ] 3.13 `P2Fold` is not worse than `DFI` on `fold_sparse` at equal unknowns (fixed seed, same explicit buffer)
- [ ] 3.14 Constant fold angle without orientation data: median angle below a fixed limit
- [ ] 3.15 Compare `gauss4` and `centroid`. Record the result in the golden README

## Integration
- [ ] 3.16 Factory test: `P2Fold` does not fall back to DFI
- [ ] 3.17 Save and load test keeps the type

## Finish
- [ ] 3.18 Whole suite runs in a time you accept. Record the time
- [ ] 3.19 Run `python -m pytest` in `LoopStructural` root. No new failures
