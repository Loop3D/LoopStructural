# Phase 2 tasks

Needs Phase 1 done.

- [x] 2.1 Search for hard-coded fold types
  - [x] `grep -rn "DFI\|DISCRETE_FOLD\|DiscreteFoldInterpolator" LoopStructural packages`
  - [x] List every place that must accept the new type
- [x] 2.2 Register the type
  - [x] Add the `InterpolatorType` member
  - [x] Add the string `"P2Fold"`
  - [x] Add `interpolator_map` and `support_interpolator_map` entries
  - [x] Export the class
  - [x] Test: `InterpolatorFactory.create_interpolator("P2Fold", bb, nelements=...)` returns the class with a P2 support
  - [x] Test: `to_dict` and `from_dict` round trip
- [x] 2.3 `_model_feature_factory.py:248-250`: allow `DFI` and `P2Fold`. Warn for other types
- [x] 2.4 `_model_feature_factory.py:346`: check the second `"DFI"` default
- [x] 2.5 `geological_model.py:784`: update the docstring. Keep the default
- [x] 2.6 `_folded_feature_builder.py:75`: take the type from the source builder
- [x] 2.7 Check `FoldedFeatureBuilder` with the new class
  - [x] `interpolator.fold` is set
  - [x] `fold_weights` pass through
  - [x] Faults and regions work (`_add_faults`)
- [x] 2.8 Log tets and unknowns when the P2 support is built. Document `nelements` for P2
- [x] 2.9 End-to-end script: small synthetic model, `create_and_add_folded_fold_frame` with `P2Fold`, evaluate on a grid, save and load
- [x] 2.10 Run the full test suite for `LoopStructural` and `packages/loop_interpolation`. Confirm no DFI change
