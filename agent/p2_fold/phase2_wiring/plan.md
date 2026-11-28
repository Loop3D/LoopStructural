# Phase 2: registration and wiring

## Goal
`create_and_add_folded_fold_frame(..., interpolatortype="P2Fold")` builds a model with the new class. `"DFI"` stays the default for folds.

## Decision
New type string `"P2Fold"`. Plain `"P2"` does not change.

## Items

### 1. Registration
File: `packages/loop_interpolation/src/loop_interpolation/__init__.py`, with `_interpolatortype.py`.
- Add `InterpolatorType.PIECEWISE_QUADRATIC_FOLD` (name can change; match the existing naming).
- Add the string `"P2Fold"` to the string map (near line 97).
- Add the class to `interpolator_map` (near line 112).
- Add `P2UnstructuredTetMesh` to `support_interpolator_map` for dimension 3.
- Export `P2FoldInterpolator` in `__all__`.
- Check `_normalise_interpolator_type` and `from_dict`. A saved model with this type must load.

### 2. Stop forcing DFI
Three places:
- `LoopStructural/modelling/core/_model_feature_factory.py:248-250`: accept `DFI` and `P2Fold`. Warn only for other types.
- `LoopStructural/modelling/core/geological_model.py:784`: the default stays `"DFI"`. Update the docstring.
- `LoopStructural/modelling/features/builders/_folded_feature_builder.py:75`: `from_feature_builder` hard-codes `'DFI'`. Take the type from the source builder if it is a fold type.
- Check `_model_feature_factory.py:346`, which also names `"DFI"`.

### 3. Fold builder
- Check that `FoldedFeatureBuilder` does not use DFI-only attributes (for example `interpolator.fold` and `fold_weights` pass-through).
- Check that the support type works with the fold frame, the region and the fault handling in `_add_faults`.

### 4. Element count
- `nelements` means tets for P2. 1e4 gives 56k tets and 86k unknowns. 1e5 gives 526k tets and 765k unknowns (solve about 2.5 min).
- Log the number of tets and unknowns when the support is built.
- Describe this in the docstring of `create_and_add_folded_fold_frame`.

### 5. Buffer
Do not change the buffer behaviour here (see `phase5_buffer_followup`). In docs and examples, say that the mesh is larger than the model box.

## Risks
- Serialisation (`to_dict` and `from_dict`) may miss the new type.
- Other code may test `interpolatortype == "DFI"`. Search for it.
- Code that checks `isinstance(interpolator, DiscreteFoldInterpolator)` will miss the new class. Search for it.

## Done when
- A model script with `interpolatortype="P2Fold"` builds, evaluates and saves/loads.
- `"DFI"` models give the same results as before.
