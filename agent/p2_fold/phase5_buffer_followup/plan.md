# Phase 5: single-point buffer fix (after Phase 4)

## Problem
The model box is expanded twice before the mesh is built. Found by reading the code. Not yet run.

1. `LoopStructural/modelling/core/_model_feature_factory.py:253` passes `model.bounding_box.with_buffer(buffer)` to the builder. Default `buffer` is 0.1.
2. `LoopStructural/modelling/features/builders/_geological_feature_builder.py:54-59` calls `InterpolatorFactory.create_interpolator(..., buffer=kwargs.get("buffer", 0.2))`.
3. `packages/loop_common/src/loop_common/supports/_support_factory.py:48-49` calls `with_buffer` again.

`buffer` is a named argument of the factory function, so it is not in `kwargs`. The user's value only sets the first expansion. The second is always 0.2.

`with_buffer` adds `buffer * (longest box side)` to every side (`_bounding_box.py:458-459`). Defaults give 0.1 + 0.2 = 0.3 per side. With the handoff `buffer=0.2`, the total is 0.4 per side. For a flat box, the mesh volume is about 12 times the box volume, and only about 7% of the elements are inside the box.

The same pattern is at `_model_feature_factory.py:164, 352, 500`.

## Fix
- Apply the buffer in one place only: the support factory.
- Remove the `with_buffer` calls in `_model_feature_factory.py` (lines 164, 253, 352, 500).
- Forward the user `buffer` to the builder, so it reaches `create_interpolator`.
- Keep the total default near the present total (0.3), or every model changes. Decide the new default with the user.
- `buffer=0.0` gives a field that is valid only inside the box. Boundary points return NaN (related to F2). Document it.

## Risks
- Every existing model changes mesh size and results. Golden files change.
- Users who pass `buffer` get a different mesh from before.
- Release note needed. This may be a breaking change for numerical results.

## Why separate
Mixing it into the P2Fold work makes the accuracy comparison hard to read. In Phase 4 all engines use an explicit buffer to avoid the issue.

## Done when
- One `with_buffer` call is applied per interpolator build.
- A test checks that the mesh extent equals the box plus the stated buffer.
- Golden files are updated and reviewed.
- Release note is written.
