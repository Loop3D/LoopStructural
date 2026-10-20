# Phase 5 tasks

Needs Phase 4 done.

- [ ] 5.1 Confirm the problem by running code
  - [ ] Build a model with `buffer=0.1` and read the support `origin` and `maximum`
  - [ ] Compare with the box plus 0.1, 0.2 and 0.3 times the longest side
  - [ ] Check the same for the code at `_model_feature_factory.py:164, 352, 500`
- [ ] 5.2 Ask the user for the new default total buffer
- [ ] 5.3 Write a failing test: mesh extent equals box plus the stated buffer
- [ ] 5.4 Remove `with_buffer` in `_model_feature_factory.py` (lines 164, 253, 352, 500)
- [ ] 5.5 Pass the user `buffer` to `GeologicalFeatureBuilder` and on to `create_interpolator`
- [ ] 5.6 Check other callers: `LoopStructural/interpolators/_api.py:24`, the fault builder (`fault_buffer`, `fault_bounding_box_buffer`), `_geological_feature_builder.py:485`
- [ ] 5.7 Review golden diffs (mesh size changes). Update only after reading them
- [ ] 5.8 Document `buffer=0.0` behaviour (NaN on the boundary)
- [ ] 5.9 Write a release note. State that mesh size and results change
- [ ] 5.10 Commit as `fix:` (or `fix!:` if the user wants a breaking-change mark)
