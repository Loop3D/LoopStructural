# Gradient output prerequisite

The prerequisite parent revision is `b22517f984971b98f6d039f3dd3cde0b31b7600d`.
Before the fix, `GeologicalFeature.evaluate_gradient` reused its N x 3 output
array `v` for the N x 4 fault-coordinate values. The resulting assignment
failed with raw shapes `(2,3) -> (2,4)`. The four regression cases in
`tests/unit/modelling/test_faulted_gradient_output.py` failed before this fix.

The fix is limited to naming the temporary fault-coordinate array
`fault_values` and using it for both sign comparisons. The existing N x 3 `v`
allocation and final masked gradient assignment remain in place. No geometry,
mask, fault ordering, refinement, step, interpolation, solver, or signature
behavior is changed.

This correction establishes valid faulted-gradient outputs before performance
work. It is the corrected pre-performance baseline prerequisite; it does not
claim parity against the affected gradient outputs from the raw parent, which
could not be produced because of the shape error.

Post-fix check: 23 tests passed in 3.12 seconds with this command:

```powershell
python -m pytest tests/unit/modelling/test_faulted_gradient_output.py tests/unit/modelling/test_geological_feature.py tests/unit/modelling/test_unconformity_feature.py tests/unit/modelling/test_faults_segment.py tests/unit/test_public_api_contract.py -q
```

The shared development Python environment was used with PYTHONPATH set to
this worktree and its loop_common/loop_interpolation source folders. Ruff
reports the same single pre-existing I001 finding in the changed module.
No lint rule was disabled; the full repository lint gate remains open.
