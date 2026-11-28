# Phase 0 tasks

Run tests from `packages/loop_interpolation`. One commit per task.

- [ ] 0.1 F1: store `origin`, `step_vector`, `nsteps_cells` on `P2UnstructuredTetMesh`
  - [ ] Write a failing test that reads the three attributes
  - [ ] Fix the constructor
  - [ ] Remove the manual attribute setting in `tests/golden/engines.py`
  - [ ] Confirm the diagnostics report `inside()` works for P2
- [ ] 0.2 Bug 5: change `neighbours > 3` to `neighbours >= 0` in `_p2interpolator.py:210`
  - [ ] Write a failing test: tets next to elements 0 to 3 are kept
  - [ ] Fix, then review the P2 golden diff
- [ ] 0.3 Bug 4: apply P2 weights after row normalisation
  - [ ] Change `test_known_p2_normal_weight_has_effect` into a normal test (remove xfail)
  - [ ] Fix `_p2interpolator.py:122,151,166`
  - [ ] Check that `npw` 1 to 100 now moves the field by more than 1e-4
  - [ ] Review the P2 golden diff
- [ ] 0.4 Bug 3: CG iteration limit returns `False`
  - [ ] Change `test_known_default_solve_is_converged` into a normal test
  - [ ] Fix `_solver_strategy.py:66-71`
  - [ ] Search for callers that use the return value. Check that `up_to_date` stays correct
  - [ ] Decide the P2 default tolerance or iteration count. Record it in the commit message
- [ ] 0.5 Bug 11: `update()` keeps `fold_weights`, `regularisation`, `directional_regularisation`
  - [ ] Change `test_known_update_keeps_fold_weights` into a normal test
  - [ ] Fix `_discrete_interpolator.py:1337`
  - [ ] Review DFI and FDFold golden diffs
- [ ] 0.6 Bug 6: filter NaN rows, not the whole family
  - [ ] Write a failing test with one NaN row in a constraint family
  - [ ] Fix `_discrete_interpolator.py:357-359`
  - [ ] Add the same filter to the shared fold setup (`_fold_setup.py`)
  - [ ] Log the count of dropped rows
- [ ] 0.7 Run the full golden suite. Update the `Known bugs` table in `tests/golden/README.md`
