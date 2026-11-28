# Phase 0: prerequisites

## Goal
Remove the P2 and shared bugs that `P2FoldInterpolator` would inherit. After this phase, P2 weights work, P2 solves converge, and fold weights survive `update()`.

## Why first
- Bug 4: P2 ignores `gpw`, `npw`, `tpw`. A P2 fold class would ignore its fold weights in the same way.
- Bug 3: The default CG solve on P2 can be wrong by up to 0.64 in a unit box and still return `True`.
- Bug 11: `update()` drops `fold_weights`. Fold classes are the main users.
- F1: The P2 support has no `origin`, `step_vector` or `nsteps_cells`. Registration in Phase 2 needs them.

## Items
All paths are in `packages/loop_interpolation/src/loop_interpolation/`.

| Item | Where | Fix |
|---|---|---|
| F1 | `P2UnstructuredTetMesh` constructor | Store `origin`, `step_vector`, `nsteps_cells`. Remove the workaround in `tests/golden/engines.py`. |
| Bug 4 | `_p2interpolator.py:122,151,166` | Apply weights after row normalisation, as P1 does. Turn `test_known_p2_normal_weight_has_effect` from strict xfail into a normal test. |
| Bug 5 | `_p2interpolator.py:210` | Change `np.all(neighbours > 3)` to `>= 0`. Border neighbours are negative. |
| Bug 3 | `_solver_strategy.py:66-71` | Return `False` when CG hits the iteration limit. Keep the warning. Turn `test_known_default_solve_is_converged` into a normal test. |
| Bug 11 | `_discrete_interpolator.py:1337` | `update()` must pass the stored `fold_weights`, `regularisation` and `directional_regularisation` to `setup_interpolator`. Turn `test_known_update_keeps_fold_weights` into a normal test. |
| Bug 6 | `_discrete_interpolator.py:357-359` and `_discrete_fold_interpolator.py` | Drop NaN rows only, not the whole family. Log how many rows you dropped. Apply the same filter in the shared fold setup. |

## Order and risk
1. F1 first. Other fixes need a working P2 test harness.
2. Bug 5, then bug 4, then bug 3. Each one changes P2 golden files. Review each diff alone.
3. Bug 11 and bug 6 change DFI and FDFold results. Review those golden diffs with care.
4. Bug 3 may make some existing solves report failure. Check callers of the return value of the solver.

## Done when
- The three `test_known_*` tests named above pass without xfail.
- Golden diffs are reviewed and accepted, each in its own commit.
- `python -m pytest tests/golden -q` is green.

## Not in this phase
Bugs 1, 2, 7, 8, 9, 10, 12, 13, 14 and the low items. They do not block P2Fold. Bug 7 has no effect while weights are uniform. Re-check it if Phase 1 uses per-sample weights.
