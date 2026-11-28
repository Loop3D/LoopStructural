# agent/

Working notes for `loop_interpolation`. Plans and reviews, not user documentation.

| File | What it is |
|---|---|
| `REVIEW_interpolation_hardening.md` | Code review of the shared interpolation path, P1, P2, FD, DFI and fold setup. Bugs, and the order of work for a matrix-free Rust port. Written from reading the code. The bugs it lists are now pinned by `packages/loop_interpolation/tests/golden` (`test_known_*`). |
| `HANDOFF_p2_fold.md` | Current P2Fold handoff with the corrected phase 4 gate table; the original prototype narrative is marked historical. The prototype scripts are in the author's Honours project (`debug_s0/scripts/p2_fold.py`), not in this repo. The historical volume-factor and per-sample-weight rules are explicitly superseded by the current section. |

Golden tests: `packages/loop_interpolation/tests/golden/README.md`.


Phase 4: [corrected Bermagui validation](HANDOFF_p2_fold.md),
[execution checklist](p2_fold/phase4_validation/tasks.md), and
[usage documentation](../docs/source/user_guide/fold_modelling.rst).
Both gauss4 gates passed at fewer unknowns. S0 retains a material error tail;
unregularised centroid was unstable. The comparison uses a shared DFI driver,
not an end-to-end all-P2 chain. Release checks after merge remain pending.
