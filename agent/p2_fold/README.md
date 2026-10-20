# P2 fold interpolation: work plan

Sources: `agent/HANDOFF_p2_fold.md`, `agent/REVIEW_interpolation_hardening.md`, `packages/loop_interpolation/tests/golden/README.md`.

## Decisions (settled)

1. Build `P2FoldInterpolator` on the assembled-matrix path. Do not wait for the matrix-free operator.
2. Add a new type string `"P2Fold"`. Plain `"P2"` does not change.
3. The double `buffer` fix is not part of this work. It is a separate `fix:` commit after Phase 4 (see `phase5_buffer_followup`).

## Phases

| Phase | Folder | Goal |
|---|---|---|
| 0 | `phase0_prereqs` | Fix the P2 and shared bugs that the fold class would inherit |
| 1 | `phase1_class` | Write `P2FoldInterpolator` and the quadrature rule |
| 2 | `phase2_wiring` | Register `"P2Fold"` and stop forcing `"DFI"` |
| 3 | `phase3_tests` | Unit, golden and accuracy tests |
| 4 | `phase4_validation` | Re-run Bermagui, write docs, release |
| 5 | `phase5_buffer_followup` | Single-point buffer fix (after Phase 4) |

Each folder has `plan.md` (why and how) and `tasks.md` (checklist).

## Rules for all phases

- One fix per commit. Write the failing test first.
- Use `fix:` and `feat:` prefixes. Release-please sets versions.
- Do not add a Co-Authored-By line to commits (standing user preference).
- Run `python -m pytest tests/golden -q` from `packages/loop_interpolation` before each commit. Review every `--update-golden` diff by hand. The golden files are the Rust parity oracle.
- Single runs vary by about 3 deg. Compare numbers only inside one run, with fixed seeds.

## Open issues outside this work

- F2: `TetMesh` returns NaN on the lower x, y, z faces.
- F3: Fold amplitude is damped (gradient ratio 0.35 to 0.66). Decide the intended energy in the Rust operator specification.
- The rotation-angle profile fit is not repeatable. `limb_wl` below 5 stays in a poor local minimum.
- Matrix-free P1/P2 operator and Rust port (review section 4, steps 3 to 5).


## Status after phase 4 validation (20 October 2026)

Phases 0-3 are on master. Phase 4 measurements and the decision gate are complete:
S2 gauss4 4.69 deg / 9,711 unknowns versus DFI 6.15 deg / 12,992; S0 gauss4
11.14 deg / 85,673 versus DFI 20.63 deg / 112,560. See the
[full table and limitations](../../validation/p2_fold/README.md).
The old prototype weighting and "nelements is tets" statements are superseded.
The API note and bundled-data example are in [PR #349](https://github.com/Loop3D/LoopStructural/pull/349). Post-merge
release scoping and Phase 5 remain pending; this phase does not change buffer
behaviour or manually bump versions.
