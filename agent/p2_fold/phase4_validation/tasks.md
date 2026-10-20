# Phase 4 tasks

Needs Phase 3 done.

## Validation
- [x] 4.1 Get the prototype scripts from the Honours project (`debug_s0/scripts/`). Port the comparison to use the library class, not the engine swap
- [x] 4.2 Run S2 (1e4) for DFI, `P2Fold` `centroid`, `P2Fold` `gauss4`. Same seed, same explicit buffer
- [x] 4.3 Run S0 (1e5, `limb_wl=None`) for the same three engines
- [x] 4.4 Record: median error, share above 30 deg, unknowns, tets, solve time, peak memory
- [x] 4.5 Compare with the old handoff table. Explain the differences
- [x] 4.6 Decision gate: does `gauss4` beat DFI at equal or fewer unknowns?
  - [x] Yes: continue
  - [ ] No: stop and ask the user (not applicable: both gates passed)

## Docs
- [x] 4.7 Write the doc note and an example script
- [x] 4.8 Update `agent/HANDOFF_p2_fold.md`: add the new table. Mark the wrong rules as fixed
- [x] 4.9 Update `agent/README.md` and `agent/p2_fold/README.md`

## Release
- [x] 4.10 Check all commits have correct prefixes and no Co-Authored-By line
- [x] 4.11 Open a PR. Describe the Phase 0 fixes and the new class
- [ ] 4.12 After merge, check the release-please PR bumps only the right packages (see the release-please path scoping memory)
- [ ] 4.13 Start Phase 5


## Execution notes - 20 October 2026

Worktree: `.claude/worktrees/p2-phase4`, branch `p2-phase4-validation`, based on
`origin/master` at b1ca2cd7. External project:
`C:/Users/lachl/scratch/bonnie/Honours_model23.09`.

All ten candidate solves completed with finite inner probes. The
[validation report](../../../validation/p2_fold/README.md) contains the complete
corrected table, protocol, costs, input hashes, historical comparison and figures.
The code uses the library factory; no engine swap or prototype class is used.
Extra smaller P2 cases make the unknown-budget decision gate meaningful.
S0 uses a fixed DFI S2 driver for every candidate; this is documented explicitly.

The existing release PRs #348 (root), #347 (loop-common), and #344
(loop-interpolation) were checked before merge and touch only their intended
package/version files plus manifest entries. This does not substitute for task
4.12 after the phase 4 PR merges. No manual version change or Phase 5 buffer fix
is included. Task 4.13 is intentionally conditional on phase 4 documentation
being merged, per the settled phase ordering.


Checks: ten external-data candidate solves; 172 interpolation tests passed
(9 skips, 2 expected failures); 9 gate/integration tests passed; bundled example
smoke-tested; Sphinx HTML built with gallery execution disabled. The 286 docs
warnings are outside the new fold page. Existing golden files are unchanged.


Release follow-up found during validation: root dependency ranges for loop-common
and loop-interpolation are `<0.1.0`, excluding the pending 0.1.0 releases. They
must be aligned before publishing the new root release. This is separate from
package path scoping, which currently passes. No release/version PR is merged
by this validation task.

PR: [#349](https://github.com/Loop3D/LoopStructural/pull/349). Merge and tasks 4.12-4.13 remain pending.
