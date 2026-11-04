# Specification index

Updated: 2026-11-04. This index records document and workflow state. It does
not certify that an implementation or a numerical acceptance gate passed.

| Spec | State | Next step or existing open point |
|---|---|---|
| [001-low-memory-interpolator](001-low-memory-interpolator/spec.md) | Spec declares ready for plan, reviewed 2026-10-27 | Follow its existing plan and task workflow. This cascade did not reassess it. |
| [002-rbf-interpolator](002-rbf-interpolator/spec.md) | Draft, not reviewed | Continue its existing review and clarification workflow. This cascade did not change it. |
| [003-fault-evaluation-performance](003-fault-evaluation-performance/spec.md) | Implementation in progress; 005 contract cascade applied | Its active plan and tasks remain in the existing LoopStructural-native worktree. OP3 and OP4 remain open; completed tasks, including concurrent T043, and measured failures are preserved. |
| [004-matrix-free-directional-constraints](004-matrix-free-directional-constraints/spec.md) | Draft; 005 contract cascade applied; task graph checked | Review the updated plan and tasks. OP1 still requires the quadrature measurements before a default is chosen. |
| [005-coordinate-frame-correctness](005-coordinate-frame-correctness/spec.md) | S1 implemented 2026-11-04 on `feat/005-coordinate-frame-correctness` (17 of 17 tasks; evidence in [verification.md](005-coordinate-frame-correctness/verification.md)); S2 to S6 not planned | Review and merge the S1 branch. OP1: plan S2 to S6 in later planner runs. OP2 closed 2026-11-04: the R2/R5 matrix stays in the default pytest run. 7 minor review findings are open (run-state.json). |

## Cascade of 005

The accepted coordinate contract keeps modelling observations, reference
geometry, masks, and public queries in world coordinates. Local assembly
and low-level core inputs use the computational frame. Each boundary
converts points and directions once.

Analytic evidence and explicit review may justify corrections to affected
golden results and fingerprints. Unaffected references, numeric bounds,
acceptance tolerances, and performance targets remain unchanged.

The cascade changes documents only. It does not update golden files,
benchmark data, solver code, or the run-state of an implementation.

The matching 003 spec copies in the main checkout and its active worktree
are kept in sync. Its plan and tasks are not copied into the main checkout.
The documentation changes create no whole-spec implementation dependency
between 005 and either target.
