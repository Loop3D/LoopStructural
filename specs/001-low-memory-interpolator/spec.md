# Feature spec: Low-memory interpolator for the browser

**Feature**: `001-low-memory-interpolator` | **Created**: 2026-10-27 | **Status**: Ready for plan (reviewed 2026-10-27)

**Input**: The user wants an implicit-function interpolator that uses very
little memory and runs in a web browser. It works on a regular grid. It
never builds the system matrix. It handles inequality constraints (for
example, "horizon A is above horizon B"). Early work exists in `loop-rs` and
`LoopAI` (see `research.md`).

<!--
A spec says WHAT the user needs and WHY. It does not say HOW.
Write in ASD-STE100 Simplified Technical English.
-->

## Context

- Today the solver code that does this work is in two Rust projects outside
  this repository. It runs as a native program or as a Python extension. It
  does not run in a browser. Nobody has measured its memory use.
- The Python package `loop_interpolation` has a matrix-based version, and a
  second version for inequality constraints. Both need much more memory.
- The goal is **not** a faster solver. The goal is a solver that **fits in a
  small, known amount of memory**, so that it can run in a browser tab.
- The existing work is described in `research.md`. It is the starting point.
  This spec does not copy its design.
- Words used in this spec:
  - **dof**: one unknown value of the field, one for each grid node.
  - **budget**: the maximum memory that the user lets the solver use.
  - **inequality constraint**: a rule such as `field(A) >= field(B)` or
    `field(P) <= 10`.
  - **stage**: one named part of a solve, for example "setup", "solve" or
    "check constraints".
  - **desktop**: the native program, or the Python call, that uses the same
    solver code outside a browser.
  - **constraint row**: one data constraint or one inequality constraint. A
    regularisation term is not a constraint row.
  - **active constraint**: an inequality constraint that holds with equality
    (within SC-004) at the answer.
  - **reference solution**: the field that the matrix-based solver of
    `loop_interpolation` gives for the same problem. For a problem that is too
    big for that solver, it is the field of the same problem on the largest
    grid that the matrix-based solver can solve, sampled at the same points.

## User scenarios and tests

### User story 1: Solve a model in a browser tab (priority: P1)

A geologist opens a web page. They load contact points and orientation data.
They choose a grid. They press "Solve". The page shows a scalar field that
honours the data. No server does the work, and no install is needed.

**Why this priority**: this is the whole purpose. Without it, nothing else in
this spec has value.

**Independent test**: run problem R1 (see Success criteria) in a headless
browser. Compare the field with the reference solution.

**Acceptance scenarios**:

1. **Given** a grid of 1 million dof and a set of value and orientation
   constraints, **when** the user solves in the browser, **then** the field is
   returned and agrees with the reference result (SC-003).
2. **Given** the same problem, **when** the solve runs, **then** peak memory
   stays inside the limit of SC-001.
3. **Given** a page with no cross-origin isolation, **when** the user solves,
   **then** the solve still works on one thread. This holds in current
   Chromium, Firefox and Safari. Automatic tests run in Chromium. Firefox and
   Safari run a smoke test on a small grid.
4. **Given** the same problem, **when** it runs on the desktop, **then** it
   gives the same field as in the browser (SC-003).
5. **Given** a problem with tangent constraints and interface constraints,
   **when** the user solves, **then** the tangent vectors are normal to the
   field gradient, and the two points of each interface have the same field
   value, within the tolerance of SC-003.
6. **Given** a problem with each regularisation type (curvature, border,
   directional), **when** the user solves with one type at a time, **then** the
   field agrees with the reference solution within the tolerance of SC-003.
7. **Given** the same problem solved with 1 thread and with 4 threads,
   **when** both solves end, **then** the two fields differ by at most the
   tolerance of SC-003.
8. **Given** a solved field, **when** the user asks for the value and the
   gradient at a point inside the grid, **then** the solver returns them. The
   values agree with the reference solution within the tolerance of SC-003.
   A point outside the grid gives an error.
9. **Given** a grid whose origin has a large coordinate value (for example
   `5 000 000`), **when** the user solves, **then** the field agrees with the
   field of the same problem moved to an origin of zero, within the tolerance
   of SC-003.
10. **Given** input with a NaN, an infinite value, an array of the wrong
    length, or a zero-length direction vector, **when** the user solves,
    **then** the solver returns an error that names the constraint. It does
    not return a field.

### User story 2: Keep horizons in the right order (priority: P1)

The geologist has two or more horizons. The field must keep them in order
everywhere, not only at the data points. The geologist adds inequality
constraints. The solver respects them.

**Why this priority**: inequality constraints are the reason this solver
exists, and they are the part that normally costs the most memory.

**Independent test**: use a problem where the unconstrained answer breaks the
order. Solve with the constraints. Check that the order holds at all constraint
points. Check peak memory.

**Acceptance scenarios**:

1. **Given** a value-bound constraint that the unconstrained answer breaks,
   **when** the user solves, **then** the bound holds within the tolerance of
   SC-004.
2. **Given** a pair constraint `field(A) - field(B) >= d` that the
   unconstrained answer breaks, **when** the user solves, **then** it holds
   within the same tolerance.
3. **Given** problem R3 (many active constraints), **when** the user solves,
   **then** peak memory stays inside the limit of SC-011. Memory does not grow
   with the **square** of the number of active constraints.
4. **Given** a constraint that is not active at the answer, **when** the user
   solves, **then** the result equals the result with no such constraint.
5. **Given** a constraint that has a lower bound and an upper bound,
   **when** the user solves, **then** both bounds hold within the tolerance of
   SC-004.
6. **Given** constraints that cannot all hold, **when** the user solves,
   **then** the solver stops and reports that the problem is infeasible. It
   gives the 20 constraints with the largest violation (FR-015). It does not
   return a field without a warning.

### User story 3: The page stays alive during a long solve (priority: P2)

A solve can take up to the time of SC-007. During it, the page must respond. The user
sees progress and can cancel.

**Why this priority**: a frozen tab looks like a crash. The browser may show a
"page not responding" message.

**Independent test**: start a long solve. Measure the delay of the page's
event loop. Cancel the solve and check the time to stop.

**Acceptance scenarios**:

1. **Given** a solve is running, **when** the user uses the page, **then**
   the page responds in less than 100 ms (SC-006).
2. **Given** a solve is running, **when** the solver reports progress, **then**
   the report gives the stage, the iteration, and the residual. The solver
   sends a report at least once each second (FR-009).
3. **Given** a solve of R1 is running, **when** the user cancels, **then** it stops
   in less than 1 s. The solver frees its memory. The page can start a new
   solve.
4. **Given** the user cancels, **when** the user asks for the result, **then**
   the solver can return the best field so far, marked as "not converged".

### User story 4: Refuse a problem that does not fit (priority: P2)

The user chooses a grid that is too big for the browser. The solver says so
before it starts. It does not crash the tab.

**Why this priority**: a crash in a browser loses the user's page state, with
no clear reason.

**Independent test**: ask for an estimate for grids of increasing size. Check
that the estimate is not lower than the real peak. Check that a too-large
request fails with a message.

**Acceptance scenarios**:

1. **Given** a grid and a set of constraints, **when** the user asks for a
   memory estimate, **then** the solver returns the expected peak memory
   without allocating the large arrays.
2. **Given** a budget, **when** the estimate is larger than the budget,
   **then** the solve does not start. The message gives the estimate, the
   budget, and the largest grid that would fit.
3. **Given** the browser refuses to grow memory during a solve, **when** this
   happens, **then** the solver returns an error. It does not trap or abort
   the module.

### User story 5: Choose precision against memory (priority: P3)

The user can trade accuracy for memory. A single-precision mode uses at most
60% of the memory of the default mode (scenario 1).

**Why this priority**: it raises the largest grid that fits by about 1.7
times. It is only useful
after stories 1 and 2 work.

**Independent test**: solve the reference problem in both modes. Compare the
fields and the peak memory.

**Acceptance scenarios**:

1. **Given** single-precision mode, **when** the user solves the reference
   problem, **then** peak memory is at most 60% of the double-precision peak.
2. **Given** single-precision mode, **when** the solve ends, **then** the
   field agrees with the double-precision field within the tolerance of
   SC-005.
3. **Given** a solve in either precision mode, **when** the solve ends,
   **then** the reported residual differs from a double-precision calculation
   of the same residual by at most 1e-6 relative error (FR-022).
4. **Given** single-precision mode cannot reach the requested tolerance,
   **when** the solve ends, **then** the result says so. It does not claim
   convergence.

### User story 6: Use the same solver from Python (priority: P3)

A Python user of `loop_interpolation` calls the same solver. This keeps one
solver, and it lets the Python test suite check the browser solver.

**Why this priority**: one solver tested in two places is cheaper to keep
correct than two solvers.

**Independent test**: run the same problem through the Python call and through
the browser build. Compare the fields.

**Acceptance scenarios**:

1. **Given** a problem built with the Python constraint classes, **when** the
   user selects the low-memory solver, **then** the field agrees with the
   matrix-based solver within the tolerance of SC-003.
2. **Given** the user does not install the new solver, **when** the user uses
   `loop_interpolation`, **then** all existing behaviour stays the same
   (FR-019). The existing test suite passes with no change.
3. **Given** a setting that the new solver does not support (for example a
   model region mask), **when** the user selects the new solver, **then** the
   call fails with an error that names the setting. It does not ignore the
   setting (FR-026).
4. **Given** the same problem, **when** it runs through the Python call and
   through the browser build, **then** the two fields agree within the
   tolerance of SC-003 (FR-017).

### Edge cases

- The grid has one node along an axis. (Error: each axis needs at least 3
  nodes.)
- The step of the grid differs by more than a factor of 100 between two axes.
  (The solver gives a warning.)
- A constraint point lies outside the grid. (It is dropped with a warning
  that counts the dropped points.)
- The user closes the page or ends the worker during a solve. (All memory is
  freed. No state is saved, see Out of scope.)
- An input value is NaN or infinite. (The solver refuses the input and names
  the constraint.)
- All constraints are inequality constraints, with no value or orientation
  data. The problem may have many answers. (The regularisation picks one.)
- A rotated grid. (Not supported. The caller rotates the data, FR-020. A
  request with a rotation gives a clear error.)
- A solve with zero constraints. (Error: no constraints.)
- An inequality constraint whose lower bound is above its upper bound. (Error
  that names the constraint.)
- A constraint weight of zero or less. (Error that names the constraint.)
- Two value constraints at one point with different values. (The solver
  accepts them and gives a least-squares compromise. It adds a warning that
  counts the conflicts.)
- A grid whose node count passes the largest count that the solver can index.
  (Error that gives the node count and the largest count. The memory estimate
  does not run first.)
- The user cancels during setup, before the first iteration. (The solve stops
  in less than 1 s, as in FR-010.)
- A second solve starts in the same page while one solve runs. (The second
  call is refused with an error that says a solve is running. The first solve
  is not changed.)
- The solver reaches the iteration limit. (The result says "not converged"
  with the final residual.)
- 2D grids. (Not in the first release. A 2D request gives a clear error.)
- A model region (a mask on the nodes), a column scaling, or a constant-norm
  iteration from the Python package. (Not in the first release. Each gives a
  clear error, FR-026.)

## Requirements

### Functional requirements

- **FR-001**: The solver must never store the system matrix or the normal
  matrix.
- **FR-023**: The solver must never store any other matrix whose size grows
  with the number of dof times the number of neighbours of a node.
- **FR-002**: The memory that the solver needs must grow in a straight line
  with the number of dof and the number of constraints. It must not grow with
  the square of either.
- **FR-003**: The solver must support these constraints: value, gradient,
  normal, tangent, and interface (equal value of two points).
- **FR-004**: The solver must support these inequality constraints: a bound on
  the field at a point, and a bound on the difference of the field at two
  points. A bound can be a lower bound, an upper bound, or both.
- **FR-005**: The solver must support the existing regularisation: curvature
  in all six directions, border regularisation, and directional
  regularisation.
- **FR-006**: The solver must run in a web browser page that has no
  cross-origin isolation. In that page it runs on one thread.
- **FR-007**: The solver may use more threads when the platform has them. The
  result must not depend on the thread count, within the tolerance of SC-003.
- **FR-008**: A solve must not block the page (story 3).
- **FR-009**: The solver must report progress (stage, iteration, residual) at
  least once each second.
- **FR-010**: The solver must stop on request in less than 1 s, free its
  memory, and be ready for a new solve.
- **FR-011**: The solver must give a memory estimate for a problem before it
  allocates the large arrays (story 4). The estimate must not be lower than
  the real peak.
- **FR-012**: The solver must accept a budget. It must refuse to start when the
  estimate is over the budget. It must return an error, not abort, when an
  allocation fails.
- **FR-013**: The solver must offer a single-precision mode and a
  double-precision mode (story 5).
- **FR-022**: In both modes, the residual that decides convergence must agree
  with a double-precision calculation of the same residual within 1e-6
  relative error.
- **FR-014**: Every result must say if the solve converged, how many
  iterations it used, the final residual, and, for inequality problems, the
  largest constraint violation.
- **FR-015**: For an infeasible inequality problem, the solver must say so and
  return the 20 constraints with the largest violation at the stopped state
  (story 2, scenario 5). The list is not a proof of a minimal conflict set.
- **FR-016**: The solver must return the field values at the grid nodes, and
  must evaluate the field and its gradient at any point inside the grid.
- **FR-017**: The solver must give the same results from Python (story 6)
  and from the browser, within the tolerance of SC-003.
- **FR-018**: Bad input (NaN, infinite, wrong array length, zero-length
  vector for a direction) must give an error that names the constraint. It
  must not give a wrong field.
- **FR-019**: The existing behaviour of `loop_interpolation` must not change
  when the user does not select the new solver.
- **FR-020**: The solver must take all points and vectors in the coordinate
  system of the grid. It must not change coordinates. The caller does all
  coordinate transformation, including any rotation of the model.
- **FR-024**: The solver must keep the full precision of the grid origin that
  the caller gives.
- **FR-021**: The user must be able to select the solver in place of the
  existing finite difference interpolator of `loop_interpolation`, with a
  solver setting (story 6).
- **FR-025**: The same constraints, weights and regularisation settings must
  work with the new solver.
- **FR-026**: A setting that the solver does not support must give a clear
  error that names the setting. The solver must not ignore it.

### Protocol and data (if the feature changes them)

- **Input**: a grid (origin, node counts, step; a rotation must be absent), a list of
  constraints, regularisation weights, solver settings (precision, tolerance,
  iteration limit, budget).
- **Output**: the field at the nodes, a convergence report, and a list of
  warnings.
- **Data format**: constraints go in as flat numeric arrays, one array for each
  field. They do not go in as text. The copies of the large arrays count in the
  peak memory of SC-001. (How this is done goes in the plan.)

### Key entities

- **Problem**: one grid, many constraints, regularisation, and settings.
- **Constraint**: a small record (a point, a target, a weight). It is not a
  matrix row.
- **Budget**: the maximum memory in bytes.
- **Memory estimate**: the expected peak bytes, split into the vectors, the
  constraint records, and the working storage.
- **Reference solution**: see the definition in Context.
- **Solve report**: convergence flag, iteration count, residual, largest
  violation, time, peak memory, and warnings.

## Success criteria

<!-- Measurable. No implementation details. -->

Reference problem **R1**: a grid of 100 x 100 x 100 nodes (1 million dof), 2 000
value constraints, 500 orientation constraints, and the default
regularisation. Reference problem **R2**: the same grid, with 200 pair
inequality constraints that are active at the answer.
Reference problem **R3**: the same grid, with 20 000 pair
inequality constraints that are active at the answer.

- **SC-001**: In double precision, R1 has a peak memory of 100 MB or less in
  the browser, including the field that the solver returns. This is about 100
  bytes for each dof.
- **SC-002**: R2 has a peak memory of at most 20% more than R1.
- **SC-003**: The browser result, the desktop result, and the reference
  result agree. The maximum difference is at most 1e-6 of the field range in
  double precision. The solver stops at a residual tolerance that gives this
  agreement. Results for different thread counts may differ in the last bits.
- **SC-004**: In R2, every inequality constraint holds with a violation of at
  most 1e-4 of the field range. (A tighter limit of 1e-8 is for a later
  phase, see Clarifications.)
- **SC-005**: In single precision, the maximum difference from double precision
  is at most 1e-3 of the field range on R1.
- **SC-006**: During a solve of R1 in the browser, the delay of the page event
  loop stays below 100 ms.
- **SC-007**: R1 solves in 60 s or less in Chromium on one thread, on the
  reference machine. The reference machine is a laptop with 4 cores and 16 GB
  of RAM. CI measures and reports the time. The release check uses the
  reference machine.
- **SC-008**: On the desktop, a grid of 250 x 250 x 250 nodes (about 15.6 million dof) in
  single precision, with a budget of 2 GB, solves with a peak memory of 1.5 GB
  or less. (A request over the budget is refused, see story 4.)
- **SC-009**: The memory estimate is never lower than the measured peak, and
  is at most 25% higher, for 10 grid sizes from 10^3 to 200^3 nodes, on the desktop. In Chromium the
  same check runs for the grid sizes that fit in browser memory.
- **SC-010**: The module for the browser has a download size of 2 MB or less
  after gzip compression.
- **SC-011**: R3 has a peak memory of at most 20% more than R1.
- **SC-012**: Peak memory for R1-type problems on 50^3, 100^3 and 150^3 nodes
  fits a straight line against the dof count with R^2 of at least 0.99
  (FR-002).

## Out of scope

- Meshes that are not regular grids (tetrahedral or other). The existing
  solvers for these stay as they are.
- Fault, fold and unconformity logic. This feature solves one scalar field. A
  model builder calls it many times.
- Radial basis function and other non-grid interpolators.
- Drawing the field or an isosurface. (The viewer does this. See `geoviewer`.)
- A server or cloud mode.
- GPU solvers (WebGPU).
- Automatic choice of the best settings.
- Saving a solve to continue it later.
- Changes to the Python matrix-based solvers.
- Coordinate transformation and grid rotation. The caller does these (FR-020).
- Fold interpolators (they stay on the matrix-based solvers).
- A constraint tolerance of 1e-8. A later spec can add it.
- 2D grids, model region masks, column scaling, constant-norm iteration and
  rotated grids. Each gives a clear error in this release.
- Browsers other than the current Chromium, Firefox and Safari.
- Phones and tablets.

## Clarifications

<!-- Record each question and its answer here, with the date. -->

- 2026-10-27: The user asked for a spec for this work. The user said that the
  work can go in `loop_interpolation` or in a new module of this library. The
  user said that existing work is in `~/dev/LoopAI` or `~/dev/loop-rs`.
  Result of the search: `loop-rs` is more complete (see `research.md`).

- 2026-10-27: The user answered: "Plan it as a drop in for the loop
  interpolation. If a web page calls it, the web page can deal with coordinate
  transformation." Result: the solver is a drop-in choice for
  `loop_interpolation` (FR-021). The caller does coordinate transformation
  (FR-020). The input is flat arrays.

- 2026-10-27: Review interview. Q: What peak memory for R1? A: 100 MB, keep
  SC-001. Q: Solve time and machine? A: 60 s on a laptop with 4 cores and
  16 GB RAM (SC-007). Q: Which browsers? A: Current Chromium, Firefox and
  Safari. Tests run in Chromium. Firefox and Safari get a smoke test (story 1,
  scenario 3). Q: Two solves at once? A: The second is refused with an error.
  Q: List of conflicting constraints? A: The 20 with the largest violation
  (FR-015). Q: Exact 1e-8 tolerance? A: A later spec (Out of scope).

- 2026-10-27: Second review interview. Q: Must the field be bit-for-bit equal
  for any thread count? A: No. Agreement within the SC-003 tolerance is
  enough (FR-007). Q: What is the SC-003 tolerance? A: 1e-6 of the field
  range. Keep it. The solver stops at a residual tolerance that gives it. Q:
  How many active constraints in R3? A: 20 000 (SC-011).
