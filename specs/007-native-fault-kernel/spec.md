# Feature spec: Native fault application

**Feature**: `007-native-fault-kernel` | **Created**: 2026-11-27 | **Status**: Draft

**Input**: Restoring points through the faults of a feature takes most of the
evaluation time of a faulted model. A Rust proof of concept restored points
6.5 to 16 times faster on one thread and 60 to 150 times faster on 20 threads.
The user wants a native fault application that LoopStructural uses when it is
available, with the Python code as the reference and the fallback.

<!--
A spec says WHAT the user needs and WHY. It does not say HOW.
Write in ASD-STE100 Simplified Technical English.
-->

## Context

- To evaluate a faulted feature, LoopStructural moves each point back through
  the faults of the feature ("restores" the point), and then evaluates the
  interpolated field at the restored point. Each fault moves a point in steps
  (10 by default) along the gradient of its fault frame.
- After spec 003 and the fixes of master #424, #425, #427, #428 and #429, each
  fault is applied once for each batch of points in one evaluation call, also
  in networks of cutting faults. The remaining cost is the fault application
  itself, in Python:
  - About 7 µs for each point and fault in the golden models.
  - About 40 ms of fixed cost for each fault frame that is cut by younger
    faults, for each call. Most of it is Python work that finds the restore
    order of the chain of cutting faults again at each call.
  - A model with 94 faults (54 faults that do not touch other faults, 40
    cutting faults, chains of cuts up to 20 faults deep) evaluates 5000 points
    in 2.1 s. The 40 cutting faults use 1.7 s of this time.
- A proof of concept (outside the repository, 2026-11-26) wrapped an existing
  Rust port of the fault application for Python, and read the model from
  exported files. Against the
  Python restoration of spec 003, at 10 000 points:
  - 6.5 to 11 times faster on one thread, 64 to 84 times faster on 20 threads,
    for the golden models with 1 to 4 faults.
  - 16 times (one thread) and 149 times (20 threads) faster for a model with
    12 faults and a chain of cuts; 58 and 505 times faster at 1000 points,
    because the Python path has a fixed cost of about 50 ms for each call
    on this model, which dominates a batch of 1000 points.
  - The restored points differ from the Python path by up to 5.4e-5 m in a
    model with a 1 km box (`f4_u0`) and 1.9e-4 m in a model with a 10 km box
    (displacements 40 to 120 m). This does not pass the golden tolerance
    (1e-8 × the box size: 1e-5 m and 1e-4 m). The crate was written for an
    older version of the fault code; the plan must find the cause (C-007).
    The NaN points are the same.
- LoopStructural already ships optional native code for its solvers (an
  optional extra of the interpolation package). The plan chooses where the
  native path lives.
- Dependencies: spec 003 (the golden tests, the benchmarks and the evaluation
  scope) must be merged before the plan of this spec starts. The fixes of
  master #425, #427, #428 and #429 are merged.
- Words used in this spec:
  - **restoration**: moving points back through the faults of a feature, in
    the restore order (each fault that cuts another fault first, and the whole
    chain of cutting faults), as the Python path does today.
  - **native path**: the compiled implementation of this spec.
  - **Python path**: the present implementation. It stays the reference.
  - **supported feature**: a feature whose faults and field the native path
    can evaluate (FR-008). Any other feature uses the Python path.

## User scenarios and tests

### User story 1: The Python path finds the restore order once (priority: P1)

A developer evaluates a model with cutting faults. The Python path computes
the restore order of the chain of cutting faults of a feature once, and not
at each call, so the fixed cost of a cut fault frame falls from about 40 ms
to at most 10 ms for each call.

**Why this priority**: it helps every user, also users without the native
package, and it makes the Python reference fair to compare with. It has no
native code, so it is the first stage.

**Independent test**: evaluate the 94-fault model at 500 points. The time of
a cutting fault falls, and the results are the same bits as before.

**Acceptance scenarios**:

1. **Given** the 94-fault model, **when** `evaluate_model` runs at 500
   points, **then** the mean time of a cutting fault is at most 10 ms
   (today about 40 ms), and the model ids and values are bit-identical to
   the results before this story.
2. **Given** a fault whose cutting faults change after the first evaluation,
   **when** the feature is evaluated again, **then** the restore order uses
   the new cutting faults.

### User story 2: Restore points natively (priority: P1)

A user evaluates a faulted model with the native package installed. The
restoration of the points of a supported feature runs in native code.

**Why this priority**: the restoration is about 90 % of the evaluation time
of a faulted model.

**Independent test**: restore the points of each golden model with the
native path and with the Python path, and compare.

**Acceptance scenarios**:

1. **Given** a golden model and its reference points, **when** the points
   are restored with the native path, **then** the restored points agree
   with the Python path within the golden tolerances (spec 003 FR-004), and
   the NaN points are the same.
2. **Given** a feature with a block (a side of a fault forced), reverse
   movement, or points that a younger fault already restored
   (`restored_through`), **when** the native path restores the points,
   **then** the result agrees with the Python path as in scenario 1.

### User story 3: Evaluate values and gradients natively (priority: P2)

A user evaluates feature values and gradients (`evaluate_value`,
`evaluate_gradient`, `evaluate_model`, `evaluate_model_gradient`) with the
native package installed. The field values and the gradients of a supported
feature are computed natively.

**Why this priority**: after story 2, the field evaluation and the gradient
(which restores the corners of a tetrahedron for each point and refines it
near faults) are the next largest costs.

**Independent test**: run the golden tests with the native path on and with
it off.

**Acceptance scenarios**:

1. **Given** each golden model, **when** the golden tests run with the
   native path on, **then** they pass with the golden tolerances, the model
   ids are the same except at a unit boundary, and the NaN points are the
   same.
2. **Given** a point near a fault where the gradient tetrahedron still
   crosses the fault after the refinement, **when** the gradient is
   evaluated natively, **then** it is NaN, as on the Python path.

### User story 4: On by default, with a way out (priority: P2)

A user installs LoopStructural with the native package (for example through
the QGIS plugin). The native path is used for supported features without any
change to the user's code. A user can turn it off.

**Why this priority**: without it, few users get the speed-up; with a switch
and a fallback, a problem in the native path never blocks a user.

**Independent test**: run the full test suite with the native path on, off,
and with the native package not installed.

**Acceptance scenarios**:

1. **Given** the native package is installed, **when** a model is evaluated,
   **then** the native path is used for each supported feature.
2. **Given** the configuration setting or the environment variable turns the
   native path off, or the native package is not installed, **when** a model
   is evaluated, **then** the Python path is used, with no error and the same
   results as today.
3. **Given** the setting turns the native path on and the variable turns it
   off, **when** a model is evaluated, **then** the native path is used (the
   setting wins).
4. **Given** the native package is installed but its version does not match,
   **when** two models are evaluated in one session, **then** the Python path
   is used and exactly one warning names both versions.
5. **Given** a clean computer with Windows, macOS or Linux, **when** a user
   installs LoopStructural with the native extra from the package index,
   **then** a prebuilt native package installs and the native path is used.

### User story 5: The native path is safe to use (priority: P2)

A developer relies on the native path in an application (for example a plugin
that evaluates in background threads). The native path keeps the public API,
falls back where it must, never uses stale data, lets other threads run, and
uses the number of threads that the developer sets.

**Why this priority**: the speed-up must not cost correctness or control.

**Independent test**: the scenarios below, with the native path on.

**Acceptance scenarios**:

1. **Given** the public classes and functions of `GeologicalModel` and of the
   features, **when** their signatures are compared before and after this
   spec, **then** they are the same.
2. **Given** a feature with a displacement profile given as a Python
   function, **when** it is evaluated, **then** the Python path is used and
   the result is the same as with the native path off.
3. **Given** a model evaluated once, **when** it is evaluated again with no
   change, **then** no file is read or written, and the second call does not
   copy the model data again (its extra memory is at most the memory of the
   points and their results).
4. **Given** a model evaluated once with the native path, **when** a fault
   displacement changes and the model is evaluated again, **then** the result
   equals the Python path result for the new state.
5. **Given** a native call on 10 million points in one thread, **when** a
   second Python thread runs at the same time, **then** the second thread
   makes progress during the native call.
6. **Given** two threads that evaluate one model at the same time, **when**
   both finish, **then** each result equals the result of that thread alone.
7. **Given** the thread count setting is 1, **when** a model is evaluated,
   **then** the native path uses one thread.

### Edge cases

- A feature has one fault or field that the native path does not support:
  the whole feature uses the Python path, so that the restore order is the
  same for all its faults.
- The model changes after a native evaluation (a fault is rebuilt, a
  displacement changes, data change): the next evaluation uses the new model
  state, never a stale copy.
- Points with a NaN coordinate, points outside the grid of a fault frame or of
  the field, and an empty point array: the same results as the Python path.
- A point on a fault that cuts a fault frame (the frame has no single
  gradient there): the step direction follows master #425.
- Two threads evaluate the same model at the same time (the QGIS plugin
  evaluates in background threads): each gets the result it gets alone.
- A very large batch (10 million points): the extra memory stays within
  spec 003 FR-010; the native path does not copy the model data for each
  call.
- A long native call: it does not hold the Python global interpreter lock,
  so other Python threads can run.
- The user interrupts a long native call (Ctrl-C): the interrupt reaches
  Python at the latest when the native call ends, and the model state does
  not change.
- The native call fails in the middle of a batch (for example an internal
  error): the error reaches Python as an exception with a message, no partial
  result is returned, and the model state does not change.

## Requirements

### Functional requirements

- **FR-001** (story 1): The Python path MUST compute the restore order of the
  faults of a feature (each cutting fault first, and the whole chain of
  cutting faults) at most once for each state of the faults, not at each
  evaluation call.
- **FR-002** (story 1): When a fault, its cutting faults or its build state
  change, the next restoration MUST use the new restore order.
- **FR-003** (story 2): The native path MUST restore points through the
  faults of a supported feature with the semantics of the Python path: the
  restore order, the faults that the points were already restored through,
  the block sides and the block extension, reverse movement, the number of
  steps, the NaN rules (master #413) and the step direction on a cutting
  fault (master #425).
- **FR-004** (story 3): The native path MUST evaluate the value and the
  gradient of a supported feature at restored points, with the gradient
  method of the Python path (a regular tetrahedron for each point, refined
  only where it crosses a fault, at most 10 times; NaN if it still crosses,
  master #421).
- **FR-005** (stories 2 and 3): Native values, gradients and restored points
  MUST agree with the Python path within the golden tolerance of spec 003
  (1e-8 × the size of the model box, C-007).
- **FR-005a**: Native results MUST have the same NaN points as the Python
  path.
- **FR-005b**: Native model ids MUST equal the Python model ids, except at
  points within the value tolerance of a unit boundary.
- **FR-005c**: The test report MUST record the largest native difference.
- **FR-006** (story 4): The native path MUST be used automatically for
  supported features when the native package is installed and its version
  matches.
- **FR-006b** (story 4): A LoopStructural configuration setting and an
  environment variable MUST each be able to turn the native path off; the
  setting wins over the variable (C-008).
- **FR-006c** (story 4): With the native path off or without the package, the
  Python path MUST give the same results as before this spec.
- **FR-006a** (story 4): When the native package is installed but its version
  does not match, LoopStructural MUST use the Python path and MUST log one
  warning in a session that names both versions (C-009).
- **FR-007** (story 5): The public `GeologicalModel` and feature API MUST NOT
  change (constitution). The native path is an internal implementation.
- **FR-008** (story 5): A feature is supported when its fault frames and its field use
  structured grids with a uniform cell lookup (an affine local frame is
  allowed), and each of its fault displacement profiles is one of the
  built-in profile types of LoopStructural (C-010). Its regions (for example
  unconformities) stay on the Python path, on features that the native path
  evaluates. A feature with a displacement profile given as a Python
  function, or with any other unsupported part, MUST use the Python path.
- **FR-009** (story 5): The native path MUST read the model data in memory (grids,
  node values, frames, profiles, constants). It MUST NOT write or read
  files, and it MUST NOT copy the model data for each evaluation call.
- **FR-010** (story 5): The native path MUST NOT keep data from a model state after
  that state changes (FR-002 applies to native data too).
- **FR-011** (story 5): The native path MUST release the Python global
  interpreter lock during native work.
- **FR-011a** (story 5): The native path MUST give the same results when two
  threads evaluate one model at the same time.
- **FR-012** (story 5): The number of threads of the native path MUST be configurable.
  The default is all available cores.
- **FR-013**: The test suite MUST run with the native path on and with it
  off, and the golden tests and the spec 003 benchmarks MUST be able to
  select either path.
- **FR-014** (story 4): The release process MUST build and publish prebuilt
  native packages for Windows, macOS and Linux (64-bit x86 on each, and 64-bit
  ARM on macOS), for the Python versions that LoopStructural supports
  (C-011).
- **FR-015**: `ROADMAP.md` MUST record the decision to add the native path,
  and each stage of this spec when it completes (constitution).

## Success criteria

- **SC-001** (story 1): For the 94-fault model at 500 points, the mean time
  of a cutting fault in one evaluation is at most 10 ms (baseline about
  40 ms), with bit-identical results.
- **SC-002** (story 2): For each of the golden models `f1_u0`, `f4_u0`,
  `f4_u3`, `cutting_faults` and `plugin` at 10 000 points, the native
  restoration on one thread takes at most 20 % of the time of the Python
  path (the proof of concept: 9 % to 15 %). Measured on the spec 003 T030
  computer (Intel i7-12700H, 20 threads, Linux).
- **SC-003** (story 2): For the same models and computer, with all 20
  threads, the native restoration at 10 000 points takes at most 3 % of the
  time of the Python path (the proof of concept: 1.2 % to 1.6 %).
- **SC-004** (stories 2 and 3): For the 94-fault model, `evaluate_model` at
  5000 points takes at most 0.5 s with the native path (baseline 2.1 s with
  the Python path of spec 003).
- **SC-005** (story 3): All golden tests pass with the native path on
  (FR-005).
- **SC-006** (story 4): The full test suite passes with the native path on,
  with it off, and without the native package.
- **SC-008** (story 4): A prebuilt native package installs on a clean
  computer for each platform of FR-014, and the golden tests pass there with
  the native path on.
- **SC-007**: For 10 million points, the extra memory of the native path is
  within spec 003 FR-010.

## Out of scope

- Precomputed displacement maps on grids (the literature review of
  2026-11-26, `research.md` in this folder). After the
  ordering fixes the cost is linear in the depth of a chain of cuts, so they
  are not needed now.
- A change of the fault displacement method (for example fewer steps, Heun
  steps, or a single step with a correction): the native path implements the
  present method. A new method needs its own spec.
- Native interpolation (the solvers), native region masks, and supports other
  than structured grids.
- The QGIS plugin and how it installs the native package (separate
  repository). This spec publishes the prebuilt packages (FR-014); the plugin
  decides how to use them.
- Fault displacement profiles given as Python functions: they use the Python
  path.

## Clarifications

Answered on 2026-11-27.

- **C-001 Exactness.** Must the native results be bit-identical to the
  Python path? **Answer**: no. They must pass the golden tests with the
  tolerances of spec 003 FR-004, with the same NaN points and model ids
  (except at a unit boundary). The proof of concept differs by at most
  1.9e-4 m.
- **C-002 Default.** Is the native path on by default? **Answer**: yes, when
  the native package is installed and its version matches, with a switch to
  turn it off.
- **C-003 Scope of the first version.** **Answer**: restoration, feature
  values and gradients (stories 2 and 3). Masks and regions stay in Python
  and use features that the native path evaluates.
- **C-004 The Python bookkeeping.** Is the per-call computation of the
  restore order part of this spec? **Answer**: yes, as stage 1 (story 1).
- **C-005 Fallback unit (default; the maintainer can change it).** If one
  fault of a feature is not supported, the whole feature uses the Python
  path, so that all its faults use one restore order.
- **C-006 Threads (default; the maintainer can change it).** The number of
  threads is configurable; the default is all available cores.
- **C-007 Native tolerance (from the review).** The proof of concept differs
  by more than the golden tolerance. Must the native path meet the golden
  tolerance of spec 003, or a looser tolerance of its own? **Answer**: spec 003 FR-004
  (1e-8 × the box size). The plan must find and remove the difference of the
  proof of concept.
- **C-008 The switch.** How does a user turn the native path off?
  **Answer**: with a LoopStructural configuration setting or an environment
  variable. The setting wins over the variable.
- **C-009 Version mismatch.** **Answer**: use the Python path and log one
  warning in a session that names both versions.
- **C-010 Supported profiles.** **Answer**: all built-in fault displacement
  profile types of LoopStructural (including the default profile). A
  profile given as a Python function uses the Python path.
- **C-011 Prebuilt packages.** Are native packages for Windows, macOS and
  Linux part of this spec? **Answer**: yes. Most QGIS users are on Windows,
  and story 4 has value only with prebuilt packages.
