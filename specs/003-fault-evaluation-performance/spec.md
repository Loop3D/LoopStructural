# Feature spec: Faster evaluation of faulted models

**Feature**: `003-fault-evaluation-performance` | **Created**: 2026-11-01 | **Status**: Clarified 2026-11-01 (C-002 waits for the refined plugin model)

**Input**: The user renders models in the QGIS plugin. A model with faults and
unconformities is slow to render, and it becomes slower with each fault and
each unconformity that the user adds. The user wants golden tests that freeze
the present results, and benchmarks that measure the cost, before anyone
changes the evaluation code. Then the user wants the evaluation to be faster.

<!--
A spec says WHAT the user needs and WHY. It does not say HOW.
Write in ASD-STE100 Simplified Technical English.
The analysis of the present code and the baseline numbers are in research.md.
-->

## Context

- To show a model, the viewer, the surface export, the block model and the
  cross sections all evaluate the features of the model at many points. They
  call `evaluate_value`, `evaluate_gradient`, `evaluate_model` or
  `evaluate_model_gradient`.
- For a faulted feature, most of the time goes into the faults, not into the
  solved field. In a real model with one group and six faults, the field of
  the group takes 0.01 s for 20 000 points. The same evaluation with the faults
  takes 8.2 s (see `research.md`).
- The same fault work is done many times in one evaluation:
  - A feature applies its faults to the points.
  - Each unconformity above the feature is a region of the feature. To make
    the mask, each region evaluates its own feature, and that feature applies
    its faults to the same points again.
  - When the caller evaluates all the groups (for example `evaluate_model`),
    each group does all of this again.
- So the number of fault applications grows with the square of the number of
  unconformities. With 4 faults and 3 unconformities, an evaluation of all the
  groups (as `evaluate_model` does) applies faults 52 times. The faults change
  the points in the same way each time, so 4 applications are sufficient.
- One fault application is also expensive: about 34 µs for each point, because
  it moves each point in 10 steps and evaluates the fault frame at each step.
- The solver does not change this. The cost is the same with `cg` and with
  `solver='native'`.
- Words used in this spec:
  - **evaluation call**: one call of a public evaluation function, for example
    one `evaluate_model(xyz)`.
  - **fault application**: one call that moves a set of points back across one
    fault (today `FaultSegment.apply_to_points`).
  - **fault set**: the ordered list of faults that a feature applies.
  - **restored points**: the points after all the faults of a fault set are
    applied.
  - **golden model**: a fixed model (data, settings and random seed) that the
    tests build again in the same way each time.
  - **golden result**: a stored output of a golden model, made with the code
    before this feature. The tests compare the new output with it.
  - **benchmark case**: one golden model, one evaluation function and one
    number of points.
  - **baseline**: the benchmark results of the code before this feature.

## User scenarios and tests

### User story 1: Keep the geology the same (priority: P1)

A developer changes the evaluation code to make it faster. The tests tell
them at once if the model output changes: a unit moves, a fault moves, or a
surface changes.

**Why this priority**: a faster model that gives a different geology is a
defect. The golden tests must exist and pass on the present code before any
other work in this spec starts.

**Independent test**: run the golden tests on the code before the change.
They pass. Change one fault displacement in the code by 1 %. They fail.

**Acceptance scenarios**:

1. **Given** the golden results made with the present code, **when** the
   golden tests run on the same code, **then** all of them pass.
2. **Given** a change that moves the points of one fault by a small amount,
   **when** the golden tests run, **then** at least one test fails and its
   message names the golden model, the output and the largest difference.
3. **Given** a change that only changes the order of floating-point
   operations, **when** the golden tests run, **then** they pass (see FR-004
   for the tolerances).

### User story 2: Render a model with many faults and unconformities (priority: P1)

A geologist renders a model with several faults and unconformities in the
QGIS plugin. The rendering time grows about linearly with the number of
faults, and stays about the same when they add an unconformity.

**Why this priority**: this is the problem that the user sees. It becomes
worse with each fault and each unconformity, so real models are the slowest.

**Independent test**: run the benchmark cases of the golden models with 4
faults and 0 to 3 unconformities. Compare the time and the number of fault
applications with the baseline.

**Acceptance scenarios**:

1. **Given** a model with F faults and U unconformities where all the groups
   have the same fault set, **when** the caller evaluates the oldest group,
   **then** the evaluation call makes at most F fault applications (today
   F × (U + 1)).
2. **Given** the same model, **when** the caller calls `evaluate_model`,
   **then** the evaluation call makes at most F fault applications for each
   different fault set in the model (today the number grows with the square
   of U: 52 for F = 4 and U = 3).
3. **Given** the plugin model of `research.md` (one group, six faults, one
   unconformity), **when** the caller evaluates the group at 20 000 points,
   **then** the time is at most half of the baseline.
4. **Given** a model that the viewer evaluated on a grid, **when** the viewer
   evaluates the same grid again and the model did not change, **then** the
   evaluation call makes no fault applications. **When** the model changed
   before the second call, **then** the result is the result of the changed
   model.

### User story 3: Measure the cost with benchmarks (priority: P2)

A developer runs one command. It builds the golden models, runs the benchmark
cases, and writes the results to a file. A second command compares two result
files and shows which cases became slower or faster.

**Why this priority**: without numbers, nobody can say if a change helps, or
if a later change makes the code slow again. The benchmarks are a tool for
developers on their own computers. CI does not run them (C-005).

**Independent test**: run the benchmarks two times on the same computer and
the same code. Compare the two files. No case shows a change of more than the
noise limit (SC-007).

**Acceptance scenarios**:

1. **Given** the benchmark command, **when** a developer runs it, **then** it
   writes one machine-readable file. For each case, the file has the time, the
   number of fault applications, the number of field evaluations and the peak
   memory, and also the versions of Python, NumPy and LoopStructural and a
   description of the computer.
2. **Given** two result files, **when** a developer compares them, **then**
   the report shows the ratio of the times and of the counts for each case,
   and marks each case that is slower than the noise limit.
3. **Given** the cases for 0 to 3 unconformities, **when** the report is made,
   **then** it shows how the time grows with the number of unconformities and
   with the number of faults.

### User story 4: Evaluate gradients of a faulted model (priority: P2)

A geologist shows the gradient (for example the dip) of a faulted feature, or
the code evaluates `evaluate_model_gradient`. This is also faster.

**Why this priority**: gradient evaluation of a faulted feature applies the
faults to more points than value evaluation (four points for each point), so
it has the same problem with a larger cost.

**Independent test**: run the gradient benchmark cases and the gradient
golden tests.

**Acceptance scenarios**:

1. **Given** a faulted feature with unconformities above it, **when** the
   caller evaluates its gradient, **then** the fault applications follow the
   same limits as in user story 2, and the golden gradient tests pass.

### User story 5: Make one fault application cheaper (priority: P3)

A developer makes one fault application faster, so that a model with many
faults is faster also when nothing is repeated.

**Why this priority**: user story 2 removes the repeated work. This story
reduces the cost of the work that stays. It needs more design work, and it
can change the results a little, so it comes last.

**Independent test**: run the benchmark case "one fault, no unconformity" and
the golden tests.

**Acceptance scenarios**:

1. **Given** the case "one fault, no unconformity, 20 000 points", **when**
   the benchmark runs, **then** the time for each point is lower than the
   baseline, and the golden tests pass.

### Edge cases

- **Different fault sets.** A fault can stop at an unconformity (the fault has
  the unconformity as a region), so a group above the unconformity and a group
  below it can have different fault sets. Each group must use its own fault
  set.
- **Faults that cut faults.** The frame of a fault can be faulted by an older
  fault. The order of the faults must stay the same as today.
- **Onlap unconformities and domain faults.** These are also regions. They
  must give the same masks as today.
- **Folded features and fold frames.** A folded feature evaluates its fold
  frame, which can be faulted. It must give the same results as today.
- **`faults_enabled = False`, `ignore_regions = True` and `reverse = True`.**
  Each of these must give the same results as today.
- **Inactive faults (displacement 0).** They must give the same results as
  today.
- **Points outside the model box, and NaN points.** They must give NaN or the
  same values as today.
- **A feature that is not built.** The first evaluation call still builds the
  feature (the builder is called), and builds it only one time.
- **The model changes between two evaluation calls.** A second call after a
  change (new data, new fault displacement, a solve) must use the new model.
  No result of the first call can stay in use.
- **Very many points.** Extra memory must stay within FR-010, also for 10
  million points.
- **Threads.** `use_threads` in `_fault_segment.py` must give the same results
  as without threads, or stay off as today.

## Requirements

### Functional requirements

**Golden tests**

- **FR-001**: The repository MUST contain golden models that cover:
  - 0, 1, 2 and 4 faults, with 0, 1, 2 and 3 unconformities;
  - a fault that stops at an unconformity (two different fault sets);
  - two faults where one fault cuts the other;
  - an onlap unconformity and a domain fault;
  - a folded feature with a faulted fold frame;
  - the plugin model of `research.md`, after the user refines it (C-002). Its
    data is stored in the repository.
- **FR-002**: For each golden model, the golden results MUST contain, at a
  fixed set of points that includes points near faults and near
  unconformities:
  - the output of `evaluate_model` (the stratigraphic ids);
  - the value and the gradient of each feature;
  - the restored points of each fault set;
  - a summary of the surfaces of each group (number of vertices, box and
    area).
- **FR-003**: The golden results MUST be made with the code before this
  feature, and MUST be stored in the repository. A developer MUST be able to
  make them again only with one explicit command, and a change of the golden
  results MUST be visible in the review of a pull request.
- **FR-004**: The golden tests MUST use these tolerances unless the
  clarification of tolerances (C-001) changes them:
  - values, gradients and restored points: an absolute difference of at most
    1e-8 × the size of the model box;
  - stratigraphic ids: equal at all points, except at points where the value
    is within the value tolerance of a unit boundary;
  - surface summaries: number of vertices equal, box and area within 1e-6
    (relative).
- **FR-005**: The golden tests MUST run in the normal test run, and MUST take
  less than 2 minutes in total.

**Evaluation**

- **FR-006**: After this feature, all public evaluation functions MUST give
  the golden results within FR-004.
- **FR-007**: In one evaluation call, each fault set MUST be applied to a set
  of points at most once, also when the feature and the regions of the
  feature use the same fault set.
- **FR-008**: In one call of `evaluate_model` or `evaluate_model_gradient`,
  the number of fault applications MUST be at most the sum, over the
  different fault sets of the model, of the number of faults in each fault
  set.
- **FR-009**: Restored points MAY be kept between evaluation calls (C-003),
  only with all of these conditions:
  - a kept result is used only for the same points and the same fault set;
  - a change of the model MUST make the kept results invalid. These changes
    are included: a solve or a build of a fault or of a feature that it
    cuts, new data, a change of a build argument, a change of a fault
    displacement, of a fault function or of `steps`, a change of a region,
    and a change of `faults_enabled`;
  - a golden test MUST evaluate the same points, change the model, evaluate
    again, and get the result of the changed model.
- **FR-010**: The extra memory of one evaluation call MUST be at most the
  restored points of each different fault set (3 floats for each point and
  each fault set). Restored points that are kept between calls MUST have a
  size limit that the user can set, and MUST be removed first when they are
  the oldest. The user MUST be able to clear them and to turn the keeping
  off.
- **FR-010a**: The evaluation MUST give the same results when two threads
  evaluate the same model at the same time (the QGIS plugin evaluates in
  background threads).
- **FR-011**: The public API MUST stay the same: the same functions, the same
  arguments and the same return values.

**Benchmarks**

- **FR-012**: The repository MUST contain a benchmark suite that uses the
  golden models. The cases MUST cover:
  - `evaluate_value` of the oldest group, all groups one by one,
    `evaluate_model`, `evaluate_gradient` and `evaluate_model_gradient`;
  - 5 000, 50 000 and 500 000 points;
  - all the fault and unconformity counts of FR-001.
- **FR-013**: The suite MUST count the fault applications and the field
  evaluations of each case. The counters MUST be off in normal use, and MUST
  cost nothing when they are off.
- **FR-014**: One command MUST run the suite and write one machine-readable
  result file (see user story 3). A second command MUST compare two files.
- **FR-015**: The benchmarks are for local development only. They MUST NOT
  run in the normal test run or in CI. A developer starts them on purpose.
  They use `pytest-benchmark` (C-004), as a development dependency only.
- **FR-016**: The repository MUST contain the baseline result file of the code
  before this feature, with the computer that made it.

### Key entities

- **Golden model**: a function that builds one model from fixed data, fixed
  settings and a fixed seed. It has a name and a list of the features it has.
- **Golden result**: the stored outputs of one golden model at its fixed
  points (FR-002), with the version of the code that made them.
- **Benchmark case**: a golden model, an evaluation function and a number of
  points.
- **Benchmark result file**: the results of all cases of one run, with the
  versions and the computer.
- **Fault set**: the ordered list of faults that a feature applies. Two
  features with the same faults in the same order have the same fault set.

## Success criteria

- **SC-001**: All golden tests pass on the code before this feature and on
  the code after it.
- **SC-002**: For the golden model with 4 faults and 3 unconformities and
  5 000 points, `evaluate_value` of the oldest group makes 4 fault
  applications (baseline: 16) and takes at most 40 % of the baseline time
  (baseline: 2.55 s on the computer in `research.md`).
- **SC-003**: For the same model, evaluating all groups makes 4 fault
  applications (baseline: 52) and takes at most 15 % of the baseline time
  (baseline: 8.65 s).
- **SC-004**: For the golden models with 4 faults, the time to evaluate all
  groups with 3 unconformities is at most 1.5 times the time with 0
  unconformities (baseline: 12.9 times).
- **SC-005**: For the plugin model of `research.md`, the evaluation of the
  group at 20 000 points takes at most 50 % of the baseline (baseline: 8.2 s).
- **SC-006**: The peak memory of each benchmark case is at most 120 % of the
  baseline.
- **SC-007**: Two runs of the benchmarks on the same computer and the same
  code differ by less than 10 % in each case. The comparison report uses this
  value as the noise limit.
- **SC-008**: The full benchmark suite takes less than 15 minutes on the
  computer in `research.md`.
- **SC-009** (user story 5): the time for each point of one fault application
  is at most 50 % of the baseline (about 34 µs for each point).
- **SC-010**: A second `evaluate_model` of the same points, with no change of
  the model, makes 0 fault applications and takes at most 10 % of the first
  call. After a change of one fault displacement, the second call gives the
  golden result of the changed model.

## Out of scope

- A change of the fault geometry, the fault displacement or the kinematics.
- A change of the default number of fault steps, unless the golden tests pass
  with the new value (this is part of user story 5).
- Changes to the solvers (`cg`, `admm`, `native`).
- The algorithm that makes the surfaces from the field (marching cubes).
- Changes in the QGIS plugin. The plugin gets the improvement through a new
  LoopStructural release.
- GPU work.

## Clarifications

Answered on 2026-11-01.

- **C-001 Tolerances.** Are the tolerances of FR-004 correct? A change that
  only reorders floating-point operations can move a point by about 1e-12 of
  the box, and this can change the stratigraphic id of a point that is almost
  on a boundary.
  **Answer**: yes, the tolerances of FR-004 are correct.
- **C-002 Real-data model.** Which dataset is the real-data golden model?
  **Answer**: the plugin model of `research.md`. Its data can be in the
  repository. The user refines the model first. The golden results of this
  model are made after that. The synthetic golden models do not wait for it.
- **C-003 Results between calls.** Is it in scope to keep restored points
  between calls when the faults did not change?
  **Answer**: yes. The risks are stale results after a model change, memory,
  the cost of the key for many points, and threads. FR-009, FR-010 and
  FR-010a control them.
- **C-004 Benchmark tool.** Is a new development dependency allowed?
  **Answer**: yes, `pytest-benchmark` or `asv`. This spec uses
  `pytest-benchmark`, because it runs with pytest, writes JSON and compares
  two runs. The plan can choose `asv` if a history across many commits is
  necessary.
- **C-005 CI.** Must a CI job run the benchmarks?
  **Answer**: no. The benchmarks are for local development only (FR-015).
