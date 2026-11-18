# Feature spec: Faster evaluation of faulted models

**Feature**: `003-fault-evaluation-performance` | **Created**: 2026-11-01 | **Status**: Clarified 2026-11-18 (scope reduced: no cache between calls; SC-006 and SC-008 changed)

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
- Related work on master (2026-11-18):
  - The structured grid code of `loop_common` moved into
    `StructuredGrid3DGeometry` (#367, #369). The faster grid lookup of user
    story 5 must use the new place of that code.
  - Spec `005-coordinate-frame-correctness` changes how the fault frames and
    the coordinate transforms are made. The golden results of this spec can
    change only under C-006.
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
  - **standard suite**: the set of benchmark cases that a developer runs
    for each measurement (FR-012a).
  - **full suite**: all the cases of FR-012, including the large ones that the
    standard suite leaves out. A developer starts it on purpose.
  - **baseline**: the benchmark results of the standard suite on the code
    before this feature.

## User scenarios and tests

### User story 1: Keep the geology the same (priority: P1)

A developer changes the evaluation code to make it faster. The tests tell
them at once if the model output changes: a unit moves, a fault moves, or a
surface changes. A correctness reference may change only for an affected result
when analytic evidence proves a coordinate defect and the review records the
expected physical effect (C-006). Geometry and frame construction repairs
belong to spec 005.

**Why this priority**: a faster model that gives a different geology is a
defect. The golden tests must exist and pass on the present code before any
other work in this spec starts.

**Independent test**: run the golden tests on the code before the change.
They pass. Change one fault displacement in the code by 1 %. They fail.

**Acceptance scenarios**:

1. **Given** the golden results for a recorded correctness epoch, **when**
   the golden tests run on that code, **then** all of them pass. A later
   reference may replace an affected result only under C-006; all unaffected
   results stay fixed.
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
4. **Given** one model, **when** two threads evaluate it at the same time,
   **then** each thread gets the result that it gets when it runs alone
   (FR-010a).
5. **Given** the public evaluation functions before and after this feature,
   **when** a caller uses them with the same arguments, **then** the
   functions, the arguments and the return types are the same (FR-011).
6. **Given** a model that the caller evaluated, **when** the model changes and
   the caller evaluates the same points again, **then** the result is the
   result of the changed model, because no restored points stay after an
   evaluation call (FR-009).

### User story 3: Measure the cost with benchmarks (priority: P2)

A developer runs one command. It builds the golden models, runs the benchmark
cases, and writes the results to a file. A second command compares two result
files and shows which cases became slower or faster.

**Why this priority**: without numbers, nobody can say if a change helps, or
if a later change makes the code slow again. The benchmarks are a tool for
developers on their own computers. CI does not run them (C-005).

**Independent test**: run the standard suite two times on the same computer
and the same code. Compare the two files. No case shows a change of more than
the noise limit (SC-007).

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
  change (new data, new fault displacement, a solve, or a changed coordinate
  transform) must use the new model. No restored points of the first call can
  stay in use (FR-009).
- **An evaluation call that fails.** The next call must not use a part of the
  work of the failed call.
- **Very many points.** Extra memory must stay within FR-010, also for 10
  million points.
- **Threads.** `use_threads` in `_fault_segment.py` must give the same results
  as without threads, or stay off as today.
- **Saved and loaded models.** A model that is saved and loaded again (for
  example with `to_file` and `from_file`) must evaluate as today.

## Requirements

### Functional requirements

**Golden tests**

- **FR-001**: The repository MUST contain golden models that cover:
  - 0, 1, 2 and 4 faults, with 0, 1, 2 and 3 unconformities;
  - a fault that stops at an unconformity (two different fault sets);
  - two faults where one fault cuts the other;
  - an onlap unconformity and a domain fault;
  - a folded feature with a faulted fold frame;
  - the plugin model of `research.md`, as it is (C-002). Its data and its
    solved fields are in `data/plugin_model/` (see FR-001a).
- **FR-001a**: A golden model MUST be able to use stored solved fields, so
  that its golden results do not depend on the solver. The solvers give
  results that differ by more than the tolerances of FR-004 (`cg` and
  `native` differ by about 1e-5 of the range).
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
  results MUST be visible in the review of a pull request. A reviewed analytic
  correctness correction MAY update an affected result (C-006). The review
  record MUST name the coordinate defect, its expected physical effect, and
  the result files that change. Unaffected results MUST stay unchanged.
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
  the golden results within FR-004. A reviewed analytic correctness correction
  may define a new correctness reference epoch for affected results (C-006).
  It MUST NOT change numeric tolerances, performance limits, or unaffected
  results. A correctness epoch MUST identify its code revision and reference
  files separately from the performance baseline epoch.
- **FR-007**: In one evaluation call, each fault set MUST be applied to a set
  of points at most once, also when the feature and the regions of the
  feature use the same fault set.
- **FR-008**: In one call of `evaluate_model` or `evaluate_model_gradient`,
  the number of fault applications MUST be at most the sum, over the
  different fault sets of the model, of the number of faults in each fault
  set.
- **FR-009**: Restored points MUST NOT be kept after the evaluation call that
  made them ends, also when the call fails (C-012). A later evaluation call
  MUST make its restored points again from the present model.
- **FR-010**: The extra memory of one evaluation call MUST be at most the
  restored points of each different fault set (3 floats for each point and
  each fault set). This memory MUST be free when the call ends.
- **FR-010a**: The evaluation MUST give the same results when two threads
  evaluate the same model at the same time (the QGIS plugin evaluates in
  background threads). This feature MUST NOT add a lock that lets only one
  evaluation call of a model run at a time (C-012).
- **FR-011**: The public API MUST stay the same: the same functions, the same
  arguments and the same return values.

**Benchmarks**

- **FR-012**: The repository MUST contain a benchmark suite that uses the
  golden models. The full suite MUST cover:
  - `evaluate_value` of the oldest group, all groups one by one,
    `evaluate_model`, `evaluate_gradient` and `evaluate_model_gradient`;
  - 5 000, 20 000 (the plugin model and the model with one fault and no
    unconformity only), 50 000, 100 000 and 500 000 points;
  - all the fault and unconformity counts of FR-001.
- **FR-012a**: The standard suite MUST be a part of the full suite. It MUST
  contain every case that a success criterion (SC-002 to SC-009) uses. The
  baseline of FR-016 MUST contain each case of the standard suite. The
  standard suite MUST contain:
  - all the 5 000-point cases (every golden model, every function);
  - the 20 000-point cases of the plugin model and of the model with one fault
    and no unconformity;
  - 9 cases at 50 000 points and 7 cases at 100 000 points (C-016). They come
    from the fixed random choice of 2026-11-02 (fixed seed), which picked
    12 cases at each size from the models with 0 or 1 faults and from the
    stopped-fault, onlap-domain and faulted-fold models. The models that cost
    the most (the models with 2 or 4 faults, the cutting-faults model and the
    plugin model) are not in this choice, because they are measured at 5 000
    and 20 000 points. The choice is the same in every run, so two result
    files have the same cases. The gradient cases at 100 000 points are not in
    the standard suite (C-014).
  The time pass runs all the cases of the standard suite. The count pass runs
  only its 5 000-point and 20 000-point cases, because SC-009 reads the count
  of a 20 000-point case. The memory pass runs only its 5 000-point cases
  (C-014). The full suite is
  not limited in time. It is the standard suite plus all the other cases of
  the full suite (FR-012), with all the passes.
- **FR-012b**: The command that runs the suite MUST run the standard suite by
  default, and MUST run the full suite only when the developer asks for it.
- **FR-013**: The suite MUST count the fault applications and the field
  evaluations of each case. The counters MUST be off in normal use, and MUST
  cost nothing when they are off.
- **FR-014**: One command MUST run the suite and write one machine-readable
  result file (see user story 3). A second command MUST compare two files.
- **FR-015**: The benchmarks are for local development only. They MUST NOT
  run in the normal test run or in CI. A developer starts them on purpose.
  They use `pytest-benchmark` (C-004), as a development dependency only.
- **FR-016**: The repository MUST contain the baseline result file of the
  standard suite on the code before this feature, with the computer that made
  it. A baseline of the full suite is not required.

### Key entities

- **Golden model**: a function that builds one model from fixed data, fixed
  settings and a fixed seed. It has a name and a list of the features it has.
- **Golden result**: the stored outputs of one golden model at its fixed
  points (FR-002), with the version of the code that made them.
- **Benchmark case**: a golden model, an evaluation function and a number of
  points.
- **Benchmark result file**: the results of all cases of one run, with the
  versions and the computer.
- **Suite**: the standard suite or the full suite (see Context).
- **Fault set**: the ordered list of faults that a feature applies. Two
  features with the same faults in the same order have the same fault set.

## Success criteria

- **SC-001**: All golden tests pass on the recorded pre-performance correctness
  reference and on the code after this feature. An affected result can change
  only under C-006. If a future correction changes a reference, its new
  correctness revision and affected files MUST be recorded separately from the
  performance baseline. Existing reference and manifest records remain valid
  historical records. A performance comparison MUST use comparable correctness
  epochs and MUST report correctness changes separately from performance.
- **SC-002**: For the golden model with 4 faults and 3 unconformities and
  5 000 points, `evaluate_value` of the oldest group makes 4 fault
  applications (baseline: 16) and takes at most 40 % of the baseline time
  (baseline: 1.96 s in the baseline file of FR-016; C-015).
- **SC-003**: For the same model, one `evaluate_model` call, which evaluates
  all groups, makes 4 fault applications (baseline: 52). Evaluating all groups
  one by one takes at most 15 % of the baseline time (baseline: 6.48 s in the
  baseline file of FR-016; C-015, C-017).
- **SC-004**: For the golden models with 4 faults and 5 000 points, the time
  of one `evaluate_model` call with 3 unconformities is at most 1.5 times the
  time with 0 unconformities (baseline: 13.3 times in the baseline file of
  FR-016; C-017).
- **SC-005**: For the plugin model of `research.md`, the evaluation of the
  group at 20 000 points takes at most 50 % of the baseline (baseline: 8.2 s).
- **SC-006**: For each case of the memory pass of the standard suite, the peak
  extra memory is at most 120 % of the baseline, or at most 16 MiB
  (16 777 216 bytes) more than the baseline (C-013). A case passes when one
  of the two limits holds. A baseline of 0 bytes uses the 16 MiB limit only.
- **SC-007**: Two runs of the standard suite on the same computer and the same
  code differ by less than 10 % in each case. The comparison report uses this
  value as the noise limit. The limit covers the cases with a baseline time of
  1 s or more. A shorter case uses more rounds, and the report shows its
  noise, but the case does not fail this criterion.
- **SC-008**: One run of the standard suite, with all its passes, takes less
  than 15 minutes on the computer in `research.md`, on the code after this
  feature (C-014).
- **SC-009** (user story 5): the time for each point of one fault application
  is at most 50 % of the baseline (about 34 µs for each point).
- **SC-010**: Removed on 2026-11-18 (C-012). The number is not used again.

## Out of scope

- A change of the fault geometry, the fault displacement or the kinematics
  as part of this performance feature. Spec 005 owns the separate frame and
  coordinate repairs. Their affected references may change only under C-006.
- A change of the default number of fault steps, unless the golden tests pass
  with the new value (this is part of user story 5).
- Keeping restored points or other evaluation results between evaluation
  calls (C-012). The QGIS plugin can keep its own output if it evaluates the
  same points again.
- Changes to the solvers (`cg`, `admm`, `native`).
- The algorithm that makes the surfaces from the field (marching cubes).
- Changes in the QGIS plugin. The plugin gets the improvement through a new
  LoopStructural release.
- GPU work.

## Clarifications

Answered on 2026-11-01. Coordinate-frame amendments from spec
`005-coordinate-frame-correctness` were added on 2026-11-03. They define the
contract at this feature's boundary; spec 005 owns repair of frame construction
and world/local conversion defects.

- **C-001 Tolerances.** Are the tolerances of FR-004 correct? A change that
  only reorders floating-point operations can move a point by about 1e-12 of
  the box, and this can change the stratigraphic id of a point that is almost
  on a boundary.
  **Answer**: yes, the tolerances of FR-004 are correct.
- **C-002 Real-data model.** Which dataset is the real-data golden model?
  **Answer**: the plugin model of `research.md`, as it is. Its data can be in
  the repository. It is in `data/plugin_model/`: the calls that the plugin
  made, the data of each feature, and the solved fields. `plugin_model.py`
  builds it without QGIS, and gives the same values as the model in QGIS at
  5000 fixed points (difference 0).
- **C-003 Results between calls.** Is it in scope to keep restored points
  between calls when the faults did not change?
  **Answer** (2026-11-01): yes. **Replaced on 2026-11-18 by C-012: no.**
- **C-004 Benchmark tool.** Is a new development dependency allowed?
  **Answer**: yes, `pytest-benchmark` or `asv`. This spec uses
  `pytest-benchmark`, because it runs with pytest, writes JSON and compares
  two runs. The plan can choose `asv` if a history across many commits is
  necessary.
- **C-005 CI.** Must a CI job run the benchmarks?
  **Answer**: no. The benchmarks are for local development only (FR-015).
- **C-006 Coordinate contract and reviewed correctness epoch** (from spec
  `005-coordinate-frame-correctness`, 2026-11-03). Which coordinate space do
  the inputs, restored points, references, and low-level arrays use, and may a
  demonstrated coordinate defect change a frozen result?
  **Answer**: public modelling inputs, fixed query points, restored points,
  and public evaluation results use world coordinates. Arrays submitted to a
  support or solver use that support's computational coordinates; each boundary
  converts once. The transform is the mapping between those spaces and is part
  of the restoration dependency. The coordinate owner must rebuild dependent
  local state after an accepted change or reject the change before mutation.
  No restored points stay after an evaluation call (FR-009, C-012), so a
  changed transform cannot give a stale restoration. A future result may change
  only when analytic evidence proves the defect and the review records its
  expected physical effect. Append new correctness provenance for that change;
  keep current references and manifests as valid historical records. Keep
  unaffected references, numeric tolerances, and performance limits fixed.
  Record each changed correctness reference separately from its performance
  baseline, with provenance that permits comparable runs.

Answered on 2026-11-02. The first full run of the 650 cases on the code
before this feature did not finish in 4 hours. The 500 000-point cases were
90 % of the time that was measured. No success criterion uses them. (These
answers were C-006 to C-010 in the text of 2026-11-02. The cascade of spec
005 removed them by mistake on 2026-11-04. They were put back with new
numbers on 2026-11-18.)

- **C-007 Suite time budget.** What is the time budget of one run of the
  standard suite?
  **Answer**: about 1 hour on the code before this feature, with all passes.
  (Changed by C-014: the count and memory passes now run fewer cases.)
- **C-008 Standard suite content.** Which cases are in the standard suite?
  **Answer**: all the 5 000-point cases; the 20 000-point cases of the plugin
  model and of the model with 1 fault and no unconformity; a fixed choice of
  50 000 and 100 000-point cases of the cheaper models (C-011). The full suite
  stays available on request (FR-012a, FR-012b).
- **C-009 Noise limit and short cases.** How does SC-007 treat short cases?
  **Answer**: the 10 % limit covers the cases with a baseline time of 1 s or
  more. Short cases use more rounds and are reported, but they do not fail the
  criterion (SC-007).
- **C-010 Time limit after the change.** What is the limit of SC-008?
  **Answer**: 15 minutes for the standard suite on the code after this
  feature (SC-008).
- **C-011 Large point counts.** A measurement on 2026-11-02 shows that the
  cutting-faults model needs 13 to 66 s for one call at 5 000 points, and the
  model with 4 faults and 3 unconformities needs 42 s for the five functions.
  At 50 000 points these models would need hours. How does the standard suite
  keep the 1 hour budget?
  **Answer**: a fixed random choice of cases at 50 000 and 100 000 points,
  taken from the cheaper models. Every time is measured; no time is inferred
  (FR-012a). The success criteria use only the 5 000 and 20 000-point cases.

Answered on 2026-11-18, after a review of the implementation and of the
failed acceptance run T030 (`t030-diagnosis.md`, `profile-d4a.md`).

- **C-012 Cache between calls.** The first evaluation call is about 29 times
  faster than the baseline without a cache (`f4_u3`, `evaluate_model`,
  5 000 points: 6.34 s to 0.22 s). A cache between calls makes a second call
  on the same points about 10 times faster again. It needs a revision counter
  in many setters, a hash of the state of each fault set in each call, a lock
  for each model and up to 64 MiB of memory. It caused most of the SC-006
  failures. Do we keep it, move it to a later spec, or remove it?
  **Answer**: remove it. This spec keeps no restored points between calls
  (FR-009, FR-010) and adds no lock for each model (FR-010a). SC-010 and the
  warm benchmark state are removed. No later spec is planned. If the plugin
  evaluates the same points again, it can keep its own output.
- **C-013 Memory limit for small cases.** The first acceptance run failed the
  memory criterion (SC-006) in 45 cases. In 34 of them the baseline is below
  10 MB, and in 16 the
  extra memory is below 1 MiB. How does SC-006 treat a small baseline?
  **Answer**: a case passes when its peak extra memory is at most 120 % of
  the baseline, or at most 16 MiB more than the baseline (SC-006).
- **C-014 Suite time.** The first acceptance run took about 3 000 s against
  the 900 s limit of SC-008. More than half of the time was the priming call
  of the warm cases. The two `onlap_domain` gradient cases at 100 000 points
  took 24 % of the time. How does the suite meet SC-008?
  **Answer**: keep the limit of 15 minutes, and make the standard suite
  shorter: no warm state (C-012), the count and memory passes run only the
  5 000-point cases, and the gradient cases at 100 000 points move to the
  full suite (FR-012a). The planner added the 20 000-point cases to the count
  pass on 2026-11-18, because SC-009 needs their counts.
- **C-015 Baseline of the time ratios.** SC-002 and SC-003 gave the baselines
  of `research.md` (2.55 s and 8.65 s, from a different script). The baseline
  file of FR-016 has 1.96 s and 6.48 s for the same golden-model cases. Which
  baseline do the time ratios use?
  **Answer**: all the time ratios use the baseline file of FR-016 (the ratios
  of SC-002, SC-003, SC-005 and SC-009). The numbers of `research.md` are for
  information only. The ratio limits do not change.
- **C-016 Large cases without the warm state.** The fixed choice of
  2026-11-02 has 24 large case ids with a cold or a warm state. How does the
  standard suite pick its large cases now?
  **Answer**: use the same 24 ids without the state. Merge the ids that
  become equal, and remove the 5 gradient cases at 100 000 points. This
  gives 9 cases at 50 000 points and 7 cases at 100 000 points. The baseline
  file has a measured value for each of them: the code before this feature
  has no cache, so its warm case was one cold call. No new baseline run is
  necessary.
- **C-017 Count and ratio of "all groups" without reuse between calls.**
  SC-003 and SC-004 measured the evaluation of all groups one by one (one
  `evaluate_value` call for each group). After C-012, each call applies its
  faults again: T051 measured 16 fault applications for `f4_u3` (4 for each
  call) and a ratio of 3.8 between `f4_u3` and `f4_u0`. One `evaluate_model`
  call gives 4 applications and a ratio of 1.03. Which call do SC-003 and
  SC-004 measure?
  **Answer** (2026-11-18): the count of SC-003 (4) and the ratio of SC-004
  (at most 1.5) use one `evaluate_model` call (benchmark case `model_ids`).
  The time limit of SC-003 (at most 15 % of the baseline) stays on the
  evaluation of all groups one by one (case `all_group_values`). The limits
  do not change.
