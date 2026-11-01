---
name: spec-planner
description: Turn a clarified geoviewer spec (specs/NNN-name/spec.md) into a plan (plan.md) and an orchestrator-ready task list (tasks.md). Each task has a model tier (simple or standard), exact files, a dependency list and a wave number, so an orchestrator agent can send tasks to sub agents in parallel. Use when the user asks to plan a spec, make tasks, or prepare a spec for implementation. Do not use for reviewing a spec or for writing production code.
---

# Spec planner

Read `docs/development/coding-tools.md` for role and model selection.

Read one spec. Write `plan.md` and `tasks.md` in the folder of the spec. The
reader of `tasks.md` is an orchestrator agent. The orchestrator sends each
task to a sub agent. So each task must be complete, alone, and machine
readable. Do not write production code.

## Input

The user gives a spec number or a path (for example `009` or
`specs/009-orientation-glyphs/spec.md`). If the user gives nothing, use the
spec that changed last in `git status`. If you cannot find one, ask.

## Steps

1. Read `.specify/memory/constitution.md`. It has priority over all other
   documents.
2. Read `.specify/templates/plan-template.md` and
   `.specify/templates/tasks-template.md`.
3. Read the spec and `specs/README.md`.
4. Stop if the spec has a `[NEEDS CLARIFICATION]` item. Tell the user to
   run the spec-reviewer skill first. Do not guess an answer.
5. Read the specs that this spec names. Read the code that the plan will
   change (crates, modules, existing tests). Search repository text and file paths. Do not
   plan from memory.
6. Write `plan.md` from the template.
7. Write `tasks.md` with the rules below.
8. Run the checks below. Fix the problems.
9. Add the plan and tasks status to the index in `specs/README.md` if the
   index has a column for it.

## Plan rules

- Follow the template. Fill the constitution check row by row.
- Name the exact types, messages, files and functions. A sub agent must
  not need to make a design decision.
- Put open points in a list at the top of `plan.md`. A task that depends on
  an open point must say so (`Blocked by: OP1`).
- Do not add a feature that the spec does not need (principle VII).

## Task format

Start `tasks.md` with the template header and a short **Orchestrator
notes** block (see below). Then list the tasks by phase. Write each task in
this form:

```
- [ ] T012 [P] [US1] [S] Test: `plane_normal` for 20 fixture pairs (SC-001).
  - Files: `geoviewer/src/orientation.rs`
  - Needs: T010
  - Wave: 2
  - Done when: `cargo test --workspace orientation::plane_normal` fails
    with the stub from T015, then passes after T016.
```

Fields:

- **Id**: `T` and three digits, in the order of the file. Never reuse an id.
- **`[P]`**: the task can run at the same time as the other tasks of its
  wave. Use it only if the `Files` list does not overlap with another task of
  the same wave.
- **`[USn]`**: the user story. Use `[Setup]`, `[Found]` or `[Finish]` for
  tasks with no story.
- **Model tier**: `[S]` or `[M]` (see below). Every task has one.
- **Files**: every file the task creates or changes. This is the lock list
  for the orchestrator. Two tasks in one wave must not share a file.
- **Needs**: the task ids that must finish first. Write `none` if there is no
  dependency. List only direct dependencies.
- **Wave**: the wave number (see below).
- **Done when**: one command or one observable result that proves the task
  is complete. For a test task, say that the test must fail first.
- **Blocked by**: only if the task waits for an open point or for another
  spec. The orchestrator must not start a blocked task.

## Model tiers

Mark a task `[S]` (simple) only if ALL of these are true:

- The task changes one file, or two files that are closely linked (a test
  and the module that it tests).
- The task names exact items: the function signature, the type with its
  fields, the message JSON, or the expected values.
- The task needs no design decision and no knowledge of other modules,
  beyond the files in `Files` and the items named in the task text.
- A test checks the result with a command in `Done when`.
- The work is of a known kind: add an enum variant, a serde field, a doc
  comment, a fixture file, a table-driven test, a re-export, a rename, a
  `Cargo.toml` entry, a user guide page, or a test that copies a scenario
  from the spec.

Mark all other tasks `[M]` (standard). Use `[M]` for:

- Algorithms and numerical code (principle IV).
- Code with threads, channels, GPU resources, shaders or error recovery.
- A change that touches 3 or more files, or that needs a decision.
- A task where the test cannot show a wrong answer easily.

When in doubt, use `[M]`. A wrong `[S]` costs a failed task. A wrong `[M]`
costs only a little more money. Do not use a tier above `[M]`. If a task
needs more, the task is too big: split it.

A task is too big if it takes more than one day, or if you cannot write its
`Done when` as one command. Split it until each part is clear.

To make more tasks simple: first write a task that adds the types and the
function stubs (bodies that fail). Then the tests and the bodies are small
and well scoped.

## Dependencies and waves

Build the dependency graph, then assign waves:

1. A task with `Needs: none` is in wave 1.
2. Any other task is in wave `1 + max(wave of its Needs)`.
3. Tasks in one wave can run at the same time if they do not share a file.
   If two tasks in a wave share a file, move one to the next wave, and add
   the dependency.

Dependency rules:

- A test task comes before the task that makes it pass (principle V). The
  implementation task `Needs` the test task.
- A task that uses a type or a function `Needs` the task that creates it.
  The stubs task is the usual root of a chain.
- Shared files (`lib.rs`, `mod.rs`, `Cargo.toml`, the protocol schema) are
  bottlenecks. Put changes to them in one early task, or in a chain. Do not
  let parallel tasks edit them.
- A generated file (for example `schema/protocol.schema.json`) has one
  generator task. Other tasks `Need` that task.
- A user story phase must be testable alone (spec rule). Do not make a
  story `Need` a task of a later story unless the task is in Foundation.
- The `Finish` tasks (docs, `cargo fmt`, `cargo clippy`, `cargo test`) are in
  the last wave. The last task `Needs` all others.
- Do not make a cycle. Do not add a dependency "to be safe". An extra
  dependency removes parallel work.

## Orchestrator notes block

Put this block after the header of `tasks.md`. Fill in the numbers.

```
## Orchestrator notes

- **Worktree**: implement this spec in its own worktree. Do not switch the
  branch of the main checkout.
- **Rules for each sub agent**: read `.specify/memory/constitution.md`; read
  `plan.md`; change only the files in the task `Files` list; for a test
  task, run the test and see it fail; write comments and docs in
  ASD-STE100 Simplified Technical English.
- **Tiers**: `[S]` = send to a small, fast model. `[M]` = send to the
  standard model. If an `[S]` task fails twice, send it to `[M]`.
- **Waves**: N waves. Start a wave when all its `Needs` are done. Run the
  tasks of a wave in parallel.
- **After each wave**: run `cargo check --workspace`. Fix a failure before
  the next wave.
- **Checkpoints**: at each checkpoint run the named command. Do not start
  the next phase until it passes.
- **Blocked tasks**: do not start a task with `Blocked by`.
- **Before done**: `cargo fmt --all --check`,
  `cargo clippy --workspace --all-targets -- -D warnings`,
  `cargo test --workspace`.
- **Commits**: Conventional Commits, `type(scope): description`.
```

After the tasks, add two summary tables:

```
## Wave table
| Wave | Tasks | Tier [S] | Tier [M] |
|---|---|---|---|

## Critical path
T001 -> T002 -> ... (the longest chain; it sets the minimum time)
```

Keep the `## Dependencies` section of the template. Write the graph in
words for a human reader, in a few lines.

## Checks before you finish

- Every FR and SC of the spec is covered by at least one task. Make the
  list, and put it in `plan.md` if a gap exists.
- Every test task comes before its implementation task.
- Every task has: a tier, `Files`, `Needs`, `Wave`, `Done when`.
- No two tasks in one wave share a file.
- No cycle. Every `Needs` id exists. Every wave number is correct.
- Each `[S]` task passes the full `[S]` test above.
- Each task text can be understood alone, with `plan.md` and the files in
  `Files`.
- The text is in ASD-STE100 Simplified Technical English.

## Finish

- Show a short summary: the number of tasks, the number of `[S]` and `[M]`
  tasks, the number of waves, the critical path, and the open points that
  block tasks.
- Follow the commit rules in `AGENTS.md`, for example
  `docs(spec): add plan and tasks for 009-orientation-glyphs`. Commit only
  if the user asks. Work in the correct worktree for this spec. Do not
  switch the branch of the main checkout.

## Clean up and disk use

- Do not make a new worktree if you can use an existing one. If you make a
  worktree for this work, tell the user which one, and remove it
  (`git worktree remove <path>`, then `git worktree prune`) only when all
  work is committed and the user agrees. Never remove a worktree with
  uncommitted or unpushed work.
- Do not run `cargo build` or other builds. This skill needs none. If you
  must run a `cargo` command, use `cargo check --workspace` through `scripts/cargo-lock.sh`. Never set
  `CARGO_TARGET_DIR`, so that no extra target directory is made.
- Put temporary files in the scratchpad and delete them when done.
