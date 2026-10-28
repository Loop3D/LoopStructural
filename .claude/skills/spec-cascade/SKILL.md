---
name: spec-cascade
description: Carry the changes of one geoviewer spec into all other specs. Reads what changed in the spec, plan and task list of a source spec, finds every other spec that the change affects (spec.md, plan.md, tasks.md), shows the impact list, asks the user, and then edits the other specs to agree. Use when a spec has changed and its tasks, "Changes to other specs" or "Handover to other specs" now require changes in other specs. Do not use to review one spec, to plan a new spec, or to write production code.
---

# Spec cascade

Read `docs/development/coding-tools.md` for tool setup and model routing.

A spec changed. Its task list (or its plan) now needs changes in other specs.
Find these changes, and make them in all affected specs. Change documents
only. Do not write production code.

The source spec is the one that changed. A target spec is any other spec that
the change affects.

## Input

The user gives the number of the source spec or a path (for example `013` or
`specs/013-object-size-limit/`). If the user gives nothing, use the spec
folder that changed last in `git status` or `git log`. If you cannot find
one, ask.

The user can also give the base for the comparison (a commit, a branch or
"working tree"). The default base is `HEAD` for uncommitted changes. If the
source spec has no uncommitted change, use the last commit that touched it
(`git log -n 1 -- specs/NNN-name`) and compare to its parent.

## Rules

- The constitution has priority over all other documents. A cascade must not
  break a rule of it.
- Keep the order spec -> clarify -> plan -> tasks. A change to a `spec.md`
  needs a new plan and task check. A change to a `plan.md` needs a task check.
- Do not guess a decision. Put `[NEEDS CLARIFICATION: ...]` in the target and
  ask the user (see "Ask the user").
- Do not change the source spec. Report a problem in it to the user.
- Do not change production code, fixtures or generated files.
- Never reuse or renumber a task id. Never change a task that is `[x]`: add a
  new task instead (see "Edit rules for tasks.md").
- Work in the correct worktree for the source spec. Do not switch the branch
  of the main checkout, because other agents can work there.
- Write all text in ASD-STE100 Simplified Technical English.

## Steps

1. Read `.specify/memory/constitution.md`, `specs/README.md` and the
   templates in `.specify/templates/`.
2. Find the change of the source spec (see "Find the change").
3. Find the targets (see "Find the targets").
4. Show the impact list and ask the user (see "Ask the user").
5. Edit the targets (see "Edit rules"). Do the targets in number order.
6. Cascade again for each target that changed in a way that affects a third
   spec (see "Second-order changes").
7. Run the checks. Fix the problems.
8. Update `specs/README.md` and show the summary.

## Find the change

Run `git diff <base> -- specs/NNN-name` and read the full diff of `spec.md`,
`plan.md` and `tasks.md`. Make a change list with these items:

- **Tasks**: added, removed and changed tasks (id, tier, `Files`, `Needs`,
  `Wave`, `Blocked by`, `Done when`). A removed or moved task is a change too.
- **Requirements**: added, removed and changed FR, SC and user stories.
- **Contract items**: protocol version, message names, error codes, object
  types, `geoh5` mappings, names of public types and functions, shared files
  (`lib.rs`, `mod.rs`, `Cargo.toml`, the protocol schema, fixtures), constants.
- **Statements about other specs**: the sections "Changes to other specs"
  (`spec.md`), "Handover to other specs" and "Open points" (`plan.md`), and
  the "Not in this task list" note (`tasks.md`).
- **Status**: the new state of the spec in the index (for example `Draft`,
  `Implemented`, or open tasks).

If the diff is empty, tell the user and stop.

## Find the targets

A target is a spec that the change affects. Use all of these sources:

1. The statements about other specs in the change list. Each one names the
   target and the change.
2. A search in every other `specs/*/` folder (repository text search) for each of: the source
   spec number and name; the changed FR, SC and task ids; the changed message
   names, error codes, type names, function names and constants; the changed
   files; the protocol version number.
3. A conflict check against the other specs:
   - the same protocol version number or the same new message name;
   - the same error code with a different meaning;
   - one file in `Files` of tasks of two specs that can run at the same time;
   - one generated file (for example `schema/protocol.schema.json` or the
     fixtures) with two generator tasks;
   - an object type, role or data type that the two specs define in two ways.
4. The specs that depend on the source spec, and the specs on which it depends
   (`Blocked by`, "Needs spec NNN", "Depends on").

Read every target in full before you decide. Do not decide from the search
result alone. A search match can be a reference that is still correct.

## Impact list

Make one row for each required change. Do not list a target that needs no
change. Use this table:

```
| # | Target | File | Place | Change | Source of the change | Type |
|---|---|---|---|---|---|---|
```

- **Place**: the section, FR/SC id or task id.
- **Source of the change**: the item of the change list that causes it.
- **Type**: `mechanical` (the new text follows from the source with no
  decision) or `decision` (a person must choose).

Add these checks for each target:

- Is the target `Draft`, in progress or `Implemented`? A change to an
  implemented target needs new tasks, and its index status changes.
- Does the change require a new plan step or a new task? Is any task now
  wrong, blocked, or not needed?
- Is the order of the specs still right? Does a `Blocked by` change?

## Ask the user

Show the impact list in the chat. Then ask with the available user-question tool (or chat):

1. A first question: apply the `mechanical` changes now? Options: apply all
   (Recommended), apply only some (the user names the rows), do not apply.
   Do not change a target before the user answers.
2. After the edits, ask the `decision` rows. Ask a small batch of questions within the tool limits.
   Give 2 to 4 options. Put your recommended option first and add
   "(Recommended)" to its label. Explain the effect of each option. Do not ask
   a question that the constitution, a spec or the source change already
   answers.
3. After each answer, update the target at once. Replace the matching
   `[NEEDS CLARIFICATION]` mark with the decision. Record the question and
   the answer, with the date, in the "Clarifications" section of the target
   `spec.md`, or in the "Answers" section of the `plan.md`.

If the user says no, stop and change nothing.

## Edit rules

### `spec.md` of a target

- Change only what the impact list names. Keep the template sections in order.
- Keep the spec about WHAT. Do not add crates, files or Rust types.
- Change requirements, scenarios, success criteria, edge cases and out-of-scope
  items together, so that each FR still has a scenario or criterion.
- Add the change to the "Clarifications" section: the date, the source spec
  and a short reason.
- Do not change the spec number or the feature id.

### `plan.md` of a target

- Change the types, messages, files, functions and open points that the source
  change touches. Keep the constitution check true.
- Update the "Requirement coverage" table for each changed or new FR and SC.
- Update the "Handover" or "Open points" entries that name the source spec.
- Do not add a feature that the spec does not need (principle VII).

### `tasks.md` of a target

- **Not done tasks (`[ ]`)**: edit the text, `Files`, `Needs` and `Done when`
  in place. Keep the id.
- **Done tasks (`[x]`)**: do not edit them. If the change makes a done task
  wrong, add a new task with the next free id that changes the result. Say in
  its text which done task it replaces.
- **Failed tasks (`[!]`)**: tell the user. Do not edit them without a
  decision.
- **Removed work**: if a not-done task is no longer needed, mark it
  `[~] ... Removed: <reason>` (keep the line, so that the id is never reused)
  and take it out of every `Needs` list. Tell the user about each one.
- **New tasks**: use the next free id (never reuse one). Follow the task
  format and the model-tier rules of the `spec-planner` skill: tier `[S]` or
  `[M]`, `Files`, `Needs`, `Wave`, `Done when`, and `Blocked by` where it
  applies. A test task comes before its implementation task (principle V).
- **Blocked by**: add `Blocked by: spec NNN` where a task waits for the
  source spec, and remove it when the source spec makes it free.
- **Recalculate the waves**: wave = `1 + max(wave of Needs)`; `Needs: none` is
  wave 1. Fix the waves of all tasks that depend on a changed task. Two tasks
  in one wave must not share a file. If they do, add a dependency and move one
  of them.
- **Update the summary parts**: the count of waves in the "Orchestrator notes",
  the "Not in this task list" note, the "Wave table", the totals, the
  "Critical path", the "Dependencies" text, and the `Needs` list of the last
  `Finish` task (it needs all other tasks).
- **Keep the other parts**: header, rules comment and orchestrator notes
  as the template gives them.

### Cross-spec conflicts

- A protocol version conflict: the spec that lands first keeps its number. Move
  the later spec to the next minor number in its `spec.md`, `plan.md` and test
  tasks. Ask the user which spec lands first if the specs do not say.
- A shared file or generated file in the same wave: make the second task
  `Need` the first, or make one chain. Do not leave two parallel tasks that
  edit one file.
- An error code or message with two meanings: do not choose. Put a
  `[NEEDS CLARIFICATION]` in both specs and ask the user.

## Second-order changes

After you edit the targets, make a change list for each target (same rules as
"Find the change") and look for new impact on other specs. Do this again until
no new row appears. Do not do more than 3 rounds without a check. After 3
rounds, show the user the rows that are left and ask if you must continue.
Never edit the same file for the same reason twice. If a target depends on the
source spec and the source spec depends on the target, that is a cycle: stop,
report it, and ask the user.

## Checks before you finish

For the source spec and for every target that changed:

- Every `Needs` id exists. There is no cycle. Each wave number is correct. No
  two tasks in one wave share a file.
- No task id is reused or renumbered. Done tasks are not edited.
- Every FR and SC of each changed spec is covered by a task (or by a handover
  entry) in its "Requirement coverage" table.
- Every test task comes before its implementation task.
- No task has a `Blocked by` that points to a finished item. No task is
  waiting for a spec that waits for it.
- No `[NEEDS CLARIFICATION]` mark is left, except those that the user chose to
  keep. List them in the summary.
- The spec numbers, the protocol versions and the error codes agree across all
  specs. Run a repository text search for the changed items to make sure.
- The text is in ASD-STE100 Simplified Technical English.

## Finish

- Update the index in `specs/README.md` (status of each changed spec).
- Show a short summary: the source spec; the targets that changed, with the
  files and the number of edits; the tasks added, changed and removed (with
  ids); the open clarifications; the conflicts that you found; and the
  specs that must run the spec-reviewer or spec-planner skill again.
- Follow the commit rules in `AGENTS.md` (Conventional Commits). For example:
  `docs(spec): carry 013-object-size-limit changes into specs 007 to 011`.
  Commit only if the user asks. Use one commit for the whole cascade. Work in
  the correct worktree for the source spec. Do not switch the branch of the
  main checkout.

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
