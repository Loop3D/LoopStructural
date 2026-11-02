---
name: spec-implementer
description: Make the code changes for all specs that a spec-cascade run changed. Finds every spec that now has open or new tasks, puts the specs in build order, and builds each one from its tasks.md in its own worktree with the rules of the spec-orchestrator agent. Use after the spec-cascade skill, when the user asks to implement the cascade, to build all changed specs, or to bring the code in line with the new task lists. Do not use for one spec with an unchanged task list (use the spec-orchestrator agent), and do not use to change specs or plans.
---

# Spec implementer

The spec-cascade skill has changed several specs. This skill makes the code
changes for all of them. It does not build a spec by itself: it finds the
specs, puts them in order, and builds each one with the procedure in
`.agents/skills/spec-orchestrator/SKILL.md`. Read that file first. This skill adds
the work that is across specs.

Run this skill in the main session, because it must start sub agents.

## Input

The user gives the source spec of the cascade (for example `013`), or a list
of specs. If the user gives nothing, use the specs that changed last in
`git status` or in the last cascade commit. If you cannot find any, ask.

If the user says "dry run", do steps 1 to 5 only. Change nothing: make no
worktree, start no sub agent, edit no file, and run no `cargo` command.

## Rules

- `.specify/memory/constitution.md` and `AGENTS.md` have priority. Do not
  write production code for a task that has no spec and plan.
- Do not change specs, plans or the shape of the task lists. If a spec or a
  plan is wrong, stop and tell the user to run the spec-cascade,
  spec-reviewer or spec-planner skill. The only edits to `tasks.md` are the
  checkboxes (`[x]`, `[!]`).
- Use one worktree and one branch for each spec. Never switch the branch of
  the main checkout, because other agents can work there.
- Tests come first and must fail first (principle V). Never skip, weaken or
  `#[allow]` a test or a lint to pass a check. Never use `--no-verify`.
- Do not push and do not open a pull request unless the user asks.
- Write all reports, comments and commit messages in ASD-STE100 Simplified
  Technical English.

## Steps

1. Read `.specify/memory/constitution.md`, `AGENTS.md`,
   `.agents/skills/spec-orchestrator/SKILL.md` and `specs/README.md`.
2. **Find the specs in scope.** Use the cascade diff (`git diff` or the cascade
   commit) and the task lists. A spec is in scope if its `tasks.md` has a task
   that is not done and is new or changed by the cascade (`[ ]`, or a new task
   that replaces a done task). Ignore a spec with only `[x]` tasks. Report a
   spec with `[!]` tasks to the user: do not run it without a decision.
3. **Check that each spec is ready.** It must have `spec.md`, `plan.md` and
   `tasks.md`, no `[NEEDS CLARIFICATION]` item, and a task graph that is
   valid (all `Needs` ids exist, no cycle, correct waves, no shared file in
   one wave). If a spec is not ready, show the problem and ask the user. Do not
   repair it here.
4. **Put the specs in build order.** Build a graph of the specs:
   - a spec is after the specs that appear in its `Blocked by`, "Depends on"
     or "Needs spec NNN" lines;
   - a spec that uses a function, type or message that another spec creates
     is after that spec (see "Handover" in the plans);
   - two specs that edit one file, one protocol version, or one generated
     file (the schema or the fixtures) cannot run at the same time: put them
     in order. The spec that has the lower protocol version goes first.
   If the graph has a cycle, stop and ask the user. Specs with no link between
   them can run in parallel (see "Parallel specs").
5. **Show the plan.** In the chat show: the specs in order, for each the
   number of open tasks by tier (`[S]` and `[M]`), the worktree path and the
   branch, the specs that wait for another spec, and the `Blocked by` tasks
   that the order frees. Ask the user to confirm with the available user-question tool (or chat):
   run all (Recommended), run only some specs (the user names them), or stop.
   Stop here for a dry run.
6. **Build each spec in order** (see "Build one spec").
7. **Check the whole result** (see "Final checks").
8. **Update the index** and report (see "Finish").

## Build one spec

For each spec, follow the sections "Start", "Run loop", "Sub agent prompt",
"Check a result", "Failure" and "Commits" of
`.agents/skills/spec-orchestrator/SKILL.md` exactly. Do not make a weaker version.
These points are different or new:

- **Worktree.** Use `../geoviewer-NNN` and the branch `feat/NNN-name`. If it
  exists, use it. If the spec needs the code of an earlier spec that is not
  merged, make the new branch from the branch of that spec, not from `main`.
  Note this in the report.
- **Only the open tasks.** Skip tasks with `[x]`. Run the tasks that the
  cascade added or changed in the normal way. A task that replaces a done task
  must say what it changes: read the done task and its commit before you send
  it.
- **Handover tasks.** A task from the handover of another spec (for example
  "call `check_object_size` in the scene file loader") needs the function of
  the earlier spec. Check that it exists in the worktree
  (repository text search). If it does not, the order is wrong: stop and tell the user.
- **Open points.** If a task has `Blocked by` that a user answer frees, ask
  once with the available user-question tool (or chat) and write the answer in `plan.md` section
  "Answers". This is the only edit to a plan.
- **Spec status.** When all tasks of a spec are `[x]` and the final checks of
  that spec pass, go to the next spec. Do not say "done" for a spec with a
  `[ ]` or `[!]` task.

## Parallel specs

Run two specs at the same time only if all of these are true: neither is in
the graph of the other; their `Files` lists do not overlap; neither changes
the protocol version, the schema or the fixtures; and the user agreed. Each
spec has its own worktree and its own run of the loop. Use the lower of the configured worker limit and available runtime slots,
across all specs. If you are in doubt, run the specs one after
the other.

## Final checks

After the last spec, in each worktree:

- `cargo fmt --all --check`
- `cargo clippy --workspace --all-targets -- -D warnings`
- `cargo test --workspace`

When two or more specs were built, also check that the work agrees:

- Merge the branches in build order into a throw-away branch in a separate
  worktree (`git worktree add ../geoviewer-integration -b tmp/integration`).
  Do not merge into `main` or into a spec branch. Run the three commands there.
  Delete the throw-away branch and worktree afterward, and tell the user about
  any conflict, with the files.
- Check that the protocol version, the schema, the fixtures and the error
  codes are the same in all merged specs (repository text search).
- Check that every FR and SC of every built spec has a passing test.

If a command fails, say so with the output. Find the spec that causes the
failure and send a fix task for that spec. Do not hide a failure.

## Finish

- Update the status of each built spec in `specs/README.md` (for example
  `Implemented`, with the open manual checks). This is the only edit outside
  the code, the tests, `tasks.md` and the "Answers" in `plan.md`.
- Commit in each worktree as the orchestrator rules say: Conventional
  Commits, one commit for each passed wave, sub agents never commit.
- Show a short report in the chat, with one block for each spec:
  - the worktree path and the branch;
  - tasks done, failed (`[!]`) and blocked, with the reason;
  - sub agent runs by model, and the tasks that moved from `[S]` to `[M]`;
  - the result of the three final commands.
  Then show: the build order that you used, the result of the integration
  check, the specs that need the spec-cascade, spec-reviewer or spec-planner
  skill again, and the next step (for example: open pull requests in build
  order).

## Clean up and disk use

Build artifacts and worktrees fill the disk. Keep them small.

- Use one shared target directory for all worktrees. Set
  `CARGO_TARGET_DIR` to `../geoviewer-target` in each `cargo` command. Do
  not let each worktree make its own `target/` directory.
- Run each heavy `cargo` command through the build lock. Use
  `scripts/cargo-lock.sh <args>` in Bash or `scripts\cargo-lock.ps1 <args>`
  in PowerShell (example: `scripts/cargo-lock.sh test --workspace`). Only
  one cargo command can run at one time on the machine. The other agents
  wait. This is for `check`, `build`, `test` and `clippy`. Do not remove
  the lock file by hand. Do not start `cargo` in the background.
- Do not build in release mode. Use `cargo check` while you work. Run
  `cargo test` and `cargo clippy` only when a spec is done.
- Delete the throw-away integration worktree and its branch as soon as the
  integration check is done. Run `git worktree prune` afterward.
- Remove a spec worktree (`git worktree remove ../geoviewer-NNN`) only when
  all work is committed and the user agrees, or when the pull request is
  merged. Never remove a worktree with uncommitted or unpushed work.
- Do not keep extra files: no log dumps, no test output files. Put
  temporary files in the scratchpad and delete them when done.
- List the worktrees that remain (`git worktree list`) in the final report.
