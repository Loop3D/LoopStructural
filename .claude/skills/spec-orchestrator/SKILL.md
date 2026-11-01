---
name: spec-orchestrator
description: Build a geoviewer spec from its tasks.md. Reads the waves, tiers and dependencies written by the spec-planner skill, sends each task to a sub agent (small model for [S], standard model for [M]), runs the tasks of one wave in parallel, checks each wave, and keeps the checkboxes in tasks.md current. Use when the user asks to implement, build or execute a spec that has a plan.md and a tasks.md. Run it as the main agent, because it must start sub agents.
---

# Spec orchestrator

Read `docs/development/coding-tools.md` for role and model selection.

You build one spec. You do not write production code yourself. You read
`tasks.md`, send each task to a sub agent, check the result, and update the
task list. Write all reports, comments and commit messages in ASD-STE100
Simplified Technical English.

## Input

The user gives a spec number or a path (for example `009` or
`specs/009-orientation-glyphs/`). If the user gives nothing, ask.

If the user says "dry run", do the "Start" steps 1 to 3 and then the "Dry
run report" below. Change nothing: make no worktree, start no sub agent,
edit no file, and run no `cargo` command.

## Dry run report

Show, in the chat:

- the result of each graph check (pass, or the exact problem);
- a table of the waves: for each task give the id, the tier, the model that
  you would use (resolved from the provider policy), the `Files`, and the `Needs`;
- the number of tasks for each model, and the critical path;
- the tasks with `Blocked by`, and the open points that you would ask about;
- the worktree path and branch that you would make;
- any task that is not clear enough for a sub agent (no exact names, no
  `Done when` command, an `[S]` task that fails the `[S]` rules). Give a
  fix for each.

## Start

1. Read `.specify/memory/constitution.md`, `AGENTS.md`, the spec,
   `plan.md` and `tasks.md`.
2. Stop and tell the user if `plan.md` or `tasks.md` is missing, or if
   `tasks.md` has no `Files`, `Needs`, `Wave` and tier fields. Tell them to
   run the spec-planner skill. Do not guess the missing fields.
3. Validate the task graph before you start any work:
   - every `Needs` id exists, and there is no cycle;
   - each wave number is `1 + max(wave of Needs)`;
   - no two tasks in one wave share a file in `Files`.
   If a check fails, show the problem and ask the user. Do not repair the
   graph yourself.
4. Make a separate worktree. Never switch the branch of the main checkout,
   because other agents can work there:
   - run `git status` in the main checkout and note the changes that are
     not yours;
   - run `git worktree add ../geoviewer-NNN -b feat/NNN-name` (use the
     existing worktree and branch if they exist);
   - do all the work in that worktree. Give its absolute path to every sub
     agent.
5. Make the task state. Tasks with `[x]` are done. Do not run them again.
   Tasks with `Blocked by` are not ready. Ask the user about each open point
   that blocks tasks, once, with the available user-question tool (or chat). Record the answer in
   `plan.md` (section "Answers") before you start the task.
6. Run `cargo check --workspace` once. If it fails before any change, tell
   the user and stop.

## Run loop

Do these steps for each wave, in order:

1. Find the ready tasks: not done, not blocked, and all `Needs` done.
2. Start one sub agent for each ready task, all in one message so that they
   run in parallel. Use at most the configured worker limit and the available runtime slots. Run the rest in a
   next batch of the same wave.
3. Wait for all sub agents to finish. Sub agents do not run `cargo`.
4. Check each result by file (see "Check a result"). This step does not
   build.
5. Run the wave check once (see "Wave check"). Attribute each failure to a
   task, and send a fix (see "Failure"). Run the wave check again after the
   fixes.
6. Mark the passed tasks `[x]` in `tasks.md`.
7. At the end of a phase, run the command of its checkpoint. Do not start
   the next phase until it passes.
8. Commit the wave (see "Commits"), then start the next wave.

Test tasks and the tasks that make them pass are in different waves. A test
task must fail first. If a test task passes with no implementation, the test
is wrong: send it back.

## Sub agent prompt

Each sub agent starts with no context. Give it everything it needs:

```
You are working on task TNNN of spec NNN-name in the worktree <ABSOLUTE PATH>.
Work only in that worktree. Do not run git switch, git checkout, git commit
or git stash. Do not change the main checkout.

Read first: .specify/memory/constitution.md and specs/NNN-name/plan.md.

Task (copy of the full task text from tasks.md, with Files, Needs, Done when):
<TASK TEXT>

Rules:
- Change only the files in "Files". If you need another file, stop and
  report it. Do not change it.
- Write the code in the style of the files around it.
- Do not run `cargo` (no check, build, test, clippy or fmt). The
  orchestrator runs one build and one test run for the whole wave. Read
  the code and the tests carefully, so that your files compile and agree
  with the signatures that you read.
- For a test task: write the test. It must fail for the right reason when
  the orchestrator runs it. Do not write the implementation.
- For an implementation task: write the code so that the named tests pass.
  Do not change the tests. If a test is wrong, stop and report it.
- Report the "Done when" command in your final report. Do not run it.
- Write comments and docs in ASD-STE100 Simplified Technical English.
- Do not add a feature, a dependency or an abstraction that the task does
  not name.

Final report (short): status (done / blocked / failed), files changed, the
"Done when" command (not run), and any problem or open question.
```

Resolve `[S]` and `[M]` through `docs/development/coding-tools.md` before
starting workers. Use the runtime's delegation tool and give each worker
fresh context. All workers of one spec share its worktree and edit only
their assigned files. Do not create a separate worktree for each worker.
If delegation or the requested model is unavailable, report the limitation.
Do not silently replace the model or claim that a worker ran.

## Check a result

Do not trust the report. For each task:

1. Run `git status` and `git diff --stat` in the worktree. Check that the
   task changed only files in its `Files` list. A change in another file is
   a failure.
2. Read the diff of every `[S]` task, and of every task with a suspect
   report. A passed test is not enough: check that the code does what the
   task says and no more.
3. Run the `Done when` commands of the wave in the wave check (see "Wave
   check"), not one by one.

Because tasks of one wave run at the same time, attribute changes by file.
If two sub agents leave files that you cannot attribute, send the files to
one fix agent, or check the tasks one at a time with a wave check for each.

## Wave check

Run these once for each wave, in the worktree, through the build lock (see
the `cargo-build` skill for the command form):

1. `cargo fmt --all`, then `git diff --stat` to see that it changed only
   files of the wave.
2. One `cargo test --workspace` command with the filters of all the
   `Done when` commands of the wave (for example `--test NAME` and `--lib
   NAME`), with `--no-fail-fast`. If the wave has no test, run
   `cargo check --workspace`.
3. Read the result for each task:
   - A test task must fail for the right reason (for example a `todo!()`
     panic or an assertion on a missing value). A compile error in a file
     that the task does not own is a failure, unless the plan says that the
     stubs fail that way.
   - An implementation task must have all its named tests pass.
4. If a build error has no clear owner, find the task by the file in the
   error text.

Do not run a separate `cargo` command for each task.

## Failure

- A task fails or breaks the build: send the error text and the diff back to
  a new sub agent for the same task. Say what failed.
- An `[S]` task that fails twice: send it to the `[M]` model from the provider policy.
- An `[M]` task that fails twice: stop that task. Mark it `[!]` in
  `tasks.md` with a one-line reason. Continue with tasks that do not depend
  on it. Tell the user at the end.
- A sub agent reports "blocked" or an open question: if the plan or the spec
  answers it, send the answer and run the task again. If not, ask the user
  with the available user-question tool (or chat), then write the answer in `plan.md`.
- A task needs a file outside its `Files` list: this is a plan error. Ask
  the user if you can add the file to the task. Do not let the sub agent
  change it by itself.
- Fix agents also do not run `cargo`. Send them the error text from the wave
  check, and run the wave check again when they finish.
- A failure in the wave check with no clear owner: find the task by the file
  in the error text. If this does not work, check the tasks of the last
  wave one at a time, with a wave check for each.
- Never skip a failing test, weaken a test, or add `#[allow]` to pass a
  check. Never use `--no-verify`.

## Commits

- Commit in the worktree after each wave that passes the wave check. You are
  the only agent that commits. Sub agents never commit.
- Use Conventional Commits in ASD-STE100: `type(scope): description`, for
  example `feat(glyphs): add orientation module`. Test tasks use `test`.
  Add a co-author line only if the session supplies one.
- Do not push, and do not open a pull request, unless the user asks.

## Finish

The `Finish` tasks of `tasks.md` run as normal tasks. Before you say that the
spec is done, run these commands in the worktree yourself:

- `cargo fmt --all --check`
- `cargo clippy --workspace --all-targets -- -D warnings`
- `cargo test --workspace`

Then check that every FR and SC of the spec has a passing test. If a command
fails, say so with the output. Do not say "done" while a task is `[ ]` or
`[!]`.

Give the user a short report:

- the worktree path and the branch;
- tasks done, failed (`[!]`) and blocked, with the reason for each;
- sub agent runs by model, and the tasks that moved from `[S]` to `[M]`
  (this tells the user if the planner tiers were right);
- the result of the three final commands;
- the next step (for example: open a pull request).

## Clean up and disk use

Build artifacts and worktrees fill the disk. Keep them small.

- Use one shared target directory for all worktrees. The cargo config of
  the user sets one `target-dir` for all worktrees. Never set
  `CARGO_TARGET_DIR`. Do not let a worktree make its
  own target directory. Read the `cargo-build` skill.
- Run each heavy `cargo` command through the build lock. Use
  `scripts/cargo-lock.sh <args>` in Bash or `scripts\cargo-lock.ps1 <args>`
  in PowerShell (example: `scripts/cargo-lock.sh test --workspace`). Only
  one cargo command can run at one time on the machine. The other agents
  wait. This is for `check`, `build`, `test` and `clippy`. Do not remove
  the lock file by hand. Do not start `cargo` in the background.
- Do not build in release mode. Only the orchestrator runs `cargo`: once
  for each wave (see "Wave check"). Run `cargo clippy` only when the spec is
  done.
- Do not keep extra files: no copies of the repository, no log dumps, no
  test output files. Put temporary files in the scratchpad and delete them
  when the task is done.
- When the spec is done, tell the user which worktree to remove. Remove it
  (`git worktree remove ../geoviewer-NNN`, then `git worktree prune`) only
  when all work is committed and the user agrees, or when the pull request
  is merged. Never remove a worktree with uncommitted or unpushed work.
- Always remove the throw-away worktrees and branches that you made
  yourself (for example `tmp/*`) before you finish.
- List the worktrees that remain (`git worktree list`) in the final report.
