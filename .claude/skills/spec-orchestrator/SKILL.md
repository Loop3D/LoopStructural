---
name: spec-orchestrator
description: Build a geoviewer spec from its tasks.md. Reads the waves, tiers and dependencies written by the spec-planner skill, sends each task to a sub agent (small model for [S], standard model for [M]), runs the tasks of one wave in parallel, checks each wave, sends each wave to a code reviewer agent, and keeps the checkboxes in tasks.md current. Use when the user asks to implement, build or execute a spec that has a plan.md and a tasks.md. Run it as the main agent, because it must start sub agents.
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

The user can also give a review mode: `review: wave` (the default),
`review: spec` or `review: off`. See "Code review".

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
- the review mode, and the models of the roles `code-reviewer-wave`,
  `code-reviewer-spec` and `review-verifier`;
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
7. Record the start commit: `git rev-parse HEAD` in the worktree. This is
   the base of the spec review.

## Run loop

Do these steps for each wave, in order:

1. Find the ready tasks: not done, not blocked, and all `Needs` done.
2. Start one sub agent for each ready task, all in one message so that they
   run in parallel. Use at most the configured worker limit and the available runtime slots. Run the rest in a
   next batch of the same wave.
3. Wait for all sub agents to finish.
4. Check each result (see "Check a result").
5. Mark the passed tasks `[x]` in `tasks.md`.
6. Run `cargo check --workspace` in the worktree. If it fails, find the
   task that caused the failure, and send a fix (see "Failure").
7. At the end of a phase, run the command of its checkpoint. Do not start
   the next phase until it passes.
8. Commit the wave (see "Commits").
9. If the review mode is `wave`, review the wave (see "Code review"). Do not
   start the next wave until the review passes or the user accepts the open
   findings.
10. Start the next wave.

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
- For a test task: write the test, run it, and show that it fails for the
  right reason. Do not write the implementation.
- For an implementation task: make the named tests pass. Do not change the
  tests. If a test is wrong, stop and report it.
- Run the "Done when" command. Run `cargo fmt --all` on your files.
- Write comments and docs in ASD-STE100 Simplified Technical English.
- Do not add a feature, a dependency or an abstraction that the task does
  not name.

Final report (short): status (done / blocked / failed), files changed, the
command that you ran with its result, and any problem or open question.
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
2. Run the `Done when` command yourself. A test task must fail for the right
   reason (not a compile error in a file that the task does not own, unless
   the plan says the stubs fail that way).
3. Read the diff of every `[S]` task that passed, and of every task with a
   suspect report. A passed test is not enough: check that the code does
   what the task says and no more.

Because tasks of one wave run at the same time, attribute changes by file.
If two sub agents leave files that you cannot attribute, check each task
alone.

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
- A failure in a shared check (`cargo check`) with no clear owner: run the
  tasks of the last wave one at a time to find it.
- Never skip a failing test, weaken a test, or add `#[allow]` to pass a
  check. Never use `--no-verify`.

## Code review

The reviewer is the `code-reviewer` agent (skill
`.agents/skills/code-reviewer/SKILL.md`). It reads and reports. It never
changes the spec worktree. You decide what to fix, and normal sub agents
make the fixes. Take the models from the provider policy:

- `code-reviewer-wave` (standard): the wave reviews;
- `code-reviewer-spec` (strong): the spec reviews and the second rounds of
  a spec review;
- `review-verifier` (strong): the verification of spec findings.

Modes:

- `wave` (default): review each wave after its commit, and review the spec
  at the end.
- `spec`: review the spec at the end only. This is faster, but fixes come
  late and can be larger.
- `off`: no review. Say this in the final report.

### Wave review

Start one reviewer sub agent with the `code-reviewer-wave` model:

```
Use the code-reviewer skill (.agents/skills/code-reviewer/SKILL.md).
Spec worktree: <ABSOLUTE PATH>
Spec: specs/NNN-name/
Job: wave N
Range: <commit before the wave>..<wave commit>
Tasks in scope: <ids>
```

Check the findings yourself (see "Act on the findings", step 1). A wave
review has no verifier.

### Spec review

Run the spec review in the `Finish` step. Go in this order:

1. Make three review worktrees, one for each focus, from the last commit:
   `git worktree add --detach ../geoviewer-NNN-review-<focus> HEAD`, for
   `correctness`, `tests` and `integration`. The reviewers write scratch
   tests there. These worktrees are throw-away.
2. Start three reviewer sub agents with the `code-reviewer-spec` model, all
   in one message so that they run in parallel. Give each one:

   ```
   Use the code-reviewer skill (.agents/skills/code-reviewer/SKILL.md).
   Spec worktree: <ABSOLUTE PATH>
   Review worktree: <ABSOLUTE PATH of the worktree of this focus>
   Spec: specs/NNN-name/
   Job: spec
   Focus: correctness | tests | integration
   Range: <start commit>..<last commit>
   ```

3. When all three finish, merge the findings. Two findings with the same
   file, line and problem are one finding. Keep the higher severity.
4. Start one verifier sub agent with the `review-verifier` model for each
   focus that has `blocker` or `major` findings, all in one message. Give it
   the review worktree of that focus, `Job: verify`, the same range, and the
   full text of those findings. Do not give it the verdict or the notes of
   the reviewer.
5. Drop a finding only if the verifier rejects it with a reason. Use the
   severity from the verifier. Then act on the remaining findings as below.
   For the second round, make one more review worktree,
   `../geoviewer-NNN-review-recheck`, from the fix commit.
6. When the review is done, remove all review worktrees:
   `git worktree remove --force ../geoviewer-NNN-review-<name>`, then
   `git worktree prune`. They hold only scratch tests, so `--force` is
   correct here. Never use `--force` on the spec worktree.

### Act on the findings

1. Do not trust the report. Read the code at the given line, and confirm
   that the problem is real. For a spec review, the verifier did this
   check. Do the check yourself only if a verifier result is not clear. If a
   finding is not real, drop it and note why.
2. If a finding goes against the spec, the plan or a user answer, do not fix
   it. Ask the user.
3. Send each confirmed `blocker` and `major` finding to a new sub agent for
   the task that owns the file. Use the sub agent prompt, with the task text
   and the finding added below it (`Review finding to fix: ...`). If the
   finding has a scratch test, give it to the sub agent. The sub agent can
   add it as a regression test, but only in a test file in the `Files` of
   the task. The `Files` rule stays. Use the model of the task tier.
   Findings for one task go to one sub agent.
4. Check each fix as in "Check a result". Run `cargo check --workspace`.
   Commit the fixes in one commit: `fix(scope): correct review findings of
   wave N` or `fix(scope): correct spec review findings`.
5. Start a second review of the fix commit only. Use one reviewer and the
   same model as the first round. For a spec review, use `Focus: all` and
   the recheck review worktree. If a
   `blocker` or `major` finding remains after this second review, stop the
   review loop. Mark the task `[!]` with a one-line reason, and ask the user
   before you continue. Do not start a third review round.

Do not fix `minor` findings during the run. Collect them for the final
report.

## Commits

- Commit in the worktree after each wave that passes `cargo check`. You are
  the only agent that commits. Sub agents never commit.
- Use Conventional Commits in ASD-STE100: `type(scope): description`, for
  example `feat(glyphs): add orientation module`. Test tasks use `test`.
  Add a co-author line only if the session supplies one.
- Do not push, and do not open a pull request, unless the user asks.

## Finish

The `Finish` tasks of `tasks.md` run as normal tasks. Then, if the review mode
is not `off`, run the spec review (see "Spec review"). Before you say that the spec is done, run these commands in
the worktree yourself:

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
- the review mode, the number of findings for each review (found,
  rejected by the verifier, dropped by you, fixed, open), the FR and SC
  coverage table of the `tests` reviewer, and the list of open `minor`
  findings;
- the result of the three final commands;
- the next step (for example: open a pull request).

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
  `cargo test` and `cargo clippy` only when a wave or the spec is done.
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
