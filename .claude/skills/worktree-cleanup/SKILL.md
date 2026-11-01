---
name: worktree-cleanup
description: Clean up the git worktrees and the build artifacts of geoviewer to free disk space. Lists all worktrees, finds those that are safe to remove (work merged or pushed, no uncommitted changes), shows the list, asks the user, and then removes the worktrees, their merged branches and the cargo target directories. Use when the user asks to clean up, free disk space or memory, remove old worktrees, or clear build artifacts, and at the end of a spec-implementer or spec-orchestrator run. Do not use to change specs or code.
---

# Worktree cleanup

Read `docs/development/coding-tools.md` for tool setup and model routing.

Remove worktrees and build artifacts that are no longer needed. Never lose
work. Do nothing that you cannot undo until the user agrees.

## Rules

- Never remove the main checkout. Never switch its branch.
- Never remove a worktree that has uncommitted changes, untracked files that
  are not build output, or commits that are not pushed or merged.
- Never use `--force` on `git worktree remove` or `git branch -D` without
  the explicit consent of the user for that one item.
- Write the report in ASD-STE100 Simplified Technical English.

## Steps

1. **List.** Run `git worktree list` and `git fetch --prune`. For each
   worktree other than the main checkout, collect:
   - the path and the branch;
   - `git -C <path> status --short` (uncommitted changes);
   - the commits that are not on `origin/main`
     (`git log origin/main..<branch> --oneline`) and if the branch is pushed
     (`git branch -r --contains <branch>`);
   - if the pull request is merged (`gh pr list --head <branch> --state
     merged`, if `gh` is available);
   - the size of its `target/` directory, if it has one.
2. **Classify.** Put each worktree in one group:
   - **Safe**: clean, and the branch is merged into `origin/main` or its
     pull request is merged.
   - **Pushed**: clean, not merged, but pushed. The work is safe on the
     remote. Remove only if the user agrees.
   - **Keep**: uncommitted changes, or commits that exist only here.
     Tell the user. Do not remove.
   - Also mark stale `tmp/*` branches and prunable worktrees (the folder is
     gone: `git worktree prune --dry-run`).
3. **Show and ask.** Show a table of the groups with the size that each
   removal frees. Ask the user which items to remove. Do not remove
   anything before the answer.
4. **Remove.** For each item that the user approved:
   - `git worktree remove <path>`;
   - `git branch -d <branch>` for merged branches (use `-D` only with
     explicit consent);
   - then `git worktree prune`.
5. **Build artifacts.** Find the large build directories:
   - `target/` in the main checkout and in the worktrees that stay;
   - the shared target directory `../geoviewer-target`, if it exists.
   Show the sizes. With the consent of the user, run `cargo clean` in each
   one (or delete the shared directory). Prefer to remove `target/debug/incremental`
   and `target/debug/deps` first if the user wants to keep a warm build.
6. **Report.** Show the worktrees removed, the worktrees kept (with the
   reason), the branches deleted, the disk space freed, and the output of
   `git worktree list` after the cleanup.
