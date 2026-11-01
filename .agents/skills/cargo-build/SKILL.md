---
name: cargo-build
description: Build, check, test and lint geoviewer with cargo in the fastest safe way. Gives the one command form to use, so that the shared target folder is not rebuilt, and tells what to do when a build is slow or the disk is full. Use before any cargo check, build, test or clippy run, when a build is slow, or when the user asks how to build. Do not use to change specs or to plan code.
---

# Cargo build

Read `AGENTS.md` for the rules on `scripts/cargo-lock.sh`. This skill adds the
rules that keep builds fast.

## Why builds are slow

- All worktrees share one target folder (set by `scripts/setup-local-build`). Cargo picks the
  Bevy features from the packages in the command. A command with `-p` and a
  command with `--workspace` use different features. Each one makes the other
  rebuild Bevy (`bevy_pbr`, `bevy_render`, `wgpu`). This costs many minutes.
- Dev dependencies build with `opt-level = 3`. A full Bevy build is slow by
  design. Do not make it happen again without a reason.
- The lock allows one cargo command at one time. A queue is normal.

## Rules

1. Always use the lock script. Never run `cargo` directly.
   - Bash (macOS, Linux, Git Bash): `scripts/cargo-lock.sh <args>`
   - PowerShell: `scripts\cargo-lock.ps1 <args>`
   Both run `scripts/cargo_lock.py` and use the same lock. If the shared
   target folder is not set, run `scripts/setup-local-build.sh` (or `.ps1`).
2. Always add `--workspace` to `check`, `build`, `test` and `clippy`. Do not
   use `-p`, `--package` or `--manifest-path`. To run one test, filter by
   name after the workspace flag:
   `scripts/cargo-lock.sh test --workspace view_tone`
   The same rule holds for `--test NAME` and `--lib`.
3. Use `check` for the inner loop. Use `build` only when you need the binary.
   Use `test` for the final check.
4. Do not change features, profiles or `RUSTFLAGS` for one command.
   Only a task that says so may use `--features`, `--no-default-features` or
   `--target` (for example the wasm checks of spec 033). Such a command
   rebuilds Bevy once. Run it as few times as you can.
5. `cargo run` does not take `--workspace`. Run
   `scripts/cargo-lock.sh build --workspace`, then start the binary from
   the `debug` folder of the shared target folder.
6. For `-- --ignored`, add a test name. A bare `--workspace -- --ignored`
   runs the ignored tests of all crates.
7. Never set `CARGO_TARGET_DIR` and never use `--target-dir`. A private
   target folder (for example `cargo-target-032`) builds Bevy a second time
   and uses 40 GB or more. If the shared folder gives stale results, run
   `test --workspace` and trust a failure over a pass.
8. Do not run `cargo clean` or `cargo update` unless the task says so.
9. Before the final report, run these three, in this order:
   - `scripts/cargo-lock.sh fmt --all --check`
   - `scripts/cargo-lock.sh clippy --workspace --all-targets -- -D warnings`
   - `scripts/cargo-lock.sh test --workspace`

## When the build is slow

1. Read the first line. If it says `another cargo command runs`, wait. Do not
   start a second command.
2. If you see `Compiling bevy_*` after a small edit, a command used a
   different form. Find it, and use the form in rule 2.
3. The lock script stops with exit code 4 when the target drive has less
   than 5 GB free (`GEOVIEWER_MIN_FREE_GB` sets the limit, 0 turns it off).
   If you see this, run the `worktree-cleanup` skill. Do not turn the check
   off. Parallel agents each add build output, so do not start a new
   worktree when the drive is low.
4. For a timing, add `--timings` to a `build --workspace` command. The report
   is in the `cargo-timings` folder of the shared target folder.

## Report

Write the report in ASD-STE100 Simplified Technical English. Give the
commands that you ran and their result. If a check fails, give the output.
