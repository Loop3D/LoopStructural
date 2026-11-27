# Spec 007 handoff: native fault application

Instructions for the person or agent who implements spec 007.

## Where to work

- Worktree: `/home/lgrose/dev/LoopStructural-wt/feat-007-native-fault-kernel`
- Branch: `feat/007-native-fault-kernel`, from master `6874364c` (spec 003 merged).
- Python environment: there is no venv in this worktree yet. Make one, or use
  `../feat-003-fault-evaluation/.venv` (Python 3.13, maturin, pytest-xdist)
  with `PYTHONPATH=$PWD:$PWD/packages/loop_common/src:$PWD/packages/loop_interpolation/src`.

## What to read, in this order

1. `.specify/memory/constitution.md`: the project rules. They win over everything else.
2. `specs/007-native-fault-kernel/spec.md`: the requirements (FR-001 to FR-015), the
   success criteria (SC-001 to SC-008) and the decisions (C-001 to C-011).
3. `specs/003-fault-evaluation-performance/spec.md`: FR-004 (the golden tolerance
   that the native path must meet) and FR-010 (the memory rule).
4. `LoopStructural/utils/_fault_block.py` (`restore_through`, `cutters_first`,
   `with_cutters`) and `LoopStructural/modelling/features/fault/_fault_segment.py`
   (`apply_to_points`, `_evaluate_frame`): the Python reference.
5. The proof of concept, outside the repository:
   `/home/lgrose/dev/LoopStructural-wt/_runs003/rust_poc/README.md`. It is a
   starting point, not code to copy as it is: it reads exported files (the spec
   asks for the in-memory model), and its results differ by up to 1.9e-4 m (the
   spec asks for 1e-8 times the box size; finding the cause is part of story 2, C-007).
6. `specs/007-native-fault-kernel/research.md`: background only. Precomputed maps are out of scope.

## Order of work

1. Run the spec-planner skill on spec 007, to write `plan.md`. Then run spec-tasks,
   to write `tasks.md`. Then build with spec-orchestrator, wave by wave.
2. Story 1 (Python only: compute the restore order once for each feature) comes
   first. It needs no native code, and it must give bit-identical results.
3. Stories 2 to 5 follow, in the order of the spec. The native code goes in
   `packages/loop_interpolation_core` (it already builds with maturin, PyO3 0.29,
   numpy 0.29 and rayon), unless the plan records another decision.

## How to measure

- Golden tests and the standard suite: `pytest tests/fault_evaluation` and
  `python -m benchmarks.fault_evaluation.run --suite standard` (see spec 003
  `measurements.md` for the full T030 commands).
- The 94-fault model (SC-001, SC-004):
  `python -m benchmarks.fault_evaluation.many_faults --output RESULT.json`.
  It builds for about 7 minutes. It prints the `evaluate_model` time and the mean
  time of an isolated and a cutting fault at 500, 2000 and 5000 points.
  - Before a change: add `--save-reference REF.npz`. After it: add `--compare REF.npz`.
    The script exits with 1 if any model id or group value differs (bit for bit).
  - The reference of master `6874364c` (Python path, 2026-11-27) is in
    `/home/lgrose/dev/LoopStructural-wt/_runs007/`: `many_faults_python.json` and
    `many_faults_python_ref.npz`. Build 408 s; `evaluate_model` 1.78 s / 1.89 s /
    2.09 s at 500 / 2000 / 5000 points (94 fault calls); a cutting fault 40-47 ms
    (mean), an isolated fault 2.6-3.5 ms. These are the SC-001 and SC-004 baselines.
- Measure on the T030 computer (Intel i7-12700H, Linux), with no other work running.
  Use `OMP_NUM_THREADS=1` for test runs with `-n`.

## Rules

- The public `GeologicalModel` API does not change (FR-011).
- The Python path stays the reference and the fallback. The full test suite passes
  with the native path on, with it off, and without the native package (SC-006).
- Before a merge: `ruff check LoopStructural packages` and `pytest` pass.
- Ask the maintainer before you open a pull request or push with force.
