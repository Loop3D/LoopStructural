# Setup decisions

- 2026-11-01: The user requested planning and implementation workflows for
  spec 003. Preserve the existing Python project commands and compatibility
  boundaries in `.dev-tools/project.json` and `ROADMAP.md`.
- 2026-11-01: The user approved `.specify/memory/constitution.md` after review
  of its rules and the corrected CI description. The active merge rule requires
  pytest and `ruff check LoopStructural packages` to pass. No CI change was made.
- The installed `dev-tools.exe` fails with a uv trampoline path error. The local
  source checkout at `C:/Users/lachl/dev/dev-tools` provides a working CLI via
  `python -m dev_tools` with that checkout's `src` on PYTHONPATH. This is a local
  execution workaround; the global launcher has not been repaired.
- The four managed setup skill files differed from their recorded hashes only
  by CRLF line endings. Their LF-normalized hashes were verified against the
  manifest before normalization and a reviewed toolkit update was applied.
- The source CLI strict doctor passed after constitution activation and update.
- Planning templates, `specs/README.md`, and
  `docs/development/coding-tools.md` are absent. Spec 003 uses existing 001/004
  document structure and reports model capabilities at implementation time.
- The spec 003 plan and tasks are in the existing `LoopStructural-native`
  worktree. That worktree must receive the approved setup and pass strict doctor
  before production work starts. Gradient count/memory interpretation and
  measured benchmark acceptance remain explicit gates in its plan.
