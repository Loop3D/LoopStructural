---
name: dev-tools-setup
description: Configure the LoopConnect development toolkit in an existing project, including project commands, documentation paths, and a reviewed constitution. Use when adopting this toolkit or completing its initial setup.
---

# Set up project workflows

Read `.dev-tools/project.json` and the project's existing agent instructions.
If the toolkit is not initialized, run `dev-tools init --dry-run` against the
target project, inspect the proposed paths, then apply the authorized setup.
If the CLI is unavailable, explain how to install it from the toolkit repository;
do not guess a package index release or recreate its managed manifest manually.

Inspect build manifests, CI workflows, test configuration, and development docs
before asking for settings. Populate `commands.check`, `commands.test`,
`commands.format`, and `commands.lint` with commands supported by that evidence.
Leave unavailable commands empty. Commands are strings run from the project root;
the installer and doctor never execute them automatically.

Set `documentation.architecture` and `documentation.specs` to existing paths when
appropriate. Ask about unresolved choices that affect the workflow. Do not invent
an architecture or specification history during installation.
Preserve existing JSON fields and project instructions, and avoid putting secrets
or machine-specific absolute paths into committed configuration.

## Constitution

Resolve `documentation.constitution` relative to the project root, defaulting to
`.specify/memory/constitution.md` for configurations without that setting.
Read an existing constitution and preserve it; initialization and updates do not
own this document. If the project already keeps its constitution elsewhere,
record that path instead of making a competing document.

If the constitution is missing, empty, or marked as a toolkit draft, follow
[the constitution workflow](references/constitution.md). Draft concrete text
from project evidence before asking for review. Do not present proposed rules
as established project policy or activate a draft without the user's review.
An existing nonempty constitution needs no fresh approval just to use it.

## Verify

Run `dev-tools doctor` against the project after editing configuration. Report
which settings were discovered, which are unset, and any verification limits.
For workflows that require a constitution, run `dev-tools doctor
--require-constitution` from the project root and resolve failures before
continuing. This checks installed files and the constitution prerequisite, not
whether the project builds or tests successfully or whether its policies are sound.

The configuration's `agents` selects installation targets, not available runtime
tools or models. Use only capabilities actually exposed by the active client.
MCP registration and authentication are not implemented by this setup workflow.
