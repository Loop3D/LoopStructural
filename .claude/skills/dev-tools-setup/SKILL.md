---
name: dev-tools-setup
description: Interview the user to configure the LoopConnect development toolkit in an existing project, including development commands, workflow decisions, documentation paths, and a reviewed constitution. Use when adopting this toolkit or completing its initial setup.
---

# Set up project workflows

Setup is a conversation with the user. Inspect the repository first, then
interview them about intent, unresolved choices, and proposed policy changes.
Do not finish setup using inferred preferences alone. These questions belong in
the coding agent conversation; the CLI only installs and checks files.

Read `.dev-tools/project.json`, the project's existing agent instructions, and
`.dev-tools/setup-notes.md` if it exists. Resume recorded decisions rather than
restarting an interview that is already in progress.
If the toolkit is not initialized, run `dev-tools init --dry-run` against the
target project, inspect the proposed paths, then apply the authorized setup.
If the CLI is unavailable, explain how to install it from the toolkit repository;
do not guess a package index release or recreate its managed manifest manually.

Inspect build manifests, CI workflows, test configuration, and development docs
before asking for settings. Populate `commands.check`, `commands.test`,
`commands.format`, and `commands.lint` with commands supported by that evidence.
Leave unavailable commands empty. Commands are strings run from the project root;
the installer and doctor never execute them automatically.

## Interview

Start by showing a brief inventory: conventions discovered with source paths,
proposed settings, and choices still open. Distinguish established project rules
from recommendations. Use the available user-question tool when supported by
the current client; otherwise ask in chat and wait for a reply.

Ask small rounds of up to three related questions. Offer two or three practical
choices when useful, put the recommended choice first, explain its tradeoff, and
allow a free-text answer. A suggested or preselected option is not an answer.
Do not ask the user to repeat facts or decisions already established by the
repository or conversation. Adapt these topics to the gaps rather than reading
them as a mandatory questionnaire:

1. **Adoption goals:** Which workflows should setup prepare for, what must be
   preserved, and are there constraints the repository does not document? Confirm
   the intended scope if the user has not already given it. Explain which skills
   are actually installed; do not promise the unmigrated spec workflows.
2. **Development decisions:** Resolve competing test/check commands, fast local
   checks versus full CI validation, and missing documentation locations. Reuse
   supported commands; ask before establishing a new requirement.
3. **Governance:** Resolve missing or conflicting compatibility, architecture,
   review, and amendment rules needed for the constitution. Existing rules are
   evidence, not an invitation to reopen every settled decision.

After each reply, record accepted decisions, their evidence, and open questions
in the project-owned `.dev-tools/setup-notes.md`, preserving existing notes.
Apply answered settings to `.dev-tools/project.json` while preserving other
fields. Keep proposed rules and pending decisions out of active policy.

If the user defers a decision, record it as unresolved rather than selecting a
default on their behalf. Continue independent setup and identify the affected
workflows as incomplete. Never treat silence as consent. Once the necessary
choices are answered, show the resulting settings and any remaining gaps, then
complete the constitution review below if one is needed.

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
