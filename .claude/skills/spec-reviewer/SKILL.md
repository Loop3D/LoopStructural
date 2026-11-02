---
name: spec-reviewer
description: Review a geoviewer spec (specs/NNN-name/spec.md) against the constitution and the spec template, ask the user if the edits can be made, then interview the user about the new spec and update it from the answers. Use when the user asks to review, check or audit a spec, before the clarify step and before the plan step. Do not use for plans, tasks or code.
---

# Spec reviewer

Read `docs/development/coding-tools.md` for role and model selection.

Review one spec. Report problems. Do not change the spec until the user
agrees. After the report, ask the user if you can make the edits. If the user
agrees, edit the spec. Then ask the user questions about the new spec, one
round at a time, and update the spec from the answers.

## Input

The user gives a spec number or a path (for example `002` or
`specs/002-protocol-and-schema/spec.md`). If the user gives nothing, review
the spec that changed last in `git status`. If you cannot find one, ask.

## Steps

1. Read `.specify/memory/constitution.md`. It has priority over all other
   documents.
2. Read `.specify/templates/spec-template.md`.
3. Read the spec. Read `specs/README.md` for the index.
4. Read the other specs that the spec names or that touch the same objects or
   messages. Look for conflicts.
5. Do the checks below. Collect findings.
6. Write the report in the format below.

## Checks

### A. Template and structure

- All sections of the template are present, in order: user scenarios,
  requirements, success criteria, out of scope, clarifications.
- The header has the feature id (`NNN-short-name`), a date and a status.
- The spec is in `specs/NNN-name/`, the number is not used twice, and the
  index in `specs/README.md` lists it.
- Each user story has a priority, a "Why this priority", an "Independent
  test" and acceptance scenarios in Given/When/Then form.
- The stories are in priority order, and each story can be tested alone.

### B. WHAT, not HOW

- The spec does not name crates, Rust types, libraries or files. Those go in
  `plan.md`. Exception: the constitution allows a name when it is part of the
  contract (for example `geoh5` names).
- Each requirement says what the viewer must do, in the form "The viewer
  must ...". Each requirement has one testable idea. Flag "and" lists that
  hide two requirements.
- No vague words without a number: "fast", "large", "easy", "robust",
  "appropriate", "etc.".

### C. Testability (Principle V)

- Each requirement (FR-NNN) is covered by at least one acceptance scenario or
  success criterion. Make a coverage list.
- Each acceptance scenario can become an automatic test with no GPU and no
  human.
- Each success criterion has a number and a unit. It is measurable.
- Edge cases include invalid input, a client that disconnects, and large
  data.

### D. Constitution (one line for each principle that applies)

- **I. One owner.** Each object has one owner. Edits to client-owned objects
  go to the owner. The spec mentions revisions, undo and redo, derived
  objects and push where they apply.
- **II. Protocol first.** New client features start with the messages in the
  spec. The spec says if the protocol version changes (major or minor). No
  QGIS or LoopStructural names or concepts in the protocol.
- **III. Object model.** Each new object type or data type gives its `geoh5`
  mapping, or says "no mapping" and how the export handles it. Names from
  `geoh5` have the same meaning as in `geoh5`.
- **IV. Exact coordinates.** The spec keeps `f64` world coordinates for all
  values that a user or a client sees.
- **VI. Never lose work.** The spec says what happens on an invalid message,
  a disconnect, a long operation (progress, cancel) and a crash during save.
- **VII. Simple.** The spec has no feature that no user story needs.
- **VIII. Plugins.** The spec has no modelling algorithm in the viewer and no
  code for one specific plugin.
- **IX. Shaders.** A spec for a new drawn object states the shader approach
  and the instance count in its performance target.
- **Platforms and performance.** A spec with large data gives a performance
  target against the reference scene.
- **Security.** The spec does not weaken the loopback and token rules, and it
  does not run code from a message or a scene file.

If a principle does not apply, write "not applicable". Do not skip it
silently.

### E. Clarity

- List each `[NEEDS CLARIFICATION]` item. The clarify step must resolve all
  of them before the plan.
- Find hidden open questions: words such as "maybe", "probably", "should
  decide", "TBD", or a requirement with two possible readings. Propose a
  `[NEEDS CLARIFICATION: ...]` mark for each.
- Terms are used in one way. Entities are defined once and used with the same
  name.
- "Out of scope" lists the things that a reader can expect, but that the spec
  does not do.

### F. Consistency with other documents

- The spec agrees with `docs/viewer-plan.md`, or the difference is on
  purpose. The constitution wins if they differ.
- The spec does not conflict with other specs (same message, same object,
  different meaning).

### G. Language

- The spec uses ASD-STE100 Simplified Technical English: short sentences
  (up to 25 words; 20 words for a procedure), one idea for each sentence,
  active voice, approved words, and the same word for the same thing. Flag
  the worst examples only. Do not rewrite the whole spec.

## Severity

- **Blocker**: breaks a constitution rule, a requirement that no test can
  check, or a conflict with another spec. The plan must not start.
- **Major**: a missing section, a missing `geoh5` mapping, a requirement with
  no scenario, an implementation detail in the spec.
- **Minor**: wording, naming, style, small gaps.

## Report format

Write the report in the chat in ASD-STE100 Simplified Technical English.
Use this structure:

```
# Review of NNN-name

**Verdict**: Ready for plan | Ready after clarify | Needs changes

## Blockers
1. [Section or FR-NNN] Problem. Why it matters (principle or rule). Fix.

## Major
...

## Minor
...

## Requirement coverage
| Requirement | Scenario or criterion | Gap |
|---|---|---|

## Open clarifications
- Item, with the question to ask.

## Constitution check
| Principle | Result | Note |
|---|---|---|
```

Rules for the report:

- Give each finding a location (section name or FR/SC number) and a
  specific fix. Do not write "improve this".
- Put only real problems in the findings. If a section is good, do not list
  it. Use "None" for an empty severity group.
- Verdict rules: any blocker gives "Needs changes". No blocker and one or
  more open clarifications gives "Ready after clarify". Otherwise give
  "Ready for plan".
- Start with the verdict. Keep the report short enough to read in one
  minute, apart from the two tables.

## After the report

Do the two phases below in order. Do not skip a phase.

### Phase 1. Apply the review fixes

Ask the user if they want you to apply the fixes from the report. Use the
available user-question tool, or ask in chat. Do not change the spec before the user answers.

If the user says no, stop. If the user says yes:

- Change only `spec.md`. Do not write plan or code.
- Apply only the fixes that need no decision from the user. Put a
  `[NEEDS CLARIFICATION: ...]` mark where a fix needs a decision.
- Do not record any answer yet. Phase 2 does that.

### Phase 2. Interview about the new spec

Read the edited `spec.md` again. Then ask the user questions about it, to
find gaps and to resolve each `[NEEDS CLARIFICATION]` item.

- Use the available user-question tool (or chat). Ask a small batch of questions within the tool limits.
- Start with the `[NEEDS CLARIFICATION]` items. Then ask about the hidden
  open questions, missing edge cases, unclear acceptance scenarios and
  success criteria with no number.
- Make each question specific. Give a few answer options within the tool limits. Put your
  recommended option first and add "(Recommended)" to its label. Explain the
  effect of each option in its description.
- Ask about WHAT the viewer does. Do not ask about HOW (crates, types,
  files). Those questions belong to the plan.
- Do not ask a question that the constitution, the spec or another spec
  already answers. Do not ask the user to approve the interview itself.
- After each answer, update `spec.md` at once:
  - Replace the matching `[NEEDS CLARIFICATION]` mark with the decision.
  - Change the requirements, scenarios, success criteria, edge cases and
    out-of-scope list that the answer affects.
  - Add the question and the answer to the "Clarifications" section, with
    the date.
- If an answer causes a new open question, ask about it in the next call.
- Stop when no `[NEEDS CLARIFICATION]` item is left and the user has no
  more answers to give. Do not ask more than 4 rounds without a check. After
  4 rounds, ask the user if they want to continue.

### Finish

- Show a short summary: the sections that changed and the new verdict (use
  the verdict rules above).
- Follow the commit rules in `AGENTS.md` (Conventional Commits, for example
  `docs(spec): fix review findings in 002-protocol-and-schema`). Commit only
  if the user asks. Work in the correct worktree for this spec. Do not
  switch the branch of the main checkout.

## Clean up and disk use

- Do not make a new worktree if you can use an existing one. If you make a
  worktree for this work, tell the user which one, and remove it
  (`git worktree remove <path>`, then `git worktree prune`) only when all
  work is committed and the user agrees. Never remove a worktree with
  uncommitted or unpushed work.
- Do not run `cargo build` or other builds. This skill needs none. If you
  must run a `cargo` command, use `cargo check` and set `CARGO_TARGET_DIR`
  to `../geoviewer-target`, so that no extra `target/` directory is made.
- Put temporary files in the scratchpad and delete them when done.
