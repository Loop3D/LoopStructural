---
name: spec-interview
description: Develop a new or rough spec by interviewing the user, with the architecture model as the anchor. Asks about the need, finds which components own the data, which boundaries the feature crosses and which specs it depends on, then writes spec.md and an arch-impact block and shows the architecture diagram of the spec. Use when the user wants to start a spec, turn an idea into a spec, or link a spec to the architecture. Do not use to review a finished spec (use spec-reviewer), to plan tasks (use spec-planner) or to write code.
---

# Spec interview

Turn an idea into a spec through a conversation. The architecture is part of
the conversation: for each answer, ask which component owns it. This
connects the spec to the system, and it finds design problems early.

Use the `spec-diagrams` skill for the model format, the `arch-impact` block
and the diagrams. Do not write production code.

## Input

The user gives an idea, a spec number, or a path to a rough `spec.md`. If the
user gives nothing, ask for the idea in one sentence. If a spec exists, read
it first and interview only about the gaps.

## Before the first question

Read these, when they exist. Do not ask the user what they already say.

1. Project rules: `AGENTS.md`, `CLAUDE.md`, the constitution
   (`.specify/memory/constitution.md`). They have priority.
2. The architecture model: the `arch` block in `docs/architecture.md`, and
   the text around it. Run the `spec-diagrams` `check` command. If the model
   is missing, tell the user and offer to make it (see that skill).
3. `.specify/templates/spec-template.md` and `specs/README.md`. If there is
   no template, use the structure in "Spec structure" below.
4. The specs that look related: same components, same words.

## Interview

Ask in rounds of up to 4 questions. Use the available user-question tool, or
ask in chat. Put your recommended option first and mark it
"(Recommended)". Give the effect of each option in its description. Do not
ask what the documents already answer. After each round, write the answers
into `spec.md` at once, then ask the next round.

### Round 1. The need

- Who has the problem, and what do they do today?
- What is the smallest result that is still useful? (This is story 1.)
- What must not change?

### Round 2. Place in the architecture

Before this round, show the user a draft: name the components that you think
the feature touches, from the model, and say why. Then ask:

- Which component owns the data or the decision? Only one may own it.
- Does the feature cross a boundary (a process, a protocol, a file format, a
  public API)? Which one?
- Is a new component needed, or does an existing one grow? Offer both, and
  say what each costs. Prefer the existing component unless it breaks the
  rules of the project.
- Which existing specs does this one depend on? Which specs does it change?
- Does the change break a rule of the constitution or the model? Name the
  rule.

### Round 3. Behaviour

- The main scenario, step by step. Each step is a `flow` line: who sends
  what to whom.
- Acceptance scenarios (Given, When, Then) that a test can check with no
  human.
- Edge cases: invalid input, a failure in the middle, a lost connection, big
  data, the old data that exists today.
- Numbers for each success criterion (a value and a unit).

### Round 4. Boundaries and limits

- What messages, fields or file changes appear at each boundary? Is the
  change compatible with old clients and old files? Does a version number
  change?
- What is out of scope? Name the things a reader can expect.
- What is the risk, and what is the way back if it fails?

Stop when each section of the spec has content and no open question is
left. After 4 rounds, ask the user if they want to continue. Mark each
question that the user cannot answer now as `[NEEDS CLARIFICATION: ...]`.
Do not guess.

## Write the spec

Use the project template. Keep the spec about WHAT and WHY. Components and
messages in plain words are fine, but types, crates, classes and files go in
the plan. Write the user text in ASD-STE100 Simplified Technical English if
the project uses it.

Add the `arch-impact` block at the end of the "Requirements" section, or in
the place the template gives. Build it from the answers of round 2 and
round 3:

- `depends`, `changes` and `adds` from round 2.
- `flow` from the main scenario of round 3.

Then:

1. Run `check`. Fix each unknown id by asking the user, not by guessing.
2. Run the `spec NNN` command and show the diagram to the user. Ask if it
   shows what they mean. A wrong diagram is the cheapest place to find a
   wrong design.
3. If the spec adds a component, ask the user if you can add it to the model
   in `docs/architecture.md` as a proposal. Do not change the model without
   an answer.
4. Offer to put the diagram into the spec with a marker pair
   (`<!-- mermaid:spec-NNN -->`) and run `inject`.
5. Add the spec to the index in `specs/README.md` if the project has one.

## Spec structure (when there is no template)

Header (id, date, status, input), user stories with priority and an
independent test, acceptance scenarios, edge cases, functional requirements
(`FR-NNN`, one testable idea each), success criteria (`SC-NNN`, with
numbers), out of scope, clarifications (question, answer, date), and the
`arch-impact` block.

## Finish

- Give a short summary: sections written, open `[NEEDS CLARIFICATION]`
  items, components touched, specs it depends on.
- Say what to do next: run the spec-reviewer skill, then plan.
- Commit only if the user asks, and follow the commit rules of the project.
  Use the worktree or branch for the spec. Do not switch the branch of the
  main checkout.
