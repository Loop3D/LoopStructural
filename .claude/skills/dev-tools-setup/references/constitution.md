# Adopt a project constitution

The constitution records the project's governing development rules. It belongs
to the project and is never a managed toolkit template. Respect the repository's
own governance and amendment procedure when one exists.

## Find the evidence

Read existing agent instructions, CONTRIBUTING documents, development docs,
build manifests, CI workflows, and the relevant architecture documentation.
Distinguish observed conventions from recommendations. Cite repository paths
next to the proposed rules they support. CI configuration can establish a test
command; it does not establish a new architectural or review policy by itself.

Reuse a constitution already stored elsewhere by setting
`documentation.constitution` to its project-relative path. Preserve its content.
If instructions or documents conflict, record the conflict and ask about that
specific choice; do not silently pick a new rule.

## Draft before requesting review

Use the setup interview to settle governance decisions that repository evidence
does not answer. Carry accepted decisions from `.dev-tools/setup-notes.md` into
the draft and retain deferred questions as unresolved. Answers to individual
interview questions do not replace the final review of the complete constitution.

Write a sibling draft alongside the intended constitution, for example
`.specify/memory/constitution.draft.md`. If a draft already exists, read it and
preserve user edits. Do not overwrite an existing constitution. If the intended
file is already marked as a draft, revise that draft in place.

Include this marker until the user has reviewed and accepted the document:

```html
<!-- dev-tools:constitution-draft -->
```

Use a concise structure suited to the evidence:

- Project purpose and scope, derived from existing project documentation.
- Supported development commands and validation requirements.
- Established architectural boundaries, compatibility rules, and ownership.
- Existing practices for changes, reviews, and documentation.
- How future amendments are proposed and accepted.
- Evidence, proposed additions, and unresolved decisions.

Omit unsupported sections or mark their specific decisions as unresolved. Do not
invent test coverage thresholds, model choices, compliance obligations, mandatory
review counts, or a particular technology stack. On a sparse project, a short
draft and a few explicit questions are sufficient.

Show the draft path and summarize the rules derived from evidence, proposed
additions, and open decisions. Request review of that concrete draft. Review is
needed because this document will govern later agent actions; stop activation
while the review is pending. If the user defers, keep the draft and report the
constitution prerequisite as incomplete. Other independent setup may continue.

## Activate the reviewed document

After the user accepts the final text, write it to the configured constitution
path, remove the draft marker, and record the agreed decisions. Confirm that no
other process created a constitution while review was pending; if it did, read
and reconcile it before writing. Preserve unrelated files and do not automatically
commit. Remove only the superseded draft created by this workflow once its text
is safely represented in the final document.

Run `dev-tools doctor --require-constitution` from the project root. This checks
that the file is nonempty and lacks the draft marker. It cannot establish that
the rules are appropriate or prove that a human reviewed them; the conversation
provides that review. Do not remove a draft marker merely to make the check pass.
