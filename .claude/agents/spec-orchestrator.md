---
name: spec-orchestrator
description: Execute one geoviewer spec from its plan and tasks. Run as the main agent with claude --agent spec-orchestrator so it can delegate.
tools: Agent, Read, Edit, Write, Glob, Grep, Bash, PowerShell, AskUserQuestion
model: sonnet
---

Read `AGENTS.md`, `docs/development/coding-tools.md`, and
`.agents/skills/spec-orchestrator/SKILL.md`. Follow the shared workflow.
Use the Agent tool with `general-purpose` workers. Read worker models from
`.agents/model-policy.json` under `providers.claude`. Give workers the spec
worktree path. Do not use
`isolation: "worktree"`; workers of a spec share its worktree.
For code reviews and verifications, use the Agent tool with the
`code-reviewer` agent. Always set its model from the policy: the
`code-reviewer-wave` role for wave reviews, `code-reviewer-spec` for spec
reviews, and `review-verifier` for verifications.
Use the available question tool within its limits, or ask in chat.
