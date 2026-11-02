---
name: spec-diagrams
description: Make Mermaid diagrams from the architecture model and from the arch-impact blocks of specs. Draws the whole architecture, the effect of one spec (changed and added components), the flow of one spec as a sequence diagram, and the dependency graph of all specs. Checks that every spec names only components that exist in the model. Keeps diagram regions in Markdown files current. Use when the user asks for an architecture diagram, a spec diagram, a dependency graph, or to link specs to the architecture. Do not use to write specs (use spec-interview) or code.
---

# Spec diagrams

Connect specs to the architecture with two small text formats and one script.
The script uses only the Python standard library (Python 3.10 or later).
Do not draw Mermaid by hand when the script can make it. Hand-drawn diagrams
go out of date.

## The two formats

**The model.** One fenced `arch` block, in `docs/architecture.md` by
default. It is the one list of components. Lines start with a word:

````
```arch
component ui "App UI" kind=app
component core "Core library"
component scene "Scene model" kind=leaf
external qgis "QGIS plugin"
store disk "Scene file"
group ws "Workspace" : ui core scene
edge ui -> core "uses"
edge qgis <-> core "live link"
```
````

- `component`, `external` and `store` need an id and a label. Ids use
  letters, digits and `_`. `kind=` only changes the colour.
- `edge A -> B "label"` is one way. `A <-> B` is two way. The label is
  optional.
- `group` draws a box around the listed ids.
- A line that starts with `#` is a comment.

**The spec impact.** One fenced `arch-impact` block in each `spec.md`. It
says how the spec touches the model:

````
```arch-impact
depends 002 004
changes core scene
adds component panel "Settings panel" kind=app
adds edge panel -> core "calls"
flow qgis -> core "SetViewSettings"
flow core -> ui "apply settings"
```
````

- `depends` lists spec numbers that must be done first.
- `changes` lists existing component ids that the spec changes.
- `adds` gives a new component or a new edge. A new component is not in the
  model until the spec is built. Then move the line into the model.
- `flow` lines are the steps of the main scenario, in order. Use ids from the
  model or from `adds`.

## Commands

Run the script from the repository root. It is in the `scripts` folder next
to this file.

```sh
python -I <skill-dir>/scripts/spec_diagrams.py model          # whole architecture
python -I <skill-dir>/scripts/spec_diagrams.py spec 32        # impact of spec 032
python -I <skill-dir>/scripts/spec_diagrams.py spec 32 --full # same, all components
python -I <skill-dir>/scripts/spec_diagrams.py flow 32        # sequence diagram
python -I <skill-dir>/scripts/spec_diagrams.py deps           # spec dependency graph
python -I <skill-dir>/scripts/spec_diagrams.py check          # validate all specs
python -I <skill-dir>/scripts/spec_diagrams.py check --strict # missing block is an error
python -I <skill-dir>/scripts/spec_diagrams.py inject docs/architecture.md
python -I <skill-dir>/scripts/spec_diagrams.py inject --check docs/architecture.md
```

Options `--root`, `--model` and `--specs` change the repository root, the
model file and the specs folder. Defaults: the current folder,
`docs/architecture.md` and `specs`.

The `spec` diagram marks changed components yellow and new components green
with a dashed edge. By default it shows only the changed and new components
and the ends of the new edges. Use `--full` to show all components.

## Keep diagrams in documents

Put a marker pair in any Markdown file. `inject` fills it with a Mermaid
block. Names: `model`, `deps`, `spec-NNN`, `spec-NNN-full`, `flow-NNN`.

```
<!-- mermaid:spec-032 -->
<!-- /mermaid:spec-032 -->
```

Run `inject` after each change to the model or to an `arch-impact` block.
Run `inject --check` in CI. It exits with 1 when a region is stale.

## Steps

1. Run `check`. Fix each error before you draw. An error names the spec and
   the unknown id. Warnings about a missing block are not errors, unless the
   user asks for `--strict`.
2. Run the command for the diagram that the user wants. Show the Mermaid
   code in a `mermaid` fenced block.
3. If the user wants the diagram in a file, add the marker pair and run
   `inject`. Do not paste the generated code by hand.
4. If the model has no `arch` block, tell the user. Offer to write it from
   the real structure: read the workspace manifest, the module list and the
   existing architecture document. Ask the user to confirm the component
   list before you write it. Do not guess components.

## Rules

- The model is the source of truth. Change the model, then draw. Never edit
  the generated Mermaid.
- Keep ids short and stable. A renamed id breaks every spec that uses it;
  `check` finds those.
- Do not put more than about 15 components in a diagram. Use `group`, or
  make a second model file, when the diagram is too big.
- Run the tests when you change the script:
  `python -I -m unittest discover <skill-dir>/scripts`.
