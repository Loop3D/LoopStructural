#!/usr/bin/env python3
"""Make Mermaid diagrams from the architecture model and from specs.

The architecture model is a fenced `arch` block in a Markdown file.
Each spec can have a fenced `arch-impact` block in its spec.md.
Only the Python standard library is used (Python 3.10 or later).

Model lines (`arch` block):
    component ID "Label" [kind=KIND]
    external  ID "Label"
    store     ID "Label"
    group     ID "Label" : ID ID ...
    edge A -> B ["label"]        (use <-> for a two-way edge)

Impact lines (`arch-impact` block):
    depends NNN NNN ...
    changes ID ID ...
    adds component|external|store ID "Label"
    adds edge A -> B ["label"]
    flow A -> B "message"
"""
from __future__ import annotations

import argparse
import re
import shlex
import sys
from dataclasses import dataclass, field
from pathlib import Path

FENCE = re.compile(r"^```(arch-impact|arch)[ \t]*\n(.*?)^```[ \t]*$", re.S | re.M)
MARKER = re.compile(
    r"(<!-- mermaid:(?P<name>\S+) -->\n)(?P<body>.*?)(<!-- /mermaid:(?P=name) -->)", re.S
)
NODE_KEYWORDS = ("component", "external", "store")
PALETTE = {
    "component": "fill:#f2f2f2,stroke:#666,color:#111",
    "leaf": "fill:#e8f1ff,stroke:#3b6fb6,color:#111",
    "app": "fill:#e9f7ec,stroke:#2e7d4f,color:#111",
    "external": "fill:#fff4e0,stroke:#b9770e,color:#111",
    "store": "fill:#fff4e0,stroke:#b9770e,color:#111",
    "changed": "fill:#fff2b3,stroke:#b38f00,stroke-width:3px,color:#111",
    "added": "fill:#d6f5d6,stroke:#2e7d4f,stroke-width:3px,stroke-dasharray:4 3,color:#111",
}


class ArchError(Exception):
    """A problem in an `arch` or `arch-impact` block."""


@dataclass
class Node:
    id: str
    label: str
    keyword: str = "component"
    kind: str = ""

    @property
    def css(self) -> str:
        return self.kind or self.keyword


@dataclass
class Edge:
    src: str
    dst: str
    label: str = ""
    both: bool = False


@dataclass
class Model:
    nodes: dict[str, Node] = field(default_factory=dict)
    edges: list[Edge] = field(default_factory=list)
    groups: dict[str, tuple[str, list[str]]] = field(default_factory=dict)


@dataclass
class Impact:
    depends: list[str] = field(default_factory=list)
    changes: list[str] = field(default_factory=list)
    nodes: list[Node] = field(default_factory=list)
    edges: list[Edge] = field(default_factory=list)
    flows: list[Edge] = field(default_factory=list)


@dataclass
class Spec:
    number: str
    slug: str
    path: Path
    impact: Impact | None


# ---------------------------------------------------------------- parsing


def _lines(block: str):
    for number, raw in enumerate(block.splitlines(), 1):
        line = raw.strip()
        if line and not line.startswith("#"):
            try:
                yield number, shlex.split(line)
            except ValueError as error:
                raise ArchError(f"line {number}: {error}") from error


def _edge(number: int, tokens: list[str]) -> Edge:
    if len(tokens) < 3 or tokens[1] not in ("->", "<->") or len(tokens) > 4:
        raise ArchError(f'line {number}: use `A -> B ["label"]` or `A <-> B ["label"]`')
    return Edge(tokens[0], tokens[2], tokens[3] if len(tokens) == 4 else "", tokens[1] == "<->")


def _node(number: int, tokens: list[str]) -> Node:
    if len(tokens) < 3:
        raise ArchError(f'line {number}: use `{tokens[0]} ID "Label"`')
    node = Node(tokens[1], tokens[2], tokens[0])
    for extra in tokens[3:]:
        key, _, value = extra.partition("=")
        if key != "kind" or not value:
            raise ArchError(f"line {number}: unknown option `{extra}`")
        node.kind = value
    return node


def parse_model(text: str) -> Model:
    model = Model()
    for _, block in _blocks(text, "arch"):
        for number, tokens in _lines(block):
            word = tokens[0]
            if word in NODE_KEYWORDS:
                node = _node(number, tokens)
                if node.id in model.nodes:
                    raise ArchError(f"line {number}: id `{node.id}` is used twice")
                model.nodes[node.id] = node
            elif word == "edge":
                model.edges.append(_edge(number, tokens[1:]))
            elif word == "group":
                if ":" not in tokens or len(tokens) < 4:
                    raise ArchError(f'line {number}: use `group ID "Label" : ID ID ...`')
                colon = tokens.index(":")
                model.groups[tokens[1]] = (tokens[2], tokens[colon + 1 :])
            else:
                raise ArchError(f"line {number}: unknown word `{word}`")
    for edge in model.edges:
        for end in (edge.src, edge.dst):
            if end not in model.nodes:
                raise ArchError(f"edge uses unknown id `{end}`")
    for gid, (_, members) in model.groups.items():
        for member in members:
            if member not in model.nodes:
                raise ArchError(f"group `{gid}` uses unknown id `{member}`")
    return model


def parse_impact(text: str) -> Impact | None:
    blocks = _blocks(text, "arch-impact")
    if not blocks:
        return None
    impact = Impact()
    for _, block in blocks:
        for number, tokens in _lines(block):
            word, rest = tokens[0], tokens[1:]
            if word == "depends":
                impact.depends += rest
            elif word == "changes":
                impact.changes += rest
            elif word == "adds" and rest and rest[0] in NODE_KEYWORDS:
                impact.nodes.append(_node(number, rest))
            elif word == "adds" and rest and rest[0] == "edge":
                impact.edges.append(_edge(number, rest[1:]))
            elif word == "flow":
                impact.flows.append(_edge(number, rest))
            else:
                raise ArchError(f"line {number}: unknown line `{' '.join(tokens)}`")
    return impact


def _normal(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n")


def _blocks(text: str, name: str) -> list[tuple[str, str]]:
    return [(n, b) for n, b in FENCE.findall(_normal(text)) if n == name]


def load_model(path: Path) -> Model:
    try:
        return parse_model(path.read_text(encoding="utf-8"))
    except ArchError as error:
        raise ArchError(f"{path}: {error}") from error
    except OSError as error:
        raise ArchError(f"cannot read the model file: {error}") from error


def load_specs(specs_dir: Path) -> list[Spec]:
    specs = []
    for folder in sorted(specs_dir.glob("[0-9][0-9][0-9]-*")):
        spec_file = folder / "spec.md"
        if not spec_file.is_file():
            continue
        number, _, slug = folder.name.partition("-")
        try:
            impact = parse_impact(spec_file.read_text(encoding="utf-8"))
        except ArchError as error:
            raise ArchError(f"{spec_file}: {error}") from error
        specs.append(Spec(number, slug, spec_file, impact))
    return specs


def find_spec(specs: list[Spec], number: str) -> Spec:
    for spec in specs:
        if spec.number == number.zfill(3):
            return spec
    raise ArchError(f"no spec with number {number}")


# -------------------------------------------------------------- rendering


def _id(raw: str) -> str:
    safe = re.sub(r"\W", "_", raw)
    return safe + "_" if safe.lower() == "end" else safe


def _text(raw: str) -> str:
    return raw.replace('"', "#quot;")


def _shape(node: Node) -> str:
    label = _text(node.label)
    if node.keyword == "store":
        return f'{_id(node.id)}[("{label}")]'
    return f'{_id(node.id)}["{label}"]'


def _merge(model: Model, impact: Impact | None) -> Model:
    """Return the model plus the nodes and edges that a spec adds."""
    merged = Model(dict(model.nodes), list(model.edges), dict(model.groups))
    if impact:
        for node in impact.nodes:
            merged.nodes[node.id] = node
    return merged


def _check_ids(model: Model, impact: Impact, where: str) -> None:
    known = set(model.nodes) | {node.id for node in impact.nodes}
    for used in (
        impact.changes
        + [e.src for e in impact.edges + impact.flows]
        + [e.dst for e in impact.edges + impact.flows]
    ):
        if used not in known:
            raise ArchError(f"{where}: unknown component id `{used}`")


def render_model(model: Model, impact: Impact | None = None, focus: bool = False) -> str:
    """Draw the model. With an impact, mark changed and added parts."""
    changed = set(impact.changes) if impact else set()
    added = {n.id for n in impact.nodes} if impact else set()
    added_edges = list(impact.edges) if impact else []
    merged = _merge(model, impact)
    all_edges = [(e, False) for e in merged.edges] + [(e, True) for e in added_edges]

    keep = set(merged.nodes)
    if focus and impact:
        core = changed | added
        keep = set(core)
        for edge in added_edges:
            keep |= {edge.src, edge.dst}

    out = ["flowchart TD"]
    grouped: set[str] = set()
    for gid, (label, members) in merged.groups.items():
        members = [m for m in members if m in keep]
        if not members:
            continue
        out.append(f'    subgraph {_id(gid)}["{_text(label)}"]')
        for member in members:
            out.append(f"        {_shape(merged.nodes[member])}")
        out.append("    end")
        grouped |= set(members)
    for nid, node in merged.nodes.items():
        if nid in keep and nid not in grouped:
            out.append(f"    {_shape(node)}")

    for edge, is_added in all_edges:
        if edge.src not in keep or edge.dst not in keep:
            continue
        label = f'|"{_text(edge.label)}"|' if edge.label else ""
        arrow = {
            (False, False): "-->",
            (True, False): "<-->",
            (False, True): "-.->",
            (True, True): "<-.->",
        }[(edge.both, is_added)]
        out.append(f"    {_id(edge.src)} {arrow}{label} {_id(edge.dst)}")

    used = {merged.nodes[n].css for n in keep}
    for nid in keep:
        if nid in changed:
            used.add("changed")
        if nid in added:
            used.add("added")
    for name in sorted(used & PALETTE.keys()):
        out.append(f"    classDef {name} {PALETTE[name]}")
    for nid in sorted(keep):
        css = "added" if nid in added else "changed" if nid in changed else merged.nodes[nid].css
        if css in PALETTE:
            out.append(f"    class {_id(nid)} {css}")
    return "\n".join(out)


def render_deps(specs: list[Spec]) -> str:
    numbers = {s.number for s in specs}
    edges = []
    for spec in specs:
        for dep in spec.impact.depends if spec.impact else []:
            edges.append((dep.zfill(3), spec.number))
    shown = {a for a, _ in edges} | {b for _, b in edges}
    shown |= {s.number for s in specs if s.impact}
    out = ["flowchart LR"]
    for spec in specs:
        if spec.number in shown:
            out.append(f'    s{spec.number}["{spec.number} {_text(spec.slug)}"]')
    for dep, number in edges:
        if dep in numbers:
            out.append(f"    s{dep} --> s{number}")
    return "\n".join(out)


def render_flow(model: Model, impact: Impact) -> str:
    merged = _merge(model, impact)
    _check_ids(model, impact, "flow")
    order: list[str] = []
    for step in impact.flows:
        for end in (step.src, step.dst):
            if end not in order:
                order.append(end)
    out = ["sequenceDiagram"]
    for nid in order:
        out.append(f"    participant {_id(nid)} as {_text(merged.nodes[nid].label)}")
    for step in impact.flows:
        out.append(f"    {_id(step.src)}->>{_id(step.dst)}: {_text(step.label) or '...'}")
    return "\n".join(out)


# ---------------------------------------------------------------- checks


def check(model: Model, specs: list[Spec], strict: bool) -> tuple[list[str], list[str]]:
    errors, warnings = [], []
    numbers = {s.number for s in specs}
    for spec in specs:
        if spec.impact is None:
            (errors if strict else warnings).append(f"{spec.number}: no `arch-impact` block")
            continue
        try:
            _check_ids(model, spec.impact, spec.number)
        except ArchError as error:
            errors.append(str(error))
        for dep in spec.impact.depends:
            if dep.zfill(3) not in numbers:
                errors.append(f"{spec.number}: depends on unknown spec `{dep}`")
            elif dep.zfill(3) == spec.number:
                errors.append(f"{spec.number}: depends on itself")
        if not spec.impact.changes and not spec.impact.nodes:
            warnings.append(f"{spec.number}: block changes and adds nothing")
    return errors, warnings


# ---------------------------------------------------------------- inject


def make(name: str, model: Model, specs: list[Spec]) -> str:
    """Make the diagram for a marker name: model, deps, spec-N, spec-N-full, flow-N."""
    if name == "model":
        return render_model(model)
    if name == "deps":
        return render_deps(specs)
    match = re.fullmatch(r"(spec|flow)-(\d+)(-full)?", name)
    if not match:
        raise ArchError(f"unknown diagram name `{name}`")
    kind, number, full = match.groups()
    spec = find_spec(specs, number)
    if spec.impact is None:
        raise ArchError(f"spec {spec.number} has no `arch-impact` block")
    _check_ids(model, spec.impact, spec.number)
    if kind == "flow":
        return render_flow(model, spec.impact)
    return render_model(model, spec.impact, focus=not full)


def inject(path: Path, model: Model, specs: list[Spec], write: bool) -> list[str]:
    """Fill the marker regions of one file. Return the names that were stale."""
    text = _normal(path.read_text(encoding="utf-8"))
    stale: list[str] = []

    def fill(match: re.Match) -> str:
        body = f"```mermaid\n{make(match['name'], model, specs)}\n```\n"
        if body != match["body"]:
            stale.append(match["name"])
        return match[1] + body + match[4]

    new_text = MARKER.sub(fill, text)
    if write and stale:
        path.write_text(new_text, encoding="utf-8", newline="\n")
    return stale


# ------------------------------------------------------------------- cli


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="repository root")
    parser.add_argument("--model", type=Path, default=Path("docs/architecture.md"))
    parser.add_argument("--specs", type=Path, default=Path("specs"))
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("model", help="draw the whole architecture")
    sub.add_parser("deps", help="draw the dependencies between specs")
    p_spec = sub.add_parser("spec", help="draw the architecture with the impact of one spec")
    p_spec.add_argument("number")
    p_spec.add_argument("--full", action="store_true", help="show all components")
    p_flow = sub.add_parser("flow", help="draw the flows of one spec as a sequence diagram")
    p_flow.add_argument("number")
    p_check = sub.add_parser("check", help="check that specs use known component ids")
    p_check.add_argument("--strict", action="store_true", help="a missing block is an error")
    p_inject = sub.add_parser("inject", help="fill <!-- mermaid:NAME --> regions in files")
    p_inject.add_argument("files", nargs="+", type=Path)
    p_inject.add_argument("--check", action="store_true", help="do not write; fail if stale")
    args = parser.parse_args(argv)

    root = args.root.resolve()
    try:
        model = load_model(root / args.model)
        specs = load_specs(root / args.specs)
        if args.command == "model":
            print(render_model(model))
        elif args.command == "deps":
            print(render_deps(specs))
        elif args.command in ("spec", "flow"):
            name = f"{args.command}-{args.number}" + ("-full" if getattr(args, "full", False) else "")
            print(make(name, model, specs))
        elif args.command == "check":
            errors, warnings = check(model, specs, args.strict)
            missing = [w.split(":")[0] for w in warnings if w.endswith("no `arch-impact` block")]
            if missing:
                print(f"warning: {len(missing)} specs have no `arch-impact` block: {' '.join(missing)}")
            for line in warnings:
                if not line.endswith("no `arch-impact` block"):
                    print(f"warning: {line}")
            for line in errors:
                print(f"error: {line}")
            return 1 if errors else 0
        elif args.command == "inject":
            code = 0
            for file in args.files:
                stale = inject(root / file, model, specs, write=not args.check)
                for name in stale:
                    print(f"{file}: {'stale' if args.check else 'updated'}: {name}")
                code = 1 if (stale and args.check) else code
            return code
    except ArchError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    raise SystemExit(main())
