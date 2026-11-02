import tempfile
import unittest
from pathlib import Path

import spec_diagrams as sd

MODEL = '''# Architecture

```arch
component ui "app ui" kind=app
component core "core"
component scene "scene" kind=leaf
component far "far away"
external qgis "QGIS"
store disk "Scene file"
group ws "Workspace" : ui core scene
edge ui -> core "uses"
edge core -> scene
edge qgis <-> core "link"
edge core -> disk
edge far -> scene
```
'''

SPEC = '''# Feature spec

```arch-impact
depends 001
changes core
adds component panel "Panel" kind=app
adds edge panel -> core "calls"
flow qgis -> core "Set"
flow core -> ui "apply"
```
'''


def make_repo(root: Path, spec_text=SPEC):
    (root / "docs").mkdir()
    (root / "docs" / "architecture.md").write_text(MODEL, encoding="utf-8")
    for folder, text in (("001-base", SPEC.replace("depends 001\n", "")), ("002-next", spec_text)):
        (root / "specs" / folder).mkdir(parents=True)
        (root / "specs" / folder / "spec.md").write_text(text, encoding="utf-8")


class ParseTests(unittest.TestCase):
    def test_model(self):
        model = sd.parse_model(MODEL)
        self.assertEqual(set(model.nodes), {"ui", "core", "scene", "far", "qgis", "disk"})
        self.assertTrue(model.edges[2].both)
        self.assertEqual(model.groups["ws"][1], ["ui", "core", "scene"])

    def test_crlf(self):
        self.assertEqual(len(sd.parse_model(MODEL.replace("\n", "\r\n")).nodes), 6)

    def test_unknown_edge_end(self):
        with self.assertRaises(sd.ArchError):
            sd.parse_model("```arch\ncomponent a \"A\"\nedge a -> b\n```\n")

    def test_duplicate_id(self):
        with self.assertRaises(sd.ArchError):
            sd.parse_model("```arch\ncomponent a \"A\"\ncomponent a \"B\"\n```\n")

    def test_impact(self):
        impact = sd.parse_impact(SPEC)
        self.assertEqual(impact.depends, ["001"])
        self.assertEqual(impact.changes, ["core"])
        self.assertEqual(impact.nodes[0].id, "panel")
        self.assertEqual(len(impact.flows), 2)

    def test_no_impact_block(self):
        self.assertIsNone(sd.parse_impact("# nothing"))

    def test_bad_impact_line(self):
        with self.assertRaises(sd.ArchError):
            sd.parse_impact("```arch-impact\nwrong line\n```\n")


class RenderTests(unittest.TestCase):
    def setUp(self):
        self.model = sd.parse_model(MODEL)
        self.impact = sd.parse_impact(SPEC)

    def test_model_diagram(self):
        out = sd.render_model(self.model)
        self.assertTrue(out.startswith("flowchart TD"))
        self.assertIn('subgraph ws["Workspace"]', out)
        self.assertIn('disk[("Scene file")]', out)
        self.assertIn('qgis <-->|"link"| core', out)

    def test_overlay_marks_changes(self):
        out = sd.render_model(self.model, self.impact, focus=False)
        self.assertIn("class core changed", out)
        self.assertIn("class panel added", out)
        self.assertIn('panel -.->|"calls"| core', out)

    def test_focus_drops_far_nodes(self):
        out = sd.render_model(self.model, self.impact, focus=True)
        self.assertIn("class core changed", out)
        self.assertNotIn("far", out)
        self.assertNotIn("scene", out)
        self.assertNotIn("qgis", out)
        self.assertIn("panel", out)

    def test_flow(self):
        out = sd.render_flow(self.model, self.impact)
        self.assertIn("participant qgis as QGIS", out)
        self.assertIn("qgis->>core: Set", out)

    def test_reserved_word_id(self):
        self.assertEqual(sd._id("end"), "end_")

    def test_quote_in_label(self):
        self.assertNotIn('"x"y', sd._shape(sd.Node("a", 'x"y')))


class RepoTests(unittest.TestCase):
    def test_check_and_inject(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_repo(root)
            model = sd.load_model(root / "docs/architecture.md")
            specs = sd.load_specs(root / "specs")
            self.assertEqual(sd.check(model, specs, strict=True), ([], []))
            self.assertIn("s001 --> s002", sd.render_deps(specs))

            doc = root / "docs" / "view.md"
            doc.write_text("a\n<!-- mermaid:spec-002 -->\nold\n<!-- /mermaid:spec-002 -->\nb\n")
            self.assertEqual(sd.inject(doc, model, specs, write=False), ["spec-002"])
            self.assertIn("old", doc.read_text())
            self.assertEqual(sd.inject(doc, model, specs, write=True), ["spec-002"])
            self.assertEqual(sd.inject(doc, model, specs, write=False), [])
            self.assertIn("```mermaid", doc.read_text())

    def test_check_finds_problems(self):
        bad = SPEC.replace("changes core", "changes nope").replace("depends 001", "depends 009")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_repo(root, bad)
            model = sd.load_model(root / "docs/architecture.md")
            errors, _ = sd.check(model, sd.load_specs(root / "specs"), strict=False)
            self.assertEqual(len(errors), 2)

    def test_missing_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_repo(root, "# no block")
            model = sd.load_model(root / "docs/architecture.md")
            specs = sd.load_specs(root / "specs")
            self.assertEqual(sd.check(model, specs, strict=False)[0], [])
            self.assertEqual(len(sd.check(model, specs, strict=True)[0]), 1)

    def test_cli(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_repo(root)
            self.assertEqual(sd.main(["--root", str(root), "check", "--strict"]), 0)
            self.assertEqual(sd.main(["--root", str(root), "spec", "2"]), 0)
            self.assertEqual(sd.main(["--root", str(root), "spec", "9"]), 2)


if __name__ == "__main__":
    unittest.main()
