"""Tests for the backend-matrix helper (FR-030, SC-018).

These tests check the helper and the real MATRIX cells (what they name, not their result).
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from tests.fixtures import backend_matrix as bm
from tests.fixtures.backend_matrix import (
    OPTIONAL_DEPENDENCIES,
    SPEC_ROWS,
    Cell,
    OptionalStatus,
)

KEYS = ("native", "surfe", "viewer", "matplotlib", "pyevtk")
MODULES = {
    "native": "loop_interpolation_core",
    "surfe": "surfepy",
    "viewer": "loopstructuralvisualisation",
    "matplotlib": "matplotlib",
    "pyevtk": "pyevtk",
}


def _cell(mandatory=True, dependency=None, refs=None, min_cases=1, stub=False, row="Python P2"):
    if refs is None:
        refs = (("a/b/test_x.py", "test_f", None),)
    return Cell(
        cell_id="c",
        spec_row=row,
        mandatory=mandatory,
        stub=stub,
        refs=refs,
        min_cases=min_cases,
        dependency=dependency,
    )


def _status(key, available):
    module = MODULES[key]
    reason = "" if available else (
        "optional dependency '%s' is not installed (%s)"
        % (module, OPTIONAL_DEPENDENCIES[key].distribution)
    )
    return OptionalStatus(key=key, module=module, available=available, reason=reason)


def _xml(cases):
    """Make a junit text from (classname, name, kind) tuples."""
    body = []
    for classname, name, kind in cases:
        inner = {
            "passed": "",
            "failed": "<failure message='x'/>",
            "skipped": "<skipped message='s'/>",
        }[kind]
        body.append('<testcase classname="%s" name="%s">%s</testcase>' % (classname, name, inner))
    return "<testsuites><testsuite name='p'>%s</testsuite></testsuites>" % "".join(body)


def _write(tmp_path, name, cases):
    path = tmp_path / name
    path.write_text(_xml(cases), encoding="utf-8")
    return path


def test_spec_rows_and_empty_matrix():
    assert SPEC_ROWS == (
        "Python FDI and PLI/P1",
        "Python P2",
        "Python FDI-fold, PLI discrete-fold, and P2-fold",
        "Python FDI and PLI constant-norm variants",
        "Common 2D and 3D geometry",
        "Existing native FDI or fold integration",
        "Optional Surfe integration",
    )
    assert "extra" not in SPEC_ROWS
    assert set(OPTIONAL_DEPENDENCIES) == set(KEYS)


@pytest.mark.parametrize("key", KEYS)
def test_optional_status_names_exact_module(key):
    absent = bm.optional_status(key, find_spec=lambda name: None)
    assert absent.available is False
    assert absent.module == MODULES[key]
    assert "'%s'" % MODULES[key] in absent.reason
    assert OPTIONAL_DEPENDENCIES[key].distribution in absent.reason
    assert absent.reason == "optional dependency '%s' is not installed (%s)" % (
        MODULES[key],
        OPTIONAL_DEPENDENCIES[key].distribution,
    )
    seen = []

    def finder(name):
        seen.append(name)
        return object()

    present = bm.optional_status(key, find_spec=finder)
    assert present.available is True
    assert present.reason == ""
    assert seen == [MODULES[key]]


@pytest.mark.parametrize("key", ["native", "surfe"])
def test_require_optional_skips_with_named_dependency(key, monkeypatch):
    real = bm.optional_status
    monkeypatch.setattr(bm, "optional_status", lambda k: real(k, find_spec=lambda name: None))
    with pytest.raises(pytest.skip.Exception) as info:
        bm.require_optional(key)
    assert MODULES[key] in str(info.value)


def test_require_optional_returns_none_when_available(monkeypatch):
    real = bm.optional_status
    monkeypatch.setattr(bm, "optional_status", lambda k: real(k, find_spec=lambda n: object()))
    assert bm.require_optional("native") is None


@pytest.mark.parametrize("outcome", ["passed", "failed", "skipped", "absent"])
@pytest.mark.parametrize("kind", ["mandatory", "optional_absent", "optional_installed"])
def test_classify_cell_rules(kind, outcome):
    if kind == "mandatory":
        cell, status = _cell(mandatory=True), None
    elif kind == "optional_absent":
        cell, status = _cell(mandatory=False, dependency="native"), _status("native", False)
    else:
        cell, status = _cell(mandatory=False, dependency="native"), _status("native", True)
    result = bm.classify_cell(cell, outcome, status)
    if outcome == "failed":
        assert result == "failed"
    elif kind == "mandatory":
        if outcome == "passed":
            assert result == "passed"
        else:
            assert result.startswith("VIOLATION: ")
    elif kind == "optional_absent":
        if outcome == "skipped":
            assert result == "excluded: " + status.reason
            assert "loop_interpolation_core" in result
        else:
            assert result.startswith("VIOLATION: ")
    else:
        if outcome == "passed":
            assert result == "passed"
        else:
            assert result.startswith("VIOLATION: ")


def test_classify_cell_tuple_one_absent_names_only_absent_module():
    cell = _cell(mandatory=False, dependency=("viewer", "matplotlib"))
    statuses = (_status("viewer", False), _status("matplotlib", True))
    result = bm.classify_cell(cell, "skipped", statuses)
    assert result.startswith("excluded: ")
    assert "loopstructuralvisualisation" in result
    assert "matplotlib" not in result
    assert bm.classify_cell(cell, "passed", statuses).startswith("VIOLATION: ")


def test_classify_cell_tuple_both_available_skipped_is_violation():
    cell = _cell(mandatory=False, dependency=("viewer", "matplotlib"))
    statuses = [_status("viewer", True), _status("matplotlib", True)]
    assert bm.classify_cell(cell, "skipped", statuses).startswith("VIOLATION: ")
    assert bm.classify_cell(cell, "passed", statuses) == "passed"


def test_classify_cell_tuple_both_absent_names_each_module():
    cell = _cell(mandatory=False, dependency=("viewer", "pyevtk"))
    statuses = (_status("viewer", False), _status("pyevtk", False))
    result = bm.classify_cell(cell, "skipped", statuses)
    assert "loopstructuralvisualisation" in result and "'pyevtk'" in result
    assert "; " in result


def test_cell_dependency_keys():
    assert bm.cell_dependency_keys(_cell(dependency=None)) == ()
    assert bm.cell_dependency_keys(_cell(dependency="surfe")) == ("surfe",)
    assert bm.cell_dependency_keys(_cell(dependency=("viewer", "matplotlib"))) == (
        "viewer",
        "matplotlib",
    )


def test_outcomes_from_junit_four_outcomes(tmp_path):
    path = _write(
        tmp_path,
        "a.xml",
        [
            ("a.b.test_p", "test_ok[x]", "passed"),
            ("a.b.test_f", "test_bad", "failed"),
            ("a.b.test_s", "test_skip", "skipped"),
        ],
    )
    cases = bm.outcomes_from_junit([path])

    def one(file, func):
        return bm.cell_outcome(
            _cell(refs=((file, func, None),)), cases
        )

    assert one("a/b/test_p.py", "test_ok") == "passed"
    assert one("a/b/test_f.py", "test_bad") == "failed"
    assert one("a/b/test_s.py", "test_skip") == "skipped"
    assert one("a/b/test_missing.py", "test_none") == "absent"


def test_outcomes_from_junit_min_cases_and_selector(tmp_path):
    path = _write(
        tmp_path,
        "a.xml",
        [("a.b.test_p", "test_ok[fdi]", "passed"), ("a.b.test_p", "test_ok[pli]", "passed")],
    )
    cases = bm.outcomes_from_junit([str(path)])
    refs = (("a/b/test_p.py", "test_ok", None),)
    assert bm.cell_outcome(_cell(refs=refs, min_cases=2), cases) == "passed"
    assert bm.cell_outcome(_cell(refs=refs, min_cases=3), cases) == "absent"
    sel = (("a/b/test_p.py", "test_ok", "fdi"),)
    assert bm.cell_outcome(_cell(refs=sel, min_cases=1), cases) == "passed"
    assert bm.cell_outcome(_cell(refs=sel, min_cases=2), cases) == "absent"


def test_outcomes_from_junit_class_test(tmp_path):
    path = _write(tmp_path, "a.xml", [("a.b.test_x.TestC", "test_m", "passed")])
    cases = bm.outcomes_from_junit([path])
    cell = _cell(refs=(("a/b/test_x.py", "test_m", None),))
    assert bm.cell_outcome(cell, cases) == "passed"
    other = _cell(refs=(("a/b/test_xy.py", "test_m", None),))
    assert bm.cell_outcome(other, cases) == "absent"


def test_outcomes_from_junit_two_files_merged(tmp_path):
    one = _write(tmp_path, "one.xml", [("a.b.test_x", "test_f", "passed")])
    two = _write(tmp_path, "two.xml", [("a.b.test_y", "test_g", "passed")])
    cases = bm.outcomes_from_junit([one, two])
    cell = _cell(
        refs=(("a/b/test_x.py", "test_f", None), ("a/b/test_y.py", "test_g", None)),
        min_cases=2,
    )
    assert bm.cell_outcome(cell, cases) == "passed"
    assert bm.cell_outcome(cell, bm.outcomes_from_junit([one])) == "absent"


def test_outcomes_from_junit_same_key_worse_outcome_wins(tmp_path):
    one = _write(tmp_path, "one.xml", [("a.b.test_x", "test_f", "passed")])
    two = _write(tmp_path, "two.xml", [("a.b.test_x", "test_f", "skipped")])
    three = _write(tmp_path, "three.xml", [("a.b.test_x", "test_f", "failed")])
    cell = _cell(refs=(("a/b/test_x.py", "test_f", None),))
    assert bm.cell_outcome(cell, bm.outcomes_from_junit([one, two])) == "skipped"
    assert bm.cell_outcome(cell, bm.outcomes_from_junit([two, three, one])) == "failed"


def test_outcomes_from_junit_module_level_skip(tmp_path):
    path = _write(tmp_path, "a.xml", [("", "a.b.test_x", "skipped")])
    cases = bm.outcomes_from_junit([path])
    refs = (("a/b/test_x.py", "test_f", None), ("a/b/test_x.py", "test_g", "sel"))
    optional = _cell(mandatory=False, dependency="native", refs=refs)
    mandatory = _cell(mandatory=True, refs=refs)
    assert bm.cell_outcome(optional, cases) == "skipped"
    assert bm.cell_outcome(mandatory, cases) == "skipped"
    result = bm.classify_cell(optional, "skipped", _status("native", False))
    assert result.startswith("excluded: ")
    assert bm.classify_cell(mandatory, "skipped", None).startswith("VIOLATION: ")


def test_rows_without_real_cell():
    def make(row, stub):
        return Cell("x", row, True, stub, (), 1, None)

    assert bm.rows_without_real_cell(()) == list(SPEC_ROWS)
    matrix = tuple(make(row, False) for row in SPEC_ROWS[:5])
    matrix += (make(SPEC_ROWS[5], True), make("extra", False))
    assert bm.rows_without_real_cell(matrix) == [SPEC_ROWS[5], SPEC_ROWS[6]]


def test_render_table_and_main(tmp_path, capsys, monkeypatch):
    cell = Cell(
        "demo", "Optional Surfe integration", False, False,
        (("a/b/test_x.py", "test_f", None),), 1, "surfe",
    )
    stub = Cell(
        "demo-stub", "Optional Surfe integration", True, True,
        (("a/b/test_x.py", "test_s", None),), 1, None,
    )
    statuses = {"surfe": _status("surfe", False)}
    table = bm.render_table(
        (cell, stub), statuses, {"demo": "skipped", "demo-stub": "passed"}
    )
    assert "demo" in table and "excluded" in table and "surfepy" in table
    assert "Python side" in table
    monkeypatch.setattr(bm, "MATRIX", ())
    path = _write(tmp_path, "a.xml", [("a.b.test_x", "test_s", "passed")])
    assert bm.main([str(path)]) == 0
    assert "Python side" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# The real MATRIX
# ---------------------------------------------------------------------------

REPO = Path(__file__).resolve().parents[2]
MANDATORY_IDS = (
    "py-fdi",
    "py-pli",
    "py-p2",
    "py-fold-fdi",
    "py-fold-pli",
    "py-fold-p2",
    "py-cn-fdi",
    "py-cn-pli",
    "common-geometry",
    "native-boundary",
    "surfe-boundary",
)
OPTIONAL_IDS = ("native-extension", "surfe-extension", "display-viewer", "export-pyevtk")
SKIP_NAMES = frozenset({"skip", "skipif", "importorskip", "xfail", "require_optional"})
STUB_HEAD = "Boundary stub cells (Python side; not a native or Surfe result)"


def _by_id(cell_id):
    return next(c for c in bm.MATRIX if c.cell_id == cell_id)


def test_matrix_has_15_cells_and_11_mandatory():
    ids = [c.cell_id for c in bm.MATRIX]
    assert len(ids) == 15
    assert len(set(ids)) == 15
    assert set(ids) == set(MANDATORY_IDS) | set(OPTIONAL_IDS)
    assert {c.cell_id for c in bm.MATRIX if c.mandatory} == set(MANDATORY_IDS)
    assert sum(1 for c in bm.MATRIX if c.mandatory) == 11
    assert sum(1 for c in bm.MATRIX if not c.mandatory) == 4
    for cell in bm.MATRIX:
        assert cell.refs
        assert cell.min_cases >= 1
        assert len(set(cell.refs)) == len(cell.refs)


def test_every_spec_row_has_a_cell():
    assert bm.rows_without_real_cell(bm.MATRIX) == []
    for cell in bm.MATRIX:
        assert cell.spec_row in SPEC_ROWS or cell.spec_row == "extra"
    # A stub cell alone does not satisfy a row: remove the real cells of rows 6 and 7.
    kept = tuple(
        c for c in bm.MATRIX if c.stub or c.spec_row not in (SPEC_ROWS[5], SPEC_ROWS[6])
    )
    assert len(kept) == len(bm.MATRIX) - 2
    assert bm.rows_without_real_cell(kept) == [SPEC_ROWS[5], SPEC_ROWS[6]]


def test_stub_cells_are_mandatory_and_listed_apart():
    stubs = [c for c in bm.MATRIX if c.stub]
    assert {c.cell_id for c in stubs} == {"native-boundary", "surfe-boundary"}
    for cell in stubs:
        assert cell.mandatory
        assert cell.dependency is None
    statuses = {key: _status(key, False) for key in KEYS}
    outcomes = {c.cell_id: "passed" for c in bm.MATRIX if c.mandatory}
    outcomes.update({c.cell_id: "skipped" for c in bm.MATRIX if not c.mandatory})
    table = bm.render_table(bm.MATRIX, statuses, outcomes)
    assert STUB_HEAD in table
    rows, _, stub_block = table.partition(STUB_HEAD)
    row6 = next(line for line in rows.splitlines() if line.startswith("| " + SPEC_ROWS[5]))
    row7 = next(line for line in rows.splitlines() if line.startswith("| " + SPEC_ROWS[6]))
    assert "native-extension" in row6 and "excluded" in row6
    assert "surfe-extension" in row7 and "excluded" in row7
    assert "boundary" not in rows
    assert "- native-boundary" in stub_block and "passed" in stub_block
    assert "- surfe-boundary" in stub_block
    assert "VIOLATION" not in table


def test_stub_passed_extension_skipped_dependency_absent_reads_excluded():
    statuses = {key: _status(key, key != "native") for key in KEYS}
    outcomes = {"native-boundary": "passed", "native-extension": "skipped"}
    outcomes.update({c.cell_id: "passed" for c in bm.MATRIX if c.cell_id not in outcomes})
    table = bm.render_table(bm.MATRIX, statuses, outcomes)
    rows, _, stub_block = table.partition(STUB_HEAD)
    row6 = next(line for line in rows.splitlines() if line.startswith("| " + SPEC_ROWS[5]))
    assert row6.endswith("| excluded: " + statuses["native"].reason + " |")
    stub_line = next(line for line in stub_block.splitlines() if "native-boundary" in line)
    assert stub_line.endswith(": passed")


def test_cell_with_two_dependencies():
    cell = _by_id("display-viewer")
    assert cell.mandatory is False
    assert cell.spec_row == "extra"
    assert bm.cell_dependency_keys(cell) == ("viewer", "matplotlib")
    statuses = {key: _status(key, True) for key in KEYS}
    statuses["matplotlib"] = _status("matplotlib", False)
    got = [statuses[k] for k in bm.cell_dependency_keys(cell)]
    result = bm.classify_cell(cell, "skipped", got)
    assert result.startswith("excluded: ")
    assert "'matplotlib'" in result
    assert "loopstructuralvisualisation" not in result
    assert bm.classify_cell(cell, "passed", got).startswith("VIOLATION: ")
    # Both installed: a pass is a pass and a skip is a violation.
    both = [_status("viewer", True), _status("matplotlib", True)]
    assert bm.classify_cell(cell, "passed", both) == "passed"
    assert bm.classify_cell(cell, "skipped", both).startswith("VIOLATION: ")
    table = bm.render_table(bm.MATRIX, statuses, {"display-viewer": "skipped"})
    line = next(row for row in table.splitlines() if "display-viewer" in row)
    assert "excluded" in line and "matplotlib" in line


def test_unknown_dependency_key_raises():
    cell = _cell(mandatory=False, dependency="nonexistent")
    statuses = {key: _status(key, True) for key in KEYS}
    with pytest.raises(KeyError):
        bm.render_table((cell,), statuses, {"c": "skipped"})


def test_dependencies_of_optional_cells():
    assert _by_id("native-extension").dependency == "native"
    assert _by_id("surfe-extension").dependency == "surfe"
    assert _by_id("export-pyevtk").dependency == "pyevtk"
    for cell_id in MANDATORY_IDS:
        assert _by_id(cell_id).dependency is None


def _parse(root, rel):
    return ast.parse((Path(root) / rel).read_text(encoding="utf-8-sig"), filename=rel)


def _top_functions(tree):
    return {
        n.name: n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
    }


def _skip_uses(node, allow=None):
    """Return (name, line) for each use of a skip name below ``node``.

    A function named ``allow`` is not entered.
    """
    hits = []
    stack = [node]
    while stack:
        cur = stack.pop()
        if allow is not None and isinstance(cur, ast.FunctionDef) and cur.name == allow:
            continue
        name = None
        if isinstance(cur, ast.Name):
            name = cur.id
        elif isinstance(cur, ast.Attribute):
            name = cur.attr
        elif isinstance(cur, ast.alias):
            name = cur.name.rsplit(".", 1)[-1]
        if name in SKIP_NAMES:
            hits.append((name, getattr(cur, "lineno", 0)))
        stack.extend(ast.iter_child_nodes(cur))
    return hits


def find_skips(root, refs):
    """List the skip uses that can touch the tests named by ``refs``.

    The scan covers the referenced functions (decorators and body), the module-level
    statements and ``pytestmark`` of each file, ``tests/conftest.py``, every
    ``tests/fixtures/*.py`` and the ``conftest.py`` of the folder of each ref. The
    definition of ``require_optional`` is the one exception.
    """
    root = Path(root)
    hits = []
    seen = set()

    def whole_file(rel):
        if rel in seen or not (root / rel).is_file():
            return
        seen.add(rel)
        found = _skip_uses(_parse(root, rel), allow="require_optional")
        hits.extend((rel, n, ln) for n, ln in found)

    whole_file("tests/conftest.py")
    for path in sorted((root / "tests" / "fixtures").glob("*.py")):
        whole_file(path.relative_to(root).as_posix())
    for rel, func, _ in refs:
        whole_file(rel.rsplit("/", 1)[0] + "/conftest.py")
        tree = _parse(root, rel)
        node = _top_functions(tree).get(func)
        if node is not None:
            hits.extend((rel + "::" + func, n, ln) for n, ln in _skip_uses(node))
        for stmt in tree.body:
            if not isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                hits.extend((rel, n, ln) for n, ln in _skip_uses(stmt))
    return hits


@pytest.mark.parametrize("cell_id", MANDATORY_IDS)
def test_mandatory_cells_name_existing_tests(cell_id):
    cell = _by_id(cell_id)
    assert cell.mandatory
    for rel, func, selector in cell.refs:
        assert "\\" not in rel
        assert (REPO / rel).is_file(), rel
        assert func in _top_functions(_parse(REPO, rel)), (rel, func)
        assert selector is None or selector


@pytest.mark.parametrize("cell_id", MANDATORY_IDS)
def test_mandatory_cells_have_no_skip(cell_id):
    cell = _by_id(cell_id)
    assert find_skips(REPO, cell.refs) == []


def test_find_skips_sees_a_skip_in_a_temporary_file(tmp_path):
    (tmp_path / "tests" / "fixtures").mkdir(parents=True)
    (tmp_path / "tests" / "sub").mkdir()
    conftest = tmp_path / "tests" / "conftest.py"
    conftest.write_text("X = 1\n", encoding="utf-8")
    test_file = tmp_path / "tests" / "sub" / "test_a.py"
    clean = "def test_a():\n    assert True\n"
    test_file.write_text(clean, encoding="utf-8")
    refs = (("tests/sub/test_a.py", "test_a", None),)
    assert find_skips(tmp_path, refs) == []
    for text in (
        "import pytest\n\n\ndef test_a():\n    pytest.skip('x')\n",
        "import pytest\n\n\n@pytest.mark.skipif(True, reason='x')\ndef test_a():\n    pass\n",
        "import pytest\npytestmark = pytest.mark.xfail\n\n\ndef test_a():\n    pass\n",
        "import pytest\npytest.importorskip('numpy')\n\n\ndef test_a():\n    pass\n",
        "def test_a():\n    require_optional('native')\n",
    ):
        test_file.write_text(text, encoding="utf-8")
        assert find_skips(tmp_path, refs), text
    test_file.write_text(clean, encoding="utf-8")
    skip_module = "import pytest\npytest.skip('x', allow_module_level=True)\n"
    conftest.write_text(skip_module, encoding="utf-8")
    assert find_skips(tmp_path, refs)
    conftest.write_text("X = 1\n", encoding="utf-8")
    helper = tmp_path / "tests" / "fixtures" / "helper.py"
    helper.write_text("def f():\n    require_optional('native')\n", encoding="utf-8")
    assert find_skips(tmp_path, refs)
    helper.unlink()
    (tmp_path / "tests" / "sub" / "conftest.py").write_text(skip_module, encoding="utf-8")
    assert find_skips(tmp_path, refs)


def test_main_prints_the_stub_note(tmp_path, capsys):
    path = _write(tmp_path, "empty.xml", [("a.b.test_x", "test_f", "passed")])
    assert bm.main([str(path)]) == 1  # the mandatory cells are absent
    out = capsys.readouterr().out
    assert STUB_HEAD in out
    assert "not a native or Surfe pass" in out
