"""Tests for the backend-matrix helper (FR-030, SC-018).

These tests check the helper only. The tests of the real MATRIX cells come later.
"""

from __future__ import annotations

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
    assert bm.MATRIX == ()
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


def test_render_table_and_main(tmp_path, capsys):
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
    path = _write(tmp_path, "a.xml", [("a.b.test_x", "test_s", "passed")])
    assert bm.main([str(path)]) == 0
    assert "Python side" in capsys.readouterr().out
