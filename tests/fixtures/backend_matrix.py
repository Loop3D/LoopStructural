"""Backend matrix: cells, optional dependencies and pure checker functions (FR-030, SC-018).

A cell is a group of tests for one backend or capability. A mandatory cell must pass.
An optional cell needs a named dependency. It is "excluded" only when the dependency is
absent and the cell was skipped.

Run ``python -m tests.fixtures.backend_matrix <junit>...`` from the repository root. The
command prints the table and exits with 1 on a violation.
"""

from __future__ import annotations

import importlib.util
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

try:  # Literal exists in typing from Python 3.8.
    from typing import Literal
except ImportError:  # pragma: no cover
    from typing_extensions import Literal  # type: ignore

Outcome = Literal["passed", "failed", "skipped", "absent"]

# The 7 rows of "Required backend coverage" in spec.md. The value "extra" is not a row.
SPEC_ROWS: Tuple[str, ...] = (
    "Python FDI and PLI/P1",
    "Python P2",
    "Python FDI-fold, PLI discrete-fold, and P2-fold",
    "Python FDI and PLI constant-norm variants",
    "Common 2D and 3D geometry",
    "Existing native FDI or fold integration",
    "Optional Surfe integration",
)

_SEVERITY = {"passed": 0, "skipped": 1, "failed": 2}


@dataclass(frozen=True)
class OptionalDependency:
    """One optional dependency.

    Parameters
    ----------
    key : str
        Key in ``OPTIONAL_DEPENDENCIES``.
    module : str
        Exact import name that ``find_spec`` looks for.
    distribution : str
        Name to install (pip name or extra).
    capability : str
        What the dependency gives.
    """

    key: str
    module: str
    distribution: str
    capability: str


OPTIONAL_DEPENDENCIES: Dict[str, OptionalDependency] = {
    d.key: d
    for d in (
        OptionalDependency(
            "native",
            "loop_interpolation_core",
            "loop-interpolation[native]",
            "native FDI and fold extension",
        ),
        # The pip name of Surfe is unverified (see T092).
        OptionalDependency("surfe", "surfepy", "surfepy", "Surfe RBF extension"),
        OptionalDependency(
            "viewer", "loopstructuralvisualisation", "loopstructuralvisualisation", "2D viewer"
        ),
        OptionalDependency("matplotlib", "matplotlib", "matplotlib", "plot output"),
        OptionalDependency("pyevtk", "pyevtk", "pyevtk", "VTK export"),
    )
}


@dataclass(frozen=True)
class OptionalStatus:
    """Status of one optional dependency.

    Parameters
    ----------
    key : str
        Key in ``OPTIONAL_DEPENDENCIES``.
    module : str
        Exact module name.
    available : bool
        True when the module can be imported.
    reason : str
        Empty when available. Else the text that names the module and the distribution.
    """

    key: str
    module: str
    available: bool
    reason: str


@dataclass(frozen=True)
class Cell:
    """One cell of the backend matrix.

    Parameters
    ----------
    cell_id : str
        Unique id of the cell.
    spec_row : str
        A value of ``SPEC_ROWS`` or ``"extra"``.
    mandatory : bool
        True when the cell must pass with no skip.
    stub : bool
        True when the cell uses a recording stub. It is not a result for a spec row.
    refs : tuple of (str, str, str or None)
        Test refs: path with "/", test function, id selector or None.
    min_cases : int
        Least number of passed cases for the outcome "passed".
    dependency : str or tuple of str or None
        Key or keys of ``OPTIONAL_DEPENDENCIES``. The cell is "dependency absent" when
        ANY key is not available.
    """

    cell_id: str
    spec_row: str
    mandatory: bool
    stub: bool
    refs: Tuple[Tuple[str, str, Optional[str]], ...]
    min_cases: int
    dependency: Union[str, Tuple[str, ...], None]


MATRIX: Tuple[Cell, ...] = ()

StatusArg = Union[OptionalStatus, Sequence[OptionalStatus], None]


def optional_status(
    key: str, find_spec: Callable[[str], object] = importlib.util.find_spec
) -> OptionalStatus:
    """Find out if an optional dependency is installed.

    Parameters
    ----------
    key : str
        Key in ``OPTIONAL_DEPENDENCIES``.
    find_spec : callable
        Function with the signature of ``importlib.util.find_spec``.

    Returns
    -------
    OptionalStatus
        The status. ``reason`` is empty when the module is available.
    """
    dep = OPTIONAL_DEPENDENCIES[key]
    try:
        found = find_spec(dep.module)
    except (ImportError, ValueError):
        found = None
    if found is not None:
        return OptionalStatus(key, dep.module, True, "")
    reason = "optional dependency '%s' is not installed (%s)" % (dep.module, dep.distribution)
    return OptionalStatus(key, dep.module, False, reason)


def require_optional(key: str) -> None:
    """Skip the running test when an optional dependency is not installed.

    Parameters
    ----------
    key : str
        Key in ``OPTIONAL_DEPENDENCIES``.
    """
    import pytest

    status = optional_status(key)
    if not status.available:
        pytest.skip(status.reason)


def cell_dependency_keys(cell: Cell) -> Tuple[str, ...]:
    """Return the dependency keys of a cell as a tuple.

    Parameters
    ----------
    cell : Cell
        The cell.

    Returns
    -------
    tuple of str
        Empty when the cell has no dependency.
    """
    dep = cell.dependency
    if dep is None:
        return ()
    if isinstance(dep, str):
        return (dep,)
    return tuple(dep)


def classify_cell(cell: Cell, outcome: str, status: StatusArg) -> str:
    """Classify the outcome of a cell.

    Parameters
    ----------
    cell : Cell
        The cell.
    outcome : str
        One of "passed", "failed", "skipped", "absent".
    status : OptionalStatus or sequence of OptionalStatus or None
        Status of the dependencies of the cell.

    Returns
    -------
    str
        "passed", "failed", "excluded: <reason>" or "VIOLATION: <text>".
    """
    if outcome == "failed":
        return "failed"
    if status is None:
        statuses: List[OptionalStatus] = []
    elif isinstance(status, OptionalStatus):
        statuses = [status]
    else:
        statuses = list(status)
    if cell.mandatory:
        if outcome == "passed":
            return "passed"
        return "VIOLATION: mandatory cell '%s' is %s" % (cell.cell_id, outcome)
    absent = [s for s in statuses if not s.available]
    if absent:
        if outcome == "skipped":
            return "excluded: " + "; ".join(s.reason for s in absent)
        return "VIOLATION: optional cell '%s' is %s while a dependency is absent (%s)" % (
            cell.cell_id,
            outcome,
            "; ".join(s.module for s in absent),
        )
    if outcome == "passed":
        return "passed"
    return "VIOLATION: optional cell '%s' is %s while its dependency is installed" % (
        cell.cell_id,
        outcome,
    )


def outcomes_from_junit(paths: Iterable[object]) -> Dict[str, Outcome]:
    """Read junit files and merge the cases.

    Parameters
    ----------
    paths : iterable of path-like
        Junit XML files.

    Returns
    -------
    dict
        Key ``"<classname>::<name>"``, value the outcome. For one key in several files,
        the worse outcome wins (failed, skipped, passed).
    """
    merged: Dict[str, Outcome] = {}
    for path in paths:
        root = ET.parse(str(path)).getroot()
        for case in root.iter("testcase"):
            if case.find("failure") is not None or case.find("error") is not None:
                outcome: Outcome = "failed"
            elif case.find("skipped") is not None:
                outcome = "skipped"
            else:
                outcome = "passed"
            key = "%s::%s" % (case.get("classname", ""), case.get("name", ""))
            old = merged.get(key)
            if old is None or _SEVERITY[outcome] > _SEVERITY[old]:
                merged[key] = outcome
    return merged


def _ref_matches(ref: Tuple[str, str, Optional[str]], classname: str, name: str) -> bool:
    path, func, selector = ref
    module = (path[:-3] if path.endswith(".py") else path).replace("/", ".")
    if classname == "":
        return name == module  # module-level skip
    if classname != module and not classname.startswith(module + "."):
        return False
    if name.split("[", 1)[0] != func:
        return False
    return selector is None or selector in name


def cell_outcome(cell: Cell, cases: Mapping[str, str]) -> Outcome:
    """Compute the outcome of a cell from merged junit cases.

    Parameters
    ----------
    cell : Cell
        The cell.
    cases : mapping
        Result of ``outcomes_from_junit``.

    Returns
    -------
    str
        "failed" if one matching case failed, else "skipped" if one was skipped, else
        "passed" if at least ``min_cases`` passed, else "absent".
    """
    matched: Dict[str, str] = {}
    for key, outcome in cases.items():
        classname, _, name = key.partition("::")
        if any(_ref_matches(ref, classname, name) for ref in cell.refs):
            matched[key] = outcome
    values = list(matched.values())
    if "failed" in values:
        return "failed"
    if "skipped" in values:
        return "skipped"
    if values.count("passed") >= max(cell.min_cases, 1):
        return "passed"
    return "absent"


def rows_without_real_cell(matrix: Iterable[Cell]) -> List[str]:
    """List the spec rows that have no real (non-stub) cell.

    Parameters
    ----------
    matrix : iterable of Cell
        The cells.

    Returns
    -------
    list of str
        Rows of ``SPEC_ROWS`` in order.
    """
    covered = {c.spec_row for c in matrix if not c.stub}
    return [row for row in SPEC_ROWS if row not in covered]


def _cell_statuses(cell: Cell, statuses: Mapping[str, OptionalStatus]) -> List[OptionalStatus]:
    return [statuses[k] for k in cell_dependency_keys(cell) if k in statuses]


def render_table(
    matrix: Iterable[Cell],
    statuses: Mapping[str, OptionalStatus],
    outcomes: Optional[Mapping[str, str]] = None,
) -> str:
    """Render the matrix as text for ``verification.md``.

    Parameters
    ----------
    matrix : iterable of Cell
        The cells.
    statuses : dict
        Status of each optional dependency, by key.
    outcomes : mapping, optional
        Outcome of each cell by ``cell_id``. A missing cell counts as "absent".

    Returns
    -------
    str
        The table. Stub cells are in a separate block.
    """
    cells = list(matrix)
    outcomes = outcomes or {}

    def result(cell: Cell) -> str:
        return classify_cell(
            cell, outcomes.get(cell.cell_id, "absent"), _cell_statuses(cell, statuses)
        )

    lines = ["| Spec row | Cell | Level | Result |", "|---|---|---|---|"]
    for row in SPEC_ROWS:
        real = [c for c in cells if c.spec_row == row and not c.stub]
        if not real:
            lines.append("| %s | (none) | - | VIOLATION: no real cell |" % row)
        for c in real:
            level = "mandatory" if c.mandatory else "optional"
            lines.append("| %s | %s | %s | %s |" % (row, c.cell_id, level, result(c)))
    for c in cells:
        if c.spec_row not in SPEC_ROWS and not c.stub:
            level = "mandatory" if c.mandatory else "optional"
            lines.append("| %s | %s | %s | %s |" % (c.spec_row, c.cell_id, level, result(c)))
    lines.append("")
    lines.append("Stub cells (they test the Python side only; not a native or Surfe pass):")
    stubs = [c for c in cells if c.stub]
    if not stubs:
        lines.append("(none)")
    for c in stubs:
        lines.append("- %s (%s): %s" % (c.cell_id, c.spec_row, result(c)))
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Print the matrix table for the given junit files.

    Parameters
    ----------
    argv : sequence of str, optional
        Paths of junit files. Default is ``sys.argv[1:]``.

    Returns
    -------
    int
        1 when a cell is a violation, else 0.
    """
    args = list(sys.argv[1:] if argv is None else argv)
    cases = outcomes_from_junit(args)
    statuses = {key: optional_status(key) for key in OPTIONAL_DEPENDENCIES}
    outcomes = {c.cell_id: cell_outcome(c, cases) for c in MATRIX}
    print(render_table(MATRIX, statuses, outcomes))
    print("Stub cells test the Python side of the boundary only. They are not a native or "
          "Surfe pass.")
    violation = any(
        classify_cell(c, outcomes[c.cell_id], _cell_statuses(c, statuses)).startswith(
            "VIOLATION"
        )
        for c in MATRIX
    )
    return 1 if violation else 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
