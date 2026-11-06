from __future__ import annotations

import ast
import subprocess
import sys

import loop_link

# Modules that the core import must not load. See SC-005.
FORBIDDEN_TOP_LEVEL = {
    "PyQt5",
    "PyQt6",
    "PySide2",
    "PySide6",
    "qgis",
    "pyvista",
    "pandas",
    "loop_common",
    "LoopStructural",
}


def _loaded_forbidden_modules(import_target: str) -> subprocess.CompletedProcess[str]:
    # Run in a fresh isolated interpreter so that modules loaded by pytest do not count.
    code = (
        f"import {import_target}\n"
        "import sys\n"
        f"forbidden = {sorted(FORBIDDEN_TOP_LEVEL)!r}\n"
        "print(sorted(m for m in sys.modules if m.split('.')[0] in set(forbidden)))\n"
    )
    return subprocess.run(
        [sys.executable, "-I", "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )


def test_core_import_loads_no_qt_or_optional_extras() -> None:
    result = _loaded_forbidden_modules("loop_link")

    assert result.returncode == 0, result.stderr
    loaded = ast.literal_eval(result.stdout.strip())
    assert loaded == [], result.stderr


def test_all_names_exist() -> None:
    missing = [name for name in loop_link.__all__ if not hasattr(loop_link, name)]

    assert missing == []


def test_adapters_import_loads_no_qt_or_optional_extras() -> None:
    result = _loaded_forbidden_modules("loop_link.adapters")

    assert result.returncode == 0, result.stderr
    loaded = ast.literal_eval(result.stdout.strip())
    assert loaded == [], result.stderr
