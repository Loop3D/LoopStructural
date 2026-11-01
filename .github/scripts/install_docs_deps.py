"""Install the documentation dependencies against the packages in this checkout.

The root pyproject pins loop-common and loop-interpolation below the versions in packages/, and the
visualisation package depends on LoopStructural, so a single ``pip install -e .[docs]`` either pulls
the old releases from PyPI or fails to resolve. Install the third-party requirements from the
pyprojects, then install the workspace packages with --no-deps.
"""

import re
import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKSPACE = ("loop-common", "loop-interpolation", "loopstructural", "loopstructuralvisualisation")
PACKAGES = ("packages/loop_common", "packages/loop_interpolation", "packages/loopstructural_visualisation")


def project(path):
    with open(ROOT / path / "pyproject.toml", "rb") as f:
        return tomllib.load(f)["project"]


def third_party(requirements):
    """Drop the requirements that are provided by this checkout."""
    return [r for r in requirements if re.split(r"[\s\[<>=!~;]", r, maxsplit=1)[0].lower().replace("_", "-") not in WORKSPACE]


root = project(".")
requirements = third_party(root["dependencies"]) + third_party(root["optional-dependencies"]["docs"])
for path in PACKAGES:
    meta = project(path)
    requirements += third_party(meta.get("dependencies", []))
    for extra in ("native", "all"):
        requirements += third_party(meta.get("optional-dependencies", {}).get(extra, []))

pip = [sys.executable, "-m", "pip", "install"]
subprocess.check_call(pip + sorted(set(requirements)))
subprocess.check_call(pip + ["--no-deps", "-e", str(ROOT)] + [a for p in PACKAGES for a in ("-e", str(ROOT / p))])
