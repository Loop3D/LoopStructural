"""See pyproject.toml for project metadata."""

import os
import runpy

from setuptools import setup

package_root = os.path.abspath(os.path.dirname(__file__))

version = runpy.run_path(os.path.join(package_root, "LoopStructural/version.py"))
version = version["__version__"]
setup()
