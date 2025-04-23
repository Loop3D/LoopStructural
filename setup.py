"""See pyproject.toml for project metadata."""

import os

from setuptools import setup

package_root = os.path.abspath(os.path.dirname(__file__))

version = {}
with open(os.path.join(package_root, "LoopStructural/version.py")) as fp:
    exec(fp.read(), version)
version = version["__version__"]
setup()
