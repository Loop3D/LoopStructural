# The version is set in pyproject.toml. Do not write a version string here:
# release-please changes each "__version__ = ..." string in every version.py
# file in the repository to the version of the release it makes.
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("loopstructuralvisualisation")
except PackageNotFoundError:
    # the package is imported from source and is not installed
    __version__ = "unknown"
