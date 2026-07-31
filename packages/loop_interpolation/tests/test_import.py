import pytest

# Import the module to test
from loop_interpolation import *

# List of classes to test for importability
classes_to_test = [
    "InterpolatorType",
    "GeologicalInterpolator",
    "DiscreteInterpolator",
    "FiniteDifferenceInterpolator",
    "PiecewiseLinearInterpolator",
    "DiscreteFoldInterpolator",
    "SurfeRBFInterpolator",
    "P1Interpolator",
    "P2Interpolator",
    "TetMesh",
    "StructuredGrid",
    "UnStructuredTetMesh",
    "P1Unstructured2d",
    "P2Unstructured2d",
    "StructuredGrid2D",
    "P2UnstructuredTetMesh",
]


@pytest.mark.parametrize("class_name", classes_to_test)
def test_import_class(class_name):
    """Test if a class can be imported from the interpolation module."""
    assert class_name in globals(), f"{class_name} is not importable from the interpolation module."
