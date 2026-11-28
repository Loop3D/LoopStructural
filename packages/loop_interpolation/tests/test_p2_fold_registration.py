"""Registration of the "P2Fold" interpolator type."""
import numpy as np
import pytest
from loop_common.geometry import BoundingBox
from loop_common.supports import P2UnstructuredTetMesh
from loop_interpolation import (
    InterpolatorFactory,
    InterpolatorType,
    P2FoldInterpolator,
    P2Interpolator,
)


@pytest.fixture
def bb():
    return BoundingBox(np.zeros(3), np.full(3, 10.0))


def test_factory_creates_p2fold(bb):
    i = InterpolatorFactory.create_interpolator("P2Fold", bb, nelements=500)
    assert type(i) is P2FoldInterpolator
    assert isinstance(i.support, P2UnstructuredTetMesh)
    assert i.type == InterpolatorType.PIECEWISE_QUADRATIC_FOLD


def test_p2_unchanged(bb):
    i = InterpolatorFactory.create_interpolator("P2", bb, nelements=500)
    assert type(i) is P2Interpolator


def test_to_dict_from_dict_round_trip(bb):
    i = InterpolatorFactory.create_interpolator("P2Fold", bb, nelements=500)
    i.set_value_constraints(np.array([[5.0, 5.0, 5.0, 1.0]]))
    d = i.to_dict()
    assert d["type"] == "PIECEWISE_QUADRATIC_FOLD"
    j = InterpolatorFactory.from_dict(d)
    assert type(j) is P2FoldInterpolator
    assert j.support.nodes.shape == i.support.nodes.shape
    assert np.allclose(j.data["value"], i.data["value"])
