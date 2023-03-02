from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius
from LoopStructural.interpolators import PiecewiseLinearInterpolator
from LoopStructural.interpolators import FiniteDifferenceInterpolator
import numpy as np

def test_create_geological_model():
    model = GeologicalModel([0, 0, 0], [5, 5, 5])
    assert (model.origin-np.array([0,0,0])).sum() == 0
    assert (model.maximum-np.array([5,5,5])).sum() == 0

def test_create_scale_factor_model():
    model = GeologicalModel([0, 0, 0], [5, 5, 5], rescale=True)
    assert model.scale_factor == 5


def test_access_feature_model():
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(data)
    s0 = model.create_and_add_foliation(
        "s0", interpolatortype="FDI", nelements=1000, solver="fake", damp=False
    )
    assert s0 == model["s0"]


def test_get_interpolator(interpolator):
    model = GeologicalModel([0, 0, 0], [5, 5, 5])
    interpolator = model.get_interpolator(
        interpolatortype="PLI", nelements=1e5, buffer=0.2
    )
    assert type(interpolator) == PiecewiseLinearInterpolator


def test_element_number_PLI():
    model = GeologicalModel([0, 0, 0], [5, 5, 5])
    interpolator = model.get_interpolator(
        interpolatortype="PLI", nelements=1e5, buffer=0.2
    )
    assert np.log10(interpolator.support.n_nodes) - 5 < 1
    interpolator = model.get_interpolator(
        interpolatortype="PLI", nelements=1e6, buffer=0.2
    )
    assert np.log10(interpolator.support.n_nodes) - 6 < 1
    interpolator = model.get_interpolator(
        interpolatortype="PLI", nelements=3e4, buffer=0.2
    )
    assert np.log10(interpolator.support.n_nodes) - 4 < 1


def test_element_number_FDI():
    model = GeologicalModel([0, 0, 0], [5, 5, 5])
    interpolator = model.get_interpolator(
        interpolatortype="FDI", nelements=1e5, buffer=0.2
    )
    assert np.log10(interpolator.support.n_nodes) - 5 < 1
    interpolator = model.get_interpolator(
        interpolatortype="FDI", nelements=1e6, buffer=0.2
    )
    assert np.log10(interpolator.support.n_nodes) - 6 < 1
    interpolator = model.get_interpolator(
        interpolatortype="FDI", nelements=3e4, buffer=0.2
    )
    assert np.log10(interpolator.support.n_nodes) - 4 < 1


def test_buffer():
    model = GeologicalModel([0, 0, 0], [5, 5, 5], rescale=True)
    interpolator = model.get_interpolator(
        interpolatortype="FDI", nelements=1e5, buffer=0.2
    )
    assert np.sum(interpolator.support.origin + 0.2) == 0
