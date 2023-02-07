import pytest
from LoopStructural.modelling.features.builders import GeologicalFeatureBuilder
from tests.fixtures.horizontal_data import horizontal_data


def test_geological_feature_builder_constructor(interpolator):
    builder = GeologicalFeatureBuilder(interpolator)
    assert builder.interpolator == interpolator


def test_get_interpolator():
    pass


def test_add_data_to_interpolator(interpolator, horizontal_data):
    builder = GeologicalFeatureBuilder(interpolator)
    builder.add_data_from_data_frame(horizontal_data)
    # assert builder.data.shape =
    pass


def test_install_gradient_constraints():
    pass


def test_get_value_constraints():
    pass


def test_get_gradient_constraints():
    pass


def test_get_tangent_constraints():
    pass


def test_norm_constraints():
    pass


def get_orientation_constraints():
    pass


def get_data_locations():
    pass


def test_test_interpolation_geometry():
    pass


def test_not_up_to_date():
    """test to make sure that the feature
    isn't interpolated when everything is set up
    """
    pass


def test_get_feature():
    pass


def test_change_up_to_date():
    pass
