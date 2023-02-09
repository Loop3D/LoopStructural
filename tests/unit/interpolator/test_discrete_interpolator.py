import pytest
import numpy as np


def test_nx(interpolator, data):
    assert interpolator.nx == 21 * 21 * 21


def test_region(interpolator, data, region_func):
    """Test to see whether restricting the interpolator to a region works"""
    # interpolator = generate_interpolator(interpolator)
    interpolator.set_value_constraints(data[["X", "Y", "Z", "val", "w"]].to_numpy())
    interpolator._setup_interpolator()
    interpolator.set_region(region_func)
    # assert np.all(interpolator.region == region_func(interpolator.support.nodes))
    print(np.sum(interpolator.region))
    interpolator.solve_system()


def test_add_constraint_to_least_squares(interpolator):
    """make sure that when incorrect sized arrays are passed it doesn't get added"""
    pass


def test_update_interpolator():
    pass
