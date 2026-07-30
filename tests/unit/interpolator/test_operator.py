import numpy as np
import pytest

from loop_interpolation._operator import Operator

ALL_MASKS = [
    "Dx_mask",
    "Dy_mask",
    "Dz_mask",
    "Dxx_mask",
    "Dyy_mask",
    "Dzz_mask",
    "Dxy_mask",
    "Dxz_mask",
    "Dyz_mask",
    "Lapacian",
]


@pytest.mark.parametrize("mask_name", ALL_MASKS)
def test_mask_shape(mask_name):
    mask = getattr(Operator, mask_name)
    assert mask.shape == (3, 3, 3)


@pytest.mark.parametrize(
    "mask_name",
    [
        "Dx_mask",
        "Dy_mask",
        "Dz_mask",
        "Dxx_mask",
        "Dyy_mask",
        "Dzz_mask",
        "Dxy_mask",
        "Dxz_mask",
        "Dyz_mask",
        "Lapacian",
    ],
)
def test_mask_sums_to_zero(mask_name):
    # All of the finite difference stencils should be translation invariant,
    # i.e. applying them to a constant field must give 0.
    mask = getattr(Operator, mask_name)
    assert np.isclose(mask.sum(), 0.0)


def test_dx_mask_values():
    # central difference coefficients along the last axis
    assert Operator.Dx_mask[1, 1, 0] == -0.5
    assert Operator.Dx_mask[1, 1, 1] == 0.0
    assert Operator.Dx_mask[1, 1, 2] == 0.5
    # everywhere else should be zero
    mask = Operator.Dx_mask.copy()
    mask[1, 1, :] = 0
    assert np.all(mask == 0)


def test_dy_mask_is_dx_swapaxes():
    assert np.array_equal(Operator.Dy_mask, Operator.Dx_mask.swapaxes(1, 2))


def test_dz_mask_is_dx_swapaxes():
    assert np.array_equal(Operator.Dz_mask, Operator.Dx_mask.swapaxes(0, 2))


def test_dxx_mask_values():
    assert Operator.Dxx_mask[1, 1, 0] == 1
    assert Operator.Dxx_mask[1, 1, 1] == -2
    assert Operator.Dxx_mask[1, 1, 2] == 1


def test_dyy_mask_is_dxx_swapaxes():
    assert np.array_equal(Operator.Dyy_mask, Operator.Dxx_mask.swapaxes(1, 2))


def test_dzz_mask_is_dxx_swapaxes():
    assert np.array_equal(Operator.Dzz_mask, Operator.Dxx_mask.swapaxes(0, 2))


def test_dxz_mask_is_dxy_swapaxes():
    assert np.array_equal(Operator.Dxz_mask, Operator.Dxy_mask.swapaxes(0, 1))


def test_dyz_mask_is_dxy_swapaxes():
    assert np.array_equal(Operator.Dyz_mask, Operator.Dxy_mask.swapaxes(0, 2))


def test_dxy_mask_scaling():
    # the mixed derivative mask is the corner differences scaled by 1/sqrt(2)
    expected_unscaled = np.array(
        [np.zeros((3, 3)), [[-0.25, 0, 0.25], [0, 0, 0], [0.25, 0, -0.25]], np.zeros((3, 3))]
    )
    assert np.allclose(Operator.Dxy_mask * np.sqrt(2), expected_unscaled)


def _grid_varying_along_axis(axis):
    """Build a 3x3x3 grid of values that increase linearly (0,1,2) along `axis`."""
    idx = np.arange(3, dtype=float)
    shape = [1, 1, 1]
    shape[axis] = 3
    return np.broadcast_to(idx.reshape(shape), (3, 3, 3))


def test_dx_mask_recovers_first_derivative_along_axis2():
    values = _grid_varying_along_axis(2)
    # central difference of f(k) = k with unit spacing gives derivative 1
    assert np.isclose(np.sum(Operator.Dx_mask * values), 1.0)


def test_dy_mask_recovers_first_derivative_along_axis1():
    values = _grid_varying_along_axis(1)
    assert np.isclose(np.sum(Operator.Dy_mask * values), 1.0)


def test_dz_mask_recovers_first_derivative_along_axis0():
    values = _grid_varying_along_axis(0)
    assert np.isclose(np.sum(Operator.Dz_mask * values), 1.0)


def test_dx_mask_zero_for_orthogonal_variation():
    # Dx_mask only touches the middle plane/row, varying values along axis 0
    # or axis 1 (rather than axis 2) should not contribute to the estimate
    # unless they fall on the row that is used (row 1 of axis1).
    values = _grid_varying_along_axis(0)
    assert np.isclose(np.sum(Operator.Dx_mask * values), 0.0)


def test_dxx_mask_recovers_second_derivative_along_axis2():
    # f(k) = (k-1)**2 -> values [1, 0, 1], f'' = 2 analytically
    values = np.zeros((3, 3, 3))
    coords = (np.arange(3) - 1) ** 2
    values[:, :, :] = coords.reshape(1, 1, 3)
    assert np.isclose(np.sum(Operator.Dxx_mask * values), 2.0)


def test_laplacian_mask_matches_discrete_laplacian_definition():
    expected = np.array(
        [
            [[0, 0, 0], [0, 1, 0], [0, 0, 0]],
            [[0, 1, 0], [1, -6, 1], [0, 1, 0]],
            [[0, 0, 0], [0, 1, 0], [0, 0, 0]],
        ]
    )
    assert np.array_equal(Operator.Lapacian, expected)


def test_laplacian_zero_for_harmonic_like_quadratic():
    # f(x,y,z) = x^2 + y^2 - 2z^2 is harmonic (Laplacian == 0)
    coords = np.arange(3) - 1
    x, y, z = np.meshgrid(coords, coords, coords, indexing="ij")
    values = x.astype(float) ** 2 + y.astype(float) ** 2 - 2 * z.astype(float) ** 2
    assert np.isclose(np.sum(Operator.Lapacian * values), 0.0)
