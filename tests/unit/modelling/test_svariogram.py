import numpy as np
import pytest

from LoopStructural.modelling.features.fold._svariogram import (
    SVariogram,
    find_peaks_and_troughs,
)


def test_find_peaks_and_troughs_alternating_series():
    x = np.arange(9)
    y = np.array([0, 1, 0, 1, 0, 1, 0, 1, 0])

    px, py = find_peaks_and_troughs(x, y)

    # every point is a local extremum in a strictly alternating series
    assert px == list(x)
    assert py == list(y)


def test_find_peaks_and_troughs_monotonic_series_only_endpoints():
    x = [0, 1, 2, 3]
    y = [0, 5, 10, 15]

    px, py = find_peaks_and_troughs(x, y)

    # for a monotonically increasing series there is no change in gradient
    # sign, so only the first and last points (always included) are
    # returned
    assert px == [0, 3]
    assert py == [0, 15]


def test_find_peaks_and_troughs_raises_on_mismatched_length():
    with pytest.raises(ValueError):
        find_peaks_and_troughs(np.arange(5), np.arange(4))


def test_svariogram_constructor_drops_nan_pairs():
    xdata = np.array([0.0, 1.0, 2.0, np.nan, 4.0, 5.0])
    ydata = np.array([0.0, 1.0, 2.0, 3.0, np.nan, 5.0])

    sv = SVariogram(xdata, ydata)

    # rows 3 and 4 (0-indexed) should be dropped because either x or y is nan
    assert np.array_equal(sv.xdata, np.array([0.0, 1.0, 2.0, 5.0]))
    assert np.array_equal(sv.ydata, np.array([0.0, 1.0, 2.0, 5.0]))


def test_svariogram_dist_and_variance_matrix_shapes_and_symmetry():
    xdata = np.array([0.0, 1.0, 3.0])
    ydata = np.array([2.0, 4.0, 8.0])

    sv = SVariogram(xdata, ydata)

    assert sv.dist.shape == (3, 3)
    assert sv.variance_matrix.shape == (3, 3)
    # distance and squared-difference matrices are symmetric
    assert np.allclose(sv.dist, sv.dist.T)
    assert np.allclose(sv.variance_matrix, sv.variance_matrix.T)
    # diagonal is always zero (distance/variance to self)
    assert np.allclose(np.diag(sv.dist), 0)
    assert np.allclose(np.diag(sv.variance_matrix), 0)
    assert np.isclose(sv.dist[0, 1], 1.0)
    assert np.isclose(sv.variance_matrix[0, 1], (2.0 - 4.0) ** 2)


def test_initialise_lags_with_explicit_step_and_nsteps():
    sv = SVariogram(np.arange(0, 10, dtype=float), np.arange(0, 10, dtype=float))
    sv.initialise_lags(step=2.0, nsteps=3)

    assert np.allclose(sv.lags, [1.0, 3.0, 5.0])


def test_initialise_lags_with_step_only_infers_nsteps():
    sv = SVariogram(np.arange(0, 10, dtype=float), np.arange(0, 10, dtype=float))
    sv.initialise_lags(step=1.0)

    # lags should cover the data range (0-9) in unit steps, offset by half
    assert np.isclose(sv.lags[0], 0.5)
    assert sv.lags[-1] < 10


def test_initialise_lags_auto_guesses_step_from_average_spacing():
    sv = SVariogram(np.arange(0, 20, dtype=float), np.arange(0, 20, dtype=float))
    sv.initialise_lags()

    assert sv.lags is not None
    assert len(sv.lags) > 0
    # nearest-neighbour spacing is 1, so guessed step should be 1 * 4 = 4
    assert np.isclose(sv.lags[1] - sv.lags[0], 4.0)


def test_initialise_lags_with_integer_dtype_input_raises_bug():
    """Documents a real bug: when xdata/ydata are integer arrays (e.g. the
    common `np.arange(0, 20)` without an explicit float dtype) and no step
    or nsteps are provided, `initialise_lags` tries to write `np.nan` into
    the (integer-dtype) copy of the distance matrix via
    `d[d == 0] = np.nan`, which raises `ValueError: cannot convert float NaN
    to integer`. The auto-guess path therefore does not work for integer
    input data, only for float input data.
    """
    sv = SVariogram(np.arange(0, 20), np.arange(0, 20))
    assert sv.xdata.dtype.kind in ("i", "u")
    with pytest.raises(ValueError):
        sv.initialise_lags()


def test_initialise_lags_cap_does_not_actually_limit_lag_count_bug():
    """Documents a real bug: when the auto-guessed number of steps exceeds
    200, the code logs "using 200" and recomputes `step` based on a local
    variable named `nstep` (200), but then builds `self.lags` using the
    *original* uncapped `nsteps` (`np.arange(step / 2.0, nsteps * step,
    step)`), not `nstep`. Because of this variable-name typo, the resulting
    number of lags is not actually capped at 200 as the log message claims.
    """
    xdata = np.linspace(0, 1000, 1000)
    ydata = np.sin(xdata)
    sv = SVariogram(xdata, ydata)
    sv.initialise_lags()

    # the log message claims lags are capped to 200, but they are not
    assert len(sv.lags) != 200


def test_calc_semivariogram_returns_expected_shapes():
    xdata = np.linspace(0, 20, 40)
    ydata = np.sin(xdata)
    sv = SVariogram(xdata, ydata)

    lags, variogram, npairs = sv.calc_semivariogram(step=1.0)

    assert lags.shape == variogram.shape == npairs.shape
    # every bin in this densely-sampled data should have at least one pair
    assert np.all(npairs > 0)
    # semivariogram values must be non-negative (they are means of squared
    # differences)
    assert np.all(variogram[~np.isnan(variogram)] >= 0)


def test_calc_semivariogram_uses_explicit_lags_when_given():
    xdata = np.linspace(0, 20, 40)
    ydata = np.sin(xdata)
    sv = SVariogram(xdata, ydata)

    custom_lags = np.array([1.0, 2.0, 3.0])
    lags, variogram, npairs = sv.calc_semivariogram(lags=custom_lags)

    assert np.array_equal(lags, custom_lags)
    assert np.array_equal(sv.lags, custom_lags)


def test_calc_semivariogram_raises_without_any_lag_information():
    # xdata/ydata with a single point can't infer a step size via nearest
    # neighbour distance (nanmin of an all-nan row), so no lags can be
    # determined and no step/nsteps/lags were provided
    sv = SVariogram(np.array([1.0]), np.array([1.0]))
    with pytest.raises(ValueError):
        sv.calc_semivariogram()


def test_find_wavelengths_detects_approximate_period_of_sinusoid():
    xdata = np.linspace(0, 100, 200)
    true_wavelength = 20.0
    ydata = 10 * np.sin(2 * np.pi * xdata / true_wavelength)

    sv = SVariogram(xdata, ydata)
    wavelengths = sv.find_wavelengths(step=1.0)

    assert len(wavelengths) == 2
    # the first (most reliable) wavelength guess should be reasonably close
    # to the true periodicity of the underlying signal
    assert abs(wavelengths[0] - true_wavelength) < 5.0


def test_find_wavelengths_falls_back_to_range_when_no_periodicity_found():
    xdata = np.linspace(0, 10, 20)
    ydata = np.linspace(0, 1, 20)  # perfectly linear, no periodicity

    sv = SVariogram(xdata, ydata)
    wavelengths = sv.find_wavelengths(step=0.5)

    assert wavelengths[0] == pytest.approx(2 * (xdata.max() - xdata.min()))
    assert wavelengths[1] == 0.0
