import numpy as np

from LoopStructural.utils import azimuthplunge2vector, strikedip2vector


def test_strikedip2vector():
    strike = [0, 45, 90]
    dip = [30, 60, 90]
    expected_result = np.array(
        [
            [0.5, 0.0, 0.8660254],
            [0.61237, -0.61237, 0.5],
            [0.0, -1.0, 0.0],
        ]
    )

    result = strikedip2vector(strike, dip)
    print(result - expected_result)
    assert np.allclose(result, expected_result, atol=1e-3)


# import numpy as np
# from LoopStructural.utils.maths import azimuthplunge2vector


def test_azimuthplunge2vector_single_values():
    plunge = 0
    plunge_dir = 90
    expected_result = np.array([[1, 0, 0]])
    result = azimuthplunge2vector(plunge, plunge_dir)
    assert np.allclose(result, expected_result)


def test_azimuthplunge2vector_array_values():
    plunge = [0, 90, 0]
    plunge_dir = [90, 90, 0]
    expected_result = np.array(
        [
            [1, 0, 0],
            [0.0, 0.0, -1.0],
            [0, 1, 0],
        ]
    )
    result = azimuthplunge2vector(plunge, plunge_dir)
    assert np.allclose(result, expected_result)


def test_azimuthplunge2vector_empty_arrays():
    plunge = []
    plunge_dir = []
    result = azimuthplunge2vector(plunge, plunge_dir)
    assert result.shape[0] == 0
