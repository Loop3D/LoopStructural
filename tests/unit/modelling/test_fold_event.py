def test_constructor():

    pass


def test_constant_fold_axis():
    pass


def test_rotation_fold_axis():
    pass


from types import SimpleNamespace

import numpy as np

from LoopStructural.modelling.features.fold import FoldEvent


class _PlaneFeature:
    """Stub fold frame coordinate with a constant gradient"""

    def __init__(self, gradient):
        self.gradient = np.array(gradient, dtype=float)

    def evaluate_value(self, points):
        return points @ self.gradient

    def evaluate_gradient(self, points):
        return np.tile(self.gradient, (points.shape[0], 1))


def _fold_event(invert_norm=False):
    frame = SimpleNamespace(
        features=[_PlaneFeature([2, 0, 0]), _PlaneFeature([0, 3, 0]), _PlaneFeature([0, 0, 1])]
    )
    return FoldEvent(
        frame,
        fold_limb_rotation=lambda gx: np.full_like(gx, 30.0),
        fold_axis=np.array([0, 1.0, 0]),
        invert_norm=invert_norm,
    )


def test_deformed_orientation_is_unit_length():
    points = np.random.default_rng(0).random((20, 3))
    fold_direction, fold_axis, dgz = _fold_event().get_deformed_orientation(points)
    assert fold_direction.shape == points.shape
    assert dgz.shape == points.shape
    assert np.allclose(np.linalg.norm(fold_direction, axis=1), 1.0)
    assert np.allclose(np.linalg.norm(dgz, axis=1), 1.0)
    assert np.allclose(np.einsum("ij,ij->i", fold_direction, fold_axis), 0.0)


def test_invert_norm_flips_dgz():
    points = np.random.default_rng(0).random((20, 3))
    _, _, dgz = _fold_event().get_deformed_orientation(points)
    _, _, inverted = _fold_event(invert_norm=True).get_deformed_orientation(points)
    assert inverted.shape == points.shape
    assert np.allclose(inverted, -dgz)
