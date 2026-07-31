import numpy as np

from LoopStructural import GeologicalModel
from LoopStructural.modelling.features._analytical_feature import (
    AnalyticalGeologicalFeature,
)


def test_evaluate_value_no_model_matches_plane_equation():
    """With model=None (the default used by the fault builder call sites),
    evaluate_value should simply be the signed distance of `pos` (already in
    world coordinates) from the plane defined by `origin`/`vector`."""
    feature = AnalyticalGeologicalFeature(
        name="plane",
        vector=np.array([1.0, 0.0, 0.0]),
        origin=np.array([15.0, 15.0, 15.0]),
    )

    pos = np.array([[16.0, 15.0, 15.0]])
    result = feature.evaluate_value(pos)

    assert np.allclose(result, [1.0])


def test_evaluate_value_with_real_model_is_not_double_transformed():
    """Regression test for a latent double-transform bug: under the
    affine-transform bounding box contract, `pos` passed to `evaluate_value`
    is already in world coordinates. `AnalyticalGeologicalFeature` must not
    additionally call `self.model.rescale` (a local->world transform) on it,
    since that would double-apply the transform whenever a real (non-None)
    `model` is supplied.

    We build a GeologicalModel whose bounding box has a non-zero local
    origin (anchored at the model's world-space `origin` via
    `set_local_transform`, as GeologicalModel.__init__ does for the
    2-argument constructor). If evaluate_value incorrectly rescaled `pos`
    through `model.rescale` (local->world, i.e. `pos + local_origin` for an
    identity rotation), the computed distance would be offset by the local
    origin instead of matching the correct, un-transformed plane distance.
    """
    origin = np.array([10.0, 10.0, 10.0])
    maximum = np.array([20.0, 20.0, 20.0])
    model = GeologicalModel(origin, maximum)
    # Sanity check: the model's local frame is anchored away from zero, so a
    # spurious rescale would actually shift the result.
    assert not np.allclose(model.bounding_box.local_origin, 0.0)

    feature = AnalyticalGeologicalFeature(
        name="plane",
        vector=np.array([1.0, 0.0, 0.0]),
        origin=np.array([15.0, 15.0, 15.0]),
        model=model,
    )

    pos = np.array([[16.0, 15.0, 15.0]])
    result = feature.evaluate_value(pos)

    # Correct (un-transformed) result: distance from x=15 plane to x=16 is 1.
    assert np.allclose(result, [1.0])
    # If the removed `model.rescale` call were still present, the result
    # would instead be offset by `local_origin[0]` (here, 11.0 not 1.0).
    assert not np.allclose(result, [1.0 + model.bounding_box.local_origin[0]])


def test_evaluate_gradient_unaffected_by_model():
    """evaluate_gradient returns the (constant) plane normal direction and is
    unaffected by whether a model is attached, consistent with pos/direction
    already being expressed in world space."""
    origin = np.array([10.0, 10.0, 10.0])
    maximum = np.array([20.0, 20.0, 20.0])
    model = GeologicalModel(origin, maximum)

    vector = np.array([0.0, 2.0, 0.0])
    feature = AnalyticalGeologicalFeature(
        name="plane",
        vector=vector,
        origin=np.array([15.0, 15.0, 15.0]),
        model=model,
    )

    pos = np.array([[16.0, 15.0, 15.0], [12.0, 11.0, 19.0]])
    gradient = feature.evaluate_gradient(pos)

    assert gradient.shape == pos.shape
    assert np.allclose(gradient, np.tile(vector, (pos.shape[0], 1)))
