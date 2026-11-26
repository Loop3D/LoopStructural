import numpy as np

from LoopStructural.modelling.features import GeologicalFeature, LambdaGeologicalFeature
from LoopStructural.modelling.features._unconformity_feature import UnconformityFeature


def test_lambda_feature_with_unconformity_and_function_regions():
    parent = GeologicalFeature("parent", builder=None)
    unconformity = UnconformityFeature(parent, 0.5)
    unconformity.evaluate_value = lambda pos, **kwargs: pos[:, 2]
    feature = LambdaGeologicalFeature(
        lambda pos: pos[:, 0], regions=[unconformity, lambda pos: pos[:, 1] < 0.5]
    )
    pos = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 1.0], [3.0, 1.0, 0.0]])
    v = feature.evaluate_value(pos)
    np.testing.assert_array_equal(v, [1.0, np.nan, np.nan])
