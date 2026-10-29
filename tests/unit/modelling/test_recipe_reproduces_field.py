"""A model recipe should rebuild a model that gives the same scalar field.

The field depends on how each feature was built: the number of elements, the interpolator type
and the regularisation weight. The recipe stores the parameters each feature was *requested*
with (not the element count that was realised, which would not reproduce the grid), and the
rebuilt model repeats the request.
"""
import numpy as np
import pytest

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius

CASES = [
    pytest.param({}, id="defaults"),
    pytest.param({"nelements": 2000}, id="nelements"),
    pytest.param({"interpolatortype": "PLI"}, id="interpolatortype"),
    pytest.param({"regularisation": 10.0}, id="regularisation"),
    pytest.param(
        {"nelements": 2000, "regularisation": 10.0, "damp": False}, id="several-parameters"
    ),
]


def build_model(**build_parameters):
    data, bb = load_claudius()
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.data = data
    model.stratigraphic_column.clear(basement=False)
    model.create_and_add_foliation("strati", **build_parameters)
    model.update()
    return model, bb


@pytest.mark.parametrize("build_parameters", CASES)
def test_rebuilt_model_gives_the_same_field(build_parameters):
    model, bb = build_model(**build_parameters)

    rebuilt = GeologicalModel.from_recipe_json(model.to_recipe_json())
    rebuilt.update()

    points = np.random.default_rng(0).uniform(bb[0], bb[1], (300, 3))
    original = model.evaluate_feature_value("strati", points)
    again = rebuilt.evaluate_feature_value("strati", points)
    np.testing.assert_allclose(again, original, atol=1e-6)

