"""A model recipe should rebuild a model that gives the same scalar field.

``GeologicalModel.to_recipe_json`` / ``from_recipe_json`` are documented as a way to describe and
rebuild a model. The existing recipe tests check the names of features, the data and the
stratigraphic column, but not the field, and the field depends on how each feature was built:
the number of elements, the interpolator type and the regularisation weight.

Those build parameters are not stored: a feature entry in the recipe has only ``name``, ``kind``,
``faults`` and ``regions``. A model built with the defaults rebuilds exactly, because the rebuild
uses the defaults. A model built with anything else rebuilds with the defaults and gives a
different field (on the Claudius data, differences of 12 to 62 on a field of range about 560).

The non-default cases are strict xfails. When the recipe stores the build parameters they will
XPASS, which fails the run: remove the markers then.
"""
import numpy as np
import pytest

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius

RECIPE_DROPS_BUILD_PARAMETERS = pytest.mark.xfail(
    strict=True,
    reason="to_recipe_json stores only name, kind, faults and regions for a feature, not the "
    "build parameters (nelements, interpolatortype, regularisation, ...), so the rebuilt "
    "model uses the defaults",
)

CASES = [
    pytest.param({}, id="defaults"),
    pytest.param(
        {"nelements": 2000}, id="nelements", marks=RECIPE_DROPS_BUILD_PARAMETERS
    ),
    pytest.param(
        {"interpolatortype": "PLI"}, id="interpolatortype", marks=RECIPE_DROPS_BUILD_PARAMETERS
    ),
    pytest.param(
        {"regularisation": 10.0}, id="regularisation", marks=RECIPE_DROPS_BUILD_PARAMETERS
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

