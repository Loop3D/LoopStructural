import warnings

import pandas as pd
import pytest

from LoopStructural.geometry import BoundingBox
from LoopStructural.modelling.features.builders import StructuralFrameBuilder


def _builder_with_data(interpolatortype):
    bounding_box = BoundingBox([0, 0, 0], [1, 1, 1])
    builder = StructuralFrameBuilder(
        interpolatortype=interpolatortype,
        bounding_box=bounding_box,
        nelements=100,
        name="frame",
        model=None,
    )
    for i in range(3):
        builder.builders[i].data = pd.DataFrame({"dummy": [1]})
    return builder


def test_feature_is_alias_of_frame(interpolatortype):
    builder = _builder_with_data(interpolatortype)
    assert builder.feature is builder.frame


def test_setup_is_deprecated_alias_of_build(interpolatortype, monkeypatch):
    builder = _builder_with_data(interpolatortype)
    calls = []
    monkeypatch.setattr(builder, "build", lambda *a, **k: calls.append((a, k)))

    with pytest.warns(DeprecationWarning):
        builder.setup(w1=2.0)

    assert calls == [((), {"w1": 2.0})]


def test_build_uses_w3_for_coordinate2_orthogonality(interpolatortype, monkeypatch):
    builder = _builder_with_data(interpolatortype)

    calls = []
    monkeypatch.setattr(
        builder.builders[1],
        "add_orthogonal_feature",
        lambda feature, w, **kwargs: calls.append((feature, w)),
    )

    builder.build(w1=2.0, w2=3.0, w3=7.0)

    # second call on coordinate 1 orthogonalises against coordinate 2's feature and must use w3, not w2
    assert calls[-1][0] is builder.builders[2].feature
    assert calls[-1][1] == 7.0
