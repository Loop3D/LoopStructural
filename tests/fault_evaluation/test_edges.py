"""Legacy edge behavior not represented by the fixed interior golden points."""

from __future__ import annotations

import numpy as np

from LoopStructural import GeologicalModel
from LoopStructural.modelling.core.stratigraphic_column import StratigraphicColumn
from LoopStructural.modelling.features import GeologicalFeature
from tests.fault_evaluation.models import build_model


def test_empty_outside_and_nan_value_queries_keep_legacy_results():
    model = build_model("f0_u0")
    feature = model["group_0"]
    bounds = model.bounding_box
    outside = bounds.maximum.copy()
    outside[0] += bounds.maximum[0] - bounds.origin[0]

    np.testing.assert_array_equal(feature.evaluate_value(np.empty((0, 3))), np.empty(0))
    values = feature.evaluate_value(np.vstack((outside, [np.nan, 250.0, -250.0])))
    assert values.shape == (2,)
    np.testing.assert_allclose(values, [np.nan, 0.0], rtol=0.0, atol=0.0, equal_nan=True)


def test_ignore_regions_bypasses_only_the_region_mask():
    class LinearInterpolator:
        def evaluate_value(self, points):
            return points[:, 0] + 2.0 * points[:, 1] - points[:, 2]

    feature = GeologicalFeature(
        "regional_field",
        builder=None,
        regions=[lambda points: points[:, 0] > 0.5],
        interpolator=LinearInterpolator(),
    )
    points = np.array([[0.25, 2.0, 1.0], [0.75, 2.0, 1.0]])

    np.testing.assert_allclose(feature.evaluate_value(points), [np.nan, 3.75], equal_nan=True)
    np.testing.assert_allclose(
        feature.evaluate_value(points, ignore_regions=True), [3.25, 3.75], rtol=0.0, atol=0.0
    )


def test_reverse_fault_and_fault_toggles_match_legacy_references():
    model = build_model("f1_u0")
    feature = model["group_0"]
    points = np.array([[450.0, 500.0, -500.0], [550.0, 500.0, -450.0]])

    # Frozen at ddc333ea; this checks the legacy reverse path on real faults.
    reverse_expected = np.array(
        [
            [450.00012670148914, 498.9067023889334, -519.6971159034398],
            [550.0000403813251, 501.0925651557694, -430.2937116988142],
        ]
    )
    reverse = feature._apply_faults(points.copy(), reverse=True)
    np.testing.assert_allclose(reverse, reverse_expected, rtol=0.0, atol=1e-5)
    assert not np.allclose(reverse, points)

    faults_enabled = feature.faults_enabled
    fault = feature.faults[0]
    displacement = fault.displacement
    try:
        feature.faults_enabled = False
        np.testing.assert_array_equal(feature._apply_faults(points.copy()), points)
        feature.faults_enabled = faults_enabled
        fault.displacement = 0.0
        np.testing.assert_array_equal(fault.apply_to_points(points), points)
    finally:
        feature.faults_enabled = faults_enabled
        fault.displacement = displacement


def test_first_public_model_evaluation_builds_a_lazy_feature_once():
    class LinearInterpolator:
        def evaluate_value(self, points):
            return points.sum(axis=1)

    class CountingBuilder:
        def __init__(self):
            self.interpolator = LinearInterpolator()
            self.faults = []
            self._built = False
            self.build_calls = 0

        def up_to_date(self):
            if not self._built:
                self.build()

        def build(self):
            self.build_calls += 1
            self._built = True

    model = GeologicalModel([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    builder = CountingBuilder()
    feature = GeologicalFeature("lazy", builder, faults=[])
    model._add_feature(feature)
    column = StratigraphicColumn()
    column.clear(basement=False)
    unit = column.add_unit("covered")
    unit.min_value, unit.max_value = 0.0, 10.0
    column.group_mapping["Group_0"] = "lazy"
    model.stratigraphic_column = column
    points = np.array([[0.25, 0.5, 0.75], [0.5, 0.5, 0.5]])

    np.testing.assert_array_equal(model.evaluate_model(points), [0, 0])
    assert builder.build_calls == 1
    np.testing.assert_array_equal(model.evaluate_model(points), [0, 0])
    assert builder.build_calls == 1
