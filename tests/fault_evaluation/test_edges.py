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
    # A point with a NaN coordinate gives NaN since master #413 (it gave 0 before).
    np.testing.assert_allclose(values, [np.nan, np.nan], rtol=0.0, atol=0.0, equal_nan=True)


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

    # Frozen from master 790036d2 (spec 005 corrected the fault frame, which changed these
    # values from the ddc333ea ones); this checks the reverse path on real faults.
    reverse_expected = np.array(
        [
            [449.99986049755563, 501.0930534680464, -480.302857700829],
            [549.9999717550921, 498.90707222676485, -469.70625323565577],
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


def test_saved_and_loaded_model_evaluates_as_before(tmp_path):
    """Spec edge case: a model saved with to_file and loaded with from_file."""
    import pytest

    pytest.importorskip("dill")
    model = build_model("f1_u0")
    bounds = model.bounding_box
    points = np.random.default_rng(3).random((500, 3)) * (bounds.maximum - bounds.origin)
    points += bounds.origin
    before = model.evaluate_model(points)
    path = tmp_path / "model.pkl"
    model.to_file(str(path))
    loaded = GeologicalModel.from_file(str(path))
    np.testing.assert_array_equal(loaded.evaluate_model(points), before)
    np.testing.assert_array_equal(
        loaded.evaluate_model_gradient(points), model.evaluate_model_gradient(points)
    )


def test_fault_frame_threads_stay_off():
    """Spec edge case: use_threads must give the same results or stay off as today."""
    from LoopStructural.modelling.features.fault import _fault_segment

    assert _fault_segment.use_threads is False


def test_points_of_dtype_object_evaluate_as_float_points():
    """Points from a DataFrame with an object column are evaluated (review round 2)."""
    from tests.fault_evaluation.models import sample_points

    model = build_model("f2_u0")
    xyz = sample_points("f2_u0", 50)
    np.testing.assert_array_equal(model.evaluate_model(xyz.astype(object)), model.evaluate_model(xyz))


def test_array_displacement_evaluates_as_a_float_displacement():
    """A displacement given as a NumPy array is accepted (review round 2)."""
    from tests.fault_evaluation.models import sample_points

    model = build_model("f2_u0")
    xyz = sample_points("f2_u0", 50)
    reference = model.evaluate_model(xyz)
    for fault in model.faults:
        fault.displacement = np.array([fault.displacement])
    np.testing.assert_array_equal(model.evaluate_model(xyz), reference)
