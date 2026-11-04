"""Reference frame checks for discrete interpolation constraints."""

import itertools

import numpy as np
import pytest
from loop_common.supports import StructuredGrid, TetMesh

from LoopStructural.modelling.features.builders import GeologicalFeatureBuilder
from tests.fixtures.coordinate_frames import (
    ANCHORS,
    ROTATIONS_DEG,
    AnalyticReference,
    model_with_frame,
    r3_gradient,
    r3_value,
    tp,
)


@pytest.fixture(params=list(itertools.product(range(len(ANCHORS)), ROTATIONS_DEG)))
def frame_case(request):
    anchor, rotation = ANCHORS[request.param[0]], request.param[1]
    model = model_with_frame(anchor, anchor + (8, 10, 12), rotation)
    return anchor, rotation, model


def _builder(engine, model):
    builder = GeologicalFeatureBuilder(
        interpolatortype=engine,
        bounding_box=model.bounding_box,
        nelements=210,
        name="r3",
        model=model,
        buffer=0,
    )
    corners = np.array(
        list(itertools.product(*zip(model.bounding_box.origin, model.bounding_box.maximum)))
    )
    local = model.bounding_box.project(corners)
    lo, hi = local.min(axis=0), local.max(axis=0)
    support_type = StructuredGrid if engine == "FDI" else TetMesh
    interp = type(builder.interpolator)(
        support_type(origin=lo, nsteps_cells=[5, 6, 7], step_vector=(hi - lo) / [5, 6, 7])
    )
    interp.bounding_box = model.bounding_box
    builder._interpolator = interp
    return builder


def _field(anchor):
    return AnalyticReference(
        lambda p: r3_value(p, anchor), lambda p: r3_gradient(p, anchor), "reference"
    )


@pytest.mark.parametrize("engine", ["FDI", "PLI"])
def test_world_local_conversion_methods(frame_case, engine):
    anchor, _, model = frame_case
    interp = _builder(engine, model).interpolator
    points = anchor + np.array([[1.2, 2.3, 3.4], [4.1, 5.2, 6.3]])
    vectors = np.array([[1.0, 2.0, 3.0], [-2.0, 1.0, 4.0]])
    for to_local, to_world, values in (
        (interp.world_to_local_points, interp.local_to_world_points, points),
        (interp.world_to_local_vectors, interp.local_to_world_vectors, vectors),
    ):
        np.testing.assert_allclose(to_world(to_local(values)), values, atol=tp(12, 12), rtol=0)


@pytest.mark.parametrize("engine", ["FDI", "PLI"])
def test_equality_region_selects_world_nodes(frame_case, engine):
    anchor, _, model = frame_case
    builder = _builder(engine, model)
    nodes = builder.interpolator.support.nodes
    world = model.bounding_box.reproject(nodes)
    expected = np.flatnonzero(world[:, 0] > anchor[0] + 4)
    field = _field(anchor)
    builder.add_equality_constraints(field, lambda p: p[:, 0] > anchor[0] + 4)
    builder.install_equality_constraints()
    entries = builder.interpolator.equal_constraints
    actual = np.concatenate([entry["col"] for entry in entries.values()])
    np.testing.assert_array_equal(np.sort(actual), expected)
    np.testing.assert_allclose(
        np.concatenate([entry["B"] for entry in entries.values()]),
        field.evaluate_value(world[expected]),
        rtol=0,
        atol=tp(50, 50),
    )
    assert "reference_equality" in entries


@pytest.mark.parametrize("engine", ["FDI", "PLI"])
@pytest.mark.parametrize("sample_rule", [None, "centroid", "gauss4"])
def test_orthogonal_rows_affine_residual(frame_case, engine, sample_rule):
    anchor, _, model = frame_case
    builder = _builder(engine, model)
    for name, world_direction, target in (
        ("x", np.array([1.0, 0.0, 0.0]), 0.0),
        ("y", np.array([0.0, 1.0, 0.0]), 1.0),
    ):
        reference = AnalyticReference(
            lambda p, offset=anchor[1]: p[:, 1] - offset,
            lambda p, direction=world_direction: np.broadcast_to(direction, (len(p), 3)).copy(),
            name,
        )
        builder.add_orthogonal_feature(reference, B=target, sample_rule=sample_rule)
    builder.install_gradient_constraint()
    support = builder.interpolator.support
    exact = model.bounding_box.reproject(support.nodes)[:, 1] - anchor[1]
    for name in ("x_orthogonal", "y_orthogonal"):
        constraint = builder.interpolator.constraints[name]
        residual = constraint["matrix"] @ exact - constraint["b"]
        assert len(residual) > 0
        np.testing.assert_allclose(residual, 0.0, rtol=0, atol=1e-10)


@pytest.mark.parametrize("engine", ["FDI", "PLI"])
@pytest.mark.parametrize("sample_rule", [None, "centroid", "gauss4"])
def test_reference_gradient_samples_use_world_points_and_local_vectors(
    frame_case, engine, sample_rule
):
    anchor, _, model = frame_case
    builder = _builder(engine, model)
    reference = _field(anchor)
    seen = {}
    original_eval = reference.evaluate_gradient

    def record_gradient(points, ignore_regions=True):
        seen["points"] = np.asarray(points).copy()
        return original_eval(points, ignore_regions=ignore_regions)

    reference.evaluate_gradient = record_gradient
    original_add = builder.interpolator.add_gradient_orthogonal_constraints

    def record_add(points, vector, **kwargs):
        seen["vectors"] = np.asarray(vector).copy()
        return original_add(points, vector, **kwargs)

    builder.interpolator.add_gradient_orthogonal_constraints = record_add
    builder.add_orthogonal_feature(reference, sample_rule=sample_rule)
    builder.install_gradient_constraint()
    assert "points" in seen
    support = builder.interpolator.support
    if sample_rule is None:
        expected_local = support.barycentre
    else:
        from loop_interpolation._quadrature import tet_quadrature_points

        expected_local, _, _ = tet_quadrature_points(support.nodes, support.elements, sample_rule)
    expected_world = model.bounding_box.reproject(expected_local)
    np.testing.assert_allclose(seen["points"], expected_world, rtol=0, atol=1e-10)
    world = expected_world
    expected = r3_gradient(world, anchor)
    expected /= np.linalg.norm(expected, axis=1)[:, None]
    expected = model.bounding_box.project_vectors(expected)
    np.testing.assert_allclose(seen["vectors"], expected, rtol=0, atol=1e-10)


@pytest.mark.parametrize("engine", ["FDI", "PLI"])
@pytest.mark.parametrize("lower", [True, False])
def test_inequality_feature_bounds_world(frame_case, engine, lower):
    anchor, _, model = frame_case
    builder = _builder(engine, model)
    builder.interpolator.add_inequality_feature(_field(anchor), lower=lower)
    bounds = builder.interpolator.ineq_constraints["inequality_feature"]["bounds"]
    world = model.bounding_box.reproject(builder.interpolator.support.nodes)
    values = r3_value(world, anchor)
    column = 0 if lower else 1
    np.testing.assert_allclose(bounds[:, column], values, rtol=0, atol=tp(200, 200))


@pytest.mark.parametrize("engine", ["FDI", "PLI"])
def test_restricted_reference_exclusion_counts(frame_case, engine):
    anchor, _, model = frame_case
    builder = _builder(engine, model)
    nodes = builder.interpolator.support.nodes
    world = model.bounding_box.reproject(nodes)
    lower, upper = anchor + (2.1, 2.1, 2.1), anchor + (5.9, 7.9, 9.9)
    ref = _field(anchor)
    ref.name = "restricted"
    ref.support_box = (lower, upper)
    selected = np.all((world >= lower) & (world <= upper), axis=1)
    builder.add_equality_constraints(ref, lambda p: np.ones(len(p), dtype=bool))
    builder.install_equality_constraints()
    equality = builder.interpolator.equal_constraints
    equality_rows = np.concatenate([entry["col"] for entry in equality.values()])
    np.testing.assert_array_equal(np.sort(equality_rows), np.flatnonzero(selected))
    assert len(equality_rows) == int(selected.sum())
    assert "restricted_equality" in equality

    builder.add_orthogonal_feature(ref, B=0.0, sample_rule=None)
    builder.install_gradient_constraint()
    sample_local = builder.interpolator.support.barycentre
    sample_world = model.bounding_box.reproject(sample_local)
    valid_samples = np.all((sample_world >= lower) & (sample_world <= upper), axis=1)
    orthogonal = builder.interpolator.constraints["restricted_orthogonal"]
    assert orthogonal["matrix"].shape[0] == int(valid_samples.sum())
    excluded = builder.excluded_rows
    assert excluded["restricted_equality"].get("invalid_reference_value", 0) == int(
        (~selected).sum()
    )
    assert excluded["restricted_orthogonal"].get("invalid_reference_direction", 0) == int(
        (~valid_samples).sum()
    )


@pytest.mark.parametrize("engine", ["FDI", "PLI"])
def test_support_expansion_contains_restored_points(frame_case, engine):
    _, _, model = frame_case
    builder = _builder(engine, model)

    class Restorer:
        def __init__(self):
            self.received = None

        def apply_to_points(self, points):
            self.received = points.copy()
            return points + np.array([20.0, -3.0, 2.0])

    restorer = Restorer()
    builder.faults.append(restorer)
    model.bounding_box.nsteps = [2, 2, 2]
    expected = model.bounding_box.regular_grid(local=False)
    builder.check_interpolation_geometry({}, buffer=0)
    np.testing.assert_allclose(restorer.received, expected, rtol=0, atol=tp(50, 50))
    restored = expected + [20.0, -3.0, 2.0]
    local = model.bounding_box.project(restored)
    support = builder.interpolator.support
    assert np.all(local >= support.origin - 1e-10)
    assert np.all(local <= support.maximum + 1e-10)


@pytest.mark.parametrize("engine", ["FDI", "PLI"])
def test_set_interpolation_geometry_encloses_rotated_box(frame_case, engine):
    anchor, rotation, model = frame_case
    builder = _builder(engine, model)
    corners = np.array(list(itertools.product(*zip(anchor, anchor + (8, 10, 12)))))
    local = model.bounding_box.project(corners)
    builder.set_interpolation_geometry(anchor, anchor + (8, 10, 12), rotation=rotation)
    support = builder.interpolator.support
    assert np.all(local >= support.origin - 1e-10)
    assert np.all(local <= support.maximum + 1e-10)
