"""Frame checks for fold constraints, directional regularisation and P2 weights.

The fold provider, the ``mask_fn`` callback, a callable regularisation
direction and ``steepness_wtfunc`` work in world coordinates. The matrix rows
work in local coordinates. These tests check that boundary.
"""

import itertools

import numpy as np
import pytest
from loop_common.supports import P2TetMesh, StructuredGrid, TetMesh
from loop_interpolation import (
    DiscreteFoldInterpolator,
    FDFoldInterpolator,
    FiniteDifferenceInterpolator,
    P1Interpolator,
    P2FoldInterpolator,
    P2Interpolator,
)
from loop_interpolation import _discrete_fold_interpolator as pli_module
from loop_interpolation import _fd_fold_interpolator as fdi_module
from loop_interpolation import _p2_fold_interpolator as p2_module
from loop_interpolation._quadrature import tet_quadrature_points

from tests.fixtures.coordinate_frames import (
    ANCHORS,
    ROTATIONS_DEG,
    AnalyticFold,
    model_with_frame,
    tp,
)

FOLD_FAMILIES = ("fdi_fold", "pli_fold", "p2_fold_centroid", "p2_fold_gauss4")
REGULARISATION_FAMILIES = ("fdi", "pli", "p2")
_FOLD_CLASSES = {
    "fdi_fold": (FDFoldInterpolator, StructuredGrid),
    "pli_fold": (DiscreteFoldInterpolator, TetMesh),
    "p2_fold_centroid": (P2FoldInterpolator, P2TetMesh),
    "p2_fold_gauss4": (P2FoldInterpolator, P2TetMesh),
}
_PLAIN_CLASSES = {
    "fdi": (FiniteDifferenceInterpolator, StructuredGrid),
    "pli": (P1Interpolator, TetMesh),
    "p2": (P2Interpolator, P2TetMesh),
}
_NORM_MODULES = {
    "fdi_fold": fdi_module,
    "pli_fold": pli_module,
    "p2_fold_centroid": p2_module,
    "p2_fold_gauss4": p2_module,
}
FOLD_DIRECTION_ROWS = ("fold orientation", "fold axis", "fold normalisation")


@pytest.fixture(params=list(itertools.product(range(len(ANCHORS)), ROTATIONS_DEG)))
def frame_case(request):
    anchor, rotation = ANCHORS[request.param[0]], request.param[1]
    model = model_with_frame(anchor, anchor + (8, 10, 12), rotation)
    return anchor, rotation, model


def _interpolator(classes, family, model):
    """Build an interpolator on the local enclosure of the model box."""
    interpolator_class, support_class = classes[family]
    corners = np.array(
        list(itertools.product(*zip(model.bounding_box.origin, model.bounding_box.maximum)))
    )
    local = model.bounding_box.project(corners)
    lo, hi = local.min(axis=0), local.max(axis=0)
    step = (hi - lo) / [5, 6, 7]
    if support_class is P2TetMesh:
        support = P2TetMesh(origin=lo, nsteps=[5, 6, 7], step_vector=step)
    else:
        support = support_class(origin=lo, nsteps_cells=[5, 6, 7], step_vector=step)
    interpolator = interpolator_class(support)
    interpolator.bounding_box = model.bounding_box
    return interpolator


def _fold_interpolator(family, model, fold):
    interpolator = _interpolator(_FOLD_CLASSES, family, model)
    interpolator.fold = fold
    return interpolator


def _sample_rule(family):
    return family.replace("p2_fold_", "")


def _local_samples(interpolator, family):
    """Return the local points where the family asks the fold for directions."""
    support = interpolator.support
    if family == "fdi_fold":
        return support.nodes
    if family == "pli_fold":
        return support.barycentre
    points, _, _ = tet_quadrature_points(support.nodes, support.elements, _sample_rule(family))
    return points


def _add_fold(interpolator, family, **kwargs):
    """Call ``add_fold_constraints`` with the arguments that the family accepts."""
    if family.startswith("p2_fold"):
        kwargs["sample_rule"] = _sample_rule(family)
    else:
        kwargs.setdefault("fold_regularisation", None)
    interpolator.add_fold_constraints(**kwargs)


def _record_rows(monkeypatch, interpolator):
    """Record the rows before normalisation of each least squares call."""
    calls = []
    original = interpolator.add_constraints_to_least_squares

    def spy(A, B, idc, w=1.0, name="undefined"):
        calls.append({"A": np.array(A, dtype=float), "w": np.array(w, dtype=float), "name": name})
        return original(A, B, idc, w=w, name=name)

    monkeypatch.setattr(interpolator, "add_constraints_to_least_squares", spy)
    return calls


def _record_orthogonal(monkeypatch, interpolator):
    """Record the points and vectors of each gradient orthogonal call."""
    calls = []
    original = interpolator.add_gradient_orthogonal_constraints

    def spy(points, vectors, *args, **kwargs):
        calls.append(
            {
                "points": np.array(points, dtype=float),
                "vectors": np.array(vectors, dtype=float),
                "name": kwargs.get("name"),
            }
        )
        return original(points, vectors, *args, **kwargs)

    monkeypatch.setattr(interpolator, "add_gradient_orthogonal_constraints", spy)
    return calls


def _set_world_normal(interpolator, world_point, world_normal):
    """Store one normal constraint from a world row. The setter projects it to local."""
    point = np.asarray(world_point, dtype=float)[None]
    normal = np.asarray(world_normal, dtype=float)[None]
    interpolator.set_normal_constraints(np.hstack([point, normal, [[1.0]]]))


@pytest.mark.parametrize("family", FOLD_FAMILIES)
def test_fold_geometry_sampled_at_world_points(frame_case, family):
    anchor, _, model = frame_case
    fold = AnalyticFold.r4(anchor)
    interpolator = _fold_interpolator(family, model, fold)
    world_point = anchor + np.array([3.0, 4.0, 5.0])
    _set_world_normal(interpolator, world_point, (0.0, -1.0, 0.0))
    _add_fold(interpolator, family, dgz_alignment="warn")
    expected = interpolator.local_to_world_points(_local_samples(interpolator, family))
    # The first call is for the samples. The second call is for the normal check.
    assert len(fold.seen_points) == 2
    np.testing.assert_allclose(fold.seen_points[0], expected, rtol=0, atol=tp(50, 50))
    np.testing.assert_allclose(fold.seen_points[1], world_point[None], rtol=0, atol=tp(50, 50))


@pytest.mark.parametrize("family", FOLD_FAMILIES)
def test_fold_rows_affine_residual(frame_case, family):
    anchor, _, model = frame_case
    fold = AnalyticFold.constant((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))
    interpolator = _fold_interpolator(family, model, fold)
    _add_fold(interpolator, family, fold_norm=1.0, dgz_alignment="none")
    bounding_box = model.bounding_box
    # Use the rigid structure to avoid cancellation at large world coordinates.
    exact = bounding_box.reproject_vectors(interpolator.support.nodes)[:, 1] + (
        bounding_box.local_origin[1] - anchor[1]
    )
    for name in FOLD_DIRECTION_ROWS:
        constraint = interpolator.constraints[name]
        residual = constraint["matrix"] @ exact - constraint["b"]
        assert len(residual) > 0
        np.testing.assert_allclose(residual, 0.0, rtol=0, atol=1e-10)


def _expected_local_directions(interpolator, fold, local_samples):
    world = interpolator.local_to_world_points(local_samples)
    return {
        "fold orientation": interpolator.world_to_local_vectors(fold.deformed_fn(world)),
        "fold axis": interpolator.world_to_local_vectors(fold.axis_fn(world)),
        "fold normalisation": interpolator.world_to_local_vectors(fold.dgz_fn(world)),
    }


@pytest.mark.parametrize("family", FOLD_FAMILIES[:1])
def test_fold_rows_follow_curved_directions_fdi(frame_case, family, monkeypatch):
    anchor, _, model = frame_case
    fold = AnalyticFold.r4(anchor)
    interpolator = _fold_interpolator(family, model, fold)
    calls = _record_orthogonal(monkeypatch, interpolator)
    _add_fold(interpolator, family, dgz_alignment="none")
    expected = _expected_local_directions(interpolator, fold, interpolator.support.nodes)
    assert {call["name"] for call in calls} == set(FOLD_DIRECTION_ROWS)
    for call in calls:
        for column in range(3):
            np.testing.assert_allclose(
                call["vectors"][:, column], expected[call["name"]][:, column], rtol=0, atol=1e-10
            )


@pytest.mark.parametrize("family", FOLD_FAMILIES[1:])
def test_fold_rows_follow_curved_directions_pli_p2(frame_case, family, monkeypatch):
    anchor, _, model = frame_case
    fold = AnalyticFold.r4(anchor)
    interpolator = _fold_interpolator(family, model, fold)
    calls = _record_rows(monkeypatch, interpolator)
    _add_fold(interpolator, family, dgz_alignment="none")
    support = interpolator.support
    local_samples = _local_samples(interpolator, family)
    direction = _expected_local_directions(interpolator, fold, local_samples)["fold orientation"]
    if family == "pli_fold":
        gradients = support.get_element_gradients(np.arange(support.n_elements))
        nodes = support.nodes[support.get_elements()[np.arange(support.n_elements)]]
        volume = np.abs(np.linalg.det(nodes[:, 1:, :] - nodes[:, 0, None, :])) / 6
        expected = volume[:, None] * np.einsum("ij,ijk->ik", direction, gradients)
    else:
        _, tet_index, _ = tet_quadrature_points(
            support.nodes, support.elements, _sample_rule(family)
        )
        shape_derivatives, _, inside = support.evaluate_shape_derivatives(
            local_samples, elements=tet_index
        )
        assert np.all(inside)
        expected = np.einsum("ij,ijk->ik", direction, shape_derivatives)
    actual = [call["A"] for call in calls if call["name"] == "fold orientation"]
    assert len(actual) == 1
    assert actual[0].shape == expected.shape
    np.testing.assert_allclose(actual[0], expected, rtol=0, atol=1e-10)


@pytest.mark.parametrize("family", FOLD_FAMILIES)
def test_fold_mask_receives_world_points(frame_case, family, monkeypatch):
    anchor, _, model = frame_case
    fold = AnalyticFold.r4(anchor)
    interpolator = _fold_interpolator(family, model, fold)
    received = []

    def mask_fn(points):
        received.append(np.array(points, dtype=float))
        return points[:, 0] - anchor[0] > 4.0

    if family == "fdi_fold":
        orthogonal = _record_orthogonal(monkeypatch, interpolator)
    else:
        rows = _record_rows(monkeypatch, interpolator)
    _add_fold(interpolator, family, dgz_alignment="none", mask_fn=mask_fn)
    expected_world = interpolator.local_to_world_points(_local_samples(interpolator, family))
    assert len(received) == 1
    np.testing.assert_allclose(received[0], expected_world, rtol=0, atol=tp(50, 50))
    analytic_excluded = int(np.count_nonzero(expected_world[:, 0] - anchor[0] > 4.0))
    assert 0 < analytic_excluded < len(expected_world)
    if family == "fdi_fold":
        # The grid interpolator drops the excluded nodes from the rows.
        for call in orthogonal:
            assert len(expected_world) - len(call["points"]) == analytic_excluded
    else:
        # The tetrahedral interpolators keep the rows with zero weight.
        for call in rows:
            assert int(np.count_nonzero(call["w"] == 0.0)) == analytic_excluded


@pytest.mark.parametrize("family", FOLD_FAMILIES)
def test_fold_norm_alignment_uses_world_normals(frame_case, family, monkeypatch):
    anchor, _, model = frame_case
    fold = AnalyticFold.constant((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))
    interpolator = _fold_interpolator(family, model, fold)
    _set_world_normal(interpolator, anchor + (3.0, 4.0, 5.0), (0.0, -1.0, 0.0))
    module = _NORM_MODULES[family]
    original = module.resolve_fold_norm_target
    targets = []

    def spy(**kwargs):
        targets.append(original(**kwargs))
        return targets[-1]

    monkeypatch.setattr(module, "resolve_fold_norm_target", spy)
    _add_fold(interpolator, family, fold_norm=1.0, dgz_alignment="correct")
    assert targets == [-1.0]


def _world_direction_function(anchor, received):
    """Return a world direction function that records its input points."""

    def direction(points):
        received.append(np.array(points, dtype=float))
        local = np.asarray(points, dtype=float) - anchor
        vectors = np.column_stack(
            (1.0 + 0.05 * local[:, 0], 0.3 + 0.02 * local[:, 1], 0.1 * local[:, 2] - 0.4)
        )
        return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)

    return direction


@pytest.mark.parametrize("use_array", [False, True], ids=["callable", "array"])
@pytest.mark.parametrize("family", REGULARISATION_FAMILIES)
def test_directional_regularisation_world_contract(frame_case, family, use_array, monkeypatch):
    anchor, _, model = frame_case
    interpolator = _interpolator(_PLAIN_CLASSES, family, model)
    received = []
    direction = _world_direction_function(anchor, received)
    if family == "p2":
        # P2 has no regularisation sample points. It rejects the setting.
        term = {"direction": direction, "weight": 1.0}
        with pytest.raises(NotImplementedError):
            interpolator.setup_interpolator(directional_regularisation=[term])
        return
    vectors_seen = []
    monkeypatch.setattr(
        interpolator,
        "_add_directional_regularisation",
        lambda weight, vectors, name="directional regularisation": vectors_seen.append(
            np.array(vectors, dtype=float)
        ),
    )
    samples_world = interpolator.local_to_world_points(
        interpolator.get_regularisation_sample_points()
    )
    reference = direction(samples_world)
    received.clear()
    term = {"direction": reference if use_array else direction, "weight": 1.0}
    interpolator.setup_interpolator(directional_regularisation=[term])
    if not use_array:
        assert len(received) == 1
        np.testing.assert_allclose(received[0], samples_world, rtol=0, atol=tp(50, 50))
    assert len(vectors_seen) == 1
    np.testing.assert_allclose(
        vectors_seen[0], interpolator.world_to_local_vectors(reference), rtol=0, atol=1e-10
    )


@pytest.mark.parametrize("use_mask", [False, True], ids=["nomask", "mask"])
@pytest.mark.parametrize("family", FOLD_FAMILIES[:2])
def test_fold_regularisation_world_contract(frame_case, family, use_mask, monkeypatch):
    """The fold regularisation gets local vectors that come from world fold vectors."""
    anchor, _, model = frame_case
    fold = AnalyticFold.r4(anchor)
    interpolator = _fold_interpolator(family, model, fold)
    recorded = {}

    def spy(weight, vectors, name="directional regularisation"):
        recorded[name] = np.array(vectors, dtype=float)

    monkeypatch.setattr(interpolator, "_add_directional_regularisation", spy)

    def mask_fn(points):
        return points[:, 0] - anchor[0] > 4.0

    kwargs = {"mask_fn": mask_fn} if use_mask else {}
    _add_fold(
        interpolator,
        family,
        fold_regularisation=[1.0, 1.0, 1.0],
        dgz_alignment="none",
        **kwargs,
    )
    # FDI samples the nodes. PLI samples the regularisation sample points.
    if family == "fdi_fold":
        local_samples = interpolator.support.nodes
    else:
        local_samples = interpolator.get_regularisation_sample_points()
    world = interpolator.local_to_world_points(local_samples)
    excluded = mask_fn(world) if use_mask else np.zeros(len(world), dtype=bool)
    if use_mask:
        assert 0 < np.count_nonzero(excluded) < len(world)
    components = {
        "fold regularisation 1": fold.dgz_fn,
        "fold regularisation 2": fold.deformed_fn,
        "fold regularisation 3": fold.axis_fn,
    }
    assert set(recorded) == set(components)
    for name, function in components.items():
        world_vectors = function(world).copy()
        world_vectors[excluded] = 0.0
        expected = interpolator.world_to_local_vectors(world_vectors)
        assert recorded[name].shape == expected.shape
        for column in range(3):
            np.testing.assert_allclose(
                recorded[name][:, column], expected[:, column], rtol=0, atol=1e-10
            )


def test_p2_steepness_wtfunc_receives_world_barycentres(frame_case):
    _, _, model = frame_case
    interpolator = _interpolator(_PLAIN_CLASSES, "p2", model)
    support = interpolator.support
    received = []

    def wtfunc(points):
        received.append(np.array(points, dtype=float))
        return np.ones(len(points))

    weighted = np.all(support.neighbours >= 0, axis=1)
    assert 0 < np.count_nonzero(weighted) < len(weighted)
    error = None
    try:
        interpolator.minimise_grad_steepness(w=0.1, wtfunc=wtfunc)
    except ValueError as exc:  # The old code gives one weight for each element.
        error = exc
    assert len(received) == 1
    assert received[0].shape == (int(np.count_nonzero(weighted)), 3)
    expected = interpolator.local_to_world_points(support.barycentre[weighted])
    np.testing.assert_allclose(received[0], expected, rtol=0, atol=tp(50, 50))
    assert error is None
