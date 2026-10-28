"""Second fold frame coordinate on a similar fold with a known answer.

See ``scripts/validation/fold_frame_coordinate1.py`` for the geometry, the
metrics and the full comparison table.
"""

from functools import lru_cache

import numpy as np
import pytest
from loop_interpolation import P2Interpolator, P2UnstructuredTetMesh

from scripts.validation.fold_frame_coordinate1 import (
    BOX,
    EVAL_LO,
    EVAL_HI,
    Ramp,
    SimilarFold,
    _Feature,
    build_coordinate1_with_builder,
    evaluate,
    sample_data,
    solve_coordinate1,
)

SEEDS = (0, 1, 2)


@pytest.fixture(scope="module")
def fold():
    return SimilarFold()


@lru_cache(maxsize=None)
def _metrics(method, constraints, n_data, seed):
    fold = SimilarFold()
    interpolator = solve_coordinate1(fold, method, sample_data(fold, n_data, seed), constraints)
    return evaluate(fold, interpolator)


def _mean(key, method, constraints, n_data):
    return float(np.mean([_metrics(method, constraints, n_data, s)[key] for s in SEEDS]))


def test_truth_is_orthogonal_and_monotonic(fold):
    points = np.random.default_rng(3).uniform(EVAL_LO, EVAL_HI, (300, 3))
    _, _, rate, normal, tangent = fold.truth(points)
    h = 1e-3
    c1 = lambda p: fold.truth(p)[1]  # noqa: E731
    along_normal = (c1(points + h * normal) - c1(points - h * normal)) / (2 * h)
    along_tangent = (c1(points + h * tangent) - c1(points - h * tangent)) / (2 * h)
    # c1 is constant along the normals of the folded surfaces and increases along them
    np.testing.assert_allclose(along_normal, 0.0, atol=1e-3)
    np.testing.assert_allclose(along_tangent, rate, rtol=1e-2)
    assert np.all(rate > 0)
    for offset in (-1.0, 0.0, 1.0):
        assert np.all(np.diff(c1(fold.isosurface(offset))) > 0)
    # no focal curves: the truth is defined everywhere in the mesh
    assert fold.valid(np.random.default_rng(4).uniform(*BOX, (300, 3))).all()


def test_truth_c0_is_the_profile_shifted_in_z(fold):
    for offset in (-1.0, 0.0, 1.5):
        np.testing.assert_allclose(fold.truth(fold.isosurface(offset))[0], offset, atol=1e-12)
    # c1 is the distance along the reference isosurface
    reference = fold.isosurface(0.0, n=2000)
    distance = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(reference, axis=0), axis=1))]
    c1 = fold.truth(reference)[1]
    np.testing.assert_allclose(c1 - c1[0], distance, atol=2e-3)


def test_ramp_has_no_skill():
    m = evaluate(SimilarFold(), Ramp())
    assert m["skill"] == pytest.approx(0.0, abs=1e-9)
    # a field that ignores the folds is monotonic, so monotonicity cannot rank methods
    assert m["reversal"] == 0.0 and m["iso_reversal"] == 0.0
    assert m["cos_med"] > 0.4


def test_gauss4_does_best_when_the_gradient_magnitude_is_fixed():
    # orthogonality fixes the isosurfaces of c1 and the norm along the tangent fixes the
    # spacing between them. The oracle target is exact, so this is the best case.
    for n_data in (4, 12, 40):
        gauss4 = _mean("skill", "P2-gauss4", "oracle norm 1", n_data)
        assert gauss4 > 0.7
        assert _mean("stretch_bad", "P2-gauss4", "oracle norm 1", n_data) < 0.03
        assert _mean("cos_med", "P2-gauss4", "oracle norm 1", n_data) < 0.1
        for other in ("P1-step10", "P2-centroid"):
            assert gauss4 > 0.1 + _mean("skill", other, "oracle norm 1", n_data)


def test_gauss4_stretches_less_erratically_than_centroid_sampling():
    # P2 has a linear gradient inside each tetrahedron, so one sample per tet leaves
    # most of the gradient free and the rate of change along the surface is erratic
    for n_data in (12, 40):
        centroid = _mean("stretch_bad", "P2-centroid", "oriented", n_data)
        gauss4 = _mean("stretch_bad", "P2-gauss4", "oriented", n_data)
        assert gauss4 < 0.5 * centroid, (n_data, gauss4, centroid)


def test_gauss4_beats_p1_with_a_flat_norm_target_and_enough_data():
    for n_data in (12, 40):
        gauss4 = _mean("skill", "P2-gauss4", "oriented", n_data)
        assert gauss4 > 0.1 + _mean("skill", "P1-step10", "oriented", n_data)


def test_gauss4_keeps_c1_constant_along_the_normals_of_the_folded_surfaces():
    for n_data in (4, 12, 40):
        gauss4 = _mean("cos_med", "P2-gauss4", "orthogonal", n_data)
        assert gauss4 < 0.05
        assert gauss4 < _mean("cos_med", "P1-step10", "orthogonal", n_data)
        assert gauss4 < _mean("cos_med", "P2-centroid", "orthogonal", n_data)


@pytest.fixture
def small_p2():
    mesh = P2UnstructuredTetMesh(
        origin=np.zeros(3), step_vector=np.ones(3), nsteps_cells=np.array([2, 2, 2])
    )
    interpolator = P2Interpolator(mesh)
    interpolator.set_value_constraints(np.array([[0.5, 0.5, 0.5, 0.0], [1.5, 1.5, 1.5, 1.0]]))
    interpolator.setup_interpolator()
    interpolator.solve_system("cg")
    return interpolator


# GeologicalFeatureBuilder needs both of these to build the second and third coordinate
# of a frame on P2.
def test_p2_orthogonal_constraint_accepts_builder_arguments(small_p2):
    small_p2.add_gradient_orthogonal_constraints(
        np.array([[1.0, 1.0, 1.0]]), np.array([[1.0, 0.0, 0.0]]), w=1.0, b=0.0, name="orthogonal"
    )
    assert "orthogonal" in small_p2.constraints


def test_p2_gradient_is_nan_outside_the_mesh(small_p2):
    gradient = small_p2.evaluate_gradient(np.array([[1.0, 1.0, 1.0], [10.0, 10.0, 10.0]]))
    assert np.isfinite(gradient[0]).all()
    assert np.isnan(gradient[1]).all()


@pytest.mark.parametrize("rule,points_per_tet", [(None, 1), ("centroid", 1), ("gauss4", 4)])
def test_frame_builder_samples_orthogonality_at_quadrature_points(rule, points_per_tet):
    import pandas as pd

    from LoopStructural import GeologicalModel
    from LoopStructural.datasets import load_noddy_single_fold

    data, bounds = load_noddy_single_fold()
    model = GeologicalModel(bounds[0], bounds[1])
    model.set_model_data(pd.concat([data[:100], data[data["feature_name"] == "s1"]]))
    frame = model.create_and_add_fold_frame(
        "s1", nelements=500, interpolatortype="P2", orthogonal_sample_rule=rule
    )
    xyz = np.random.default_rng(0).uniform(bounds[0] + 1, bounds[1] - 1, (100, 3))
    assert np.isfinite(frame[1].evaluate_value(xyz)).all()
    interpolator = frame.builder.builders[1].interpolator
    constraint = interpolator.constraints["s1__0_orthogonal"]
    n_tets = len(interpolator.support.elements[::10])
    # points outside the first coordinate's mesh have no gradient and are dropped
    assert 0.3 * n_tets * points_per_tet < constraint["matrix"].shape[0] <= n_tets * points_per_tet
    expected = 0.5 if rule == "gauss4" else 1.0
    np.testing.assert_allclose(constraint["w"], expected * np.ones_like(constraint["w"]))


@lru_cache(maxsize=None)
def _builder_metrics(method, n_data, seed):
    fold = SimilarFold()
    data = sample_data(fold, n_data, seed)
    return evaluate(fold, _Feature(build_coordinate1_with_builder(fold, method, data)))


def _builder_mean(key, method, n_data):
    return float(np.mean([_builder_metrics(method, n_data, s)[key] for s in SEEDS]))


def test_builder_default_second_coordinate_is_smooth_and_monotonic():
    # the reference the P2 options are compared with
    for n_data in (12, 40):
        assert _builder_mean("reversal", "FDI step10 (default)", n_data) == 0.0
        assert _builder_mean("stretch_bad", "FDI step10 (default)", n_data) < 0.1
        assert _builder_mean("skill", "FDI step10 (default)", n_data) > 0.1


def test_builder_gauss4_improves_p2_over_centroid_sampling():
    for n_data in (12, 40):
        gauss4 = _builder_mean("skill", "P2 gauss4 step1", n_data)
        centroid = _builder_mean("skill", "P2 centroid step1", n_data)
        assert gauss4 > 0.1 + centroid, (n_data, gauss4, centroid)


def test_builder_gauss4_improves_p2_over_barycentre_sampling_at_the_default_step():
    for n_data in (12, 40):
        gauss4 = _builder_mean("stretch_bad", "P2 gauss4 step10", n_data)
        barycentre = _builder_mean("stretch_bad", "P2 barycentre step10", n_data)
        assert gauss4 < barycentre, (n_data, gauss4, barycentre)


def test_builder_gauss4_follows_the_folded_surfaces_best():
    for n_data in (4, 12, 40):
        gauss4 = _builder_mean("cos_med", "P2 gauss4 step1", n_data)
        assert gauss4 < 0.5 * _builder_mean("cos_med", "FDI step10 (default)", n_data)
        assert gauss4 < _builder_mean("cos_med", "P2 barycentre step10", n_data)
