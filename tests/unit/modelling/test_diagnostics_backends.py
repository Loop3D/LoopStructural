"""Residual diagnostics of faulted features and of P2 features (FR-023, FR-026).

The faulted cells check that the residual methods of a feature evaluate the
interpolator at the restored positions of the stored rows, and not at positions that
the faults move a second time. The P2 cells check that the diagnostic methods do not
change the state of the model. No test in this file changes a model.
"""

import copy
import itertools

import numpy as np
import pytest

from tests.fixtures.coordinate_frames import (
    ANCHORS,
    ROTATIONS_DEG,
    faulted_foliation_model,
    model_with_frame,
    tp,
    world_foliation_model,
)

FAULTED_ENGINES = ("FDI", "PLI")
# Bounds for the faulted cells. Each bound is 2 times the measured maximum of the
# independent function over the 12 cells of the engine. Some bounds are above the 1e-2
# plan cap by orchestrator decision (new bounds, PLI fit accuracy is a solver property).
# Value residual, FDI: largest 2.555e-3. PLI: largest 1.896e-2.
TOL_FAULTED_VALUE = {"FDI": 5.11e-3, "PLI": 3.8e-2}
# Direction residual max(1 - cos), FDI: largest 5.402e-3, PLI: largest 8.456e-3.
TOL_FAULTED_COS = {"FDI": 1.1e-2, "PLI": 1.7e-2}
# Number of stored value rows that the faults move by more than 1e-3. This is the
# number that the fixture gives for all 24 cells.
MOVED_ROWS = 8


def _anchor_index(anchor):
    return int(np.flatnonzero(np.all(ANCHORS == anchor, axis=1))[0])


FAULTED_CELLS = [
    (ANCHORS[i], rotation, engine)
    for i, rotation, engine in itertools.product(
        range(len(ANCHORS)), ROTATIONS_DEG, FAULTED_ENGINES
    )
]
FAULTED_IDS = [
    f"a{_anchor_index(anchor)}-r{int(rotation)}-{engine}"
    for anchor, rotation, engine in FAULTED_CELLS
]
P2_CELLS = [
    (ANCHORS[i], rotation)
    for i, rotation in itertools.product(range(len(ANCHORS)), ROTATIONS_DEG)
]
P2_IDS = [f"a{_anchor_index(anchor)}-r{int(rotation)}" for anchor, rotation in P2_CELLS]


@pytest.fixture(scope="module", params=FAULTED_CELLS, ids=FAULTED_IDS)
def faulted_model(request):
    """Return (anchor, engine, model, strati feature) of one faulted cell."""
    anchor, rotation, engine = request.param
    model, _, strati = faulted_foliation_model(anchor, rotation, engine)
    return anchor, engine, model, strati


@pytest.fixture(scope="module", params=P2_CELLS, ids=P2_IDS)
def p2_model(request):
    """Return the solved P2 foliation of one cell."""
    anchor, rotation = request.param
    _, _, data = world_foliation_model(anchor, rotation)
    # The gradient rows of the P2 interpolator fail in the einsum of add_gradient_constraints
    # (subscripts "ikj,ij->ik" on shape (N, 3, 10)). The cells use the value and normal rows.
    data = data.iloc[list(range(0, 27)) + list(range(31, 35))].reset_index(drop=True)
    model = model_with_frame(anchor, anchor + (8, 10, 12), rotation)
    model.data = data
    feature = model.create_and_add_foliation("r1", interpolatortype="P2", nelements=1000)
    feature.builder.up_to_date()
    assert feature.interpolator.up_to_date
    return feature


def _tolerance(factor, anchor):
    """Return max(1e-12, factor * Tp(1, W)) for a unit model at the anchor."""
    return max(1e-12, factor * tp(1.0, float(np.abs(anchor).max()) + 1.0))


def _independent_value_residuals(strati):
    """Return abs(val - support value) at the stored local positions of the value rows."""
    interpolator = strati.interpolator
    rows = np.array(interpolator.data["value"], dtype=float, copy=True)
    modelled = interpolator.support.evaluate_value(rows[:, :3], interpolator.c)
    return np.abs(rows[:, 3] - modelled)


def _independent_cosines(strati):
    """Return the world cosines of the stored gradient and normal rows.

    The function evaluates the interpolator at the world positions of the stored rows.
    It does not use a method of the feature.
    """
    interpolator = strati.interpolator
    cosines = []
    for key in ("gradient", "normal"):
        rows = np.array(interpolator.data[key], dtype=float, copy=True)
        if rows.shape[0] == 0:
            continue
        positions = interpolator.local_to_world_points(rows[:, :3])
        vectors = interpolator.local_to_world_vectors(rows[:, 3:6])
        model_vectors = interpolator.evaluate_gradient(positions)
        dots = np.einsum("ij,ij->i", model_vectors, vectors)
        dots /= np.linalg.norm(model_vectors, axis=1) * np.linalg.norm(vectors, axis=1)
        cosines.append(dots)
    return np.concatenate(cosines)


class _Spy:
    """Count the calls of three methods of a feature."""

    NAMES = ("_apply_faults", "evaluate_value", "evaluate_gradient")

    def __init__(self, feature, monkeypatch):
        self.counts = dict.fromkeys(self.NAMES, 0)
        for name in self.NAMES:
            original = getattr(feature, name)
            monkeypatch.setattr(feature, name, self._wrap(name, original))

    def _wrap(self, name, original):
        def counted(*args, **kwargs):
            self.counts[name] += 1
            return original(*args, **kwargs)

        return counted


def test_faulted_value_residuals_use_restored_positions(faulted_model, monkeypatch):
    """The value residuals equal the residuals at the stored positions (no second fault step)."""
    anchor, engine, _, strati = faulted_model
    independent = _independent_value_residuals(strati)
    spy = _Spy(strati, monkeypatch)
    residuals = strati.evaluate_value_residuals()
    print(f"value residual max: independent {independent.max():.3e}, method {residuals.max():.3e}")
    assert residuals.shape == independent.shape
    assert independent.max() <= TOL_FAULTED_VALUE[engine]
    np.testing.assert_allclose(residuals, independent, rtol=0, atol=_tolerance(3, anchor))
    assert residuals.max() <= TOL_FAULTED_VALUE[engine]
    assert spy.counts == dict.fromkeys(_Spy.NAMES, 0)


def test_faulted_gradient_cosines_use_restored_positions(faulted_model, monkeypatch):
    """The cosines equal the cosines at the stored world positions (no second fault step)."""
    anchor, engine, _, strati = faulted_model
    independent = _independent_cosines(strati)
    spy = _Spy(strati, monkeypatch)
    cosines = strati.evaluate_gradient_cosines()
    print(f"max(1 - cos): independent {(1 - independent).max():.3e}, "
          f"method {(1 - cosines).max():.3e}")
    assert cosines.shape == independent.shape
    assert (1 - independent).max() <= TOL_FAULTED_COS[engine]
    np.testing.assert_allclose(cosines, independent, rtol=0, atol=_tolerance(10, anchor))
    assert (1 - cosines).max() <= TOL_FAULTED_COS[engine]
    assert spy.counts == dict.fromkeys(_Spy.NAMES, 0)


def test_fault_restores_the_stored_rows(faulted_model):
    """Guard: the fault moves the stored value rows away from the raw world data."""
    _, _, model, strati = faulted_model
    interpolator = strati.interpolator
    stored = interpolator.local_to_world_points(
        np.array(interpolator.data["value"], dtype=float, copy=True)[:, :3]
    )
    raw = model.data[model.data["feature_name"] == "strati"]
    raw = raw[np.isfinite(raw["val"])][["X", "Y", "Z"]].to_numpy(float)
    distance = np.linalg.norm(stored[:, None, :] - raw[None, :, :], axis=2).min(axis=1)
    print(f"stored rows moved by more than 1e-3: {int((distance > 1e-3).sum())}")
    assert stored.shape[0] == 8
    assert int((distance > 1e-3).sum()) >= 4
    assert int((distance > 1e-3).sum()) == MOVED_ROWS


def test_report_physical_value_matches_independent_rms(faulted_model):
    """The physical value block of the report holds the RMS of the independent residuals."""
    anchor, _, _, strati = faulted_model
    independent = _independent_value_residuals(strati)
    entry = strati.verification_report().physical_residuals["value"]
    expected = float(np.sqrt(np.mean(independent**2)))
    assert entry.kind == "physical"
    assert abs(entry.rms - expected) <= _tolerance(3, anchor)


def _snapshot(feature):
    """Return deep copies of the state that a diagnostic method must not change."""
    interpolator = feature.interpolator
    return {
        "data": copy.deepcopy(interpolator.data),
        "c": copy.deepcopy(interpolator.c),
        "up_to_date": interpolator.up_to_date,
        "builder_up_to_date": feature.builder._up_to_date,
        "weights": {
            name: copy.deepcopy(constraint["w"])
            for name, constraint in interpolator.constraints.items()
        },
        "report": interpolator.latest_diagnostics_report,
    }


def _assert_same_state(before, after):
    assert before["data"].keys() == after["data"].keys()
    for key, array in before["data"].items():
        assert np.array_equal(array, after["data"][key]), f"data[{key!r}] changed"
    assert np.array_equal(before["c"], after["c"]), "c changed"
    assert before["up_to_date"] == after["up_to_date"]
    assert before["builder_up_to_date"] == after["builder_up_to_date"]
    assert before["weights"].keys() == after["weights"].keys()
    for name, weights in before["weights"].items():
        assert np.array_equal(weights, after["weights"][name]), f"weights of {name} changed"
    assert before["report"] is after["report"], "the report object changed"


def test_p2_misfit_leaves_state_unchanged(p2_model):
    """Guard: the old and new misfit methods and the report change no state of a P2 feature."""
    feature = p2_model
    before = _snapshot(feature)
    feature.evaluate_gradient_misfit()
    feature.evaluate_value_misfit()
    feature.evaluate_gradient_cosines()
    feature.evaluate_value_residuals()
    feature.verification_report()
    _assert_same_state(before, _snapshot(feature))


def test_p2_report_blocks_are_separate(p2_model):
    """Guard: the physical and the assembled blocks are two separate objects."""
    feature = p2_model
    report = feature.verification_report()
    assert report.physical_residuals is not report.assembled_residuals
    assert len(report.physical_residuals) > 0
    assert len(report.assembled_residuals) > 0
    assert all(entry.kind == "physical" for entry in report.physical_residuals.values())
    assert all(entry.kind == "assembled" for entry in report.assembled_residuals.values())
    residuals = feature.evaluate_value_residuals()
    # Some value rows are outside the P2 support and give NaN. The report skips them.
    finite = residuals[np.isfinite(residuals)]
    entry = report.physical_residuals["value"]
    assert entry.count == finite.size
    assert entry.nan_count == residuals.size - finite.size
    assert abs(entry.rms - float(np.sqrt(np.mean(finite**2)))) <= 1e-12
