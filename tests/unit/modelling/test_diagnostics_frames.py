"""Frame checks for the misfit methods of a geological feature (U4.4; FR-023; SC-009).

The data of the foliation are world coordinates. The interpolator stores them in
its local frame. The misfit methods must compare world vectors with world
vectors, and they must not change any state of the model.
"""

import copy
import itertools
import warnings

import numpy as np
import pandas as pd
import pytest

from tests.fixtures.coordinate_frames import (
    ANCHORS,
    R1_GRADIENT,
    ROTATIONS_DEG,
    model_with_frame,
    world_foliation_model,
)

ENGINES = ("FDI", "PLI")
R1_UNIT = R1_GRADIENT / np.linalg.norm(R1_GRADIENT)
# Bound of |1 - cosine| for the real solved model. The independent measurement
# over the 24 cells gave a maximum of 4.41e-4. The bound is 2 times this value
# (8.8e-4), rounded up to one significant digit.
TOL_SOLVED = 9e-4


def _cell_id(cell):
    anchor, rotation, engine = cell
    index = int(np.flatnonzero(np.all(ANCHORS == anchor, axis=1))[0])
    return f"a{index}-r{int(rotation)}-{engine}"


CELLS = [
    (ANCHORS[i], rotation, engine)
    for i, rotation, engine in itertools.product(range(len(ANCHORS)), ROTATIONS_DEG, ENGINES)
]
CELL_IDS = [_cell_id(cell) for cell in CELLS]
VALUE_CELLS = [
    (ANCHORS[i], rotation, "FDI")
    for i, rotation in itertools.product(range(len(ANCHORS)), ROTATIONS_DEG)
]
VALUE_IDS = [_cell_id(cell) for cell in VALUE_CELLS]

# Row ranges of the data frame of the fixture.
VALUE_ROWS = slice(0, 27)
GRADIENT_ROWS = slice(27, 31)
NORMAL_ROWS = slice(31, 35)


def _build(anchor, rotation, engine="FDI", select=None):
    """Return a solved model, its foliation and the world data frame."""
    model, feature, data = world_foliation_model(anchor, rotation)
    if engine != "FDI" or select is not None:
        if select is not None:
            data = select(data)
        model = model_with_frame(anchor, anchor + (8, 10, 12), rotation)
        model.data = data
        feature = model.create_and_add_foliation("r1", interpolatortype=engine, nelements=1000)
    feature.builder.up_to_date()
    assert feature.interpolator.up_to_date
    return model, feature, data


def _snapshot(feature):
    """Return deep copies of the state that a misfit method must not change."""
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
        "report_id": id(interpolator.latest_diagnostics_report),
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
    assert before["report_id"] == after["report_id"]


def _independent_cosines(feature):
    """Return the world cosines of the stored rows, with no change to the state.

    The function copies the stored rows and uses only the public frame methods.
    """
    interpolator = feature.interpolator
    cosines = []
    for key in ("gradient", "normal"):
        rows = np.array(interpolator.data[key], dtype=float, copy=True)
        if rows.shape[0] == 0:
            continue
        positions = interpolator.local_to_world_points(rows[:, :3])
        vectors = interpolator.local_to_world_vectors(rows[:, 3:6])
        model_vectors = feature.evaluate_gradient(positions)
        dots = np.einsum("ij,ij->i", model_vectors, vectors)
        dots /= np.linalg.norm(model_vectors, axis=1) * np.linalg.norm(vectors, axis=1)
        cosines.append(dots)
    return np.concatenate(cosines)


def _world_value_rows(feature, data):
    """Return the world points and values of the value rows."""
    rows = data.iloc[VALUE_ROWS]
    return rows[["X", "Y", "Z"]].to_numpy(float), rows["val"].to_numpy(float)


@pytest.mark.parametrize("cell", CELLS, ids=CELL_IDS)
def test_misfit_leaves_state_unchanged(cell):
    """The old misfit methods change no array, flag, weight or report."""
    anchor, rotation, engine = cell
    _, feature, _ = _build(anchor, rotation, engine)
    before = _snapshot(feature)
    feature.evaluate_gradient_misfit()
    feature.evaluate_value_misfit()
    _assert_same_state(before, _snapshot(feature))


@pytest.mark.parametrize("cell", CELLS, ids=CELL_IDS)
def test_new_misfit_methods_leaves_state_unchanged(cell):
    """The new misfit methods and the report change no array, flag, weight or report."""
    anchor, rotation, engine = cell
    _, feature, _ = _build(anchor, rotation, engine)
    before = _snapshot(feature)
    feature.evaluate_gradient_cosines()
    feature.evaluate_value_residuals()
    feature.verification_report()
    _assert_same_state(before, _snapshot(feature))


@pytest.mark.parametrize("cell", CELLS, ids=CELL_IDS)
def test_gradient_cosine_matches_analytic_direction_exact_field(cell, monkeypatch):
    """With an exact field the world cosine of each stored row is 1."""
    anchor, rotation, engine = cell
    _, feature, _ = _build(anchor, rotation, engine)

    def analytic_gradient(points, *args, **kwargs):
        return np.tile(R1_UNIT, (len(np.atleast_2d(points)), 1))

    monkeypatch.setattr(feature, "evaluate_gradient", analytic_gradient)
    old = np.asarray(feature.evaluate_gradient_misfit(), dtype=float)
    assert old.shape == (2, 4)
    assert np.all(np.abs(1.0 - old) <= 1e-10)
    new = feature.evaluate_gradient_cosines()
    assert new.shape == (8,)
    assert np.all(np.abs(1.0 - new) <= 1e-10)


@pytest.mark.parametrize("cell", CELLS, ids=CELL_IDS)
def test_gradient_cosine_matches_analytic_direction_exact_field_solved_model(cell):
    """The real solved field has a world cosine near 1 at each stored row."""
    anchor, rotation, engine = cell
    _, feature, _ = _build(anchor, rotation, engine)
    independent = _independent_cosines(feature)
    print(f"max |1 - cos| independent {_cell_id(cell)}: {np.max(np.abs(1.0 - independent)):.3e}")
    new = feature.evaluate_gradient_cosines()
    assert new.shape == (8,)
    assert np.max(np.abs(1.0 - new)) <= TOL_SOLVED


@pytest.mark.parametrize("cell", VALUE_CELLS, ids=VALUE_IDS)
def test_value_residuals_are_absolute_difference(cell):
    """The residual is the absolute difference at the 27 value rows."""
    anchor, rotation, engine = cell
    _, feature, data = _build(anchor, rotation, engine)
    points, values = _world_value_rows(feature, data)
    expected = np.abs(values - feature.evaluate_value(points))
    residuals = feature.evaluate_value_residuals()
    assert residuals.shape == (27,)
    assert not np.any(np.isnan(residuals))
    np.testing.assert_allclose(residuals, expected, rtol=0, atol=1e-8)


@pytest.mark.parametrize("cell", VALUE_CELLS, ids=VALUE_IDS)
def test_value_misfit_is_normalised_by_field_range(cell):
    """The old value misfit is the residual divided by the field range."""
    anchor, rotation, engine = cell
    _, feature, _ = _build(anchor, rotation, engine)
    expected = feature.evaluate_value_residuals() / (feature.max() - feature.min())
    np.testing.assert_allclose(feature.evaluate_value_misfit(), expected, rtol=0, atol=1e-12)


@pytest.mark.parametrize("mode", ["equal", "nan"])
def test_value_misfit_with_zero_range_does_not_divide(mode, monkeypatch):
    """A zero or NaN field range gives the unnormalised difference, with no warning."""
    _, feature, data = _build(ANCHORS[1], 30.0, "FDI")
    points, values = _world_value_rows(feature, data)
    expected = np.abs(values - feature.evaluate_value(points))
    low, high = (2.0, 2.0) if mode == "equal" else (np.nan, np.nan)
    monkeypatch.setattr(feature, "min", lambda: low)
    monkeypatch.setattr(feature, "max", lambda: high)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = feature.evaluate_value_misfit()
    assert not np.any(np.isnan(result))
    np.testing.assert_allclose(result, expected, rtol=0, atol=1e-8)


def test_misfit_with_no_rows_is_empty():
    """A model with no gradient rows or no value rows gives empty arrays."""
    _, feature, _ = _build(ANCHORS[1], 30.0, "FDI", select=lambda d: d.iloc[VALUE_ROWS])
    cosines = feature.evaluate_gradient_cosines()
    assert cosines.shape == (0,)
    _, feature, _ = _build(
        ANCHORS[1], 30.0, "FDI", select=lambda d: d.iloc[27:35].reset_index(drop=True)
    )
    residuals = feature.evaluate_value_residuals()
    assert residuals.shape == (0,)


def _rows(gradient, normal):
    """Return a selector for the value rows and the first rows of each family."""

    def select(data):
        parts = [data.iloc[VALUE_ROWS]]
        parts.append(data.iloc[GRADIENT_ROWS].iloc[:gradient])
        parts.append(data.iloc[NORMAL_ROWS].iloc[:normal])
        return pd.concat(parts, ignore_index=True)

    return select


@pytest.mark.parametrize(
    ("gradient", "normal", "shape"),
    [(4, 0, (1, 4)), (0, 4, (1, 4)), (4, 4, (2, 4)), (4, 2, (1, 6))],
    ids=["gradient-only", "normal-only", "equal-families", "different-families"],
)
def test_old_misfit_shapes_are_kept(gradient, normal, shape):
    """The old method keeps its shapes and holds the values of the new method."""
    _, feature, _ = _build(ANCHORS[1], 30.0, "FDI", select=_rows(gradient, normal))
    old = np.asarray(feature.evaluate_gradient_misfit(), dtype=float)
    assert old.shape == shape
    np.testing.assert_array_equal(old.ravel(), feature.evaluate_gradient_cosines())
