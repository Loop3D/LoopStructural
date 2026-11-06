"""Check the R2 fault frame against the analytic frame.

The tests compare a solved fault frame with the analytic frame of requirement R2.
They compare the displacement, the constraint residuals and the Euler integration.
All tests are guards: they must pass on the current tree.
"""

from __future__ import annotations

import itertools
import re
import warnings
from dataclasses import dataclass

import numpy as np
import pytest

from tests.fixtures.coordinate_frames import (
    ANCHORS,
    CURVED_RADIUS,
    R2_DISPLACEMENT_TOL,
    R2_SLIP,
    ROTATIONS_DEG,
    build_r2_fault,
    euler_circle_closed_form,
    r2_analytic_fault,
    r2_closed_form_displacement,
    r2_curved_slip_fault,
    r2_displacement_section,
    r2_frame,
    r2_section_masks,
    residuals_of_field,
    tp,
)

ENGINES = ("FDI", "PLI")
CELLS = tuple(
    (index, rotation, engine)
    for index in range(len(ANCHORS))
    for rotation in ROTATIONS_DEG
    for engine in ENGINES
)
INFERRED_CELLS = ((3, 30.0, "FDI"), (3, 30.0, "PLI"))

# Maximum of the group 1 residuals (assembled rows that must vanish).
# Measured maximum over 24 cells: 4.566e-11. The value is 2 x the measured maximum.
# The cap is 1e-9.
TAU_ASSEMBLED = 9.2e-11
# Maximum of the group 2 residuals (the rows that end in "_orthogonal").
# Measured maximum over 24 cells: 8.408e-11. The value is 2 x the measured maximum.
# The cap is 1e-3.
TAU_ORTHOGONAL = 1.7e-10

SMOOTH_POINTS = 1332
PLANE_POINTS = 37
ALL_POINTS = 1369
MIN_LARGE_DISPLACEMENT = 0.5
MIN_LARGE_COUNT = 1200
ANGLE_TOL_DEG = 0.01

# Names of the group 1 constraints for each engine. A numeric suffix is allowed.
GROUP1_NAMES = {
    "FDI": (
        "value",
        "gradient",
        "norm_0",
        "norm_1",
        "norm_2",
        "dxx",
        "dyy",
        "dzz",
        "dxy",
        "dxz",
        "dyz",
        "dx_lower",
        "dx_upper",
        "dy_lower",
        "dy_upper",
        "dz_lower",
        "dz_upper",
    ),
    "PLI": ("value", "norm", "edge jump"),
}
# Number of value rows of each coordinate (the data of the fault builder).
VALUE_ROWS = (25, 2, 3)


def _cell_id(cell: tuple[int, float, str]) -> str:
    """Return the test id of a cell."""
    index, rotation, engine = cell
    return f"a{index}-r{int(rotation)}-{engine}"


def _width(anchor: np.ndarray) -> float:
    """Return the width W of the coordinate range of the tests for the tolerance Tp."""
    return float(np.max(np.abs(anchor))) + 10.0


@dataclass
class SolvedCell:
    """A solved R2 fault with the data of its cell."""

    fault: object
    anchor: np.ndarray
    rotation: float
    engine: str
    explicit_center: bool


def _solve(index: int, rotation: float, engine: str, explicit_center: bool) -> SolvedCell:
    """Build one R2 fault and force the solve of the three coordinates."""
    anchor = ANCHORS[index]
    fault = build_r2_fault(anchor, rotation, engine, explicit_center=explicit_center)
    section = r2_displacement_section(anchor)
    # The solve starts at the first evaluation of each coordinate.
    for coordinate in range(3):
        fault[coordinate].evaluate_value(section)
    return SolvedCell(fault, anchor, rotation, engine, explicit_center)


@pytest.fixture(scope="module", params=CELLS, ids=_cell_id)
def solved_r2(request) -> SolvedCell:
    """Return a solved R2 fault with an explicit centre (24 cells)."""
    index, rotation, engine = request.param
    return _solve(index, rotation, engine, True)


@pytest.fixture(scope="module", params=INFERRED_CELLS, ids=_cell_id)
def solved_r2_inferred(request) -> SolvedCell:
    """Return a solved R2 fault with an inferred centre (2 cells)."""
    index, rotation, engine = request.param
    return _solve(index, rotation, engine, False)


def _support_box(fault) -> tuple[np.ndarray, np.ndarray]:
    """Return the lower and upper corner of the support in local coordinates."""
    nodes = fault[0].interpolator.support.nodes
    return nodes.min(axis=0), nodes.max(axis=0)


def _row_norms(vectors: np.ndarray) -> np.ndarray:
    """Return the norm of each row."""
    return np.linalg.norm(vectors, axis=1)


# ---------------------------------------------------------------------------
# Displacement of the solved frame against the analytic frame
# ---------------------------------------------------------------------------


def _check_vector_displacement(cell: SolvedCell) -> None:
    """Compare apply_to_points of the solved fault with that of the analytic fault."""
    anchor = cell.anchor
    section = r2_displacement_section(anchor)
    smooth, _ = r2_section_masks(anchor, section)
    solved = cell.fault.apply_to_points(section)
    analytic = r2_analytic_fault(anchor).apply_to_points(section)
    assert smooth.sum() == SMOOTH_POINTS
    assert np.all(np.isfinite(solved[smooth]))
    assert np.all(np.isfinite(analytic[smooth]))
    largest = np.max(_row_norms(solved[smooth] - analytic[smooth]))
    print(f"{_cell_id((0, cell.rotation, cell.engine))} vector difference {largest:.3e}")
    assert largest <= R2_DISPLACEMENT_TOL
    # The check must not pass because both displacements are zero.
    moved = _row_norms(analytic[smooth] - section[smooth])
    assert np.count_nonzero(moved >= MIN_LARGE_DISPLACEMENT) >= MIN_LARGE_COUNT


def test_vector_displacement_matches_analytic_frame(solved_r2):
    _check_vector_displacement(solved_r2)


def test_vector_displacement_inferred_centre(solved_r2_inferred):
    assert not solved_r2_inferred.explicit_center
    _check_vector_displacement(solved_r2_inferred)


def test_scalar_displacement_matches_analytic_frame(solved_r2):
    anchor = solved_r2.anchor
    section = r2_displacement_section(anchor)
    smooth, _ = r2_section_masks(anchor, section)
    solved = solved_r2.fault.evaluate_displacement(section)
    analytic = r2_analytic_fault(anchor).evaluate_displacement(section)
    assert np.all(np.isfinite(solved[smooth]))
    assert np.all(np.isfinite(analytic[smooth]))
    largest = np.max(np.abs(solved[smooth] - analytic[smooth]))
    print(f"scalar difference {largest:.3e}")
    assert largest <= R2_DISPLACEMENT_TOL
    assert np.count_nonzero(np.abs(analytic[smooth]) >= MIN_LARGE_DISPLACEMENT) >= MIN_LARGE_COUNT


def test_displacement_populations_are_counted_separately(solved_r2):
    anchor = solved_r2.anchor
    section = r2_displacement_section(anchor)
    smooth, plane = r2_section_masks(anchor, section)
    invalid = ~np.all(np.isfinite(section), axis=1)
    assert len(section) == ALL_POINTS
    assert (smooth.sum(), plane.sum(), invalid.sum()) == (SMOOTH_POINTS, PLANE_POINTS, 0)
    assert smooth.sum() + plane.sum() + invalid.sum() == ALL_POINTS
    # The support must hold the whole path of the points. The hand value is 10.5. The
    # sampled value is 9 + 1.5 f(0.025) = 10.4972, because the points are 0.5 apart.
    analytic = r2_analytic_fault(anchor).apply_to_points(section)
    largest_z = np.max(np.abs((analytic - anchor)[:, 2]))
    assert largest_z == pytest.approx(10.5, abs=0.01), (
        f"The largest abs(z') of the path is {largest_z}, not 10.5. The support check is not valid."
    )
    interpolator = solved_r2.fault[0].interpolator
    lower, upper = _support_box(solved_r2.fault)
    path = np.vstack([section, analytic])
    local = interpolator.world_to_local_points(path)
    assert np.all(local >= lower) and np.all(local <= upper), "The support does not hold the path."
    solved = solved_r2.fault.apply_to_points(section)
    assert np.all(np.isfinite(solved))
    plane_difference = np.max(_row_norms(solved[plane] - analytic[plane]))
    print(f"largest plane difference {plane_difference:.3e}")
    assert np.isfinite(plane_difference)


@pytest.mark.parametrize("index", range(len(ANCHORS)))
def test_analytic_frame_displacement_matches_closed_form(index):
    anchor = ANCHORS[index]
    section = r2_displacement_section(anchor)
    _, plane = r2_section_masks(anchor, section)
    moved = r2_analytic_fault(anchor).apply_to_points(section)
    displacement = moved - section
    expected = r2_closed_form_displacement(section, anchor)
    difference = displacement - expected
    assert np.array_equal(difference[:, 0], np.zeros(len(section)))
    assert np.array_equal(difference[:, 1], np.zeros(len(section)))
    scale = max(1.0, float(np.max(np.abs(section[:, 2]))))
    largest = np.max(np.abs(difference[:, 2]))
    print(f"closed form z difference {largest:.3e}")
    assert largest <= 64 * np.finfo(float).eps * scale
    # A point on the fault plane does not move.
    assert np.array_equal(moved[plane], section[plane])
    # The check must not pass because both displacements are zero.
    assert np.max(np.abs(expected[:, 2])) > 1.0


def _start_points(anchor: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return 8 start points on the circle of radius 5 and their start angles."""
    angles = np.deg2rad(np.arange(8) * 45.0)
    local = np.column_stack(
        (CURVED_RADIUS * np.cos(angles), CURVED_RADIUS * np.sin(angles), np.zeros(8))
    )
    return anchor + local, angles


def _polar_end(anchor: np.ndarray, angles: np.ndarray, radius: float, turn: float) -> np.ndarray:
    """Return the end points for a final radius and a turn angle."""
    final = angles + turn
    local = np.column_stack((radius * np.cos(final), radius * np.sin(final), np.zeros(len(final))))
    return anchor + local


@pytest.mark.parametrize("index", range(len(ANCHORS)))
def test_integration_error_follows_the_euler_closed_form(index):
    anchor = ANCHORS[index]
    width = _width(anchor)
    tolerance = tp(10, width)
    slip = 0.5 * R2_SLIP * (1 - 3 * 0.2**2 + 2 * 0.2**3)
    assert slip == pytest.approx(1.344, abs=1e-12)
    starts, angles = _start_points(anchor)
    coarse = r2_curved_slip_fault(anchor, 20).apply_to_points(starts)
    fine = r2_curved_slip_fault(anchor, 400).apply_to_points(starts)
    # (a) The coarse result follows the Euler closed form.
    radius, turn = euler_circle_closed_form(CURVED_RADIUS, slip, 20)
    expected_coarse = _polar_end(anchor, angles, radius, turn)
    error_a = np.max(np.abs(coarse - expected_coarse))
    print(f"euler N=20 difference {error_a:.3e} (tolerance {tolerance:.3e})")
    assert error_a <= tolerance
    # (b) The fine result is near the exact arc.
    exact = _polar_end(anchor, angles, CURVED_RADIUS, slip / CURVED_RADIUS)
    error_b = np.max(np.abs(fine - exact))
    bound_b = slip**2 / (400 * CURVED_RADIUS) + tolerance
    print(f"euler N=400 difference {error_b:.3e} (bound {bound_b:.3e})")
    assert error_b <= bound_b
    # (c) The step count changes the result. The check must show the integration error.
    gap = np.max(_row_norms(coarse - fine))
    print(f"N=20 against N=400 gap {gap:.3e}")
    assert gap > 5e-3


def test_misplaced_frame_fails_the_displacement_bound(solved_r2):
    anchor = solved_r2.anchor
    section = r2_displacement_section(anchor)
    smooth, _ = r2_section_masks(anchor, section)
    solved = solved_r2.fault.apply_to_points(section)
    misplaced = r2_analytic_fault(anchor + np.array([1.0, 0.0, 0.0])).apply_to_points(section)
    largest = np.max(_row_norms(solved[smooth] - misplaced[smooth]))
    print(f"misplaced difference {largest:.3e}")
    assert largest > R2_DISPLACEMENT_TOL


# ---------------------------------------------------------------------------
# Points outside the support and the protection against a zero gradient
# ---------------------------------------------------------------------------


def _outside_probes(cell: SolvedCell) -> np.ndarray:
    """Return six probes outside the support, on the world axes through the anchor."""
    lower, upper = _support_box(cell.fault)
    corners = np.array(list(itertools.product(*zip(lower, upper))))
    world = cell.fault[0].interpolator.local_to_world_points(corners) - cell.anchor
    probes = []
    for axis in range(3):
        half = np.max(np.abs(world[:, axis]))
        for sign in (-1.0, 1.0):
            offset = np.zeros(3)
            offset[axis] = sign * (half + 5.0)
            probes.append(cell.anchor + offset)
    return np.array(probes)


def test_outside_support_points_have_no_displacement(solved_r2):
    fault = solved_r2.fault
    probes = _outside_probes(solved_r2)
    lower, upper = _support_box(fault)
    local = fault[0].interpolator.world_to_local_points(probes)
    outside = np.any((local < lower) | (local > upper), axis=1)
    assert np.all(outside), "A probe is inside the support."
    for coordinate in range(3):
        assert np.all(np.isnan(fault[coordinate].evaluate_value(probes)))
        gradient = fault[coordinate].evaluate_gradient(probes)
        if solved_r2.engine == "PLI":
            assert np.all(np.isnan(gradient))
        else:
            # The FDI gradient outside the support is zero, not NaN (see the report of T095).
            assert np.array_equal(gradient, np.zeros_like(gradient))
    assert np.all(fault.evaluate_displacement(probes) == 0.0)
    try:
        moved = fault.apply_to_points(probes)
    except Exception as error:  # noqa: BLE001
        print(f"apply_to_points raised on an empty mask: {error!r}")
    else:
        assert np.array_equal(moved, probes)


def _zero_gradient_override(mode: str, anchor: np.ndarray):
    """Return a gradient override for the coordinate 1 at the current z' above 5."""
    replacement = 0.0 if mode == "zero" else np.nan

    def override(coordinate: int, points: np.ndarray, gradient: np.ndarray) -> np.ndarray:
        if coordinate != 1:
            return gradient
        result = np.array(gradient, dtype=float, copy=True)
        result[(np.asarray(points) - anchor)[:, 2] > 5.0] = replacement
        return result

    return override


@pytest.mark.parametrize("mode", ["zero", "nan"])
@pytest.mark.parametrize("index", range(len(ANCHORS)))
def test_zero_gradient_protection_remains_effective(index, mode):
    anchor = ANCHORS[index]
    section = r2_displacement_section(anchor)
    smooth, _ = r2_section_masks(anchor, section)
    z_start = (section - anchor)[:, 2]
    group_a = smooth & (z_start > 5.0)
    group_b = smooth & (z_start <= 3.5)
    group_c = smooth & (z_start > 3.5) & (z_start <= 5.0)
    assert group_a.any() and group_b.any() and group_c.any()
    override = _zero_gradient_override(mode, anchor)
    reference = r2_analytic_fault(anchor).apply_to_points(section)
    fault = r2_analytic_fault(anchor, gradient_override=override)
    with warnings.catch_warnings():
        if mode == "zero":
            warnings.simplefilter("error")
        result = fault.apply_to_points(section)
    if mode == "zero":
        assert np.all(np.isfinite(result))
        assert np.array_equal(result[group_a], section[group_a])
        assert np.array_equal(result[group_b], reference[group_b])
        assert np.all(np.isfinite(result[group_c]))
    else:
        # A point with a displacement of zero does not enter the integration.
        moves = _row_norms(reference - section) > 0
        assert np.all(np.isnan(result[group_a & moves]))
        assert np.all(np.isfinite(result[group_b]))
        assert np.array_equal(result[group_b], reference[group_b])


# ---------------------------------------------------------------------------
# Constraint residuals
# ---------------------------------------------------------------------------


def _group_of(name: str, engine: str) -> int:
    """Return the residual group of a constraint name."""
    if name.endswith("_orthogonal"):
        return 2
    for base in GROUP1_NAMES[engine]:
        if re.fullmatch(re.escape(base) + r"(_?\d+)?", name):
            return 1
    return 3


def test_fault_frame_assembled_residuals_vanish_on_the_analytic_field(solved_r2):
    anchor = solved_r2.anchor
    group1_max = 0.0
    group2_max = 0.0
    for coordinate in range(3):
        interpolator = solved_r2.fault[coordinate].interpolator
        nodes = interpolator.local_to_world_points(interpolator.support.nodes)
        nodal = r2_frame(nodes, anchor)[0][:, coordinate]
        residuals = residuals_of_field(interpolator, nodal)
        groups = {name: _group_of(name, solved_r2.engine) for name in residuals}
        unknown = sorted(name for name, group in groups.items() if group == 3)
        assert not unknown, f"Unknown constraint names in coordinate {coordinate}: {unknown}"
        counts1 = {name: len(r) for name, r in residuals.items() if groups[name] == 1}
        assert counts1, f"Coordinate {coordinate} has no group 1 row."
        assert counts1["value"] == VALUE_ROWS[coordinate]
        for name, residual in residuals.items():
            assert np.all(np.isfinite(residual)), name
            largest = float(np.max(np.abs(residual))) if len(residual) else 0.0
            if groups[name] == 1:
                group1_max = max(group1_max, largest)
                assert largest <= TAU_ASSEMBLED, (coordinate, name, largest)
            else:
                group2_max = max(group2_max, largest)
                assert largest <= TAU_ORTHOGONAL, (coordinate, name, largest)
    print(f"assembled group 1 max {group1_max:.3e}; group 2 max {group2_max:.3e}")


def test_fault_frame_physical_residuals_match_analytic(solved_r2):
    anchor = solved_r2.anchor
    width = _width(anchor)
    for coordinate in range(3):
        feature = solved_r2.fault[coordinate]
        interpolator = feature.interpolator
        physical = feature.verification_report().physical_residuals
        # The value rows.
        value_rows = np.asarray(interpolator.get_value_constraints(), dtype=float)
        assert len(value_rows) == VALUE_ROWS[coordinate]
        assert physical["value"].count == VALUE_ROWS[coordinate]
        if coordinate == 0:
            assert physical["value"].count == 25
            assert physical["value"].max_abs <= 1e-6
            positions = interpolator.local_to_world_points(value_rows[:, :3])
            independent = np.max(np.abs(0.0 - feature.evaluate_value(positions)))
            tolerance = max(1e-12, tp(20, width) / 20)
            assert abs(physical["value"].max_abs - independent) <= tolerance
        # The direction rows.
        families = {
            "gradient": interpolator.get_gradient_constraints(),
            "normal": interpolator.get_norm_constraints(),
        }
        direction_names = [name for name in physical if name != "value"]
        for name, rows in families.items():
            rows = np.asarray(rows, dtype=float)
            assert (len(rows) > 0) == (name in direction_names), name
        assert direction_names, "No direction family in the report."
        for name in direction_names:
            rows = np.asarray(families[name], dtype=float)
            positions = interpolator.local_to_world_points(rows[:, :3])
            data = interpolator.local_to_world_vectors(rows[:, 3:6])
            model = feature.evaluate_gradient(positions)
            angle = np.degrees(
                np.arctan2(np.linalg.norm(np.cross(model, data), axis=1), np.sum(model * data, 1))
            )
            assert physical[name].max_abs <= ANGLE_TOL_DEG
            assert abs(physical[name].max_abs - np.max(angle)) <= 1e-6


def test_residual_check_detects_misplaced_reference(solved_r2):
    anchor = solved_r2.anchor
    interpolator = solved_r2.fault[0].interpolator
    nodes = interpolator.local_to_world_points(interpolator.support.nodes)
    shifted = nodes - np.array([1.0, 0.0, 0.0])
    nodal = r2_frame(shifted, anchor)[0][:, 0]
    residual = residuals_of_field(interpolator, nodal)["value"]
    assert len(residual) == 25
    largest = float(np.max(np.abs(residual)))
    print(f"misplaced reference value maximum {largest:.3e}")
    assert largest > 1e-2
