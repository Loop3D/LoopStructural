"""Contract tests for rigid coordinate transforms and grid resize stability."""

import numpy as np
import pytest
from loop_common.geometry import (
    BoundingBox,
    RigidRotationError,
    StructuredGrid3DGeometry,
    validate_rigid_rotation,
)


@pytest.mark.parametrize("dimensions", [2, 3])
@pytest.mark.parametrize("angle", [0.0, 30.0, 90.0])
def test_validate_rigid_rotation_accepts_proper_rotations(dimensions, angle):
    """Accept proper rotations in two and three dimensions."""
    radians = np.deg2rad(angle)
    rotation = np.eye(dimensions)
    rotation[:2, :2] = [[np.cos(radians), -np.sin(radians)], [np.sin(radians), np.cos(radians)]]

    result = validate_rigid_rotation(rotation, dimensions)

    np.testing.assert_allclose(result, rotation)
    assert result is not rotation


@pytest.mark.parametrize(
    "matrix",
    [
        np.diag([2.0, 1.0, 1.0]),
        np.array([[1.0, 0.2, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        np.diag([1.0, 1.0, 0.0]),
        np.diag([-1.0, 1.0, 1.0]),
        np.array([[np.nan, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        np.array([[np.inf, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        np.eye(2),
    ],
    ids=["scale", "shear", "singular", "reflection", "nan", "infinity", "wrong-shape"],
)
def test_validate_rigid_rotation_rejects_invalid(matrix):
    """Reject transforms that are not finite proper rotations."""
    assert issubclass(RigidRotationError, ValueError)
    with pytest.raises(RigidRotationError, match=""):
        validate_rigid_rotation(matrix, 3)


@pytest.mark.parametrize(
    "anchor",
    [
        (0.0, 0.0, 0.0),
        (100.0, 200.0, 300.0),
        (-100.0, -200.0, -300.0),
        (500000.0, 6500000.0, 100.0),
    ],
)
@pytest.mark.parametrize("dimensions", [2, 3])
@pytest.mark.parametrize("angle", [0.0, 30.0, 90.0])
def test_point_and_direction_round_trips(anchor, dimensions, angle):
    """Round-trip points and vectors without translation leaking into vectors."""
    anchor = np.asarray(anchor[:dimensions], dtype=float)
    box = BoundingBox(
        origin=anchor - np.array([100.0, 200.0, 50.0][:dimensions]),
        maximum=anchor + np.array([100.0, 200.0, 50.0][:dimensions]),
        dimensions=dimensions,
    )
    radians = np.deg2rad(angle)
    rotation = np.array([[np.cos(radians), -np.sin(radians)], [np.sin(radians), np.cos(radians)]])
    if dimensions == 3:
        rotation = np.pad(rotation, ((0, 1), (0, 1)))
        rotation[2, 2] = 1.0
    box.set_local_transform(local_origin=np.array(anchor), rotation_matrix=rotation)
    offset = np.array([3.0, -7.0, 11.0][:dimensions])
    points = np.array([anchor, anchor + offset])
    vectors = np.array([[1.0, 2.0, 3.0][:dimensions], [-4.0, 5.0, -6.0][:dimensions]])
    length = np.linalg.norm(box.maximum - box.origin)
    world_magnitude = max(np.linalg.norm(box.origin), np.linalg.norm(box.maximum))
    point_tolerance = 1e-10 * length + 32 * np.finfo(float).eps * max(world_magnitude, length)

    np.testing.assert_allclose(
        box.reproject(box.project(points)), points, atol=point_tolerance, rtol=0.0
    )
    np.testing.assert_allclose(
        box.reproject(box.project(points[0])), points[0], atol=point_tolerance, rtol=0.0
    )
    np.testing.assert_allclose(
        box.reproject_vectors(box.project_vectors(vectors)), vectors, atol=1e-12, rtol=0.0
    )
    np.testing.assert_allclose(
        box.reproject_vectors(box.project_vectors(vectors[0])),
        vectors[0],
        atol=1e-12,
        rtol=0.0,
    )


@pytest.mark.parametrize("angle", [0.0, 30.0, 90.0])
def test_local_enclosure_contains_projected_corners(angle):
    """Enclose every projected world corner within the returned local bounds."""
    origin = np.array([-10.0, 20.0, -30.0])
    maximum = np.array([40.0, 80.0, 60.0])
    box = BoundingBox(origin=origin, maximum=maximum)
    radians = np.deg2rad(angle)
    rotation = np.array(
        [
            [np.cos(radians), -np.sin(radians), 0.0],
            [np.sin(radians), np.cos(radians), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    box.set_local_transform(local_origin=origin, rotation_matrix=rotation)
    local_minimum, local_maximum = box.local_enclosure()
    tolerance = 1e-10 * np.linalg.norm(maximum - origin) + 32 * np.finfo(float).eps * np.max(
        maximum - origin
    )

    assert np.all(box.project(box.corners) >= local_minimum - tolerance)
    assert np.all(box.project(box.corners) <= local_maximum + tolerance)


def test_invalid_rotation_keeps_prior_state():
    """Keep the prior transform when a rotation update fails validation."""
    box = BoundingBox(origin=np.zeros(3), maximum=np.ones(3))
    old_transform = box.world_to_local_matrix
    with pytest.raises(RigidRotationError):
        box.set_local_transform(rotation_matrix=np.diag([2.0, 1.0, 1.0]))
    np.testing.assert_array_equal(box.world_to_local_matrix, old_transform)

    grid = StructuredGrid3DGeometry(rotation_xy=30.0)
    old_rotation = grid.rotation_xy.copy()
    with pytest.raises(RigidRotationError):
        grid.rotation_xy = np.diag([-1.0, 1.0, 1.0])
    np.testing.assert_array_equal(grid.rotation_xy, old_rotation)


def test_reapply_unchanged_bounds_is_stable():
    """Keep geometry unchanged when setting its current bounds repeatedly."""
    grid = StructuredGrid3DGeometry(
        origin=np.array([-10.0, 20.0, 30.0]),
        nsteps=np.array([11, 21, 31]),
        step_vector=np.array([2.0, 3.0, 4.0]),
    )
    origin = grid.origin.copy()
    maximum = grid.maximum.copy()
    step_vector = grid.step_vector.copy()
    nsteps = grid.nsteps.copy()

    for _ in range(10):
        grid.origin = origin.copy()
        grid.maximum = maximum.copy()

    np.testing.assert_array_equal(grid.origin, origin)
    np.testing.assert_array_equal(grid.maximum, maximum)
    np.testing.assert_array_equal(grid.step_vector, step_vector)
    np.testing.assert_array_equal(grid.nsteps, nsteps)
