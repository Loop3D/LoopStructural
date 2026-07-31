import numpy as np
import pytest
from loop_common.geometry import BoundingBox
from loop_interpolation import InterpolatorBuilder, InterpolatorType


@pytest.mark.parametrize("interpolator_type", ["PLI", "FDI"])
@pytest.mark.parametrize("magnitude", [0.1, 0.5, 1.0, 2.0, 5.0])
@pytest.mark.parametrize(
    "normal_direction",
    [
        [1, 0, 0],  # x-axis
        [0, 1, 0],  # y-axis
        [0, 0, 1],  # z-axis
        [1, 1, 0],  # xy diagonal
        [0, 1, 1],  # yz diagonal
        [1, 0, 1],  # xz diagonal
        [1, 1, 1],  # xyz diagonal
        [-1, 1, 0],  # negative x
        [0, -1, 1],  # negative y
        [1, 0, -1],  # negative z
        [2, 1, 3],  # arbitrary non-axis
        [-2, 2, 1],  # arbitrary non-axis
        [0.5, -1.5, 2],  # arbitrary non-axis
        [1, 2, -2],  # arbitrary non-axis
        [-1, -1, 2],  # arbitrary non-axis
    ],
)
def test_gradient_magnitude_with_normal_constraint(interpolator_type, magnitude, normal_direction):
    # Create a bounding box and builder
    bounding_box = BoundingBox(np.array([0, 0, 0]), np.array([1, 1, 1]))
    interpolatortype = (
        InterpolatorType.PIECEWISE_LINEAR
        if interpolator_type == "PLI"
        else InterpolatorType.FINITE_DIFFERENCE
    )
    builder = InterpolatorBuilder(
        interpolatortype=interpolatortype,
        bounding_box=bounding_box,
        nelements=1000,
        buffer=0.2,
    )

    # Set up a single normal constraint at the center
    center = np.array([[0.5, 0.5, 0.5]])
    normal = np.array([normal_direction], dtype=float)
    normal = normal / np.linalg.norm(normal) * magnitude
    normal_constraints = np.hstack([center, normal, [[np.nan]]])

    # Add constraints and build the interpolator
    builder.add_normal_constraints(normal_constraints)
    interpolator = builder.build()

    # Evaluate the gradient at the constraint location
    grad = interpolator.evaluate_gradient(center)[0]
    grad_mag = np.linalg.norm(grad)

    # The direction should match, and the magnitude should be close to the input magnitude
    assert np.allclose(grad / grad_mag, normal[0] / magnitude, atol=1e-2)
    assert np.isclose(grad_mag, magnitude, atol=0.2)
