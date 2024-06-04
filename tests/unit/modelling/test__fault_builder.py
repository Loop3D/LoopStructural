import numpy as np
import pandas as pd
import pytest
from LoopStructural.modelling.features.builders._fault_builder import FaultBuilder
from LoopStructural.datatypes import BoundingBox
from LoopStructural import GeologicalModel


def test_fault_builder_update_geometry(interpolatortype):
    # Create a FaultBuilder instance
    bounding_box = BoundingBox([0, 0, 0], [1, 1, 1])
    nelements = 1000
    model = None
    fault_bounding_box_buffer = 0.2
    fault_builder = FaultBuilder(
        interpolatortype, bounding_box, nelements, model, fault_bounding_box_buffer
    )

    # Update the geometry with some points
    points = np.array([[0.5, 0.5, 0.5], [0.6, 0.6, 0.6], [0.7, 0.7, 0.7]])
    fault_builder.update_geometry(points)

    # Check if the origin and maximum are updated correctly
    expected_origin = np.array([0.5, 0.5, 0.5])
    expected_maximum = np.array([0.7, 0.7, 0.7])
    assert np.allclose(fault_builder.origin, expected_origin)
    assert np.allclose(fault_builder.maximum, expected_maximum)


def test_fault_builder_create_data_from_geometry(interpolatortype):
    # Create a FaultBuilder instance
    bounding_box = BoundingBox([0, 0, 0], [1, 1, 1])
    nelements = 1000
    model = GeologicalModel([0, 0, 0], [1, 1, 1])
    fault_bounding_box_buffer = 0.2
    fault_builder = FaultBuilder(
        interpolatortype, bounding_box, nelements, model, fault_bounding_box_buffer
    )

    # Create some test data
    fault_frame_data = pd.DataFrame(
        {
            "coord": [0.0, 0.0, 1.0, 2.0],
            "val": [0.0, 0.0, 1.0, 1.0],
            "gx": [np.nan, np.nan, np.nan, np.nan],
            "gy": [np.nan, np.nan, np.nan, np.nan],
            "gz": [np.nan, np.nan, np.nan, np.nan],
            "nx": [0.0, 0.0, 0.0, 1.0],
            "ny": [1.0, 1.0, 0.0, 0.0],
            "nz": [0.0, 0.0, 1.0, 0.0],
            "X": [0.1, 0.2, 0.3, 0.4],
            "Y": [0.5, 0.6, 0.7, 0.8],
            "Z": [0.9, 1.0, 1.1, 1.2],
        }
    )

    # Call the create_data_from_geometry method
    fault_builder.create_data_from_geometry(fault_frame_data)

    # Check if the fault attributes are updated correctly
    expected_fault_normal_vector = np.array([0, 1.0, 0])
    expected_fault_slip_vector = np.array([0.0, 0.0, 1.0])
    # expected_fault_centre = np.array([0.15, 0.6, 1.05])
    expected_fault_normal_vector /= np.linalg.norm(expected_fault_normal_vector)
    expected_fault_slip_vector /= np.linalg.norm(expected_fault_slip_vector)
    # expected_fault_minor_axis = 0.5
    # expected_fault_major_axis = 1.0
    # expected_fault_intermediate_axis = 1.0
    calculated_fault_normal_vector = fault_builder.fault_normal_vector
    calculated_fault_slip_vector = fault_builder.fault_slip_vector
    calculated_fault_normal_vector /= np.linalg.norm(calculated_fault_normal_vector)
    calculated_fault_slip_vector /= np.linalg.norm(calculated_fault_slip_vector)
    assert np.allclose(calculated_fault_normal_vector, expected_fault_normal_vector)
    assert np.allclose(calculated_fault_slip_vector, expected_fault_slip_vector)
    # assert np.allclose(fault_builder.fault_centre, expected_fault_centre)
    # assert np.isclose(fault_builder.fault_minor_axis, expected_fault_minor_axis)
    # assert np.isclose(fault_builder.fault_major_axis, expected_fault_major_axis)
    # assert np.isclose(fault_builder.fault_intermediate_axis, expected_fault_intermediate_axis)


# Run the tests
if __name__ == "__main__":
    pytest.main([__file__])
