"""Regression test for PointSet JSON/YAML serialization (P0 fix)."""

import pytest
import numpy as np
import json

from loop_common.observations import PointSet


def test_pointset_to_json_roundtrip():
    """Test that PointSet can be serialized to JSON and back without losing data."""
    points = np.array([
        [0.0, 1.0, 2.0],
        [3.0, 4.0, 5.0],
        [6.0, 7.0, 8.0],
    ])
    original = PointSet(name="test_points", coords=points)

    # Serialize to JSON string
    json_str = original.to_json()

    # Verify it's valid JSON
    json_data = json.loads(json_str)
    assert isinstance(json_data, dict)

    # coords should be serialized as a list
    assert "coords" in json_data
    assert isinstance(json_data["coords"], list)
    assert len(json_data["coords"]) == 3

    # Deserialize back
    restored = PointSet.from_json(json_str)

    # Verify the data matches
    assert restored.name == original.name
    assert np.allclose(restored.coords, original.coords)


def test_pointset_model_dump_json():
    """Test that PointSet.model_dump(mode='json') properly serializes numpy arrays."""
    points = np.array([
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
    ])
    pointset = PointSet(name="test", coords=points)

    # This should not raise an error about unserializable numpy arrays
    dumped = pointset.model_dump(mode="json")

    # coords should be a list in the dumped dict
    assert isinstance(dumped["coords"], list)
    assert dumped["coords"] == [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]


def test_pointset_list_input_coerced_to_array():
    """Test that PointSet accepts list input and coerces it to numpy array."""
    points_list = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
    pointset = PointSet(name="test", coords=points_list)

    # Should be coerced to numpy array
    assert isinstance(pointset.coords, np.ndarray)
    assert pointset.coords.shape == (2, 3)


def test_pointset_yaml_serialization():
    """Test that PointSet can be serialized to YAML (if yaml support is available)."""
    points = np.array([
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
    ])
    pointset = PointSet(name="yaml_test", coords=points)

    # This should not raise (to_yaml calls model_dump(mode='json') internally)
    yaml_str = pointset.to_yaml()
    assert isinstance(yaml_str, str)
    assert "coords" in yaml_str
