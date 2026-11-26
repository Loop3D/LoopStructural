import numpy as np
import pytest

pyevtk = pytest.importorskip("pyevtk")

from LoopStructural.export import exporters
from LoopStructural.export.file_formats import FileFormat
from LoopStructural.geometry import BoundingBox

# ---------------------------------------------------------------------------
# Lightweight fakes standing in for a GeologicalModel/BoundingBox, so these
# tests can exercise the dispatch functions in exporters.py without needing to
# build a full geological model.
# ---------------------------------------------------------------------------


class _FakeModelForVolume:
    """Fake model exposing what `write_cubeface`/`write_vol` touch. Uses a
    real `BoundingBox` since that is what `GeologicalModel.bounding_box`
    actually is, and the volume writers call real `BoundingBox` methods
    (`regular_grid`) as well as raw indexing (`[:]`)."""

    def __init__(self):
        self.bounding_box = BoundingBox(np.array([0.0, 0.0, 0.0]), np.array([1.0, 1.0, 1.0]))

    def rescale(self, points):
        return points

    def evaluate_model(self, points, scale=True):
        centre = np.array([0.5, 0.5, 0.5])
        distance = np.linalg.norm(points - centre, axis=1) - 0.3
        return (distance > 0).astype(np.int64)


# ---------------------------------------------------------------------------
# write_cubeface
# ---------------------------------------------------------------------------


def test_write_cubeface_vtk_writes_file(tmp_path):
    model = _FakeModelForVolume()
    file_name = tmp_path / "cube"
    result = exporters.write_cubeface(
        model, str(file_name), "label", np.array([5, 5, 5]), FileFormat.VTK
    )
    assert result is True
    assert (tmp_path / "cube.vtu").exists()


def test_write_cubeface_unsupported_format_returns_false(tmp_path):
    model = _FakeModelForVolume()
    file_name = tmp_path / "cube_gocad"
    result = exporters.write_cubeface(
        model, str(file_name), "label", np.array([5, 5, 5]), FileFormat.GOCAD
    )
    assert result is False


# ---------------------------------------------------------------------------
# write_vol
# ---------------------------------------------------------------------------


def test_write_vol_vtk_writes_file(tmp_path):
    model = _FakeModelForVolume()
    file_name = tmp_path / "vol"
    result = exporters.write_vol(
        model, str(file_name), "label", np.array([5, 5, 5]), FileFormat.VTK
    )
    assert result is True
    assert (tmp_path / "vol.vtu").exists()


def test_write_vol_unsupported_format_returns_false(tmp_path):
    model = _FakeModelForVolume()
    file_name = tmp_path / "vol_obj"
    result = exporters.write_vol(
        model, str(file_name), "label", np.array([5, 5, 5]), FileFormat.OBJ
    )
    assert result is False
