import numpy as np
import pytest

pyevtk = pytest.importorskip("pyevtk")

from LoopStructural.export import exporters
from LoopStructural.export.file_formats import FileFormat
from LoopStructural.geometry import BoundingBox, Surface
from LoopStructural.utils.exceptions import LoopValueError

# ---------------------------------------------------------------------------
# Lightweight fakes standing in for a GeologicalModel/BoundingBox, so these
# tests can exercise the dispatch functions in exporters.py without needing to
# build a full geological model.
# ---------------------------------------------------------------------------


class _SphericalFeature:
    """A fake geological feature whose scalar field is a signed distance to a
    sphere centred in the unit cube - guarantees a clean isosurface at 0."""

    def evaluate_value(self, points):
        centre = np.array([0.5, 0.5, 0.5])
        return np.linalg.norm(points - centre, axis=1) - 0.3


class _FakeBoundingBoxArrayLike:
    """Stands in for the parts of `model.bounding_box` that `write_feat_surfs`
    accesses directly through numpy-style indexing (`.bb[...]`)."""

    def __init__(self, bb):
        self._bb = np.asarray(bb, dtype=float)

    @property
    def bb(self):
        return self._bb


class _FakeModelForSurfaces:
    """Fake model exposing only what `write_feat_surfs` touches."""

    def __init__(self, feature_name="strati"):
        self.bounding_box = _FakeBoundingBoxArrayLike([[0, 0, 0], [1, 1, 1]])
        self.nsteps = np.array([12, 12, 12])
        self._features = {feature_name: _SphericalFeature()}

    def __contains__(self, name):
        return name in self._features

    def __getitem__(self, name):
        return self._features[name]

    def rescale(self, points):
        # identity rescale, in-place like the real implementation would allow
        return points


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
# write_feat_surfs
# ---------------------------------------------------------------------------


def test_write_feat_surfs_feature_not_in_model_returns_false_empty():
    model = _FakeModelForSurfaces()
    result = exporters.write_feat_surfs(model, "not_a_feature", file_format=FileFormat.NUMPY)
    assert result == (False, [])


def test_write_feat_surfs_isovalue_outside_range_returns_false_empty():
    model = _FakeModelForSurfaces()
    result = exporters.write_feat_surfs(
        model, "strati", file_format=FileFormat.NUMPY, isovalue=999.0
    )
    assert result == (False, [])


def test_write_feat_surfs_unsupported_format_returns_false_empty():
    model = _FakeModelForSurfaces()
    result = exporters.write_feat_surfs(model, "strati", file_format=FileFormat.OBJ)
    assert result == (False, [])


def test_write_feat_surfs_numpy_format_success():
    # NOTE: this documents a real bug - the docstring for write_feat_surfs
    # promises a `(bool, [Surface, ...])` tuple return, and every early-exit
    # path in the function does return such a tuple, but the final
    # success-path `return result` (LoopStructural/export/exporters.py) only
    # returns the bare boolean, breaking the documented contract. Calling code
    # written against the docstring (`ok, surfaces = write_feat_surfs(...)`)
    # would raise `TypeError: cannot unpack non-iterable bool object`.
    model = _FakeModelForSurfaces()
    result = exporters.write_feat_surfs(model, "strati", file_format=FileFormat.NUMPY)
    assert result is True


def test_write_feat_surfs_vtk_format_is_broken_for_real_surface(tmp_path):
    # NOTE: this documents a second bug - `_write_feat_surfs_evtk` reads
    # `surf.verts` / `surf.faces` but the `Surface` dataclass that
    # `write_feat_surfs` builds via marching_cubes only exposes `.vertices`
    # and `.triangles`. Any call to write_feat_surfs with FileFormat.VTK
    # therefore always raises AttributeError.
    model = _FakeModelForSurfaces()
    file_name = tmp_path / "iso"
    with pytest.raises(AttributeError, match="verts"):
        exporters.write_feat_surfs(
            model, "strati", file_format=FileFormat.VTK, file_name=str(file_name)
        )


def test_write_feat_surfs_gocad_format_writes_ts_file(tmp_path):
    model = _FakeModelForSurfaces()
    file_name = tmp_path / "iso_gocad"
    result = exporters.write_feat_surfs(
        model, "strati", file_format=FileFormat.GOCAD, file_name=str(file_name)
    )
    assert result is True
    ts_file = tmp_path / "iso_gocad.ts"
    assert ts_file.exists()
    content = ts_file.read_text()
    assert "GOCAD TSurf 1" in content
    assert "name: strati" in content


# ---------------------------------------------------------------------------
# _write_feat_surfs_evtk / _write_feat_surfs_gocad (called with correctly
# shaped objects, to isolate the writer logic from the attribute-name bug
# above)
# ---------------------------------------------------------------------------


class _EvtkCompatibleSurf:
    def __init__(self):
        self.verts = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        self.faces = np.array([[0, 1, 2]])
        self.values = np.array([1.0, 2.0, 3.0], dtype=np.float32)
        self.normals = np.array([[0.0, 0.0, 1.0]] * 3)
        self.name = "evtk_surf"


def test_write_feat_surfs_evtk_writes_file_when_attribute_names_match(tmp_path):
    surf = _EvtkCompatibleSurf()
    file_name = tmp_path / "compatible"
    result = exporters._write_feat_surfs_evtk(surf, str(file_name))
    assert result is True
    assert (tmp_path / "compatible.vtu").exists()


def test_write_feat_surfs_gocad_direct_call(tmp_path):
    surf = Surface(
        vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        triangles=np.array([[0, 1, 2]]),
        name="direct_surf",
    )
    file_name = tmp_path / "direct"
    result = exporters._write_feat_surfs_gocad(surf, str(file_name))
    assert result is True
    content = (tmp_path / "direct.ts").read_text()
    assert "name: direct_surf" in content
    assert "TRGL 1 2 3" in content


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


def test_write_vol_gocad_is_broken_for_real_bounding_box(tmp_path):
    # NOTE: this documents a third bug - `_write_vol_gocad` does
    # `bbox = model.bounding_box[:]`, treating `model.bounding_box` as a raw
    # numpy array. In practice (both here and in GeologicalModel) it is a
    # `BoundingBox` object whose `__getitem__` only supports string names (or
    # falls through to raising `LoopValueError` for anything else, including
    # a bare slice). So `write_vol(..., file_format=FileFormat.GOCAD)` always
    # raises instead of writing a VOXET file.
    model = _FakeModelForVolume()
    file_name = tmp_path / "vol_gocad"
    with pytest.raises(LoopValueError):
        exporters.write_vol(
            model, str(file_name), "label", np.array([5, 5, 5]), FileFormat.GOCAD
        )
