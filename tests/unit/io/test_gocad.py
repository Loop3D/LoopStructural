import logging

import numpy as np
import pytest

from LoopStructural.datatypes import StructuredGrid, Surface
from LoopStructural.export.gocad import (
    _normalise_voxet_property,
    _write_feat_surfs_gocad,
    _write_structured_grid_gocad,
)


def _read(path):
    with open(path) as fd:
        return fd.read()


# ---------------------------------------------------------------------------
# _normalise_voxet_property
# ---------------------------------------------------------------------------


def test_normalise_voxet_property_small_int_uses_octet():
    info = _normalise_voxet_property(np.array([1, 2, 3]), "prop", np.array([3]))
    assert info["storage_type"] == "Octet"
    assert info["element_size"] == 1
    assert info["values"].dtype == np.int8
    assert info["no_data_value"] is None


def test_normalise_voxet_property_large_int_uses_integer():
    info = _normalise_voxet_property(np.array([1000, -2000, 3000]), "prop", np.array([3]))
    assert info["storage_type"] == "Integer"
    assert info["element_size"] == 4
    assert info["values"].dtype == np.dtype(">i4")


def test_normalise_voxet_property_float_uses_float_and_nan_fill():
    values = np.array([1.0, np.nan, 3.0])
    info = _normalise_voxet_property(values, "prop", np.array([3]))
    assert info["storage_type"] == "Float"
    assert info["element_size"] == 4
    assert info["no_data_value"] == -999999.0
    assert info["values"][1] == np.float32(-999999.0)


def test_normalise_voxet_property_unsupported_dtype_raises():
    with pytest.raises(ValueError):
        _normalise_voxet_property(np.array(["a", "b"]), "prop", np.array([2]))


def test_normalise_voxet_property_wrong_size_raises():
    with pytest.raises(ValueError):
        _normalise_voxet_property(np.array([1.0, 2.0]), "prop", np.array([3]))


def test_normalise_voxet_property_reshapes_matching_grid_shape():
    values = np.arange(8, dtype=float).reshape((2, 2, 2))
    info = _normalise_voxet_property(values, "prop", np.array([2, 2, 2]))
    assert info["values"].shape == (8,)


# ---------------------------------------------------------------------------
# _write_structured_grid_gocad
# ---------------------------------------------------------------------------


def test_write_structured_grid_gocad_point_properties(tmp_path):
    grid = StructuredGrid(
        origin=np.array([0.0, 0.0, 0.0]),
        step_vector=np.array([1.0, 1.0, 1.0]),
        nsteps=np.array([3, 3, 3]),
        name="mygrid",
    )
    grid.properties["val"] = np.arange(27).astype(float)

    file_name = tmp_path / "grid"
    result = _write_structured_grid_gocad(grid, str(file_name))

    assert result is True
    vo_file = tmp_path / "grid.vo"
    data_file = tmp_path / "grid_val@@"
    assert vo_file.exists()
    assert data_file.exists()

    content = _read(vo_file)
    assert "GOCAD Voxet 1" in content
    assert "name: mygrid" in content
    assert "AXIS_N 3 3 3" in content
    assert "PROPERTY 1 val" in content
    assert "PROP_FILE 1 grid_val@@" in content

    # exported data should round-trip as big-endian float32
    raw = np.fromfile(data_file, dtype=np.dtype(">f4"))
    assert raw.shape[0] == 27


def test_write_structured_grid_gocad_cell_properties(tmp_path):
    grid = StructuredGrid(
        origin=np.array([0.0, 0.0, 0.0]),
        step_vector=np.array([1.0, 1.0, 1.0]),
        nsteps=np.array([3, 3, 3]),
    )
    grid.cell_properties["rock"] = np.arange(8).astype(np.int64)

    file_name = tmp_path / "cellgrid"
    result = _write_structured_grid_gocad(grid, str(file_name))

    assert result is True
    content = _read(tmp_path / "cellgrid.vo")
    # cell properties are exported on the (nsteps - 1) grid of cell centres
    assert "AXIS_N 2 2 2" in content
    raw = np.fromfile(tmp_path / "cellgrid_rock@@", dtype=np.int8)
    assert raw.shape[0] == 8


def test_write_structured_grid_gocad_prefers_point_properties_and_warns(tmp_path, caplog):
    grid = StructuredGrid(
        origin=np.array([0.0, 0.0, 0.0]),
        step_vector=np.array([1.0, 1.0, 1.0]),
        nsteps=np.array([3, 3, 3]),
    )
    grid.properties["val"] = np.arange(27).astype(float)
    grid.cell_properties["ignored"] = np.arange(8).astype(float)

    # LoopStructural's getLogger() sets `propagate = False` on every logger it
    # creates, so records never reach the root logger that caplog listens on
    # by default. Attach caplog's handler directly to the module logger to
    # work around that.
    module_logger = logging.getLogger("LoopStructural.export.gocad")
    module_logger.addHandler(caplog.handler)
    previous_level = module_logger.level
    module_logger.setLevel(logging.WARNING)
    try:
        with caplog.at_level("WARNING"):
            file_name = tmp_path / "bothgrid"
            result = _write_structured_grid_gocad(grid, str(file_name))
    finally:
        module_logger.removeHandler(caplog.handler)
        module_logger.setLevel(previous_level)

    assert result is True
    assert not (tmp_path / "bothgrid_ignored@@").exists()
    assert (tmp_path / "bothgrid_val@@").exists()
    assert any("cell_properties were not exported" in message for message in caplog.messages)


def test_write_structured_grid_gocad_no_properties_raises(tmp_path):
    grid = StructuredGrid(
        origin=np.array([0.0, 0.0, 0.0]),
        step_vector=np.array([1.0, 1.0, 1.0]),
        nsteps=np.array([2, 2, 2]),
    )
    with pytest.raises(ValueError, match="no properties to export"):
        _write_structured_grid_gocad(grid, str(tmp_path / "empty"))


def test_write_structured_grid_gocad_sanitises_property_names(tmp_path):
    grid = StructuredGrid(
        origin=np.array([0.0, 0.0, 0.0]),
        step_vector=np.array([1.0, 1.0, 1.0]),
        nsteps=np.array([2, 2, 2]),
    )
    grid.properties["weird name!"] = np.arange(8).astype(float)

    file_name = tmp_path / "weird"
    _write_structured_grid_gocad(grid, str(file_name))

    assert (tmp_path / "weird_weird_name@@").exists()


# ---------------------------------------------------------------------------
# _write_feat_surfs_gocad
# ---------------------------------------------------------------------------


def test_write_feat_surfs_gocad_basic(tmp_path):
    surf = Surface(
        vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        triangles=np.array([[0, 1, 2]]),
        name="TestSurf",
    )
    file_name = tmp_path / "surf"
    result = _write_feat_surfs_gocad(surf, str(file_name))

    assert result is True
    content = _read(tmp_path / "surf.ts")
    assert "GOCAD TSurf 1" in content
    assert "name: TestSurf" in content
    assert "VRTX 1 0.0 0.0 0.0" in content
    assert "TRGL 1 2 3" in content
    assert "PROPERTIES" not in content


def test_write_feat_surfs_gocad_with_properties(tmp_path):
    surf = Surface(
        vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        triangles=np.array([[0, 1, 2]]),
        name="TestSurf",
        properties={"myprop": np.array([1.0, 2.0, 3.0])},
    )
    file_name = tmp_path / "surf_with_props"
    result = _write_feat_surfs_gocad(surf, str(file_name))

    assert result is True
    content = _read(tmp_path / "surf_with_props.ts")
    assert "PROPERTIES myprop" in content
    assert "PROPERTY_CLASSES myprop" in content
    # each VRTX line should have the property value appended
    assert "VRTX 1 0.0 0.0 0.0 1.0" in content
    assert "VRTX 2 1.0 0.0 0.0 2.0" in content
    assert "VRTX 3 0.0 1.0 0.0 3.0" in content


def test_write_feat_surfs_gocad_skips_nan_vertices_and_touching_triangles(tmp_path):
    # Surface.__post_init__ removes NaN vertices (and any triangles that
    # reference them) automatically, so build the surface with only valid
    # vertices/triangles to test the file writer's own NaN handling logic in
    # isolation by constructing the vertices array by hand after the fact.
    surf = Surface(
        vertices=np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 1.0]]
        ),
        triangles=np.array([[0, 1, 2], [1, 2, 3]]),
        name="PartialSurf",
    )
    # Manually reintroduce a NaN vertex bypassing Surface's own cleanup, to
    # directly exercise _write_feat_surfs_gocad's own NaN-skip behaviour.
    surf.vertices[3] = [np.nan, np.nan, np.nan]

    file_name = tmp_path / "partial"
    result = _write_feat_surfs_gocad(surf, str(file_name))

    assert result is True
    content = _read(tmp_path / "partial.ts")
    # only 3 VRTX lines since the 4th vertex was NaN
    assert content.count("VRTX") == 3
    # triangle referencing the NaN vertex should be skipped, only the first
    # remains
    assert content.count("TRGL") == 1
    assert "TRGL 1 2 3" in content
