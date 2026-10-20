"""Exercise real VTK pipelines against the workspace LoopStructural API."""

from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
import pyvista as pv
from loopstructuralvisualisation import Loop2DView, Loop3DView, StratigraphicColumnView

from LoopStructural import GeologicalModel
from LoopStructural.geometry import ValuePoints, VectorPoints
from LoopStructural.modelling.features import AnalyticalGeologicalFeature, StructuralFrame
from LoopStructural.modelling.features.fault import FaultSegment


@pytest.fixture
def model():
    model = GeologicalModel([10, 20, 30], [14, 26, 38])
    model.bounding_box.nsteps = np.array([5, 6, 7])
    feature = AnalyticalGeologicalFeature("strati", [0, 0, 1], [10, 20, 30], model=model)
    model["strati"] = feature
    column = model.stratigraphic_column
    column.clear(basement=False)
    column.add_unit("lower", thickness=4, colour="red")
    column.add_unit("upper", thickness=4, colour="blue")
    column.group_mapping["Group_0"] = "strati"
    return model


@pytest.fixture
def view(model):
    plotter = Loop3DView(model, off_screen=True, window_size=(160, 120))
    yield plotter
    plotter.close()


@pytest.mark.parametrize("slicer", [False, True])
def test_scalar_field_limits_and_scalar_bars(view, model, slicer):
    options = {"clim": (0, 8)}
    actor = view.plot_scalar_field(
        model["strati"], vmin=1, vmax=7, slicer=slicer, pyvista_kwargs=options
    )
    assert isinstance(actor, pv.Actor)
    assert actor.mapper.scalar_range == (1, 7)
    assert options == {"clim": (0, 8)}
    assert len(view.scalar_bars) == 0
    view.plot_scalar_field(model["strati"], show_scalar_bar=True)
    assert len(view.scalar_bars) == 1
    view.plot_scalar_field(model["strati"])
    assert len(view.scalar_bars) == 1  # Preserve an existing bar.


def test_scalar_field_limits_without_options(view, model):
    actor = view.plot_scalar_field(model["strati"], vmin=1, vmax=7)
    assert actor.mapper.scalar_range == (1, 7)


@pytest.mark.parametrize("slicer", [False, True])
def test_surface_paint_uses_world_coordinates(view, model, slicer):
    paint = AnalyticalGeologicalFeature("paint", [1, 0, 0], [10, 20, 30], model=model)
    options = {}
    actor = view.plot_surface(
        model["strati"], value=[4.0], paint_with=paint, slicer=slicer, pyvista_kwargs=options
    )
    assert isinstance(actor, pv.Actor)
    mesh = actor.mapper.dataset
    np.testing.assert_allclose(mesh["values"], mesh.points[:, 0] - 10, atol=1e-6)
    assert options == {}
    assert len(view.scalar_bars) == 0


def test_empty_surface(view, model):
    assert view.plot_surface(model["strati"], value=[100.0]) is None


@pytest.mark.parametrize("slicer", [False, True])
def test_block_model(view, model, slicer):
    options = {}
    actor = view.plot_block_model(slicer=slicer, pyvista_kwargs=options)
    assert isinstance(actor, pv.Actor)
    assert options == {}
    assert len(view.scalar_bars) == 0
    assert {0, 1}.issubset(np.unique(actor.mapper.dataset["stratigraphy"]))


def test_model_surfaces_and_vectors(view, model):
    actors = view.plot_model_surfaces(faults=False)
    assert actors and all(isinstance(a, pv.Actor) for a in actors)
    assert len(view.scalar_bars) == 0
    actor = view.plot_vector_field(model["strati"], normalise=True)
    assert actor.mapper.dataset.n_points > 0


def test_plot_data_names_and_frame(view, model, monkeypatch):
    feature = model["strati"]
    points = np.array([[11.0, 21.0, 31.0], [12.0, 22.0, 32.0]])
    data = [
        ValuePoints(points, np.array([1.0, 2.0]), name="observations"),
        VectorPoints(points, np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]), name="normals"),
    ]
    monkeypatch.setattr(feature, "get_data", lambda: data)
    frame = StructuralFrame("frame", [feature], model=model)
    actors = view.plot_data(frame)
    assert len(actors) == 2
    assert "observations_values" in view.actors
    assert "normals_vectors" in view.actors


def test_offscreen_render(view, model):
    view.plot_surface(model["strati"], value=[4.0])
    image = view.screenshot(return_img=True)
    assert image.shape[:2] == (120, 160)
    assert np.any(image != image[0, 0])


def test_block_model_api_screenshot(model, tmp_path, monkeypatch):
    from loopstructuralvisualisation.api import plot_block_model

    monkeypatch.setattr(pv, 'OFF_SCREEN', True)
    filename = tmp_path / 'block-model.png'
    plot_block_model(model, filename=filename)
    assert filename.is_file()
    assert filename.stat().st_size > 0


def test_fault_plotting(view, model):
    coordinates = [
        AnalyticalGeologicalFeature(f"fault_{i}", direction, [12, 23, 34], model=model)
        for i, direction in enumerate(np.eye(3))
    ]
    fault = FaultSegment("fault", coordinates, model=model)
    actors = view.plot_fault(fault)
    assert len(actors) == 3
    assert all(isinstance(actor, pv.Actor) for actor in actors)
    actor = view.plot_fault_displacements(fault_list=[fault], name="test_displacement")
    assert "test_displacement" in view.actors
    mesh = actor.mapper.dataset
    np.testing.assert_allclose(
        mesh["displacement"], fault.displacementfeature.evaluate_value(mesh.points)
    )


def test_map_uses_world_coordinates_and_column(model):
    view = Loop2DView(model)
    feature = AnalyticalGeologicalFeature("x", [1, 0, 0], [10, 20, 30], model=model)
    try:
        image = view.add_scalar_field(feature, z=34)
        np.testing.assert_allclose(image.get_array(), (view.xx - 10).reshape(view.nsteps).T)
        assert view.add_contour(feature, [2.0], z=34) is not None
        assert view.add_model(z=32) is not None
    finally:
        plt.close(view.ax.figure)


def test_stratigraphic_column_custom_cmap(model):
    view = StratigraphicColumnView(model, cmap="viridis")
    try:
        assert view.plot() is not None
        assert len(view.ax.collections) == 1
    finally:
        plt.close(view.ax.figure)


def test_map_data_and_fault_ellipses(model):
    feature = model["strati"]
    empty = np.empty((0, 7))
    feature.builder = SimpleNamespace(
        get_value_constraints=lambda: np.array([[11.0, 22.0, 33.0, 3.0]]),
        get_gradient_constraints=lambda: np.array([[11.0, 22.0, 33.0, 1.0, 1.0, 0.0, 1.0]]),
        get_norm_constraints=lambda: empty,
    )
    fault = SimpleNamespace(
        fault_centre=np.array([12.0, 23.0, 34.0]),
        fault_strike_vector=np.array([0.0, 1.0, 0.0]),
        fault_major_axis=2.0,
        fault_minor_axis=1.0,
    )
    view = Loop2DView(model)
    try:
        view.add_data(feature, dip=False)
        np.testing.assert_allclose(view.ax.collections[0].get_offsets(), [[11.0, 22.0]])
        view.add_fault_ellipse(faults=[fault, fault], edgecolor="red")
        assert len(view.ax.patches) == 2
        assert view.ax.patches[0].center == (12.0, 23.0)
        assert view.ax.patches[0].angle == 90
    finally:
        plt.close(view.ax.figure)
