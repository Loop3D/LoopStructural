"""Optional notebook integration, without launching a web server."""

import pytest
import pyvista as pv

pytest.importorskip("trame")
pytest.importorskip("trame.widgets.vuetify3")

from loopstructuralvisualisation import Loop3DView
from loopstructuralvisualisation.trame import initialize
from loopstructuralvisualisation.trame.ui.vuetify3 import LoopViewer
from trame.app import get_server


@pytest.mark.parametrize("plotter_class", [Loop3DView, pv.Plotter])
def test_notebook_initialization(plotter_class):
    plotter = plotter_class(off_screen=True)
    try:
        plotter.add_mesh(pv.Sphere(), name="sphere")
        server = get_server(f"compat_{plotter._id_name}", client_type="vue3")
        viewer = initialize(server, plotter, mode="client")
        assert viewer.layout is not None
        if plotter_class is Loop3DView:
            assert isinstance(viewer, LoopViewer)
            name = next(k for k, a in plotter.actors.items() if isinstance(a, pv.Actor))
            viewer.toggle_visibility(**{f"{name}__visibility": False})
            assert not plotter.actors[name].visibility
            viewer.set_opacity(**{f"{name}__opacity": 0.5})
            assert plotter.actors[name].prop.opacity == 0.5
    finally:
        plotter.close()
