from pyvista.trame.ui import UI_TITLE
from pyvista.trame.ui import get_viewer
from .ui.vuetify3 import LoopViewer as Viewer
from .. import Loop3DView


def initialize(
    server,
    plotter,
    mode=None,
    default_server_rendering=True,
    collapse_menu=False,
    **kwargs,
):  # numpydoc ignore=PR01,RT01
    """Generate the UI for a given plotter."""
    state = server.state
    state.trame__title = UI_TITLE
    if issubclass(type(plotter), Loop3DView):
        # only use the loopviewer if the plotter is a Loop3DView
        viewer = Viewer(plotter, server=server)

    else:
        # if pyvista use trame.ui.Viewer
        viewer = get_viewer(
            plotter,
            server=server,
            suppress_rendering=mode == "client",
        )
    with viewer.make_layout(server, template_name=plotter._id_name) as layout:
        viewer.layout = layout
        viewer.ui(
            mode=mode,
            default_server_rendering=default_server_rendering,
            collapse_menu=collapse_menu,
            **kwargs,
        )
        if issubclass(type(plotter), Loop3DView):
            viewer.object_menu()

    return viewer
