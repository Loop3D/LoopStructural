# ruff: noqa: D102
"""PyVista Trame Viewer class for a Vue 3 client.

This class, derived from `pyvista.trame.ui.base_viewer`,
is intended for use with a trame application where the client type is "vue3".
Therefore, the `ui` method implemented by this class utilizes the API of Vuetify 3.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from trame.ui.vuetify3 import SinglePageWithDrawerLayout
from trame.widgets import vuetify3 as vuetify

import pyvista
from pyvista.trame.ui.vuetify3 import Viewer


if TYPE_CHECKING:  # pragma: no cover
    from trame_client.ui.core import AbstractLayout


class LoopViewer(Viewer):
    def __init__(self, *args, **kwargs):
        """Overwrite the pyvista trame layout to use a singlepage layout
        and add an object visibility menu to the drawer
        """
        super().__init__(*args, **kwargs)

    def make_layout(self, *args, **kwargs) -> AbstractLayout:

        return SinglePageWithDrawerLayout(*args, **kwargs)

    def ui(self, *args, **kwargs):
        with self.layout as layout:
            layout.title.set_text("LoopStructural Viewer")
        with self.layout.content:

            return super().ui(*args, **kwargs)

    def toggle_visibility(self, **kwargs):
        """Toggle the visibility of an object in the plotter.
        this is a slot called by the state change, the kwargs are the current state
        so we need to check the keys and update accordingly
        """
        for k in kwargs.keys():
            object_name = k.split("__visibility")[0]
            if object_name in self.plotter.actors:
                self.plotter.actors[object_name].visibility = kwargs[k]
        self.update()
        # self.actors[k].visibility = not self.actors[k].visibility

    def set_opacity(self, **kwargs):
        """Set the opacity of an object in the plotter.
        this is a slot called by the state change, the kwargs are the current state
        so we need to check the keys and update accordingly
        """
        for k in kwargs.keys():
            object_name = k.split("__opacity")[0]
            if object_name in self.plotter.actors:
                self.plotter.actors[object_name].prop.opacity = kwargs[k]
        self.update()
        # self.actors[k].visibility = not self.actors[k].visibility

    def object_menu(self):
        with self.layout.drawer as drawer:
            with vuetify.VCard():

                for k, a in self.plotter.actors.items():
                    if type(a) is not pyvista.plotting.actor.Actor:
                        continue
                    drawer.server.state[f"{k}__visibility"] = True
                    drawer.server.state[f"{k}__control_visibility"] = False
                    drawer.server.state[f"{k}__opacity"] = a.prop.opacity
                    drawer.server.state.change(f"{k}__visibility")(self.toggle_visibility)
                    drawer.server.state.change(f"{k}__opacity")(self.set_opacity)
                    with vuetify.VRow(
                        classes='pa-0 ma-0 align-center fill-height',
                        style='flex-wrap: nowrap',
                    ):
                        with vuetify.VCol():

                            vuetify.VCheckbox(
                                label=k,
                                classes="ma-0 pa-0",
                                v_model=(f"{k}__visibility"),
                                # click=(self.toggle_visibility("test")),
                            )
                        with vuetify.VCol():
                            vuetify.VBtn(
                                icon='mdi-dots-horizontal',
                                click=(f'{k}__control_visibility=!{k}__control_visibility'),
                                # click=(self.toggle_visibility("test")),
                            )
                    with vuetify.VCard(classes="ma-0 pa-0", v_show=f'{k}__control_visibility'):

                        with vuetify.VRow(
                            classes="ma-0 pa-0 d-flex align-center",
                        ):

                            vuetify.VSlider(
                                v_model=(f"{k}__opacity"),
                                label="Opacity",
                                min=0,
                                max=1,
                                step=0.1,
                                thumb_label=True,
                            )
