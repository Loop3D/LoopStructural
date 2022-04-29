from LoopStructural.utils import LoopImportError, LoopTypeError, LoopValueError

try:
    from LoopProjectFile import ProjectFile
except ImportError:
    raise LoopImportError("LoopProjectFile cannot be imported")

from .process_data import ProcessInputData
import numpy as np
import pandas as pd
from matplotlib.colors import to_hex
from LoopStructural.utils import getLogger

logger = getLogger(__name__)


class LoopProjectfileProcessor(ProcessInputData):
    def __init__(self, projectfile, use_thickness=None):
        if isinstance(projectfile, ProjectFile) == False:
            raise LoopTypeError("projectife must be of type ProjectFile")
        column_map = {"easting": "X", "northing": "Y", "altitude": "Z"}
        self.projectfile = projectfile
        contacts = self.projectfile.stratigraphyLocations
        orientations = self.projectfile.stratigraphyOrientations
        fault_orientations = self.projectfile.faultOrientations
        fault_locations = self.projectfile.faultLocations

        orientations.rename(columns=column_map, inplace=True)
        contacts.rename(columns=column_map, inplace=True)
        fault_locations.rename(columns=column_map, inplace=True)
        fault_orientations.rename(columns=column_map, inplace=True)
        fault_locations["fault_name"] = [
            f"Fault_{eventid}" for eventid in fault_locations["eventId"]
        ]
        fault_orientations["fault_name"] = [
            f"Fault_{eventid}" for eventid in fault_orientations["eventId"]
        ]
        thicknesses = dict(
            zip(
                projectfile["stratigraphicLog"].name,
                projectfile["stratigraphicLog"].thickness,
            )
        )
        fault_properties = self.projectfile.faultLog
        fault_properties.rename(
            columns={
                "avgDisplacement": "displacement",
                "influenceDistance": "minor_axis",
                "verticalRadius": "intermediate_axis",
                "horizontalRadius": "major_axis",
                # "name": "fault_name",
            },
            inplace=True,
        )
        colours = dict(
            zip(
                self.projectfile.stratigraphicLog.name,
                [
                    to_hex(c)
                    for c in self.projectfile.stratigraphicLog.loc[
                        :, ["colour1Red", "colour1Green", "colour1Blue"]
                    ].to_numpy()
                    / 255
                ],
            )
        )
        super().__init__(
            contacts=contacts,
            contact_orientations=orientations,
            stratigraphic_order=[
                ("sg", list(self.projectfile.stratigraphicLog.name))
            ],  # needs to be updated,
            thicknesses=thicknesses,
            fault_orientations=fault_orientations,
            fault_locations=fault_locations,
            fault_properties=fault_properties,
            fault_edges=[],  # list(fault_graph.edges),
            colours=colours,
            fault_stratigraphy=None,
            intrusions=None,
            use_thickness=use_thickness,
            origin=self.projectfile.origin,
            maximum=self.projectfile.maximum
            #                     fault_edge_properties=fault_edge_properties
        )
