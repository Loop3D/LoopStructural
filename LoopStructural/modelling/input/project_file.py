from ...utils import LoopImportError, LoopTypeError

try:
    from LoopProjectFile import ProjectFile
except ImportError:
    raise LoopImportError("LoopProjectFile cannot be imported")

from .process_data import ProcessInputData
from matplotlib.colors import to_hex
from ...utils import getLogger

logger = getLogger(__name__)


class LoopProjectfileProcessor(ProcessInputData):
    def __init__(self, projectfile, use_thickness=None):
        if not isinstance(projectfile, ProjectFile):
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
                "name": "fault_name",
            },
            inplace=True,
        )
        fault_locations = fault_properties.reset_index()[["fault_name", "eventId"]].merge(
            fault_locations, on="eventId"
        )
        fault_orientations = fault_properties.reset_index()[["fault_name", "eventId"]].merge(
            fault_orientations, on="eventId"
        )
        fault_properties.set_index("fault_name", inplace=True)
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
            maximum=self.projectfile.maximum,
            #                     fault_edge_properties=fault_edge_properties
        )
