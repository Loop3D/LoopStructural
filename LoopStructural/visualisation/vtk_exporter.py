from .model_plotter import BaseModelPlotter
from LoopStructural.utils import getLogger
from LoopStructural.utils import LoopImportError
from LoopStructural.modelling.features import GeologicalFeature

logger = getLogger(__name__)

import numpy as np

try:
    import meshio
except ImportError:
    raise ImportError("Cannot use VtkExporter without meshio")


class VtkExporter(BaseModelPlotter):
    def __init__(self, model, file_name, bounding_box=None, nsteps=None, **kwargs):
        logger.info("VtkExporter")
        super().__init__(model)
        self.bounding_box = bounding_box
        self.nsteps = nsteps
        if model is not None:
            self.bounding_box = model.bounding_box
            self.nsteps = model.nsteps
            logger.debug("Using bounding box from model")
        if self.bounding_box is None or self.nsteps is None:
            logger.error("Plot area has not been defined.")
        self.bounding_box = np.array(self.bounding_box)
        self.file_name = file_name

    def _add_surface(
        self,
        vertices,
        faces,
        name,
        colour="red",
        paint_with=None,
        paint_with_value=None,
        **kwargs
    ):
        """Virtual function to be overwritten by subclasses for adding surfaces to the viewer

        Parameters
        ----------
        faces : numpy array
            indices of the triangles
        vertices : numy array
            vertices of the surface
        name : string
            name of the surface in the viewer
        colour : str, optional
            matplotlib colour, by default 'red'
        paint_with : GeologicalFeature, optional
            geological feature to evaluate on vertices, by default None
        """
        point_data = {}
        if paint_with_value is not None:
            paint_with = paint_with_value
        if paint_with is not None:
            # add a property to the surface nodes for visualisation
            # calculate the mode value, just to get the most common value
            surfaceval = np.zeros(vertices.shape[0])
            if isinstance(paint_with, GeologicalFeature):
                surfaceval[:] = paint_with.evaluate_value(
                    self.model.scale(vertices, inplace=False)
                )
            if callable(paint_with):
                surfaceval[:] = paint_with(self.model.scale(vertices))
            if isinstance(paint_with, (float, int)):
                surfaceval[:] = paint_with
            point_data = {"surfaceval": surfaceval}
        # write points to paraview, set the point data to the property being painted.
        # TODO allow multiple properties
        meshio.write_points_cells(
            "{}_{}.vtk".format(self.file_name, name.split(".")[0]),
            vertices,
            [("triangle", faces)],
            point_data=point_data,
        )

    def _add_points(self, points, name, value=None, **kwargs):
        """Virtual function to be overwritten by subclasses for adding points to the viewer

        Parameters
        ----------
        points : np.array
            location of the points
        name : str
            name of the points in the viewer
        value : np.array, optional
            value to assign to the points
        """
        logger.warning("Cannot export point data to vtk")

        # if points.shape[0] < 1:
        #     raise ValueError("Points array must have at least one element")
        # if name is None:
        #     name = 'Unnamed points'
        # point_data = {}
        # if value is not None:
        #     point_data = {"value":value}
        # meshio.write_points_cells('{}_{}.vtk'.format(self.file_name,name.split('.')[0]),
        # points,
        # point_data=point_data
        # )

    def _add_vector_marker(self, location, vector, name, symbol_type="arrow", **kwargs):
        """Virtual function to be overwritten by subclasses for adding vectors to the viewer

        Parameters
        ----------
        location : numpy array
            location array
        vector : numpy array
            vector component array
        symbol_type : str, optional
            type of glyph to display the vector by, by default 'arrow'
        name : string
            name of the object in the visualisation
        """
        logger.warning("Cannot export vector field to vtk")
        pass
        # if location.shape[0] != vector.shape[0]:
        #     raise ValueError("Location and vector arrays must be the same length")
        # if location.shape[0] < 1:
        #     raise ValueError("Location array must have at least one element")
        # if name is None:
        #     name = 'Unnamed points'
        # if symbol_type == 'arrow':
        #     vectorfield = self.lv.vectors(name, **kwargs)
        #     vectorfield.vertices(location)
        #     vectorfield.vectors(vector)
        # elif symbol_type == 'disk':
        #     scaleshapes = kwargs.get('scaleshapes',np.max(self.model.maximum-self.model.origin)*0.014)
        #     vector /= np.linalg.norm(vector, axis=1)[:, None]
        #     vectorfield = self.lv.shapes(name, scaleshapes=scaleshapes,shapelength=0,**kwargs)
        #     vectorfield.vertices(location)
        #     vectorfield.vectors(vector)
