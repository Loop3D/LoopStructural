from .model_plotter import BaseModelPlotter
from ..utils import getLogger
from ..modelling.features import GeologicalFeature
from matplotlib.colors import to_rgb

logger = getLogger(__name__)

import numpy as np


class DashView(BaseModelPlotter):
    def __init__(self, model, bounding_box=None, nsteps=None, **kwargs):
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

    def _add_surface(
        self, vertices, faces, name, colour="red", paint_with=None, paint_with_value=None, **kwargs
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
        surfaceval = np.zeros(vertices.shape[0])
        if paint_with_value is not None:
            paint_with = paint_with_value
        if paint_with is not None:
            # add a property to the surface nodes for visualisation
            # calculate the mode value, just to get the most common value
            if isinstance(paint_with, GeologicalFeature):
                surfaceval[:] = paint_with.evaluate_value(self.model.scale(vertices, inplace=False))
            if callable(paint_with):
                surfaceval[:] = paint_with(self.model.scale(vertices))
            if isinstance(paint_with, (float, int)):
                surfaceval[:] = paint_with
        # write points to paraview, set the point data to the property being painted.
        # TODO allow multiple properties
        faces = np.hstack([np.ones((faces.shape[0], 1), dtype=int) + 2, faces]).flatten()
        vertices -= self.model.origin
        return {
            "vertices": vertices.flatten(),
            "elements": faces,
            "point_data": surfaceval,
            "colour": to_rgb(colour),
        }

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
