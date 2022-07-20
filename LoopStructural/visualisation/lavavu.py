from pyrsistent import v
from .model_plotter import BaseModelPlotter
from LoopStructural.utils import getLogger
from LoopStructural.utils import LoopImportError
from LoopStructural.modelling.features import GeologicalFeature

logger = getLogger(__name__)

import numpy as np

try:
    import lavavu
    from lavavu.vutils import is_notebook
# catch the import lavavu error and provide more information
except ImportError:
    raise LoopImportError(
        "lavavu", additional_information="Please install lavavu: pip install lavavu"
    )

_OPEN_VIEWERS = {}


def close_all():
    _OPEN_VIEWERS.clear()
    return True


class LavaVuModelViewer(BaseModelPlotter):
    def __init__(self, model=None, bounding_box=None, nsteps=None, **kwargs):
        if lavavu is None:
            logger.error(
                "Cannot use LavaVuModelViewer: Lavavu isn't installed \n pip install lavavu"
            )
            return
        self._id_name = "{}-{}".format(str(hex(id(self))), len(_OPEN_VIEWERS))
        _OPEN_VIEWERS[self._id_name] = self
        self.lv = lavavu.Viewer(**kwargs)
        self.lv["orthographic"] = True
        self.objects = {}

        super().__init__(model)
        self.bounding_box = bounding_box
        self.nsteps = nsteps
        if model is not None:
            self.bounding_box = model.bounding_box
            self.nsteps = model.nsteps
            logger.debug("Using bounding box from model")
        if self.bounding_box is None or self.nsteps is None:
            logger.warning(
                "Limited functionality for plot, bounding box is not defined."
            )
        self.bounding_box = np.array(self.bounding_box)

    def _parse_kwargs(self, kwargs):
        """
        remove any None kwargs from the list
        """
        return {k: v for k, v in kwargs.items() if v is not None}

    def _add_surface(
        self,
        vertices: np.ndarray,
        faces: np.ndarray,
        name: str,
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
        kwargs = self._parse_kwargs(kwargs)
        surf = self.lv.triangles(name)
        surf.vertices(vertices)
        surf.indices(faces)
        if paint_with is None:
            surf.colours(colour)
        surf["opacity"] = kwargs.get("opacity", 1)
        if paint_with_value is not None:
            paint_with = paint_with_value
        if paint_with is not None:
            # add a property to the surface nodes for visualisation
            # calculate the mode value, just to get the most common value
            surfaceval = np.zeros(vertices.shape[0])
            if isinstance(paint_with, GeologicalFeature):
                # paint with a geological feature
                # TODO make sure everything that could be
                # a feature derives from the same base
                # class, currently structuralframes and
                # faults don't work here..
                # or just change to using __call__
                surfaceval[:] = paint_with.evaluate_value(
                    self.model.scale(vertices, inplace=False)
                )
                surf.values(surfaceval, "paint_with")
            if callable(paint_with):
                # paint with a callable function e.g. (xyz)->value
                surfaceval[:] = paint_with(self.model.scale(vertices))
                surf.values(surfaceval, "paint_with")
            if isinstance(paint_with, (float, int)):
                # paint with array a constant value
                surfaceval[:] = paint_with
                surf.values(surfaceval, "paint_with")
            surf["colourby"] = "paint_with"
            cmap = kwargs.get("cmap", self.default_cmap)
            vmin = kwargs.get("vmin", np.nanmin(surfaceval))
            vmax = kwargs.get("vmax", np.nanmax(surfaceval))
            surf.colourmap(cmap, range=(vmin, vmax))

    def _add_points(self, points: np.ndarray, name: str, value=None, c=None, **kwargs):
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
        kwargs = self._parse_kwargs(kwargs)
        if points.shape[0] < 1:
            raise ValueError("Points array must have at least one element")
        if name is None:
            name = "Unnamed points"
        p = self.lv.points(name, **kwargs)
        p.vertices(points)
        if value is None and c is not None:
            value = c
        if value is not None:
            p.values(value, "v")
            p["colourby"] = "v"

            vmin = kwargs.get("vmin", np.nanmin(value))
            vmax = kwargs.get("vmax", np.nanmax(value))

            logger.info("vmin {} and vmax {}".format(vmin, vmax))
            cmap = kwargs.get("cmap", self.default_cmap)
            p.colourmap(cmap, range=(vmin, vmax))

    def _add_vector_marker(
        self,
        location: np.ndarray,
        vector: np.ndarray,
        name: str,
        symbol_type="arrow",
        **kwargs
    ):
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
        kwargs = self._parse_kwargs(kwargs)
        if location.shape[0] != vector.shape[0]:
            raise ValueError("Location and vector arrays must be the same length")
        if location.shape[0] < 1:
            raise ValueError("Location array must have at least one element")
        if name is None:
            name = "Unnamed points"
        if symbol_type == "arrow":
            vectorfield = self.lv.vectors(name, **kwargs)
            vectorfield.vertices(location)
            vectorfield.vectors(vector)
        elif symbol_type == "disk":
            scaleshapes = kwargs.get(
                "scaleshapes", np.max(self.model.maximum - self.model.origin) * 0.014
            )
            vector /= np.linalg.norm(vector, axis=1)[:, None]
            vectorfield = self.lv.shapes(
                name, scaleshapes=scaleshapes, shapelength=0, **kwargs
            )
            vectorfield.vertices(location)
            vectorfield.vectors(vector)

    def interactive(self, popout=False):
        """
        Runs the lavavu viewer as either a jupyter notebook
        inline interactive viewer or as a separate window

        Returns
        -------

        """
        if is_notebook() and popout is False:
            self.lv.control.Panel()
            self.lv.control.ObjectList()
            self.lv.control.show()
        if not is_notebook() or popout:
            self.lv.control.Panel()
            self.lv.control.ObjectList()
            self.lv.interactive()

    def set_zscale(self, zscale):
        """Set the vertical scale for lavavu

        just a simple wrapper for lavavu modelscale([xscale,yscale,zscale])

        Parameters
        ----------
        zscale : float
            vertical scale
        """
        self.lv.modelscale([1, 1, zscale])

    def set_viewer_rotation(self, rotation):
        """
        Set the viewer rotation given a list of rotations x,y,z

        Parameters
        ----------
        rotation numpy array of 3 rotation

        Returns
        -------

        """
        self.lv.rotate(rotation)

    def save(self, fname, **kwargs):
        """
        Calls lavavu.Viewer.image to save the viewer current state as an image

        Parameters
        ----------
        fname - file name string including relative path
        kwargs - optional kwargs to give to lavavu e.g. transparent, resolution

        Returns
        -------

        """
        self.lv.image(fname, **kwargs)

    def export_to_webgl(self, fname, **kwargs):

        self.lv.webgl(fname, **kwargs)

    def display(self, fname=None, **kwargs):
        """
        Calls the lv object display function. Shows a static image of the viewer inline.

        Returns
        -------

        """
        if fname:
            self.lv.image(fname, **kwargs)

        self.lv.display()

    def image(self, name, **kwargs):
        """
        Calls the lv object image function to save the display state

        Parameters
        ----------
        name : string
            name of the image file to save
        kwargs

        Returns
        -------

        """
        self.lv.image(name)

    def image_array(self, **kwargs):
        """Return the current viewer image image data as a numpy array

        Returns
        -------
        image : np.array
            image as a numpy array
        """
        return self.lv.rawimage(**kwargs).data

    def rotatex(self, r):
        """
        Rotate the viewer in the x plane

        Parameters
        ----------
        r : double
            degrees to rotate, can be +ve or -ve

        Returns
        -------

        """
        self.lv.rotatex(r)

    def rotatey(self, r):
        """
        Rotate the viewer in the Y plane

        Parameters
        ----------
        r : double
            degrees to rotate, can be +ve or -ve

        Returns
        -------

        """
        self.lv.rotatey(r)

    def rotatez(self, r):
        """
        Rotate the viewer in the z plane

        Parameters
        ----------
        r : double
            degrees to rotate, can be +ve or -ve

        Returns
        -------

        """
        self.lv.rotatez(r)

    def rotate(self, r):
        """
        Rotate by a vector of rotation angles

        Parameters
        ----------
        r : list/numpy array
            a vector of rotations

        Returns
        -------

        """
        self.lv.rotate(r)

    @property
    def rotation(self):
        """Accessor for the viewer rotation
        Returns
        -------
        list
            x,y,z rotations
        """
        return self.lv["xyzrotate"]

    @rotation.setter
    def rotation(self, xyz):
        """Set the rotation of the viewer

        Parameters
        ----------
        xyz : list like
            x y z rotations
        """
        self.lv.rotation(xyz)

    @property
    def border(self):
        """The width of the border around the model area

        Returns
        -------
        border : double
            [description]
        """
        return self.lv["border"]

    @border.setter
    def border(self, border):
        """Setter for the border

        Parameters
        ----------
        border : double
            set the thickness of the border around objects
        """
        self.lv["border"] = border

    def clear(self):
        """Remove all objects from the viewer"""
        self.lv.clear()

    @property
    def camera(self):
        return self.lv.camera()

    @camera.setter
    def camera(self, camera):
        # for some reason lavavu complains if the keys
        # attributes aren't already in the dict, so add them and then
        # call camera
        for key, value in camera.items():
            self.lv[key] = value

        self.lv.camera(camera)

    @property
    def xmin(self):
        return self.lv["xmin"]

    @xmin.setter
    def xmin(self, xmin):
        self.lv["xmin"] = xmin

    @property
    def xmax(self):
        return self.lv["xmax"]

    @xmax.setter
    def xmax(self, xmax):
        self.lv["xmax"] = xmax

    @property
    def ymin(self):
        return self.lv["ymin"]

    @ymin.setter
    def ymin(self, ymin):
        self.lv["ymin"] = ymin

    @property
    def ymax(self):
        return self.lv["ymax"]

    @ymax.setter
    def ymax(self, ymax):
        self.lv["ymax"] = ymax

    @property
    def zmin(self):
        return self.lv["zmax"]

    @zmin.setter
    def zmin(self, zmin):
        self.lv["zmin"] = zmin

    @property
    def zmax(self):
        return self.lv["zmax"]

    @zmax.setter
    def zmax(self, zmax):
        self.lv["zmax"] = zmax
