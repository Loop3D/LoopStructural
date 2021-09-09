from .model_plotter import BaseModelPlotter
from LoopStructural.utils import getLogger
from LoopStructural.utils import LoopImportError
from LoopStructural.modelling.features import GeologicalFeature
logger = getLogger(__name__)

import numpy as np
try:
    import lavavu
    from lavavu.vutils import is_notebook
#catch the import lavavu error and provide more information
except ImportError:
    raise LoopImportError('lavavu',additional_information="Please install lavavu: pip install lavavu")

_OPEN_VIEWERS = {}

def close_all():
    _OPEN_VIEWERS.clear()
    return True


class LavaVuModelViewer(BaseModelPlotter):
    def __init__(self,model, bounding_box=None, nsteps=None, **kwargs):
        if lavavu is None:
            logger.error("Lavavu isn't installed: pip install lavavu")
            return
        self._id_name = "{}-{}".format(str(hex(id(self))), len(_OPEN_VIEWERS))
        _OPEN_VIEWERS[self._id_name] = self
        self.lv = lavavu.Viewer(**kwargs)
        self.lv['orthographic'] = True
        self.objects = {}
        
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
        self.default_cmap = 'rainbow'

    def _add_surface(self,
                    vertices, 
                    faces, 
                    name,
                    colour='red', 
                    paint_with=None, 
                    paint_with_value=None,
                    min_property_val=None, 
                    max_property_val=None,
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
        surf = self.lv.triangles(name)
        surf.vertices(vertices)
        surf.indices(faces)
        if paint_with is None:
            surf.colours(colour)
        surf["opacity"] = kwargs.get('opacity',1)
        if paint_with_value is not None:
            paint_with = paint_with_value
        if paint_with is not None: 
            # add a property to the surface nodes for visualisation
            # calculate the mode value, just to get the most common value
            surfaceval = np.zeros(vertices.shape[0])
            if isinstance(paint_with,GeologicalFeature):
                surfaceval[:] = paint_with.evaluate_value(self.model.scale(vertices))
                surf.values(surfaceval, 'paint_with')
            if callable(paint_with):
                surfaceval[:] = paint_with(self.model.scale(vertices))
                surf.values(surfaceval, 'paint_with')
            if isinstance(paint_with,(float,int)):
                surfaceval[:] = paint_with
                surf.values(surfaceval, 'paint_with')
             
            surf["colourby"] = 'paint_with'     
            cmap = kwargs.get('cmap', self.default_cmap)          
            vmin = kwargs.get('vmin', min_property_val)
            vmax = kwargs.get('vmax', max_property_val)
            surf.colourmap(cmap, range=(vmin, vmax)) 

    def _add_points(self, points, name, value= None, **kwargs):
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
        if points.shape[0] < 1:
            raise ValueError("Points array must have at least one element")
        if name is None:
            name = 'Unnamed points'
        p = self.lv.points(name, **kwargs)
        p.vertices(points)
        if value is not None:
            p.values(value,'v')
            p['colourby'] = "v"
        
            vmin = kwargs.get('vmin',np.nanmin(value))
            vmax = kwargs.get('vmax',np.nanmax(value))

            logger.info('vmin {} and vmax {}'.format(vmin,vmax))
            p.colourmap(cmap, range=(vmin, vmax))

    def _add_vector_marker(self, location, vector, name, symbol_type='arrow',**kwargs):
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
        if location.shape[0] != vectors.shape[0]:
            raise ValueError("Location and vector arrays must be the same length")
        if location.shape[0] < 1:
            raise ValueError("Location array must have at least one element")
        if name is None:
            name = 'Unnamed points'
        if symbol_type == 'arrow':
            vectorfield = self.lv.vectors(name, **kwargs)
            vectorfield.vertices(location)
            vectorfield.vectors(vector)
        elif symbol_type == 'disk':
            scaleshapes = kwargs.get('scaleshapes',np.max(self.model.maximum-self.model.origin)*0.014)
            vector /= np.linalg.norm(vector, axis=1)[:, None]
            vectorfield = self.lv.shapes(name, scaleshapes=scaleshapes,shapelength=0,**kwargs)
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

