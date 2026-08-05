import pyvista as pv
import numpy as np
import re

from LoopStructural.datatypes import VectorPoints, ValuePoints
from LoopStructural.modelling.features import BaseFeature, StructuralFrame

from LoopStructural.modelling.features.fault import FaultSegment
from LoopStructural.datatypes import BoundingBox
from LoopStructural import GeologicalModel
from LoopStructural.utils import getLogger
from typing import Callable, Union, Optional, List

logger = getLogger(__name__)


class Loop3DView(pv.Plotter):
    def __init__(self, model=None, background='white', *args, **kwargs):
        """Loop3DView is a subclass of pyvista. Plotter that is designed to
        interface with the LoopStructural geological modelling package.

        Parameters
        ----------
        model : GeologicalModel, optional
            A loopstructural model used as reference for some methods, by default None
        background : str, optional
            colour for the background, by default 'white'
        """
        if 'shape' in kwargs:
            logger.warning('shape argument is not used in Loop3DView')
            kwargs.pop('shape')
        super().__init__(*args, **kwargs)
        self.set_background(background)
        self.model = model
        self.objects = {}

    def subplot(self, *args, **kwargs):
        logger.warning('subplot is not supported in Loop3DView')
        return self

    def add_mesh(self, *args, **kwargs):
        if 'name' not in kwargs or kwargs['name'] is None:
            name = 'unnamed_object'
            kwargs['name'] = name
            logger.warning(
                f'No name provided, using {name}. Pass name argument to add_mesh to remove this error'
            )
        kwargs['name'] = kwargs['name'].replace(' ', '_')
        kwargs['name'] = re.sub(r'[^a-zA-Z0-9_$]', '_', kwargs['name'])
        if kwargs['name'][0].isdigit():
            kwargs['name'] = 'ls_' + kwargs['name']
        if kwargs['name'][0] == '_':
            kwargs['name'] = 'ls' + kwargs['name']
        kwargs['name'] = self.increment_name(kwargs['name'])
        if '__opacity' in kwargs['name']:
            raise ValueError('Cannot use __opacity in name')
        if '__visibility' in kwargs['name']:
            raise ValueError('Cannot use __visibility in name')
        if '__control_visibility' in kwargs['name']:
            raise ValueError('Cannot use __control_visibility in name')
        return super().add_mesh(*args, **kwargs)

    def increment_name(self, name):
        parts = name.split('_')
        if len(parts) == 1:
            name = name + '_1'
        while name in self.actors:
            parts = name.split('_')
            try:
                parts[-1] = str(int(parts[-1]) + 1)
            except ValueError:
                parts.append('1')
            name = '_'.join(parts)
        return name

    def _check_model(self, model: GeologicalModel) -> GeologicalModel:
        """helper method to assign a geological model"""
        if model is None:
            model = self.model
        if model is None:
            raise ValueError("No model provided")
        return model

    def _get_vector_scale(self, scale: Optional[Union[float, int]]) -> float:
        autoscale = 1.0
        if self.model is not None:
            # automatically scale vector data to be 5% of the bounding box length
            autoscale = self.model.bounding_box.length.max() * 0.05
        if scale is None:
            scale = autoscale
        else:
            scale = scale * autoscale
        if scale > 10 * autoscale:
            logger.warning(
                "Vector scale magnitude greater than half of the model bounding box length, is this correct?"
            )

        return scale

    def plot_surface(
        self,
        geological_feature: BaseFeature,
        value: Optional[Union[float, int]] = None,
        paint_with: Optional[BaseFeature] = None,
        colour: Optional[str] = "red",
        cmap: Optional[str] = None,
        opacity: Optional[float] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        pyvista_kwargs: dict = {},
        show_scalar_bar: bool = False,
        slicer: bool = False,
        name: Optional[str] = None,
        bounding_box: Optional[BoundingBox] = None,
    ):
        """Add an isosurface of a geological feature to the model

        Parameters
        ----------
        geological_feature : BaseFeature
            The geological feature to plot
        value : Optional[Union[float, int, List[float]]], optional
            isosurface value, or list of values, by default average value of feature
        paint_with : Optional[BaseFeature], optional
            Paint the surface with the value of another geological feature, by default None
        colour : Optional[str], optional
            colour of the surface, by default "red"
        cmap : Optional[str], optional
            matplotlib colourmap, by default None
        opacity : Optional[float], optional
            opacity of the surface, by default None
        vmin : Optional[float], optional
            minimum value of the colourmap, by default None
        vmax : Optional[float], optional
            maximum value of the colourmap, by default None
        pyvista_kwargs : dict, optional
            other parameters passed to Plotter.add_mesh, by default {}
        name : Optional[str], optional
            name of the object, by default None
        slicer : bool, optional
            If an interactive plane slicing tool should be added, by default False
        show_scalar_bar : bool, optional
            Whether to show the scalar bar, by default False
        """

        if name is None:
            name = geological_feature.name + '_surfaces'
        name = self.increment_name(name)  # , 'surface')

        surfaces = geological_feature.surfaces(value, bounding_box=bounding_box)
        meshes = []
        for surface in surfaces:
            s = surface.vtk()
            if paint_with is not None:
                clim = [paint_with.min(), paint_with.max()]
                if vmin is not None:
                    clim[0] = vmin
                if vmax is not None:
                    clim[1] = vmax
                pyvista_kwargs["clim"] = clim
                pts = np.copy(surface.vertices)
                if self.model is not None:
                    pts = self.model.scale(pts)
                scalars = paint_with(pts)
                s["values"] = scalars
                s.set_active_scalars("values")
                colour = None
            meshes.append(s)
        mesh = pv.MultiBlock(meshes).combine()
        actor = None
        try:

            if slicer:
                actor = self.add_mesh_clip_plane(
                    mesh,
                    color=colour,
                    cmap=cmap,
                    opacity=opacity,
                    name=name,
                    **pyvista_kwargs,
                )
            else:
                actor = self.add_mesh(
                    mesh,
                    color=colour,
                    cmap=cmap,
                    opacity=opacity,
                    name=name,
                    **pyvista_kwargs,
                )

        except ValueError:
            logger.warning("No surfaces to plot")
        if paint_with is not None and not show_scalar_bar:
            self.remove_scalar_bar('values')
        return actor

    def plot_scalar_field(
        self,
        geological_feature: BaseFeature,
        cmap: str = "viridis",
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        opacity: Optional[float] = None,
        pyvista_kwargs: dict = {},
        show_scalar_bar: bool = False,
        slicer: bool = False,
        name: Optional[str] = None,
        bounding_box: Optional[BoundingBox] = None,
    ):
        """Plot a volume with the scalar field as the property
        calls feature.scalar_field() to get the scalar field and
        then pyvista add_mesh(feature.scalar_field().vtk())

        Parameters
        ----------
        geological_feature : BaseFeature
            The geological feature to plot the scalar field of
        cmap : str, optional
            matplotlib colourmap to use, by default "viridis"
        vmin : Optional[float], optional
            minimum value for cmap, by default None
        vmax : Optional[float], optional
            max value for cmap, by default None
        opacity : Optional[float], optional
            opacity of the object, by default None
        pyvista_kwargs : dict, optional
            additional kwargs sent to add_mesh, by default {}
        show_scalar_bar : bool, optional
            whether to show or hide the scalar bar, by default False
        slicer : bool, optional
            whether to plot using a plane slicer widget, by default False
        name : Optional[str], optional
            name for the object to appear in the object list, by default None

        Returns
        -------
        pv.Actor
            a reference to the actor that is added to the mesh
        """

        if name is None:
            name = geological_feature.name + '_scalar_field'
        name = self.increment_name(name)  # , 'scalar_field')

        volume = geological_feature.scalar_field(bounding_box=bounding_box).vtk()
        if vmin is not None:
            pyvista_kwargs["clim"][0] = vmin
        if vmax is not None:
            pyvista_kwargs["clim"][1] = vmax
        if slicer:
            actor = self.add_mesh_clip_plane(
                volume, cmap=cmap, opacity=opacity, name=name, **pyvista_kwargs
            )
        else:
            actor = self.add_mesh(volume, cmap=cmap, opacity=opacity, name=name, **pyvista_kwargs)
        if not show_scalar_bar:
            self.remove_scalar_bar(geological_feature.name)
        return actor

    def plot_block_model(
        self,
        cmap=None,
        model=None,
        pyvista_kwargs={},
        show_scalar_bar: bool = False,
        slicer: bool = False,
        threshold: Optional[Union[float, List[float]]] = None,
        name: Optional[str] = None,
    ):
        """Plot a voxel model where the stratigraphic id is the active scalar.
        It will use the colours defined in the stratigraphic column of the model
        unless a cmap is provided.
        Min/max range of cmap are defined by the min/max values of the stratigraphic ids or if
        clim is provided in pyvista_kwargs

        Parameters
        ----------
        cmap : str, optional
            matplotlib cmap string, by default None
        model : GeologicalModel, optional
            the model to pass if it is not the active geologicalmodel, by default None
        pyvista_kwargs : dict, optional
            additional arguments to be passed to pyvista add_mesh, by default {}
        show_scalar_bar : bool, optional
            whether show/hide the scalar bar, by default False
        slicer : bool, optional
            If an interactive plane slicing tool should be added, by default False
        threshold : Optional[Union[float, List[float]]], optional
            Whether to threshold values of the stratigraphy. Uses same syntax as pyvista threshold., by default None
        """
        model = self._check_model(model)
        if name is None:
            name = 'block_model'
        name = self.increment_name(name)  # , 'block_model')
        block, codes = model.get_block_model()
        block = block.vtk()
        block.set_active_scalars('stratigraphy')
        actor = None
        if cmap is None:
            cmap = self._build_stratigraphic_cmap(model)
        if "clim" not in pyvista_kwargs:
            pyvista_kwargs["clim"] = (np.min(block['stratigraphy']), np.max(block['stratigraphy']))
        if threshold is not None:
            if isinstance(threshold, float):
                block = block.threshold(threshold)
            elif isinstance(threshold, (list, tuple, np.ndarray)) and len(threshold) == 2:
                block = block.threshold((threshold[0], threshold[1]))
        if slicer:
            actor = self.add_mesh_clip_plane(block, cmap=cmap, name=name, **pyvista_kwargs)
        else:
            actor = self.add_mesh(block, cmap=cmap, name=name, **pyvista_kwargs)

        if not show_scalar_bar:
            self.remove_scalar_bar('stratigraphy')
        return actor

    def plot_fault_displacements(
        self,
        fault_list: Optional[List[FaultSegment]] = None,
        bounding_box: Optional[BoundingBox] = None,
        model=None,
        cmap="rainbow",
        pyvista_kwargs={},
        show_scalar_bar: bool = False,
        name: Optional[str] = None,
    ):
        """Plot the dispalcement magnitude for faults in the model
        on a voxel block

        Parameters
        ----------
        fault_list : _type_, optional
            list of faults to plot the model, by default None
        bounding_box : _type_, optional
            _description_, by default None
        model : _type_, optional
            _description_, by default None
        cmap : str, optional
            _description_, by default "rainbow"
        pyvista_kwargs : dict, optional
            _description_, by default {}
        show_scalar_bar : bool, optional
            _description_, by default False
        """
        if name is None:
            name = 'fault_displacement'
        name = self.increment_name(name)  # , 'fault_displacement_map')
        if fault_list is None:
            model = self._check_model(model)
            fault_list = model.faults
        if bounding_box is None:
            model = self._check_model(model)
            bounding_box = model.bounding_box
        pts = bounding_box.regular_grid()
        displacement_value = np.zeros(pts.shape[0])
        for f in fault_list:
            disp = f.displacementfeature.evaluate_value(bounding_box.vtk().points)
            displacement_value[~np.isnan(disp)] += disp[~np.isnan(disp)]
        volume = bounding_box.vtk()
        volume['displacement'] = displacement_value
        actor = self.add_mesh(volume, cmap=cmap, **pyvista_kwargs)
        if not show_scalar_bar:
            self.remove_scalar_bar('displacement')
        return actor

    def plot_model_surfaces(
        self,
        strati: bool = True,
        faults: bool = True,
        cmap: Optional[str] = None,
        model: Optional[GeologicalModel] = None,
        fault_colour: str = "black",
        pyvista_kwargs: dict = {},
        show_scalar_bar: bool = False,
        name: Optional[str] = None,
    ):
        """Plot the surfaces of the model

        Parameters
        ----------
        strati : bool, optional
            should stratigraphy surfaces be plotted, by default True
        faults : bool, optional
            should faults be plotted, by default True
        cmap : Optional[str], optional
            What cmap to use for the stratigraphy ids, by default None
        model : Optional[GeologicalModel], optional
            a GeologicalModel, if not provided will use self.model, by default None
        fault_colour : str, optional
            colour for the fault surfaces, by default "black"
        pyvista_kwargs : dict, optional
            Additional kwargs to send to add_mesh, by default {}
        show_scalar_bar : bool, optional
            whether to add the scalar bar, by default False
        name : Optional[str], optional
            name to add objects to object list with, by default None

        Returns
        -------
        pv.Actor
            The actor that is added to the scene
        """
        model = self._check_model(model)

        actors = []
        if strati:
            strati_surfaces = []
            surfaces = model.get_stratigraphic_surfaces()
            if cmap is None:
                cmap = model.stratigraphic_column.cmap().colors
            for s in surfaces:
                strati_surfaces.append(s.vtk())
            if name is None:
                object_name = 'model_surfaces'
            else:
                object_name = f'{name}_model_surfaces'
            object_name = self.increment_name(object_name)  # , 'model_surfaces')
            actors.append(
                self.add_mesh(
                    pv.MultiBlock(strati_surfaces).combine(),
                    cmap=cmap,
                    name=object_name,
                    **pyvista_kwargs,
                )
            )
            if not show_scalar_bar:
                self.remove_scalar_bar()
        if faults:
            fault_list = model.get_fault_surfaces()
            for f in fault_list:
                if name is None:
                    object_name = f'{f.name}_surface'
                if name is not None:
                    object_name = f'{name}_{f.name}_surface'
                object_name = self.increment_name(object_name)  # , 'fault_surfaces')
                actors.append(
                    self.add_mesh(f.vtk(), color=fault_colour, name=object_name, **pyvista_kwargs)
                )
        return actors

    def plot_vector_field(
        self,
        geological_feature: BaseFeature,
        scale: Optional[float] = None,
        name: Optional[str] = None,
        geom='arrow',
        scalars: Optional[np.ndarray] = None,
        normalise: bool = False,
        scale_function: Optional[Callable[[np.ndarray], np.ndarray]] = None,
        pyvista_kwargs: dict = {},
        bounding_box: Optional[BoundingBox] = None,
    ) -> pv.Actor:
        """Plot a vector field

        Parameters
        ----------
        geological_feature : BaseFeature
            Geological feature to plot the vector field of
        scale : float, optional
            magnitude scale for the glyphs, by default 1.0
        name : Optional[str], optional
            name for the viewer object list, by default None
        pyvista_kwargs : dict, optional
            additional kwargs to pass to add_mesh, by default {}

        Returns
        -------
        pv.Actor
            actor that is added to the scene
        """
        if name is None:
            name = geological_feature.name + '_vector_field'
        name = self.increment_name(name)  # , 'vector_field')
        vectorfield = geological_feature.vector_field(bounding_box=bounding_box)
        scale = self._get_vector_scale(scale)
        return self.add_mesh(
            vectorfield.vtk(
                scale=scale,
                geom=geom,
                normalise=normalise,
                scalars=scalars,
                scale_function=scale_function,
            ),
            name=name,
            **pyvista_kwargs,
        )

    def plot_data(
        self,
        feature: Union[BaseFeature, StructuralFrame],
        value: bool = True,
        vector: bool = True,
        scale: Optional[Union[float, int]] = None,
        geom: str = "arrow",
        name: Optional[str] = None,
        scalars: Optional[np.ndarray] = None,
        normalise: bool = True,
        pyvista_kwargs: dict = {},
    ) -> List[pv.Actor]:
        """Add the data associated with a feature to the plotter

        Parameters
        ----------
        feature : Union[BaseFeature, StructuralFrame]
            feature to add data from
        value : bool, optional
            whether to add value data, by default True
        vector : bool, optional
            whether to plot vector data, by default True
        scale : Union[float, int], optional
            vector scale, by default 1
        geom : str, optional
            vector glyph, by default "arrow"
        name : Optional[str], optional
            name to use in object list, by default None
        normalise: bool, optional
            normalise the vectors to be unit norm, by default True
        pyvista_kwargs : dict, optional
            additional kwargs to pass to pyvista add_mesh, by default {}

        Returns
        -------
        List[pv.Actor]
            list of actors added to the pv plotter

        Notes
        ------
        When plotting a vector the bounding box is used to scale the vectors. By default
        the length of the arrows will be 5% of the bounding box. The scale parameter is a
        multiplier for this value. If you sent normalise to False the vectors will not be normalised

        """
        if issubclass(type(feature), BaseFeature):
            feature = [feature]
        logger.info(f"Scale vectors by {scale}")
        scale = self._get_vector_scale(scale)
        logger.info(f"Vector scale is {scale}")
        actors = []
        bb = self.model.bounding_box if self.model is not None else None
        for f in feature:
            for d in f.get_data():
                if isinstance(d, ValuePoints):
                    if value:
                        if name is None:
                            object_name = d.name + '_values'
                        else:
                            object_name = f'{d.name}_values_{name}'
                        object_name = self.increment_name(object_name)  # , 'values')
                        actors.append(
                            self.add_mesh(
                                d.vtk(scalars=scalars), name=object_name, **pyvista_kwargs
                            )
                        )
                if isinstance(d, VectorPoints):
                    if vector:
                        if name is None:
                            object_name = d.name + '_vectors'
                        else:
                            object_name = f'{d.name}_vectors_{name}'
                        object_name = self.increment_name(object_name)  # , 'vectors')
                        actors.append(
                            self.add_mesh(
                                d.vtk(
                                    geom=geom,
                                    scale=scale,
                                    scalars=scalars,
                                    bb=bb,
                                    tolerance=None,
                                    normalise=normalise,
                                ),
                                name=name,
                                **pyvista_kwargs,
                            )
                        )
        return actors

    def plot_fold(self, folded_feature: BaseFeature, pyvista_kwargs={}):

        # folded_feature.
        pass

    def plot_fault(
        self,
        fault: FaultSegment,
        surface: bool = True,
        slip_vector: bool = True,
        displacement_scale_vector: bool = True,
        fault_volume: bool = True,
        vector_scale: Optional[Union[float, int]] = None,
        name: Optional[str] = None,
        geom: str = "arrow",
        pyvista_kwargs: dict = {},
        bounding_box: Optional[BoundingBox] = None,
    ) -> List[pv.Actor]:
        """Plot a fault including the surface, slip vector and displacement volume

        Parameters
        ----------
        fault : FaultSegment
            the fault to plot
        surface : bool, optional
            flag for the 0.0 surface, by default True
        slip_vector : bool, optional
            flag for scaled vector field, by default True
        displacement_scale_vector : bool, optional
            _description_, by default True
        fault_volume : bool, optional
            fault displacement scalar field, by default True
        vector_scale : Union[float, int], optional
            scale factor for vectors, by default 200
        name : Optional[str], optional
            name of the object for pyvista, by default None
        pyvista_kwargs : dict, optional
            additional kwargs for the pyvista plotter, by default {}

        Returns
        -------
        List[pv.Actor]
            list of actors added to the plot
        """
        actors = []
        if surface:
            if name is None:
                surface_name = fault.name + '_surface'
            else:
                surface_name = f'{fault.name}_surface_{name}'
            surface_name = self.increment_name(surface_name)
            surf = fault.surfaces([0], bounding_box=bounding_box)[0]
            actors.append(self.add_mesh(surf.vtk(), name=surface_name, **pyvista_kwargs))
        if slip_vector:
            if name is None:
                vector_name = fault.name + '_vector'
            else:
                vector_name = f'{fault.name}_vector_{name}'
            vector_name = self.increment_name(vector_name)

            vectorfield = fault.vector_field(bounding_box=bounding_box)
            vector_scale = self._get_vector_scale(vector_scale)
            actors.append(
                self.add_mesh(
                    vectorfield.vtk(scale=vector_scale, normalise=False),
                    name=vector_name,
                    **pyvista_kwargs,
                )
            )
        if fault_volume:
            if name is None:
                volume_name = fault.name + '_volume'
            else:
                volume_name = f'{fault.name}_volume_{name}'

            volume = fault.displacementfeature.scalar_field(bounding_box=bounding_box)

            volume = volume.vtk().threshold([-1.0, 1.0])
            if geom == "arrow":
                geom = pv.Arrow()
            elif geom == "disc":
                geom = pv.Disc()
                geom = geom.rotate_y(90)
            else:
                raise ValueError(f"Unknown glyph type {geom}")
            actors.append(self.add_mesh(volume, name=volume_name, **pyvista_kwargs))
        if len(actors) == 0:
            logger.warning(f"Nothing added to plot for {fault.name}")
        return actors

    def plot_fault_ellipsoid(
        self, fault: FaultSegment, name: Optional[str] = None, pyvista_kwargs: dict = {}
    ) -> pv.Actor:
        """Plot the fault ellipsoid

        Parameters
        ----------
        fault : FaultSegment
            the fault to plot
        name : Optional[str], optional
            name of the object for pyvista, by default None
        pyvista_kwargs : dict, optional
            additional kwargs for the pyvista plotter, by default {}

        Returns
        -------
        pv.Actor
            actor added to the plot
        """
        if name is None:
            name = fault.name + '_ellipsoid'
        name = self.increment_name(name)
        ellipsoid = fault.fault_ellipsoid()
        return self.add_mesh(ellipsoid, name=name, **pyvista_kwargs)

    def rotate(self, angles: np.ndarray):
        """Rotate the camera by the given angles
        order is roll, azimuth, elevation as defined by
        pyvista

        Parameters
        ----------
        angles : np.ndarray
            roll, azimuth, elevation
        """
        self.camera.roll += angles[0]
        self.camera.azimuth += angles[1]
        self.camera.elevation += angles[2]

    def display(self):
        self.show(interactive=False)
