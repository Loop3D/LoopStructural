from LoopStructural.utils import getLogger
from LoopStructural.utils import LoopImportError
from LoopStructural.modelling.features import FeatureType

logger = getLogger(__name__)


import numpy as np

try:
    from skimage.measure import marching_cubes
except ImportError:
    logger.warning("Using deprecated version of scikit-image")
    from skimage.measure import marching_cubes_lewiner as marching_cubes


from LoopStructural.modelling.features import (
    FeatureType,
    GeologicalFeature,
    BaseFeature,
)
from LoopStructural.utils.helper import create_surface, get_vectors, create_box


class BaseModelPlotter:
    def __init__(self, model=None):
        """

        Parameters
        ----------
        model
        """
        self.model = model
        self.default_vector_symbol = "disk"
        self.default_cmap = "rainbow"

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        if model is not None:
            self.bounding_box = np.array(model.bounding_box)
            self.nsteps = np.array(model.nsteps)
            self._model = model
            self._nelements = self.nsteps[0] * self.nsteps[1] * self.nsteps[2]
            logger.debug("Using bounding box from model")

    @property
    def nelements(self) -> int:
        """The number of elements to use for evaluating the isosurface

        Returns
        -------
        nelements : int
            number of elements to use for isosurfacing
        """
        return self._nelements

    @nelements.setter
    def nelements(self, nelements: int):
        """Setter for nelements, automatically caculates the number of equally sized elements
        to isosurface. Better than specifying step distance manually

        Parameters
        ----------
        nelements : int
            [description]
        """
        box_vol = (
            (self.bounding_box[1, 0] - self.bounding_box[0, 0])
            * (self.bounding_box[1, 1] - self.bounding_box[0, 1])
            * (self.bounding_box[1, 2] - self.bounding_box[0, 2])
        )
        ele_vol = box_vol / nelements
        # calculate the step vector of a regular cube
        step_vector = np.zeros(3)
        step_vector[:] = ele_vol ** (1.0 / 3.0)
        # step_vector /= np.array([1,1,2])
        # number of steps is the length of the box / step vector
        nsteps = np.ceil(
            (self.bounding_box[1, :] - self.bounding_box[0, :]) / step_vector
        ).astype(int)
        self.nsteps = nsteps
        logger.info(
            "Using grid with dimensions {} {} {}".format(
                nsteps[0], nsteps[1], nsteps[2]
            )
        )

    @property
    def nsteps(self) -> np.ndarray:
        return self._nsteps

    @nsteps.setter
    def nsteps(self, nsteps: np.ndarray):
        self._nsteps = np.array(nsteps)

    def _add_surface(
        self, tri, vertices, name, colour="red", paint_with=None, **kwargs
    ):
        """Virtual function to be overwritten by subclasses for adding surfaces to the viewer

        Parameters
        ----------
        tri : numpy array
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
        pass

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
        pass

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
        pass

    def add_section(
        self,
        geological_feature=None,
        axis="x",
        value=None,
        paint_with=None,
        name=None,
        **kwargs,
    ):
        """

        Plot a section/map thru the model and paint with a geological feature

        Parameters
        ----------
        geological_feature : Geological feature
            The feature to paint the section with
        axis : string
            which axis, x,y,z
        value : float
            Where to make the section
        kwargs
            additional kwargs passes to lavavu for colourmaps etc

        Returns
        -------

        """
        if axis == "x":
            tri, yy, zz = create_surface(
                self.bounding_box[:, [1, 2]], self.nsteps[[1, 2]]
            )
            xx = np.zeros(zz.shape)
            if value is None:
                value = np.nanmean(self.bounding_box[:, 0])
            xx[:] = value
        if axis == "y":
            tri, xx, zz = create_surface(
                self.bounding_box[:, [0, 2]], self.nsteps[[0, 2]]
            )
            yy = np.zeros(xx.shape)
            if value is None:
                value = np.nanmean(self.bounding_box[:, 1])
            yy[:] = value
        if axis == "z":
            tri, xx, yy = create_surface(self.bounding_box[:, 0:2], self.nsteps[0:2])
            zz = np.zeros(xx.shape)
            if value is None:
                value = np.nanmean(self.bounding_box[:, 2])
            zz[:] = value
        if geological_feature == "model" and self.model is not None:
            name = kwargs.get("name", "model_section")
            if paint_with == None:
                paint_with = lambda xyz: self.model.evaluate_model(xyz, scale=False)
        elif geological_feature is not None and isinstance(
            geological_feature, BaseFeature
        ):
            name = kwargs.get("name", geological_feature.name)
            paint_with = geological_feature
        elif callable(geological_feature):
            if name is None:
                name = "unnamed_callable"
            paint_with = geological_feature

        name = "{}_section_at_{}_of_{}".format(axis, value, name)
        colour = kwargs.get("colour", "red")

        # create an array to evaluate the feature on for the section
        points = np.zeros((len(xx), 3))  #
        points[:, 0] = xx
        points[:, 1] = yy
        points[:, 2] = zz

        # set the surface to be painted with the geological feature, but if a painter is specified, use that instead
        # if 'paint_with' not in kwargs:
        #     kwargs['paint_with'] = geological_feature
        return self._add_surface(
            self.model.rescale(points, inplace=False),
            tri,
            name,
            colour=colour,
            paint_with=paint_with,
            **kwargs,
        )

    def add_isosurface(
        self,
        geological_feature,
        value=None,
        isovalue=None,
        paint_with=None,
        slices=None,
        colour="red",
        nslices=None,
        cmap=None,
        filename=None,
        names=None,
        colours=None,
        opacity=None,
        function=None,
        **kwargs,
    ):
        """Plot the surface of a geological feature


        Parameters
        ----------
        geological_feature : GeologicalFeature
            [description]
        value : float, optional
            value of the scalar field to isosurface
        isovalue : float, optional
            value of the scalar field to isosurface, by default None
        paint_with : GeologicalFeature, optional
            a geological feature to paint the surface with its evaluate_value results, by default None
        slices : list, optional
            values to isosurface, by default None
        colour : string, optional
            matplotlib color, by default None
        nslices : int, optional
            number of slices to evenly distribute in the model, by default None
        cmap : string, optional
            matplotlib colormap, by default None
        names: list, optional
            list of names same length as slices
        colours: list, optional
            list of colours same length as slices
        opacity: double, optional
            change the opacity of the surface(s)
        callback_function:
            called with verts, tri and surface name - e.g.
            callback_function(verts,tri,name)

        Returns
        -------
        [type]
            [description]
        """
        return_container = {}
        if geological_feature is None:
            logger.error("Cannot add isosurface GeologicalFeature does not exist")
        # update the feature to make sure its current

        # do isosurfacing of support using marching tetras/cubes
        x = np.linspace(
            self.bounding_box[0, 0], self.bounding_box[1, 0], self.nsteps[0]
        )
        y = np.linspace(
            self.bounding_box[0, 1], self.bounding_box[1, 1], self.nsteps[1]
        )
        z = np.linspace(
            self.bounding_box[1, 2], self.bounding_box[0, 2], self.nsteps[2]
        )
        xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
        points = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        val = geological_feature.evaluate_value(points)
        mean_val = np.nanmean(val)  # geological_feature.mean()
        max_val = np.nanmax(val)  # geological_feature.max()
        min_val = np.nanmin(val)  # geological_feature.min()
        if paint_with is not None and "vmin" not in kwargs and "vmax" not in kwargs:
            paint_val = np.zeros(points.shape[0])
            if isinstance(paint_with, GeologicalFeature):
                paint_val = paint_with.evaluate_value(points)
            if callable(paint_with):
                paint_val = paint_with(points)
            # get the stats to check what we are plotting
            kwargs["vmin"] = np.nanmin(paint_val)  # geological_feature.min()
            kwargs["vmax"] = np.nanmax(paint_val)  # geological_feature.max()
        # set default parameters
        slices_ = [mean_val]
        painter = None
        voxet = None
        tris = None
        nodes = None
        # parse kwargs for parameters
        if isovalue is not None:
            slices_ = [isovalue]
        if value is not None:
            slices_ = [value]
        if slices is not None:
            slices_ = slices
        if nslices is not None:
            var = max_val - min_val
            # buffer slices by 5%
            slices_ = np.linspace(min_val + var * 0.05, max_val - var * 0.05, nslices)
        base_name = kwargs.pop("name", geological_feature.name)

        region = kwargs.get("region", None)
        if region is not None:
            val[
                ~region(np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T)
            ] = np.nan
        if self.model.dtm is not None:
            xyz = self.model.rescale(
                np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T, inplace=False
            )
            dtmv = self.model.dtm(xyz[:, :2])
            val[xyz[:, 2] > dtmv] = np.nan
        step_vector = np.array([x[1] - x[0], y[1] - y[0], z[1] - z[0]])
        for i, isovalue in enumerate(slices_):
            logger.info(
                "Creating isosurface of %s at %f" % (geological_feature.name, isovalue)
            )

            if isovalue > np.nanmax(val) or isovalue < np.nanmin(val):
                logger.warning(
                    f"{geological_feature.name}: Isovalue doesn't exist inside bounding box"
                )
                continue
            try:
                verts, faces, normals, values = marching_cubes(
                    val.reshape(self.nsteps, order="C"), isovalue, spacing=step_vector
                )
                verts += np.array(
                    [
                        self.bounding_box[0, 0],
                        self.bounding_box[0, 1],
                        self.bounding_box[1, 2],
                    ]
                )
                self.model.rescale(verts)

            except (ValueError, RuntimeError) as e:
                print(e)
                logger.warning(
                    "Cannot isosurface {} at {}, skipping".format(
                        geological_feature.name, isovalue
                    )
                )
                continue

            name = "{}_{}".format(base_name, isovalue)
            if names is not None and len(names) == len(slices_):
                name = names[i]

            if colours is not None and len(colours) == len(slices_):
                colour = colours[i]

            if function is not None:
                function(verts, faces, name)

            paint_with_value = None
            if paint_with == geological_feature:
                paint_with_value = isovalue
            return_container[name] = self._add_surface(
                verts,
                faces,
                name,
                colour=colour,
                opacity=opacity,
                paint_with=paint_with,
                paint_with_value=paint_with_value,
                cmap=cmap,
                **kwargs,
            )

        return return_container

    def add_scalar_field(
        self,
        geological_feature,
        name=None,
        cmap="rainbow",
        vmin=None,
        vmax=None,
        opacity=None,
        paint_with=None,
        **kwargs,
    ):
        """Add a block the size of the model area painted with the scalar field value

        Parameters
        ----------
        geological_feature : GeologicalFeature
            the geological feature to colour the scalar field by
        name : string, optional
            Name of the object for lavavu, needs to be unique for the viewer object, by default uses feature name
        cmap : str, optional
            mpl colourmap reference, by default 'rainbow'
        vmin : double, optional
            minimum value of the colourmap, by default None
        vmax : double, optional
            maximum value of the colourmap, by default None
        opacity : double, optional
            change the opacity of the block
        """
        if name == None:
            if geological_feature is None:
                name = "unnamed scalar field"
            else:
                name = geological_feature.name + "_scalar_field"

        points, tri = create_box(self.bounding_box, self.nsteps)
        if self.model is None:
            raise ValueError("Model not set")
        pts = self.model.rescale(points, inplace=False)
        if paint_with is None:
            paint_with = geological_feature
        return self._add_surface(
            pts,
            tri,
            name,
            paint_with=paint_with,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            opacity=opacity,
            **kwargs,
        )

    def add_box(self, bounding_box, name, colour="red", **kwargs):
        """Adds a box to the viewer defined by a bounding box

        Parameters
        ----------
        bounding_box : np.array
            [[xmin,ymin,zmin], [xmax,ymax,zmax]]
        name : string
            name of object in viewer
        colour : mpl readable colour format
            [description], by default 'red'
        """
        points, tri = create_box(bounding_box, self.nsteps)
        return self._add_surface(points, tri, name, colour, **kwargs)

    def add_model(self, cmap=None, **kwargs):
        """Add a block model painted by stratigraphic id to the viewer

        Calls self.model.evaluate_model() for a cube surrounding the model.

        Parameters
        ----------
        cmap : matplotlib cmap, optional
            colourmap name or object from mpl

        Notes
        ------
        It is sensible to increase the viewer step sizes before running this function to
        increase the resolution of the model as its not possible to interpolate a discrete
        colourmap and this causes the model to look like a lego block.
        You can update the model resolution by changing the attribute nsteps
        >>> viewer.nsteps = np.array([100,100,100])

        """
        name = kwargs.get("name", "geological_model")
        points, tri = create_box(self.bounding_box, self.nsteps)

        if cmap is None:
            try:
                import matplotlib.colors as colors
            except ImportError:
                logger.warning(
                    "Cannot use predefined colours as I can't import matplotlib"
                )
                cmap = "tab20"
            colours = []
            boundaries = []
            data = []
            for g in self.model.stratigraphic_column.keys():
                if g == "faults":
                    continue
                for u, v in self.model.stratigraphic_column[g].items():
                    data.append((v["id"], v["colour"]))
                    colours.append(v["colour"])
                    boundaries.append(v["id"])  # print(u,v)
            cmap = colors.ListedColormap(colours).colors

        return self._add_surface(
            self.model.rescale(points, inplace=False),
            tri,
            name,
            paint_with=lambda xyz: self.model.evaluate_model(xyz, scale=False),
            cmap=cmap,
            **kwargs,
        )

    def add_fault_displacements(self, cmap="rainbow", **kwargs):
        """Add a block model painted by the fault displacement magnitude

        Calls fault.displacementfeature.evaluate_value(points) for all faults

        Parameters
        ----------
        cmap : matplotlib cmap, optional
            colourmap name or object from mpl

        Notes
        ------
        It is sensible to increase the viewer step sizes before running this function to
        increase the resolution of the model as its not possible to interpolate a discrete
        colourmap and this causes the model to look like a lego block.
        You can update the model resolution by changing the attribute nsteps
        >>> viewer.nsteps = np.array([100,100,100])

        """

        name = kwargs.get("name", "fault_displacements")
        points, tri = create_box(self.bounding_box, self.nsteps)
        self._add_surface(
            tri,
            points,
            name,
            paint_with=self.model.evaluate_fault_displacements,
            cmap=cmap,
        )  # need to pass colour somehow

    def add_fault(self, fault, step=100):
        self.add_isosurface(fault, value=0, name=fault.name)
        self.add_vector_field(fault, locations=self.model.regular_grid()[::step])

    def add_structural_frame(self, frame, step=100, data=True, **kwargs):
        for i in range(3):
            self.add_isosurface(frame[i], slices=[-1, 0, 1], **kwargs)
            if data:
                self.add_data(frame[i])

    def unfault_grid(self, feature, grid=None):
        if grid is None:
            grid = self.model.regular_grid()
        # apply all faults associated with a feature to a regular grid
        self.add_value_data(
            self.model.rescale(grid, inplace=False),
            grid[:, 2],
            name="Regular grid before faults",
            pointsize=10,
        )

        for f in feature.faults:
            grid = f.apply_to_points(grid)
        self.add_value_data(
            self.model.rescale(grid, inplace=False),
            grid[:, 2],
            name="Regular grid after faults",
            pointsize=10,
        )

    def add_model_surfaces(
        self,
        strati=True,
        faults=True,
        cmap=None,
        fault_colour="black",
        displacement_cmap=None,
        **kwargs,
    ):
        """Add surfaces for all of the interfaces in the model


        Parameters
        ----------
        strati : bool, optional
            whether to draw stratigraphy
        faults : bool, optional
            whether to draw faults, by default True
        cmap : string
            matplotlib cmap
        fault_colour : string
            colour string for faults
        displacement_cmap : string/None
            if string is specified uses this cmap to colour
            faults by displacement
        Notes
        ------
        Other parameters are passed to self.add_isosurface()

        """
        try:
            from matplotlib import cm
            from matplotlib import colors
        except ImportError:
            logger.warning("Cannot add model surfaces without matplotlib \n")
            return
        from ..modelling.features import LambdaGeologicalFeature
        import time
        from tqdm.auto import tqdm

        return_container = {}
        start = time.time()
        logger.info("Updating model")
        self.model.update()
        logger.info("Model update took: {} seconds".format(time.time() - start))
        start = time.time()
        logger.info("Isosurfacing")
        n_units = 0  # count how many discrete colours
        name_suffix = kwargs.pop("name", "")
        if strati and self.model.stratigraphic_column:
            for g in self.model.stratigraphic_column.keys():
                if g in self.model.feature_name_index:
                    for u in self.model.stratigraphic_column[g].keys():
                        n_units += 1
        n_faults = 0
        for f in self.model.features:
            if f.type == FeatureType.FAULT:
                n_faults += 1

        if self.model.stratigraphic_column and cmap is None:

            colours = []
            boundaries = []
            data = []
            for g in self.model.stratigraphic_column.keys():
                if g == "faults":
                    # skip anything saved in faults here
                    continue
                for u, v in self.model.stratigraphic_column[g].items():
                    data.append((v["id"], v["colour"]))
                    colours.append(v["colour"])
                    boundaries.append(v["id"])
            cmap = colors.ListedColormap(colours)
        else:
            cmap = cm.get_cmap("tab20", n_units)
        ci = 0
        cmap_colours = colors.to_rgba_array(cmap.colors)
        n_surfaces = 0
        if strati:
            n_surfaces += n_units
        if faults:
            n_surfaces += n_faults
        with tqdm(total=n_surfaces) as pbar:

            if strati and self.model.stratigraphic_column:
                for g in self.model.stratigraphic_column.keys():
                    if g in self.model.feature_name_index:
                        feature = self.model.features[self.model.feature_name_index[g]]
                        names = []
                        values = []
                        colours = []
                        for u, vals in self.model.stratigraphic_column[g].items():
                            names.append(u + name_suffix)
                            values.append(vals["min"])
                            colours.append(cmap_colours[ci, :])
                            ci += 1
                        pbar.set_description("Isosurfacing {}".format(feature.name))
                        return_container.update(
                            self.add_isosurface(
                                feature,
                                slices=values,
                                names=names,
                                colours=colours,
                                **kwargs,
                            )
                        )
                        pbar.update(len(values))

            if faults:
                for f in self.model.features:
                    if f.type == FeatureType.FAULT:

                        def mask(x):
                            val = f.displacementfeature.evaluate_value(x)
                            val[np.isnan(val)] = 0
                            maskv = np.zeros(val.shape).astype(bool)
                            maskv[np.abs(val) > 0.001] = 1

                            return maskv

                        if (
                            self.model.stratigraphic_column
                            and f.name in self.model.stratigraphic_column["faults"]
                        ):
                            fault_colour = self.model.stratigraphic_column["faults"][
                                f.name
                            ].get("colour", ["red"])
                        pbar.set_description("Isosurfacing {}".format(f.name))
                        if displacement_cmap is not None:
                            fault_colour = [None]
                            kwargs["cmap"] = displacement_cmap
                            kwargs["vmin"] = np.min(
                                self.model.faults_displacement_magnitude
                            )
                            kwargs["vmax"] = np.max(
                                self.model.faults_displacement_magnitude
                            )
                            kwargs["paint_with"] = LambdaGeologicalFeature(
                                lambda xyz: np.zeros(xyz.shape[0]) + f.displacement
                            )
                            #  = feature
                        region = kwargs.pop("region", None)
                        return_container.update(
                            self.add_isosurface(
                                f,
                                isovalue=0,
                                region=mask,
                                colour=fault_colour,
                                name=f.name + name_suffix,
                                **kwargs,
                            )
                        )
                        pbar.update(1)
            logger.info("Adding surfaces took {} seconds".format(time.time() - start))
        return return_container

    def add_vector_field(self, geological_feature, **kwargs):
        """

        Plot the gradient of a geological feature at given locations

        Parameters
        ----------
        geological_feature : Geological Feature to evaluate gradient
        locations : ((N,3)) array of evaluation locations
        kwargs : kwargs for lavavu vector

        Returns
        -------

        """
        # if isinstance(geological_feature,GeologicalFeature):
        #     raise ValueError("{} is not a GeologicalFeature".format(type(geological_feature)))
        logger.info("Adding vector field for %s " % (geological_feature.name))
        locations = kwargs.get("locations", None)
        name = kwargs.get("name", geological_feature.name)
        if locations is None:
            x = np.linspace(
                self.bounding_box[0, 0], self.bounding_box[1, 0], self.nsteps[0]
            )
            y = np.linspace(
                self.bounding_box[0, 1], self.bounding_box[1, 1], self.nsteps[1]
            )
            z = np.linspace(
                self.bounding_box[1, 2], self.bounding_box[0, 2], self.nsteps[2]
            )
            xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
            locations = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        vector = geological_feature.evaluate_gradient(locations)
        # normalise
        mask = ~np.any(np.isnan(vector), axis=1)
        vector[mask, :] /= np.linalg.norm(vector[mask, :], axis=1)[:, None]
        self._add_vector_marker(
            self.model.rescale(locations, inplace=False), vector, name, **kwargs
        )

        return

    def add_data(self, feature, disks=False, vectors=False, **kwargs):
        """

        Plot the data linked to the feature, can choose whether to plot all data types
        using value and grad kwargs

        Parameters
        ----------
        feature
        kwargs

        Returns
        -------

        """
        name = feature.name
        add_grad = True
        add_value = True
        add_tang = True
        add_interface = True
        if "name" in kwargs:
            name = kwargs["name"]
            del kwargs["name"]
        if "grad" in kwargs:
            add_grad = kwargs["grad"]
        if "value" in kwargs:
            add_value = kwargs["value"]
        if "tang" in kwargs:
            add_tang = kwargs["tang"]
        if "interface" in kwargs:
            add_interface = kwargs["interface"]
        grad = feature.builder.get_gradient_constraints()
        norm = feature.builder.get_norm_constraints()
        value = feature.builder.get_value_constraints()
        tang = feature.builder.get_tangent_constraints()
        interface = feature.builder.get_interface_constraints()
        symbol_type = self.default_vector_symbol
        if disks:
            symbol_type = "disk"
        if vectors:
            symbol_type = "arrow"
        if vectors and disks:
            logger.warning("Cannot use both disks and arrows, using disks")
            symbol_type = "disk"
        if grad.shape[0] > 0 and add_grad:
            self.add_vector_data(
                self.model.rescale(grad[:, :3], inplace=False),
                grad[:, 3:6],
                name + "_grad_cp",
                symbol_type=symbol_type,
                **kwargs,
            )

        if norm.shape[0] > 0 and add_grad:
            self.add_vector_data(
                self.model.rescale(norm[:, :3], inplace=False),
                norm[:, 3:6],
                name + "_norm_cp",
                symbol_type=symbol_type,
                **kwargs,
            )
        if value.shape[0] > 0 and add_value:
            kwargs["range"] = [feature.min(), feature.max()]
            self.add_value_data(
                self.model.rescale(value[:, :3], inplace=False),
                value[:, 3],
                name + "_value_cp",
                **kwargs,
            )
        if tang.shape[0] > 0 and add_tang:
            self.add_vector_data(
                self.model.rescale(tang[:, :3], inplace=False),
                tang[:, 3:6],
                name + "_tang_cp",
                **kwargs,
            )
        if interface.shape[0] > 0 and add_interface:
            self.add_points(
                self.model.rescale(interface[:, :3], inplace=False),
                name + "_interface_cp",
            )

    def add_intersection_lineation(self, feature, **kwargs):
        name = feature.name
        if "name" in kwargs:
            name = kwargs["name"]
            del kwargs["name"]
        intersection = feature.builder.fold.foldframe.calculate_intersection_lineation(
            feature.builder
        )
        gpoints = feature.builder.interpolator.get_gradient_constraints()[:, :6]
        npoints = feature.builder.interpolator.get_norm_constraints()[:, :6]
        points = []
        if gpoints.shape[0] > 0:
            points.append(gpoints)
        if npoints.shape[0] > 0:
            points.append(npoints)
        points = np.vstack(points)
        if intersection.shape[0] > 0:
            self.add_vector_data(
                self.model.rescale(points[:, :3], inplace=False),
                intersection,
                name + "_intersection",
            )

    def add_points(self, points, name, **kwargs):
        """

        Plot points location in the lavavu viewer

        Parameters
        ----------
        points : numpy array of the points locations
        name :  string name of the object for lavavu
        **kwargs : lavavu points kwargs

        Returns
        -------

        """
        self._add_points(points, name, **kwargs)

    def add_dtm(self, name="dtm", colour="brown", paint_with=None, **kwargs):
        if self.model == None:
            raise BaseException("No model cannot add dtm")
        if self.model.dtm == None:
            raise BaseException("No dtm for model add one")

        # generate a triangular mesh
        tri_mask_even = np.array([[0, 3, 2], [0, 3, 1]])
        tri_mask_odd = np.array([[2, 0, 1], [2, 1, 3]])
        c_xi, c_yi = np.meshgrid(
            np.arange(0, self.model.nsteps[0] - 1),
            np.arange(0, self.model.nsteps[1] - 1),
            indexing="ij",
        )
        c_xi = c_xi.flatten(order="F")
        c_yi = c_yi.flatten(order="F")
        # get cell corners
        xcorner = np.array([0, 1, 0, 1])
        ycorner = np.array([0, 0, 1, 1])
        xi = c_xi[:, None] + xcorner[None, :]
        yi = (
            c_yi[:, None] + ycorner[None, :]
        )  # , yi, zi = self.cell_corner_indexes(c_xi, c_yi, c_zi)
        even_mask = (c_xi + c_yi) % 2 == 0
        gi = xi + yi * self.model.nsteps[0]
        triangles = np.zeros((c_xi.shape[0], 2, 3)).astype("int64")
        triangles[even_mask, :, :] = gi[even_mask, :][:, tri_mask_even]
        triangles[~even_mask, :, :] = gi[~even_mask, :][:, tri_mask_odd]

        triangles = triangles.reshape(
            (triangles.shape[0] * triangles.shape[1], triangles.shape[2])
        )
        points = np.array(
            np.meshgrid(
                np.linspace(
                    self.model.origin[0], self.model.maximum[0], self.model.nsteps[0]
                ),
                np.linspace(
                    self.model.origin[1], self.model.maximum[1], self.model.nsteps[1]
                ),
            )
        ).T.reshape(-1, 2)
        points = np.hstack([points, np.zeros((points.shape[0], 1))])
        points[:, 2] = self.model.dtm(points[:, :2])
        self._add_surface(
            points, triangles, name, colour=colour, paint_with=paint_with, **kwargs
        )

    def add_vector_data(
        self, location, vector, name, normalise=True, symbol_type="arrow", **kwargs
    ):
        """

        Plot point data with a vector component into the lavavu viewer

        Parameters
        ----------
        position : numpy array N,3 for xyz locations
        vector : numpy array of vector N,3
        name :  string name for the object in lavavu
        kwargs to pass to lavavu

        Returns
        -------

        """
        if "colour" not in kwargs:
            kwargs["colour"] = "black"
        # normalise
        if location.shape[0] > 0:
            if normalise:
                vector /= np.linalg.norm(vector, axis=1)[:, None]
            self._add_vector_marker(
                location, vector, name, symbol_type=symbol_type, **kwargs
            )
            return

    def add_value_data(self, position, value, name, **kwargs):
        """

        Plot points data with a value component

        Parameters
        ----------
        position : numpy array N,3 for xyz locations
        value : N array of values
        name :  string name of the object for lavavu
        kwargs : kwargs to pass to lavavu

        Returns
        -------

        """
        if "pointtype" not in kwargs:
            kwargs["pointtype"] = "sphere"
        if "pointsize" not in kwargs:
            kwargs["pointsize"] = 4
        # set the colour map to diverge unless user decides otherwise
        self._add_points(position, name, value=value, **kwargs)

    def add_fold(self, fold, **kwargs):
        """
        Draw the vector components of the fold at the locations

        Parameters
        ----------
        fold - fold object
        locations - numpy array of xyz

        Returns
        -------

        """
        locations = kwargs.get("locations", None)
        if locations is None:
            x = np.linspace(
                self.bounding_box[0, 0], self.bounding_box[1, 0], self.nsteps[0]
            )
            y = np.linspace(
                self.bounding_box[0, 1], self.bounding_box[1, 1], self.nsteps[1]
            )
            z = np.linspace(
                self.bounding_box[1, 2], self.bounding_box[0, 2], self.nsteps[2]
            )
            xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
            locations = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        r2r, fold_axis, dgz = fold.get_deformed_orientation(locations)
        locations = self.model.rescale(locations, inplace=False)
        self.add_vector_data(locations, r2r, fold.name + "_direction", colour="red")
        self.add_vector_data(locations, fold_axis, fold.name + "_axis", colour="black")
        self.add_vector_data(locations, dgz, fold.name + "_norm", colour="green")

    def add_support_box(self, geological_feature, paint=False, **kwargs):
        name = kwargs.get("name", geological_feature.name + "_support")
        box = np.vstack(
            [
                geological_feature.interpolator.support.origin,
                geological_feature.interpolator.support.maximum,
            ]
        )
        points, tri = create_box(box, self.nsteps)
        paint_with = None
        if paint:
            paint_with = geological_feature
        self._add_surface(
            self.model.rescale(points, inplace=False),
            tri,
            name,
            paint_with=paint_with,
            **kwargs,
        )
