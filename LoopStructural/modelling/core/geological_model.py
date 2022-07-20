"""
Main entry point for creating a geological model
"""
import logging

import numpy as np
import pandas as pd
import LoopStructural
from LoopStructural.datasets import normal_vector_headers

try:
    from LoopStructural.interpolators import DiscreteFoldInterpolator as DFI

    dfi = True
except ImportError:
    dfi = False
from LoopStructural.interpolators import FiniteDifferenceInterpolator as FDI

try:
    from LoopStructural.interpolators import PiecewiseLinearInterpolator as PLI

    pli = True
except ImportError:
    pli = False

# if LoopStructural.experimental:
from LoopStructural.interpolators import P2Interpolator

try:
    from LoopStructural.interpolators import SurfeRBFInterpolator as Surfe

    surfe = True

except ImportError:
    surfe = False

from LoopStructural.interpolators import StructuredGrid
from LoopStructural.interpolators import TetMesh
from LoopStructural.modelling.features.fault import FaultSegment
from LoopStructural.interpolators import DiscreteInterpolator

from LoopStructural.interpolators import StructuredGrid
from LoopStructural.interpolators import TetMesh
from LoopStructural.modelling.features.builders import (
    FaultBuilder,
    GeologicalFeatureBuilder,
    StructuralFrameBuilder,
    FoldedFeatureBuilder,
)
from LoopStructural.modelling.features import (
    UnconformityFeature,
    StructuralFrame,
    GeologicalFeature,
    FeatureType,
)
from LoopStructural.modelling.features.fold import (
    FoldRotationAngle,
    FoldEvent,
    FoldFrame,
)

from LoopStructural.utils.exceptions import LoopException, InterpolatorError
from LoopStructural.utils.helper import (
    all_heading,
    gradient_vec_names,
    strike_dip_vector,
    get_vectors,
)

intrusions = True
try:
    from LoopStructural.modelling.intrusions import IntrusionBuilder

    from LoopStructural.modelling.intrusions import IntrusionFrameBuilder
except ImportError as e:
    print(e)
    intrusions = False
from LoopStructural.utils import getLogger, log_to_file

logger = getLogger(__name__)


class GeologicalModel:
    """
    A geological model is the recipe for building a 3D model and  can include
    the rescaling of the model between 0 and 1.

    Attributes
    ----------
    features : list
        Contains all features youngest to oldest
    feature_name_index : dict
        maps feature name to the list index of the features
    data : pandas dataframe
        the dataframe used for building the geological model
    nsteps : tuple/np.array(3,dtype=int)
        the number of steps x,y,z to evaluate the model
    origin : tuple/np.array(3,dtype=doubles)
        the origin of the model box
    parameters : dict
        a dictionary tracking the parameters used to build the model
    scale_factor : double
        the scale factor used to rescale the model


    """

    def __init__(
        self,
        origin: np.ndarray,
        maximum: np.ndarray,
        data=None,
        rescale=False,
        nsteps=(50, 50, 25),
        reuse_supports=False,
        logfile=None,
        loglevel="info",
    ):
        """
        Parameters
        ----------
        origin : numpy array
            specifying the origin of the model
        maximum : numpy array
            specifying the maximum extent of the model
        rescale : bool
            whether to rescale the model to between 0/1

        Examples
        --------
        Demo data

        >>> from LoopStructural.datasets import load_claudius
        >>> from LoopStructural import GeologicalModel

        >>> data, bb = load_claudius()

        >>> model = GeologicalModel(bb[:,0],bb[:,1]
        >>> model.set_model_data(data)
        >>> model.create_and_add_foliation('strati')

        >>> y = np.linspace(model.bounding_box[0, 1], model.bounding_box[1, 1],
                        nsteps[1])
        >>> z = np.linspace(model.bounding_box[1, 2], model.bounding_box[0, 2],
                        nsteps[2])
        >>> xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        >>> xyz = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        >>> model.evaluate_feature_value('strati',xyz,scale=False)


        """
        # print('tet')
        if logfile:
            self.logfile = logfile
            log_to_file(logfile, loglevel)

        logger.info("Initialising geological model")
        self.features = []
        self.feature_name_index = {}
        self._data = None
        self.data = data
        self.nsteps = nsteps

        # we want to rescale the model area so that the maximum length is
        # 1
        self.origin = np.array(origin).astype(float)
        originstr = f"Model origin: {self.origin[0]} {self.origin[1]} {self.origin[2]}"
        logger.info(originstr)
        self.maximum = np.array(maximum).astype(float)
        maximumstr = "Model maximum: {} {} {}".format(
            self.maximum[0], self.maximum[1], self.maximum[2]
        )
        logger.info(maximumstr)

        lengths = self.maximum - self.origin
        self.scale_factor = 1.0
        self.bounding_box = np.zeros((2, 3))
        self.bounding_box[1, :] = self.maximum - self.origin
        self.bounding_box[1, :] = self.maximum - self.origin
        if rescale:
            self.scale_factor = float(np.max(lengths))
            logger.info(
                "Rescaling model using scale factor {}".format(self.scale_factor)
            )

        self.bounding_box /= self.scale_factor
        self.support = {}
        self.reuse_supports = reuse_supports
        if self.reuse_supports:
            logger.warning(
                "Supports are shared between geological features \n"
                "this may cause unexpected behaviour and should only\n"
                "be use by advanced users"
            )
        logger.info("Reusing interpolation supports: {}".format(self.reuse_supports))
        self.stratigraphic_column = None

        self.tol = 1e-10 * np.max(self.bounding_box[1, :] - self.bounding_box[0, :])
        self._dtm = None

    def __str__(self):
        lengths = self.maximum - self.origin
        _str = "GeologicalModel - {} x {} x {}\n".format(*lengths)
        _str += "------------------------------------------ \n"
        _str += "The model contains {} GeologicalFeatures \n".format(len(self.features))
        _str += ""
        _str += "------------------------------------------ \n"
        _str += ""
        _str += "Model origin: {} {} {}\n".format(
            self.origin[0], self.origin[1], self.origin[2]
        )
        _str += "Model maximum: {} {} {}\n".format(
            self.maximum[0], self.maximum[1], self.maximum[2]
        )
        _str += "Model rescale factor: {} \n".format(self.scale_factor)
        _str += "------------------------------------------ \n"
        _str += "Feature list: \n"
        for feature in self.features:
            _str += "  {} \n".format(feature)
        return _str

    def _ipython_key_completions_(self):
        return self.feature_name_index.keys()

    @classmethod
    def from_map2loop_directory(
        cls,
        m2l_directory,
        foliation_params={},
        fault_params={},
        use_thickness=True,
        vector_scale=1,
        gradient=False,
        **kwargs,
    ):
        """Alternate constructor for a geological model using m2l output

        Uses the information saved in the map2loop files to build a geological model.
        You can specify kwargs for building foliation using foliation_params and for
        faults using fault_params.  faults is a flag that allows for the faults to be skipped.

        Parameters
        ----------
        m2l_directory : string
            path to map2loop directory

        Returns
        -------
        (GeologicalModel, dict)
            the created geological model and a dictionary of the map2loop data
        """
        from LoopStructural.modelling.input.map2loop_processor import Map2LoopProcessor

        log_to_file(f"{m2l_directory}/loopstructural_log.txt")
        logger.info("Creating model from m2l directory")
        processor = Map2LoopProcessor(m2l_directory, use_thickness)
        processor._gradient = gradient
        processor.vector_scale = vector_scale
        for foliation_name in processor.stratigraphic_column.keys():
            if foliation_name != "faults":
                if foliation_name in foliation_params.keys():
                    processor.foliation_properties[foliation_name] = foliation_params[
                        foliation_name
                    ]
                else:
                    processor.foliation_properties[foliation_name] = foliation_params

        for fault_name in processor.fault_names:
            if fault_name in fault_params.keys():
                for param_name, value in fault_params[fault_name].items():
                    processor.fault_properties.loc[fault_name, param_name] = value
            else:
                for param_name, value in fault_params.items():
                    processor.fault_properties.loc[fault_name, param_name] = value

        model = GeologicalModel.from_processor(processor)
        return model, processor

    @classmethod
    def from_processor(cls, processor):
        logger.info("Creating model from processor")
        model = GeologicalModel(processor.origin, processor.maximum)
        model.data = processor.data
        if processor.fault_properties is not None:
            for i in processor.fault_network.faults:
                logger.info(f"Adding fault {i}")
                model.create_and_add_fault(
                    i,
                    **processor.fault_properties.to_dict("index")[i],
                    faultfunction="BaseFault",
                )
            for (
                edge,
                properties,
            ) in processor.fault_network.fault_edge_properties.items():
                if model[edge[1]] is None or model[edge[0]] is None:
                    logger.warning(
                        f"Cannot add splay {edge[1]} or {edge[0]} are not in the model"
                    )
                    continue
                splay = False
                if "angle" in properties:
                    if float(properties["angle"]) < 30 and (
                        "dip_dir"
                        not in processor.stratigraphic_column["faults"][edge[0]]
                        or np.abs(
                            processor.stratigraphic_column["faults"][edge[0]]["dip_dir"]
                            - processor.stratigraphic_column["faults"][edge[1]][
                                "dip_dir"
                            ]
                        )
                        < 90
                    ):
                        # splay
                        region = model[edge[1]].builder.add_splay(model[edge[0]])

                        model[edge[1]].splay[model[edge[0]].name] = region
                        splay = True
                if splay == False:
                    positive = None
                    if (
                        "downthrow_dir"
                        in processor.stratigraphic_column["faults"][edge[0]]
                    ):
                        positive = (
                            np.abs(
                                processor.stratigraphic_column["faults"][edge[0]][
                                    "downthrow_dir"
                                ]
                                - processor.stratigraphic_column["faults"][edge[1]][
                                    "downthrow_dir"
                                ]
                            )
                            < 90
                        )
                    model[edge[1]].add_abutting_fault(
                        model[edge[0]],
                        positive=positive,
                    )
        for s in processor.stratigraphic_column.keys():
            if s != "faults":
                faults = None
                if processor.fault_stratigraphy is not None:
                    faults = processor.fault_stratigraphy[s]
                logger.info(f"Adding foliation {s}")
                f = model.create_and_add_foliation(
                    s, **processor.foliation_properties[s], faults=faults
                )
                # check feature was built, and is an interpolated feature.
                if f is not None and f.type == FeatureType.INTERPOLATED:
                    model.add_unconformity(f, 0)
        model.stratigraphic_column = processor.stratigraphic_column
        return model

    @classmethod
    def from_file(cls, file):
        """Load a geological model from file

        Parameters
        ----------
        file : string
            path to the file

        Returns
        -------
        GeologicalModel
            the geological model object
        """
        try:
            import dill as pickle
        except ImportError:
            logger.error("Cannot import from file, dill not installed")
            return None
        model = pickle.load(open(file, "rb"))
        if type(model) == GeologicalModel:
            logger.info("GeologicalModel initialised from file")
            return model
        else:
            logger.error(f"{file} does not contain a geological model")
            return None

    def __getitem__(self, feature_name):
        """Accessor for feature in features using feature_name_index

        Parameters
        ----------
        feature_name : string
            name of the feature to return
        """
        return self.get_feature_by_name(feature_name)

    def __contains__(self, feature_name):
        return feature_name in self.feature_name_index

    @property
    def dtm(self):
        return self._dtm

    @dtm.setter
    def dtm(self, dtm):
        if not callable(dtm):
            raise BaseException(
                "DTM must be a callable function \n"
                "use LoopStructural.utils.dtm_creator to build one"
            )
        else:
            self._dtm = dtm

    @property
    def faults(self):
        faults = []
        for f in self.features:
            if type(f) == FaultSegment:
                faults.append(f)

        return faults

    @property
    def series(self):
        series = []
        for f in self.features:
            if f.type == FeatureType.INTERPOLATED:
                series.append(f)
        return series

    @property
    def faults_displacement_magnitude(self):
        displacements = []
        for f in self.faults:
            displacements.append(f.displacement)
        return np.array(displacements)

    def feature_names(self):
        return self.feature_name_index.keys()

    def fault_names(self):
        pass

    def check_inialisation(self):
        if self.data is None:
            logger.error("Data not associated with GeologicalModel. Run set_data")
            return False

    def to_file(self, file):
        """Save a model to a pickle file requires dill

        Parameters
        ----------
        file : string
            path to file location
        """
        try:
            import dill as pickle
        except ImportError:
            logger.error(
                "Cannot write to file, dill not installed \n" "pip install dill"
            )
            return
        try:
            logger.info(f"Writing GeologicalModel to: {file}")
            pickle.dump(self, open(file, "wb"))
        except pickle.PicklingError:
            logger.error("Error saving file")

    def _add_feature(self, feature):
        """
        Add a feature to the model stack

        Parameters
        ----------
        feature : GeologicalFeature
            the geological feature to add

        """

        if feature.name in self.feature_name_index:
            logger.info(
                f"Feature {feature.name} already exists at {self.feature_name_index[feature.name]}, overwriting"
            )
            self.features[self.feature_name_index[feature.name]] = feature
        else:
            self.features.append(feature)
            self.feature_name_index[feature.name] = len(self.features) - 1
            logger.info(
                f"Adding {feature.name} to model at location {len(self.features)}"
            )
        self._add_domain_fault_above(feature)
        self._add_unconformity_above(feature)
        feature.model = self

    def data_for_feature(self, feature_name):
        return self.data.loc[self.data["feature_name"] == feature_name, :]

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        """
        Set the data array for the model

        Parameters
        ----------
        data : pandas data frame
            with column headers corresponding to the
         type, X, Y, Z, nx, ny, nz, val, strike, dip, dip_dir, plunge,
         plunge_dir, azimuth

        Returns
        -------
        Note
        ----
        Type can be any unique identifier for the feature the data point
        'eg' 'S0', 'S2', 'F1_axis'
        it is then used by the create functions to get the correct data
        """
        if data is None:
            return
        if type(data) != pd.DataFrame:
            logger.warning(
                "Data is not a pandas data frame, trying to read data frame " "from csv"
            )
            try:
                data = pd.read_csv(data)
            except:
                logger.error("Could not load pandas data frame from data")
                raise BaseException("Cannot load data")
        logger.info(f"Adding data to GeologicalModel with {len(data)} data points")
        self._data = data.copy()
        self._data["X"] -= self.origin[0]
        self._data["Y"] -= self.origin[1]
        self._data["Z"] -= self.origin[2]
        self._data["X"] /= self.scale_factor
        self._data["Y"] /= self.scale_factor
        self._data["Z"] /= self.scale_factor
        if "type" in self._data:
            logger.warning("'type' is deprecated replace with 'feature_name' \n")
            self._data.rename(columns={"type": "feature_name"}, inplace=True)
        if "feature_name" not in self._data:
            logger.error("Data does not contain 'feature_name' column")
            raise BaseException("Cannot load data")
        for h in all_heading():
            if h not in self._data:
                self._data[h] = np.nan
                if h == "w":
                    self._data[h] = 1.0
                if h == "coord":
                    self._data[h] = 0
                if h == "polarity":
                    self._data[h] = 1.0
        # LS wants polarity as -1 or 1, change 0 to -1
        self._data.loc[self._data["polarity"] == 0, "polarity"] = -1.0
        self.data.loc[np.isnan(self.data["w"]), "w"] = 1.0
        if "strike" in self._data and "dip" in self._data:
            logger.info("Converting strike and dip to vectors")
            mask = np.all(~np.isnan(self._data.loc[:, ["strike", "dip"]]), axis=1)
            self._data.loc[mask, gradient_vec_names()] = (
                strike_dip_vector(
                    self._data.loc[mask, "strike"], self._data.loc[mask, "dip"]
                )
                * self._data.loc[mask, "polarity"].to_numpy()[:, None]
            )
            self._data.drop(["strike", "dip"], axis=1, inplace=True)

    def set_model_data(self, data):
        logger.warning(
            "deprecated method. Model data can now be set using the data attribute"
        )
        self.data = data

    def extend_model_data(self, newdata):
        """
        Extends the data frame

        Parameters
        ----------
        newdata : pandas data frame
            data to add to the existing dataframe
        Returns
        -------
        """
        logger.warning("Extend data is untested and may have unexpected consequences")
        data_temp = newdata.copy()
        data_temp["X"] -= self.origin[0]
        data_temp["Y"] -= self.origin[1]
        data_temp["Z"] -= self.origin[2]
        data_temp["X"] /= self.scale_factor
        data_temp["Y"] /= self.scale_factor
        data_temp["Z"] /= self.scale_factor
        self.data.concat([self.data, data_temp], sort=True)

    def set_stratigraphic_column(self, stratigraphic_column, cmap="tab20"):
        """
        Adds a stratigraphic column to the model

        Parameters
        ----------
        stratigraphic_column : dictionary
        cmap : matplotlib.cmap
        Returns
        -------

        Notes
        -----
        stratigraphic_column is a nested dictionary with the format
        {'group':
                {'series1':
                            {'min':0., 'max':10.,'id':0,'colour':}
                }
        }

        """
        # if the colour for a unit hasn't been specified we can just sample from
        # a colour map e.g. tab20
        logger.info("Adding stratigraphic column to model")
        random_colour = True
        n_units = 0
        for g in stratigraphic_column.keys():
            for u in stratigraphic_column[g].keys():
                if "colour" in stratigraphic_column[g][u]:
                    random_colour = False
                    break
                n_units += 1
        if random_colour:
            import matplotlib.cm as cm

            cmap = cm.get_cmap(cmap, n_units)
            cmap_colours = cmap.colors
            ci = 0
            for g in stratigraphic_column.keys():
                for u in stratigraphic_column[g].keys():
                    stratigraphic_column[g][u]["colour"] = cmap_colours[ci, :]

        self.stratigraphic_column = stratigraphic_column

    def get_interpolator(
        self,
        interpolatortype="FDI",
        nelements=1e4,
        buffer=0.2,
        element_volume=None,
        **kwargs,
    ):
        """
        Returns an interpolator given the arguments, also constructs a
        support for a discrete interpolator

        Parameters
        ----------
        interpolatortype : string
            define the interpolator type
        nelements : int
            number of elements in the interpolator
        buffer : double or numpy array 3x1
            value(s) between 0,1 specifying the buffer around the bounding box
        data_bb : bool
            whether to use the model boundary or the boundary around
        kwargs : no kwargs used, this just catches any additional arguments

        Returns
        -------
        """
        # get an interpolator for
        interpolator = None
        bb = np.copy(self.bounding_box)
        # add a buffer to the interpolation domain, this is necessary for
        # faults but also generally a good
        # idea to avoid boundary problems
        # buffer = bb[1, :]
        buffer = (np.min(bb[1, :] - bb[0, :])) * buffer
        bb[0, :] -= buffer  # *(bb[1,:]-bb[0,:])
        bb[1, :] += buffer  # *(bb[1,:]-bb[0,:])
        box_vol = (bb[1, 0] - bb[0, 0]) * (bb[1, 1] - bb[0, 1]) * (bb[1, 2] - bb[0, 2])
        if interpolatortype == "PLI" and pli:
            if element_volume is None:
                # nelements /= 5
                element_volume = box_vol / nelements
            # calculate the step vector of a regular cube
            step_vector = np.zeros(3)
            step_vector[:] = element_volume ** (1.0 / 3.0)
            # step_vector /= np.array([1,1,2])
            # number of steps is the length of the box / step vector
            nsteps = np.ceil((bb[1, :] - bb[0, :]) / step_vector).astype(int)
            if np.any(np.less(nsteps, 3)):
                axis_labels = ["x", "y", "z"]
                for i in range(3):
                    if nsteps[i] < 3:
                        logger.error(
                            f"Number of steps in direction {axis_labels[i]} is too small, try increasing nelements"
                        )
                logger.error("Cannot create interpolator: number of steps is too small")
                raise ValueError("Number of steps too small cannot create interpolator")
            # create a structured grid using the origin and number of steps
            if self.reuse_supports:
                mesh_id = f"mesh_{nelements}"
                mesh = self.support.get(
                    mesh_id,
                    TetMesh(origin=bb[0, :], nsteps=nsteps, step_vector=step_vector),
                )
                if mesh_id not in self.support:
                    self.support[mesh_id] = mesh
            else:
                if "meshbuilder" in kwargs:
                    mesh = kwargs["meshbuilder"](bb, nelements)
                else:
                    mesh = TetMesh(
                        origin=bb[0, :], nsteps=nsteps, step_vector=step_vector
                    )
            logger.info(
                "Creating regular tetrahedron mesh with %i elements \n"
                "for modelling using PLI" % (mesh.ntetra)
            )

            return PLI(mesh)
        if interpolatortype == "P2":
            if element_volume is None:
                # nelements /= 5
                element_volume = box_vol / nelements
            # calculate the step vector of a regular cube
            step_vector = np.zeros(3)
            step_vector[:] = element_volume ** (1.0 / 3.0)
            # step_vector /= np.array([1,1,2])
            # number of steps is the length of the box / step vector
            nsteps = np.ceil((bb[1, :] - bb[0, :]) / step_vector).astype(int)
            if "meshbuilder" in kwargs:
                mesh = kwargs["meshbuilder"](bb, nelements)
            else:
                raise NotImplementedError(
                    "Cannot use P2 interpolator without external mesh"
                )
            logger.info(
                "Creating regular tetrahedron mesh with %i elements \n"
                "for modelling using P2" % (mesh.ntetra)
            )
            return P2Interpolator(mesh)
        if interpolatortype == "FDI":

            # find the volume of one element
            if element_volume is None:
                element_volume = box_vol / nelements
            # calculate the step vector of a regular cube
            step_vector = np.zeros(3)
            step_vector[:] = element_volume ** (1.0 / 3.0)
            # number of steps is the length of the box / step vector
            nsteps = np.ceil((bb[1, :] - bb[0, :]) / step_vector).astype(int)
            if np.any(np.less(nsteps, 3)):
                logger.error("Cannot create interpolator: number of steps is too small")
                axis_labels = ["x", "y", "z"]
                for i in range(3):
                    if nsteps[i] < 3:
                        logger.error(
                            f"Number of steps in direction {axis_labels[i]} is too small, try increasing nelements"
                        )
                raise ValueError("Number of steps too small cannot create interpolator")
            # create a structured grid using the origin and number of steps
            if self.reuse_supports:
                grid_id = "grid_{}".format(nelements)
                grid = self.support.get(
                    grid_id,
                    StructuredGrid(
                        origin=bb[0, :], nsteps=nsteps, step_vector=step_vector
                    ),
                )
                if grid_id not in self.support:
                    self.support[grid_id] = grid
            else:
                grid = StructuredGrid(
                    origin=bb[0, :], nsteps=nsteps, step_vector=step_vector
                )
            logger.info(
                f"Creating regular grid with {grid.n_elements} elements \n"
                "for modelling using FDI"
            )
            return FDI(grid)

        if interpolatortype == "DFI" and dfi == True:  # "fold" in kwargs:
            if element_volume is None:
                nelements /= 5
                element_volume = box_vol / nelements
            # calculate the step vector of a regular cube
            step_vector = np.zeros(3)
            step_vector[:] = element_volume ** (1.0 / 3.0)
            # number of steps is the length of the box / step vector
            nsteps = np.ceil((bb[1, :] - bb[0, :]) / step_vector).astype(int)
            # create a structured grid using the origin and number of steps
            if "meshbuilder" in kwargs:
                mesh = kwargs["meshbuilder"].build(bb, nelements)
            else:
                mesh = kwargs.get(
                    "mesh",
                    TetMesh(origin=bb[0, :], nsteps=nsteps, step_vector=step_vector),
                )
            logger.info(
                f"Creating regular tetrahedron mesh with {mesh.ntetra} elements \n"
                "for modelling using DFI"
            )
            return DFI(mesh, kwargs["fold"])
        if interpolatortype == "Surfe" or interpolatortype == "surfe":
            # move import of surfe to where we actually try and use it
            if not surfe:
                logger.warning("Cannot import Surfe, try another interpolator")
                raise ImportError("Cannot import surfepy, try pip install surfe")
            method = kwargs.get("method", "single_surface")
            logger.info("Using surfe interpolator")
            return Surfe(method)
        logger.warning("No interpolator")
        raise InterpolatorError("Could not create interpolator")

    def create_and_add_foliation(
        self, series_surface_data, tol=None, faults=None, **kwargs
    ):
        """
        Parameters
        ----------
        series_surface_data : string
            corresponding to the feature_name in the data
        kwargs

        Returns
        -------
        feature : GeologicalFeature
            the created geological feature
        """
        if self.check_inialisation() == False:
            return
        # if tol is not specified use the model default
        if tol is None:
            tol = self.tol

        interpolator = self.get_interpolator(**kwargs)
        series_builder = GeologicalFeatureBuilder(
            interpolator, name=series_surface_data, **kwargs
        )
        # add data
        series_data = self.data[self.data["feature_name"] == series_surface_data]
        if series_data.shape[0] == 0:
            logger.warning("No data for {series_surface_data}, skipping")
            return
        series_builder.add_data_from_data_frame(series_data)
        self._add_faults(series_builder, features=faults)

        # build feature
        # series_feature = series_builder.build(**kwargs)
        series_feature = series_builder.feature
        series_builder.build_arguments = kwargs
        series_builder.build_arguments["tol"] = tol
        series_feature.type = FeatureType.INTERPOLATED
        self._add_feature(series_feature)
        return series_feature

    def create_and_add_dtm(self, series_surface_data, **kwargs):
        """
        Parameters
        ----------
        series_surface_data : string
            corresponding to the feature_name in the data
        kwargs

        Returns
        -------
        feature : GeologicalFeature
            the created geological feature
        """
        if self.check_inialisation() == False:
            return False

        interpolator = self.get_interpolator(**kwargs)
        series_builder = GeologicalFeatureBuilder(
            interpolator, name=series_surface_data, **kwargs
        )
        # add data
        series_data = self.data[self.data["feature_name"] == series_surface_data]
        if series_data.shape[0] == 0:
            logger.warning("No data for {series_surface_data}, skipping")
            return
        series_builder.add_data_from_data_frame(series_data)
        # self._add_faults(series_builder)

        # build feature
        # series_feature = series_builder.build(**kwargs)
        series_feature = series_builder.feature
        series_builder.build_arguments = kwargs
        series_feature.type = "dtm"
        # see if any unconformities are above this feature if so add region
        # self._add_unconformity_above(series_feature)self._add_feature(series_feature)
        self._add_feature(series_feature)
        return series_feature

    def create_and_add_fold_frame(self, foldframe_data, tol=None, **kwargs):
        """
        Parameters
        ----------
        foldframe_data : string
            unique string in feature_name column

        kwargs

        Returns
        -------
        fold_frame : FoldFrame
            the created fold frame
        """
        if self.check_inialisation() == False:
            return False
        if tol is None:
            tol = self.tol

        # create fault frame
        interpolator = self.get_interpolator(**kwargs)
        #
        fold_frame_builder = StructuralFrameBuilder(
            interpolator, name=foldframe_data, frame=FoldFrame, **kwargs
        )
        # add data
        fold_frame_data = self.data[self.data["feature_name"] == foldframe_data]
        fold_frame_builder.add_data_from_data_frame(fold_frame_data)
        self._add_faults(fold_frame_builder[0])
        self._add_faults(fold_frame_builder[1])
        self._add_faults(fold_frame_builder[2])
        kwargs["tol"] = tol
        fold_frame_builder.setup(**kwargs)
        fold_frame = fold_frame_builder.frame
        # for i in range(3):
        #     self._add_unconformity_above(fold_frame[i])
        fold_frame.type = FeatureType.STRUCTURALFRAME
        fold_frame.builder = fold_frame_builder
        self._add_feature(fold_frame)

        return fold_frame

    def create_and_add_folded_foliation(
        self, foliation_data, fold_frame=None, svario=True, tol=None, **kwargs
    ):
        """
        Create a folded foliation field from data and a fold frame

        Parameters
        ----------
        foliation_data : str
            unique string in type column of data frame
        fold_frame :  FoldFrame
        svario  : Boolean
            whether to calculate svariograms, saves time if avoided
        kwargs
            additional kwargs to be passed through to other functions

        Returns
        -------
        feature : GeologicalFeature
            created geological feature
        """
        if self.check_inialisation() == False:
            return False
        if tol is None:
            tol = self.tol

        if fold_frame is None:
            logger.info("Using last feature as fold frame")
            fold_frame = self.features[-1]
        assert type(fold_frame) == FoldFrame, "Please specify a FoldFrame"
        fold = FoldEvent(fold_frame, name=f"Fold_{foliation_data}")
        fold_interpolator = self.get_interpolator("DFI", fold=fold, **kwargs)
        if "fold_weights" not in kwargs:
            kwargs["fold_weights"] = {}

        series_builder = FoldedFeatureBuilder(
            interpolator=fold_interpolator, fold=fold, name=foliation_data, **kwargs
        )

        series_builder.add_data_from_data_frame(
            self.data[self.data["feature_name"] == foliation_data]
        )
        self._add_faults(series_builder)
        # series_builder.add_data_to_interpolator(True)
        self._add_faults(series_builder)
        # build feature

        kwargs["tol"] = tol

        # series_feature = series_builder.build(**kwargs)
        series_feature = series_builder.feature
        series_builder.build_arguments = kwargs
        series_feature.type = FeatureType.INTERPOLATED
        series_feature.fold = fold
        # see if any unconformities are above this feature if so add region
        # self._add_unconformity_above(series_feature)self._add_feature(series_feature)
        # result['support'] = series_feature.get_interpolator().support
        self._add_feature(series_feature)
        return series_feature

    def create_and_add_folded_fold_frame(
        self, fold_frame_data, fold_frame=None, tol=None, **kwargs
    ):
        """

        Parameters
        ----------
        fold_frame_data : string
            name of the feature to be added

        fold_frame : StructuralFrame, optional
            the fold frame for the fold if not specified uses last feature added

        kwargs

        Returns
        -------
        fold_frame : FoldFrame
            created fold frame
        """
        if self.check_inialisation() == False:
            return False
        if tol is None:
            tol = self.tol

        if fold_frame is None:
            logger.info("Using last feature as fold frame")
            fold_frame = self.features[-1]
        assert type(fold_frame) == FoldFrame, "Please specify a FoldFrame"
        fold = FoldEvent(fold_frame, name=f"Fold_{fold_frame_data}")
        fold_interpolator = self.get_interpolator("DFI", fold=fold, **kwargs)
        gy_fold_interpolator = self.get_interpolator("DFI", fold=fold, **kwargs)
        frame_interpolator = self.get_interpolator(**kwargs)
        interpolators = [
            fold_interpolator,
            gy_fold_interpolator,
            frame_interpolator.copy(),
        ]
        fold_frame_builder = StructuralFrameBuilder(
            interpolators=interpolators,
            name=fold_frame_data,
            fold=fold,
            frame=FoldFrame,
            **kwargs,
        )
        fold_frame_builder.add_data_from_data_frame(
            self.data[self.data["feature_name"] == fold_frame_data]
        )

        for i in range(3):
            self._add_faults(fold_frame_builder[i])
        # build feature
        kwargs["frame"] = FoldFrame
        kwargs["tol"] = tol
        fold_frame_builder.setup(**kwargs)
        # fold_frame_builder.build_arguments = kwargs
        folded_fold_frame = fold_frame_builder.frame
        folded_fold_frame.builder = fold_frame_builder

        folded_fold_frame.type = "structuralframe"
        # see if any unconformities are above this feature if so add region
        # for i in range(3):
        #     self._add_unconformity_above(fold_frame[i])

        self._add_feature(folded_fold_frame)

        return folded_fold_frame

    def create_and_add_intrusion(
        self,
        intrusion_name,
        intrusion_frame_name,
        intrusion_lateral_extent_model=None,
        intrusion_vertical_extent_model=None,
        intrusion_network_parameters={},
        lateral_extent_sgs_parameters={},
        vertical_extent_sgs_parameters={},
        **kwargs,
    ):
        """
        An intrusion in built in two main steps:
        (1) Intrusion builder: intrusion builder creates the intrusion structural frame.
            This object is curvilinear coordinate system of the intrusion constrained with intrusion network points,
            and flow and inflation measurements (provided by the user).
            The intrusion network is a representation of the approximated location of roof or floor contact of the intrusion.
            This object might be constrained using the anisotropies of the host rock if the roof (or floor) contact is not well constrained.

        (2) Intrusion feature: simulation of lateral and vertical extent of intrusion within the model volume.
            The simulations outcome consist in thresholds distances along the structural frame coordinates
            that are used to constrained the extent of the intrusion.

        Parameters
        ----------
        intrusion_name :  string,
            name of intrusion feature in model data
        intrusion_frame_name :  string,
            name of intrusion frame in model data
        intrusion_lateral_extent_model = function,
            geometrical conceptual model for simulation of lateral extent
        intrusion_vertical_extent_model = function,
            geometrical conceptual model for simulation of vertical extent
        intrusion_network_parameters : dictionary, optional
            contact : string, contact of the intrusion to be used to create the network (roof or floor)
            type : string, type of algorithm to create the intrusion network (interpolated or shortest path).
                    Shortest path is recommended when intrusion contact is not well constrained
            contacts_anisotropies : list of series-type features involved in intrusion emplacement
            structures_anisotropies : list of fault-type features involved in intrusion emplacement
            sequence_anisotropies : list of anisotropies to look for the shortest path. It could be only starting and end point.
        lateral_extent_sgs_parameters = dictionary, optional
            parameters for sequential gaussian simulation of lateral extent
        vertical_extent_sgs_parameters = dictionary, optional
            parameters for sequential gaussian simulation of vertical extent

        kwargs

        Returns
        -------
        intrusion feature

        """
        if intrusions == False:
            logger.error("Libraries not installed")
            raise Exception("Libraries not installed")

        intrusion_data = self.data[self.data["feature_name"] == intrusion_name].copy()
        intrusion_frame_data = self.data[
            self.data["feature_name"] == intrusion_frame_name
        ].copy()

        # -- get variables for intrusion frame interpolation
        gxxgz = kwargs.get("gxxgz", 0)
        gxxgy = kwargs.get("gxxgy", 0)
        gyxgz = kwargs.get("gyxgz", 0)

        interpolatortype = kwargs.get("interpolatortype", "FDI")
        nelements = kwargs.get("nelements", 1e2)

        weights = [gxxgz, gxxgy, gyxgz]
        interpolator = self.get_interpolator(interpolatortype=interpolatortype)

        intrusion_frame_builder = IntrusionFrameBuilder(
            interpolator, name=intrusion_frame_name, model=self, **kwargs
        )

        # -- create intrusion network
        intrusion_frame_builder.set_intrusion_network_parameters(
            intrusion_data, intrusion_network_parameters
        )
        intrusion_network_geometry = intrusion_frame_builder.create_intrusion_network()

        # -- create intrusion frame using intrusion network points and flow/inflation measurements
        intrusion_frame_builder.set_intrusion_frame_data(
            intrusion_frame_data, intrusion_network_geometry
        )

        ## -- create intrusion frame
        intrusion_frame_builder.setup(
            nelements=nelements,
            w2=weights[0],
            w1=weights[1],
            gxygz=weights[2],
        )

        intrusion_frame = intrusion_frame_builder.frame

        # -- create intrusion builder to simulate distance thresholds along frame coordinates
        intrusion_builder = IntrusionBuilder(
            intrusion_frame, model=self, name=f"{intrusion_name}_feature"
        )
        intrusion_builder.lateral_extent_model = intrusion_lateral_extent_model
        intrusion_builder.vertical_extent_model = intrusion_vertical_extent_model

        # logger.info("setting data for thresholds simulation")
        intrusion_builder.set_data_for_extent_simulation(intrusion_data)
        intrusion_builder.build_arguments = {
            "lateral_extent_sgs_parameters": lateral_extent_sgs_parameters,
            "vertical_extent_sgs_parameters": vertical_extent_sgs_parameters,
        }
        intrusion_feature = intrusion_builder.feature
        self._add_feature(intrusion_feature)

        return intrusion_feature

    def _add_faults(self, feature_builder, features=None):
        """Adds all existing faults to a geological feature builder

        Parameters
        ----------
        feature_builder : GeologicalFeatureBuilder/StructuralFrameBuilder
            The feature buider to add the faults to
        features : list, optional
            A specific list of features rather than all features in the model
        Returns
        -------

        """
        if features is None:
            features = self.features
        for f in reversed(features):
            if isinstance(f, str):
                f = self.__getitem__(f)
            if f.type == "fault":
                feature_builder.add_fault(f)
            # if f.type == 'unconformity':
            #     break

    def _add_domain_fault_above(self, feature):
        """
        Looks through the feature list and adds any domain faults to the feature. The domain fault masks everything
        where the fault scalar field is < 0 as being active when added to feature.

        Parameters
        ----------
        feature : GeologicalFeatureBuilder
            the feature being added to the model where domain faults should be added

        Returns
        -------

        """
        for f in reversed(self.features):
            if f.name == feature.name:
                continue
            if f.type == "domain_fault":
                feature.add_region(lambda pos: f.evaluate_value(pos) < 0)
                break

    def _add_domain_fault_below(self, domain_fault):
        """
        Looks through the feature list and adds any the domain_fault to the features that already exist in the stack
        until an unconformity is reached. domain faults to the feature. The domain fault masks everything
        where the fault scalar field is < 0 as being active when added to feature.

        Parameters
        ----------
        feature : GeologicalFeatureBuilder
            the feature being added to the model where domain faults should be added

        Returns
        -------

        """
        for f in reversed(self.features):
            if f.name == domain_fault.name:
                continue
            f.add_region(lambda pos: domain_fault.evaluate_value(pos) > 0)
            if f.type == FeatureType.UNCONFORMITY:
                break

    def _add_unconformity_above(self, feature):
        """

        Adds a region to the feature to prevent the value from being
        interpolated where the unconformities exists above e.g.
        if there is another feature above and the unconformity is at 0
        then the features added below (after) will only be visible where the
        uncomformity is <0

        Parameters
        ----------
        feature - GeologicalFeature

        Returns
        -------

        """
        for f in reversed(self.features):
            if f.type == FeatureType.UNCONFORMITY and f.name != feature.name:
                feature.add_region(lambda pos: f.evaluate(pos))
                break

    def create_and_add_unconformity(self, unconformity_surface_data, **kwargs):
        """
        Parameters
        ----------
        unconformity_surface_data : string
            name of the unconformity data in the data frame

        Returns
        -------
        """
        if not self.check_initialisation():
            return False
        interpolator = self.get_interpolator(**kwargs)
        unconformity_feature_builder = GeologicalFeatureBuilder(
            interpolator, name=unconformity_surface_data
        )
        # add data
        unconformity_data = self.data[
            self.data["feature_name"] == unconformity_surface_data
        ]

        unconformity_feature_builder.add_data_from_data_frame(unconformity_data)
        # look through existing features if there is a fault before an
        # unconformity
        # then add to the feature, once we get to an unconformity stop
        self._add_faults(unconformity_feature_builder)

        # build feature
        # uc_feature_base = unconformity_feature_builder.build(**kwargs)
        uc_feature_base = unconformity_feature_builder.feature
        unconformity_feature_builder.build_arguments = kwargs
        uc_feature_base.type = "unconformity_base"
        # uc_feature = UnconformityFeature(uc_feature_base,0)
        # iterate over existing features and add the unconformity as a
        # region so the feature is only
        # evaluated where the unconformity is positive
        return self.add_unconformity(uc_feature_base, 0)

    def add_unconformity(
        self, feature: GeologicalFeature, value: float
    ) -> UnconformityFeature:
        """
        Use an existing feature to add an unconformity to the model.

        Parameters
        ----------
        feature : GeologicalFeature
            existing geological feature
        value : float
            scalar value of isosurface that represents

        Returns
        -------
        unconformity : GeologicalFeature
            unconformity feature

        """
        logger.debug(f"Adding {feature.name} as unconformity at {value}")
        if feature is None:
            logger.warning(f"Cannot add unconformtiy, base feature is None")
            return
        # look backwards through features and add the unconformity as a region until
        # we get to an unconformity
        uc_feature = UnconformityFeature(feature, value)

        for f in reversed(self.features):
            if f.type == FeatureType.UNCONFORMITY:
                logger.debug(f"Reached unconformity {f.name}")
                break
            logger.debug(f"Adding {uc_feature.name} as unconformity to {f.name}")
            f.add_region(lambda pos: ~uc_feature.evaluate(pos))
        # now add the unconformity to the feature list
        self._add_feature(uc_feature)
        return uc_feature

    def add_onlap_unconformity(
        self, feature: GeologicalFeature, value: float
    ) -> GeologicalFeature:
        """
        Use an existing feature to add an unconformity to the model.

        Parameters
        ----------
        feature : GeologicalFeature
            existing geological feature
        value : float
            scalar value of isosurface that represents

        Returns
        -------
        unconformity_feature : GeologicalFeature
            the created unconformity

        """

        uc_feature = UnconformityFeature(feature, value, True)
        feature.add_region(lambda pos: ~uc_feature.evaluate(pos))
        for f in reversed(self.features):
            if f.type == FeatureType.UNCONFORMITY:
                break
            if f != feature:
                f.add_region(lambda pos: uc_feature.evaluate(pos))
        # for f in self.features:
        #     f.add_region(lambda pos: uc_feature.evaluate(pos))

        # see if any unconformities are above this feature if so add region
        # self._add_unconformity_above(uc_feature)
        # self._add_unconformity_below(uc_feature)  # , uc_feature)
        self._add_feature(uc_feature)

        return uc_feature

    def create_and_add_domain_fault(self, fault_surface_data, **kwargs):
        """
        Parameters
        ----------
        fault_surface_data : string
            name of the domain fault data in the data frame

        Returns
        -------
        domain_Fault : GeologicalFeature
            the created domain fault

        """
        interpolator = self.get_interpolator(**kwargs)
        domain_fault_feature_builder = GeologicalFeatureBuilder(
            interpolator, name=fault_surface_data
        )
        # add data
        unconformity_data = self.data[self.data["feature_name"] == fault_surface_data]

        domain_fault_feature_builder.add_data_from_data_frame(unconformity_data)
        # look through existing features if there is a fault before an
        # unconformity
        # then add to the feature, once we get to an unconformity stop
        self._add_faults(domain_fault_feature_builder)

        # build feature
        # domain_fault = domain_fault_feature_builder.build(**kwargs)
        domain_fault = domain_fault_feature_builder.feature
        domain_fault_feature_builder.build_arguments = kwargs
        domain_fault.type = "domain_fault"
        self._add_feature(domain_fault)
        self._add_domain_fault_below(domain_fault)

        domain_fault_uc = UnconformityFeature(domain_fault, 0)
        # iterate over existing features and add the unconformity as a
        # region so the feature is only
        # evaluated where the unconformity is positive
        return domain_fault_uc

    def create_and_add_fault(
        self,
        fault_surface_data,
        displacement,
        tol=None,
        fault_slip_vector=None,
        fault_center=None,
        major_axis=None,
        minor_axis=None,
        intermediate_axis=None,
        faultfunction="BaseFault",
        **kwargs,
    ):
        """
        Parameters
        ----------
        fault_surface_data : string
            name of the fault surface data in the dataframe
        displacement : displacement magnitude
        major_axis : [type], optional
            [description], by default None
        minor_axis : [type], optional
            [description], by default None
        intermediate_axis : [type], optional
            [description], by default None
        kwargs : additional kwargs for Fault and interpolators

        Returns
        -------
        fault : FaultSegment
            created fault
        """
        if "fault_extent" in kwargs and major_axis == None:
            major_axis = kwargs["fault_extent"]
        if "fault_influence" in kwargs and minor_axis == None:
            minor_axis = kwargs["fault_influence"]
        if "fault_vectical_radius" in kwargs and intermediate_axis == None:
            intermediate_axis = kwargs["fault_vectical_radius"]

        logger.info(f'Creating fault "{fault_surface_data}"')
        logger.info(f"Displacement: {displacement}")
        logger.info(f"Tolerance: {tol}")
        logger.info(f"Fault function: {faultfunction}")
        logger.info(f"Fault slip vector: {fault_slip_vector}")
        logger.info(f"Fault center: {fault_center}")
        logger.info(f"Major axis: {major_axis}")
        logger.info(f"Minor axis: {minor_axis}")
        logger.info(f"Intermediate axis: {intermediate_axis}")
        fault_slip_vector = np.array(fault_slip_vector, dtype="float")
        fault_center = np.array(fault_center, dtype="float")

        for k, v in kwargs.items():
            logger.info(f"{k}: {v}")

        if tol is None:
            tol = self.tol

        if displacement == 0:
            logger.warning(f"{fault_surface_data} displacement is 0")

        if "data_region" in kwargs:
            kwargs.pop("data_region")
            logger.error("kwarg data_region currently not supported, disabling")
        displacement_scaled = displacement / self.scale_factor
        # create fault frame
        interpolator = self.get_interpolator(**kwargs)
        # faults arent supported for surfe
        if isinstance(interpolator, DiscreteInterpolator) == False:
            logger.error(
                "Change interpolator to a discrete interpolation algorithm FDI/PLI"
            )
            interpolatortype = kwargs["interpolatortype"]
            raise InterpolatorError(f"Faults not supported for {interpolatortype}")
        fault_frame_builder = FaultBuilder(
            interpolator, name=fault_surface_data, model=self, **kwargs
        )
        # add data
        fault_frame_data = self.data[
            self.data["feature_name"] == fault_surface_data
        ].copy()
        trace_mask = np.logical_and(
            fault_frame_data["coord"] == 0, fault_frame_data["val"] == 0
        )
        logger.info(f"There are {np.sum(trace_mask)} points on the fault trace")
        if np.sum(trace_mask) == 0:
            logger.error(
                "You cannot model a fault without defining the location of the fault"
            )
            raise ValueError(f"There are no points on the fault trace")

        mask = np.logical_and(
            fault_frame_data["coord"] == 0, ~np.isnan(fault_frame_data["gz"])
        )
        vector_data = fault_frame_data.loc[mask, ["gx", "gy", "gz"]].to_numpy()
        mask2 = np.logical_and(
            fault_frame_data["coord"] == 0, ~np.isnan(fault_frame_data["nz"])
        )
        vector_data = np.vstack(
            [vector_data, fault_frame_data.loc[mask2, ["nx", "ny", "nz"]].to_numpy()]
        )
        fault_normal_vector = np.mean(vector_data, axis=0)
        logger.info(f"Fault normal vector: {fault_normal_vector}")

        mask = np.logical_and(
            fault_frame_data["coord"] == 1, ~np.isnan(fault_frame_data["gz"])
        )
        if fault_slip_vector is None:
            if (
                "avgSlipDirEasting" in kwargs
                and "avgSlipDirNorthing" in kwargs
                and "avgSlipDirAltitude" in kwargs
            ):
                fault_slip_vector = np.array(
                    [
                        kwargs["avgSlipDirEasting"],
                        kwargs["avgSlipDirNorthing"],
                        kwargs["avgSlipDirAltitude"],
                    ],
                    dtype=float,
                )
            else:
                fault_slip_vector = (
                    fault_frame_data.loc[mask, ["gx", "gy", "gz"]]
                    .mean(axis=0)
                    .to_numpy()
                )
        if np.any(np.isnan(fault_slip_vector)):
            logger.info("Fault slip vector is nan, estimating from fault normal")
            strike_vector, dip_vector = get_vectors(fault_normal_vector[None, :])
            fault_slip_vector = dip_vector[:, 0]
            logger.info(f"Estimated fault slip vector: {fault_slip_vector}")

        if fault_center is not None and ~np.isnan(fault_center).any():
            fault_center = self.scale(fault_center, inplace=False)
        else:
            # if we haven't defined a fault centre take the center of mass for lines assocaited with
            # the fault trace
            if (
                ~np.isnan(kwargs.get("centreEasting", np.nan))
                and ~np.isnan(kwargs.get("centreNorthing", np.nan))
                and ~np.isnan(kwargs.get("centreAltitude", np.nan))
            ):
                fault_center = self.scale(
                    np.array(
                        [
                            kwargs["centreEasting"],
                            kwargs["centreNorthing"],
                            kwargs["centreAltitude"],
                        ],
                        dtype=float,
                    ),
                    inplace=False,
                )
            else:
                mask = np.logical_and(
                    fault_frame_data["coord"] == 0, fault_frame_data["val"] == 0
                )
                fault_center = (
                    fault_frame_data.loc[mask, ["X", "Y", "Z"]].mean(axis=0).to_numpy()
                )
        if minor_axis:
            minor_axis = minor_axis / self.scale_factor
        if major_axis:
            major_axis = major_axis / self.scale_factor
        if intermediate_axis:
            intermediate_axis = intermediate_axis / self.scale_factor
        fault_frame_builder.create_data_from_geometry(
            fault_frame_data,
            fault_center,
            fault_normal_vector,
            fault_slip_vector,
            minor_axis=minor_axis,
            major_axis=major_axis,
            intermediate_axis=intermediate_axis,
            points=kwargs.get("points", False),
        )
        if "force_mesh_geometry" not in kwargs:

            fault_frame_builder.set_mesh_geometry(kwargs.get("fault_buffer", 0.2), 0)
        if "splay" in kwargs and "splayregion" in kwargs:
            fault_frame_builder.add_splay(kwargs["splay"], kwargs["splayregion"])

        kwargs["tol"] = tol
        fault_frame_builder.setup(**kwargs)
        fault = fault_frame_builder.frame
        fault.displacement = displacement_scaled
        fault.faultfunction = faultfunction

        for f in reversed(self.features):
            if f.type == FeatureType.UNCONFORMITY:
                fault.add_region(lambda pos: f.evaluate_value(pos) <= 0)
                break
        if displacement == 0:
            fault.type = "fault_inactive"
        self._add_feature(fault)

        return fault

    def rescale(self, points, inplace=True):
        """
        Convert from model scale to real world scale - in the future this
        should also do transformations?

        Parameters
        ----------
        points : np.array((N,3),dtype=double)
        inplace : boolean
            whether to return a modified copy or modify the original array

        Returns
        -------
        points : np.array((N,3),dtype=double)

        """
        if inplace == False:
            points = points.copy()
        points *= self.scale_factor
        points += self.origin
        return points

    def scale(self, points, inplace=True):
        """Take points in UTM coordinates and reproject
        into scaled model space

        Parameters
        ----------
        points : np.array((N,3),dtype=float)
            points to
        inplace : bool, optional default = True
            whether to copy the points array or update the passed array
        Returns
        -------
        points : np.array((N,3),dtype=double)

        """
        points = np.array(points).astype(float)
        if inplace == False:
            points = points.copy()
        # if len(points.shape) == 1:
        #     points = points[None,:]
        # if len(points.shape) != 2:
        #     logger.error("cannot scale array of dimensions".format(len(points.shape)))
        points -= self.origin
        points /= self.scale_factor
        return points

    def regular_grid(self, nsteps=None, shuffle=True, rescale=False, order="C"):
        """
        Return a regular grid within the model bounding box

        Parameters
        ----------
        nsteps : tuple
            number of cells in x,y,z

        Returns
        -------
        xyz : np.array((N,3),dtype=float)
            locations of points in regular grid
        """
        if nsteps is None:
            nsteps = self.nsteps
        x = np.linspace(self.bounding_box[0, 0], self.bounding_box[1, 0], nsteps[0])
        y = np.linspace(self.bounding_box[0, 1], self.bounding_box[1, 1], nsteps[1])
        z = np.linspace(self.bounding_box[1, 2], self.bounding_box[0, 2], nsteps[2])
        xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
        locs = np.array(
            [xx.flatten(order=order), yy.flatten(order=order), zz.flatten(order=order)]
        ).T
        if shuffle:
            logger.info("Shuffling points")
            np.random.shuffle(locs)
        if rescale:
            locs = self.rescale(locs)
        return locs

    def evaluate_model(self, xyz, scale=True):
        """Evaluate the stratigraphic id at each location


        Parameters
        ----------
        xyz : np.array((N,3),dtype=float)
            locations
        scale : bool
            whether to rescale the xyz before evaluating model

        Returns
        -------
        stratigraphic_id : np.array(N,dtype=int)
            the stratigraphic index for locations

        Examples
        --------
        Evaluate on a voxet

        >>> x = np.linspace(model.bounding_box[0, 0], model.bounding_box[1, 0],
                        nsteps[0])
        >>> y = np.linspace(model.bounding_box[0, 1], model.bounding_box[1, 1],
                        nsteps[1])
        >>> z = np.linspace(model.bounding_box[1, 2], model.bounding_box[0, 2],
                        nsteps[2])
        >>> xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        >>> xyz = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        >>> model.evaluate_model(xyz,scale=False)

        Evaluate on points defined by regular grid function

        >>> model.evaluate_model(model.regular_grid(shuffle=False),scale=False)


        Evaluate on a map

        >>> x = np.linspace(self.bounding_box[0, 0], self.bounding_box[1, 0],
                        nsteps[0])
        >>> y = np.linspace(self.bounding_box[0, 1], self.bounding_box[1, 1],
                        nsteps[1])
        >>> xx, yy = np.meshgrid(x, y, indexing='ij')
        >>> zz = np.zeros_like(yy)
        >>> xyz = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        >>> model.evaluate_model(model.regular_grid(shuffle=False),scale=False)

        Evaluate on points in reference coordinate system
        >>> model.evaluate_model(xyz,scale=True)

        """
        xyz = np.array(xyz)
        if scale:
            xyz = self.scale(xyz, inplace=False)
        strat_id = np.zeros(xyz.shape[0], dtype=int)
        for group in self.stratigraphic_column.keys():
            if group == "faults":
                continue
            feature_id = self.feature_name_index.get(group, -1)
            if feature_id >= 0:
                feature = self.features[feature_id]
                vals = feature.evaluate_value(xyz)
                for series in self.stratigraphic_column[group].values():
                    strat_id[
                        np.logical_and(
                            vals < series.get("max", feature.max()),
                            vals > series.get("min", feature.min()),
                        )
                    ] = series["id"]
            if feature_id == -1:
                logger.error(f"Model does not contain {group}")
        return strat_id

    def evaluate_fault_displacements(self, points, scale=True):
        """Evaluate the fault displacement magnitude at each location


        Parameters
        ----------
        xyz : np.array((N,3),dtype=float)
            locations
        scale : bool
            whether to rescale the xyz before evaluating model

        Returns
        -------
        fault_displacement : np.array(N,dtype=float)
            the fault displacement magnitude
        """
        if scale:
            points = self.scale(points, inplace=False)
        vals = np.zeros(points.shape[0])
        for f in self.features:
            if f.type == FeatureType.FAULT:
                disp = f.displacementfeature.evaluate_value(points)
                vals[~np.isnan(disp)] += disp[~np.isnan(disp)]
        return (
            vals * -self.scale_factor
        )  # convert from restoration magnutude to displacement

    def get_feature_by_name(self, feature_name):
        """Returns a feature from the mode given a name


        Parameters
        ----------
        feature_name : string
            the name of the feature

        Returns
        -------
        feature : GeologicalFeature
            the geological feature with the specified name, or none if no feature



        """
        feature_index = self.feature_name_index.get(feature_name, -1)
        if feature_index > -1:
            return self.features[feature_index]
        else:
            logger.error(f"{feature_name} does not exist!")
            return None

    def evaluate_feature_value(self, feature_name, xyz, scale=True):
        """Evaluate the scalar value of the geological feature given the name at locations
        xyz

        Parameters
        ----------
        feature_name : string
            name of the feature
        xyz : np.array((N,3))
            locations to evaluate
        scale : bool, optional
            whether to scale real world points into model scale, by default True

        Returns
        -------
        np.array((N))
            vector of scalar values

        Examples
        --------
        Evaluate on a voxet using model boundaries

        >>> x = np.linspace(model.bounding_box[0, 0], model.bounding_box[1, 0],
                        nsteps[0])
        >>> y = np.linspace(model.bounding_box[0, 1], model.bounding_box[1, 1],
                        nsteps[1])
        >>> z = np.linspace(model.bounding_box[1, 2], model.bounding_box[0, 2],
                        nsteps[2])
        >>> xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        >>> xyz = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        >>> model.evaluate_feature_vaue('feature',xyz,scale=False)

        Evaluate on points in UTM coordinates

        >>> model.evaluate_feature_vaue('feature',utm_xyz)

        """
        feature = self.get_feature_by_name(feature_name)
        if feature:
            scaled_xyz = xyz
            if scale:
                scaled_xyz = self.scale(xyz, inplace=False)
            return feature.evaluate_value(scaled_xyz)
        else:
            return np.zeros(xyz.shape[0])

    def evaluate_feature_gradient(self, feature_name, xyz, scale=True):
        """Evaluate the gradient of the geological feature at a location

        Parameters
        ----------
        feature_name : string
            name of the geological feature
        xyz : np.array((N,3))
            locations to evaluate
        scale : bool, optional
            whether to scale real world points into model scale, by default True

        Returns
        -------
        results : np.array((N,3))
            gradient of the scalar field at the locations specified
        """
        feature = self.get_feature_by_name(feature_name)
        if feature:
            scaled_xyz = xyz
            if scale:
                scaled_xyz = self.scale(xyz, inplace=False)
            return feature.evaluate_gradient(scaled_xyz)
        else:
            return np.zeros(xyz.shape[0])

    def update(self, verbose=False, progressbar=True):
        total_dof = 0
        nfeatures = 0
        for f in self.features:
            if f.type == FeatureType.FAULT:
                nfeatures += 3
                total_dof += f[0].interpolator.nx * 3
                continue
            if isinstance(f, StructuralFrame):
                nfeatures += 3
                total_dof += f[0].interpolator.nx * 3
                continue
            if f.type == FeatureType.INTERPOLATED:
                nfeatures += 1
                total_dof += f.interpolator.nx
                continue
        if verbose == True:
            print(
                f"Updating geological model. There are: \n {nfeatures} geological features that need to be interpolated\n"
            )

        from tqdm.auto import tqdm
        import time

        start = time.time()
        sizecounter = 0

        # Load tqdm with size counter instead of file counter
        with tqdm(total=nfeatures) as pbar:
            buf = 0
            for f in self.features:
                pbar.set_description(f"Interpolating {f.name}")
                f.builder.up_to_date(callback=pbar.update)

        if verbose:
            print(f"Model update took: {time.time()-start} seconds")
