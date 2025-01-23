"""
Main entry point for creating a geological model
"""

from ...utils import getLogger, log_to_file

import numpy as np
import pandas as pd
from typing import List
import pathlib
from ...modelling.features.fault import FaultSegment

from ...modelling.features.builders import (
    FaultBuilder,
    GeologicalFeatureBuilder,
    StructuralFrameBuilder,
    FoldedFeatureBuilder,
)
from ...modelling.features import (
    UnconformityFeature,
    StructuralFrame,
    GeologicalFeature,
    FeatureType,
)
from ...modelling.features.fold import (
    FoldEvent,
    FoldFrame,
)

from ...utils.helper import (
    all_heading,
    gradient_vec_names,
)
from ...utils import strikedip2vector
from ...datatypes import BoundingBox

from ...modelling.intrusions import IntrusionBuilder

from ...modelling.intrusions import IntrusionFrameBuilder


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
        epsion : float
            a fudge factor for isosurfacing, used to make sure surfaces appear
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
        if logfile:
            self.logfile = logfile
            log_to_file(logfile, level=loglevel)

        logger.info("Initialising geological model")
        self.features = []
        self.feature_name_index = {}
        self._data = pd.DataFrame()  # None
        if data is not None:
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

        self.scale_factor = 1.0

        self.bounding_box = BoundingBox(
            dimensions=3,
            origin=np.zeros(3),
            maximum=self.maximum - self.origin,
            global_origin=self.origin,
        )

        self.stratigraphic_column = None

        self.tol = 1e-10 * np.max(self.bounding_box.maximum - self.bounding_box.origin)
        self._dtm = None

    def to_dict(self):
        """
        Convert the geological model to a json string

        Returns
        -------
        json : str
            json string of the geological model
        """
        json = {}
        json["model"] = {}
        json["model"]["features"] = [f.name for f in self.features]
        # json["model"]["data"] = self.data.to_json()
        # json["model"]["origin"] = self.origin.tolist()
        # json["model"]["maximum"] = self.maximum.tolist()
        # json["model"]["nsteps"] = self.nsteps
        json["model"]["stratigraphic_column"] = self.stratigraphic_column
        # json["features"] = [f.to_json() for f in self.features]
        return json

    # @classmethod
    # def from_json(cls,json):
    #     """
    #     Create a geological model from a json string

    #     Parameters
    #     ----------
    #     json : str
    #         json string of the geological model

    #     Returns
    #     -------
    #     model : GeologicalModel
    #         a geological model
    #     """
    #     model = cls(json["model"]["origin"],json["model"]["maximum"],data=None)
    #     model.stratigraphic_column = json["model"]["stratigraphic_column"]
    #     model.nsteps = json["model"]["nsteps"]
    #     model.data = pd.read_json(json["model"]["data"])
    #     model.features = []
    #     for feature in json["features"]:
    #         model.features.append(GeologicalFeature.from_json(feature,model))
    #     return model
    def __str__(self):
        lengths = self.maximum - self.origin
        _str = "GeologicalModel - {} x {} x {}\n".format(*lengths)
        _str += "------------------------------------------ \n"
        _str += "The model contains {} GeologicalFeatures \n".format(len(self.features))
        _str += ""
        _str += "------------------------------------------ \n"
        _str += ""
        _str += "Model origin: {} {} {}\n".format(self.origin[0], self.origin[1], self.origin[2])
        _str += "Model maximum: {} {} {}\n".format(
            self.maximum[0], self.maximum[1], self.maximum[2]
        )
        _str += "Model rescale factor: {} \n".format(self.scale_factor)
        _str += "------------------------------------------ \n"
        _str += "Feature list: \n"
        for feature in self.features:
            _str += "  {} \n".format(feature.name)
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
        faults using fault_params.  faults is a flag that allows for the faults to be
        skipped.

        Parameters
        ----------
        m2l_directory : string
            path to map2loop directory

        Returns
        -------
        (GeologicalModel, dict)
            the created geological model and a dictionary of the map2loop data

        Notes
        ------
        For additional information see :class:`LoopStructural.modelling.input.Map2LoopProcessor`
        and :meth:`LoopStructural.GeologicalModel.from_processor`
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
        """Builds a model from a :class:`LoopStructural.modelling.input.ProcessInputData` object
        This object stores the observations and order of the geological features

        Parameters
        ----------
        processor : ProcessInputData
            any type of ProcessInputData

        Returns
        -------
        GeologicalModel
            a model with all of the features, need to call model.update() to run interpolation
        """
        logger.info("Creating model from processor")
        model = GeologicalModel(processor.origin, processor.maximum)
        model.data = processor.data
        if processor.fault_properties is not None:
            for i in processor.fault_network.faults:
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
                    logger.warning(f"Cannot add splay {edge[1]} or {edge[0]} are not in the model")
                    continue
                splay = False
                if "angle" in properties:
                    if float(properties["angle"]) < 30 and (
                        "dip_dir" not in processor.stratigraphic_column["faults"][edge[0]]
                        or np.abs(
                            processor.stratigraphic_column["faults"][edge[0]]["dip_dir"]
                            - processor.stratigraphic_column["faults"][edge[1]]["dip_dir"]
                        )
                        < 90
                    ):
                        # splay
                        region = model[edge[1]].builder.add_splay(model[edge[0]])

                        model[edge[1]].splay[model[edge[0]].name] = region
                        splay = True
                if splay is False:
                    positive = None
                    if "downthrow_dir" in processor.stratigraphic_column["faults"][edge[0]]:
                        positive = (
                            np.abs(
                                processor.stratigraphic_column["faults"][edge[0]]["downthrow_dir"]
                                - processor.stratigraphic_column["faults"][edge[1]]["downthrow_dir"]
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
                if not f:
                    logger.warning(f"Foliation {s} not added")
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
        if isinstance(model, GeologicalModel):
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
        """Set a dtm to the model.
        The dtm is a function that can be called for dtm(xy) where xy is
        a numpy array of xy locations. The function will return an array of
        z values corresponding to the elevation at xy.

        Parameters
        ----------
        dtm : callable

        """
        if not callable(dtm):
            raise BaseException("DTM must be a callable function \n")
        else:
            self._dtm = dtm

    @property
    def faults(self):
        """Get all of the fault features in the model

        Returns
        -------
        list
            a list of :class:`LoopStructural.modelling.features.FaultSegment`
        """
        faults = []
        for f in self.features:
            if isinstance(f, FaultSegment):
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
    def intrusions(self):
        intrusions = []
        for f in self.features:
            if f.type == "intrusion":
                intrusions.append(f)
        return intrusions

    @property
    def faults_displacement_magnitude(self):
        displacements = []
        for f in self.faults:
            displacements.append(f.displacement)
        return np.array(displacements)

    def feature_names(self):
        return self.feature_name_index.keys()

    def fault_names(self):
        """Get name of all faults in the model

        Returns
        -------
        list
            list of the names of the faults in the model
        """
        return [f.name for f in self.faults]

    def check_inialisation(self):
        if self.data is None:
            logger.error("Data not associated with GeologicalModel. Run set_data")
            return False
        if self.data.shape[0] > 0:
            return True

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
            logger.error("Cannot write to file, dill not installed \n" "pip install dill")
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
            logger.info(f"Adding {feature.name} to model at location {len(self.features)}")
        self._add_domain_fault_above(feature)
        if feature.type == FeatureType.INTERPOLATED:
            self._add_unconformity_above(feature)
        feature.model = self

    def data_for_feature(self, feature_name: str) -> pd.DataFrame:
        """Get all of the data associated with a geological feature

        Parameters
        ----------
        feature_name : str
            the unique identifying name of the feature

        Returns
        -------
        pd.DataFrame
            data frame containing all of the data in the model associated with this feature
        """
        return self.data.loc[self.data["feature_name"] == feature_name, :]

    @property
    def data(self) -> pd.DataFrame:
        return self._data

    @data.setter
    def data(self, data: pd.DataFrame):
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
        if not issubclass(type(data), pd.DataFrame):
            logger.warning("Data is not a pandas data frame, trying to read data frame " "from csv")
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
        # LS wants polarity as -1 or 1, change 0  to -1
        self._data.loc[self._data["polarity"] == 0, "polarity"] = -1.0
        self._data.loc[np.isnan(self._data["w"]), "w"] = 1.0
        if "strike" in self._data and "dip" in self._data:
            logger.info("Converting strike and dip to vectors")
            mask = np.all(~np.isnan(self._data.loc[:, ["strike", "dip"]]), axis=1)
            self._data.loc[mask, gradient_vec_names()] = (
                strikedip2vector(self._data.loc[mask, "strike"], self._data.loc[mask, "dip"])
                * self._data.loc[mask, "polarity"].to_numpy()[:, None]
            )
            self._data.drop(["strike", "dip"], axis=1, inplace=True)
        self._data[['X', 'Y', 'Z', 'val', 'nx', 'ny', 'nz', 'gx', 'gy', 'gz', 'tx', 'ty', 'tz']] = (
            self._data[
                ['X', 'Y', 'Z', 'val', 'nx', 'ny', 'nz', 'gx', 'gy', 'gz', 'tx', 'ty', 'tz']
            ].astype(float)
        )

    def set_model_data(self, data):
        logger.warning("deprecated method. Model data can now be set using the data attribute")
        self.data = data

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
                    ci += 1
        self.stratigraphic_column = stratigraphic_column

    def create_and_add_foliation(
        self,
        series_surface_data: str,
        interpolatortype: str = "FDI",
        nelements: int = 1000,
        tol=None,
        faults=None,
        **kwargs,
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

        Notes
        ------
        This function creates an instance of a
        :class:`LoopStructural.modelling.features.builders.GeologicalFeatureBuilder` and will return
        a :class:`LoopStructural.modelling.features.builders.GeologicalFeature`
        The feature is not interpolated until either
        :meth:`LoopStructural.modelling.features.builders.GeologicalFeature.evaluate_value` is called or
        :meth:`LoopStructural.modelling.core.GeologicalModel.update`

        An interpolator will be chosen by calling :meth:`LoopStructural.GeologicalModel.get_interpolator`

        """
        if not self.check_inialisation():
            logger.warning(f"{series_surface_data} not added, model not initialised")
            return
        # if tol is not specified use the model default
        if tol is None:
            tol = self.tol

        series_builder = GeologicalFeatureBuilder(
            bounding_box=self.bounding_box,
            interpolatortype=interpolatortype,
            nelements=nelements,
            name=series_surface_data,
            model=self,
            **kwargs,
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
        # this support is built for the entire model domain? Possibly would
        # could just pass a regular grid of points - mask by any above unconformities??
        series_builder.build_arguments['domain'] = True
        series_builder.build_arguments["tol"] = tol
        series_feature.type = FeatureType.INTERPOLATED
        self._add_feature(series_feature)
        return series_feature

    def create_and_add_fold_frame(
        self,
        foldframe_data,
        interpolatortype="FDI",
        nelements=1000,
        tol=None,
        buffer=0.1,
        **kwargs,
    ):
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
        if not self.check_inialisation():
            return False
        if tol is None:
            tol = self.tol

        # create fault frame
        #
        fold_frame_builder = StructuralFrameBuilder(
            interpolatortype=interpolatortype,
            bounding_box=self.bounding_box.with_buffer(buffer),
            name=foldframe_data,
            frame=FoldFrame,
            nelements=nelements,
            model=self,
            **kwargs,
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

        fold_frame.type = FeatureType.STRUCTURALFRAME
        fold_frame.builder = fold_frame_builder
        self._add_feature(fold_frame)

        return fold_frame

    def create_and_add_folded_foliation(
        self,
        foliation_data,
        interpolatortype="DFI",
        nelements=10000,
        buffer=0.1,
        fold_frame=None,
        svario=True,
        tol=None,
        invert_fold_norm=False,
        **kwargs,
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

        Notes
        -----

        - Building a folded foliation uses the fold interpolation code from Laurent et al., 2016
        and fold profile fitting from Grose et al., 2017. For more information about the fold modelling
        see :class:`LoopStructural.modelling.features.fold.FoldEvent`,
        :class:`LoopStructural.modelling.features.builders.FoldedFeatureBuilder`

        """
        if not self.check_inialisation():
            return False
        if tol is None:
            tol = self.tol

        if fold_frame is None:
            logger.info("Using last feature as fold frame")
            fold_frame = self.features[-1]
        assert isinstance(fold_frame, FoldFrame), "Please specify a FoldFrame"

        fold = FoldEvent(fold_frame, name=f"Fold_{foliation_data}", invert_norm=invert_fold_norm)

        if interpolatortype != "DFI":
            logger.warning("Folded foliation only supports DFI interpolator, changing to DFI")
            interpolatortype = "DFI"
        series_builder = FoldedFeatureBuilder(
            interpolatortype=interpolatortype,
            bounding_box=self.bounding_box.with_buffer(buffer),
            nelements=nelements,
            fold=fold,
            name=foliation_data,
            svario=svario,
            model=self,
            **kwargs,
        )

        series_builder.add_data_from_data_frame(
            self.data[self.data["feature_name"] == foliation_data]
        )
        self._add_faults(series_builder)
        # series_builder.add_data_to_interpolator(True)
        # build feature

        kwargs["tol"] = tol

        # series_feature = series_builder.build(**kwargs)
        series_feature = series_builder.feature
        series_builder.build_arguments = kwargs
        series_feature.type = FeatureType.INTERPOLATED
        series_feature.fold = fold

        self._add_feature(series_feature)
        return series_feature

    def create_and_add_folded_fold_frame(
        self,
        fold_frame_data,
        interpolatortype="FDI",
        nelements=10000,
        fold_frame=None,
        tol=None,
        **kwargs,
    ):
        """

        Parameters
        ----------
        fold_frame_data : string
            name of the feature to be added

        fold_frame : StructuralFrame, optional
            the fold frame for the fold if not specified uses last feature added

        kwargs : dict
            parameters passed to child functions

        Returns
        -------
        fold_frame : FoldFrame
            created fold frame

        Notes
        -----
        This function build a structural frame where the first coordinate is constrained
        with a fold interpolator.
        Keyword arguments can be included to constrain

        - :meth:`LoopStructural.GeologicalModel.get_interpolator`
        - :class:`LoopStructural.StructuralFrameBuilder`
        - :meth:`LoopStructural.StructuralFrameBuilder.setup`
         - Building a folded foliation uses the fold interpolation code from Laurent et al., 2016
        and fold profile fitting from Grose et al., 2017. For more information about the fold modelling
        see :class:`LoopStructural.modelling.features.fold.FoldEvent`,
        :class:`LoopStructural.modelling.features.builders.FoldedFeatureBuilder`
        """
        if not self.check_inialisation():
            return False
        if tol is None:
            tol = self.tol

        if fold_frame is None:
            logger.info("Using last feature as fold frame")
            fold_frame = self.features[-1]
        assert isinstance(fold_frame, FoldFrame), "Please specify a FoldFrame"
        fold = FoldEvent(fold_frame, name=f"Fold_{fold_frame_data}")

        interpolatortypes = [
            "DFI",
            "FDI",
            "FDI",
        ]
        fold_frame_builder = StructuralFrameBuilder(
            interpolatortype=interpolatortypes,
            bounding_box=self.bounding_box.with_buffer(kwargs.get("buffer", 0.1)),
            nelements=[nelements, nelements, nelements],
            name=fold_frame_data,
            fold=fold,
            frame=FoldFrame,
            model=self,
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

        folded_fold_frame.type = FeatureType.STRUCTURALFRAME

        self._add_feature(folded_fold_frame)

        return folded_fold_frame

    def create_and_add_intrusion(
        self,
        intrusion_name,
        intrusion_frame_name,
        intrusion_frame_parameters={},
        intrusion_lateral_extent_model=None,
        intrusion_vertical_extent_model=None,
        geometric_scaling_parameters={},
        **kwargs,
    ):
        """

        Note
        -----
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
        intrusion_frame_parameters = dictionary

        kwargs

        Returns
        -------
        intrusion feature

        """
        # if intrusions is False:
        #     logger.error("Libraries not installed")
        #     raise Exception("Libraries not installed")

        intrusion_data = self.data[self.data["feature_name"] == intrusion_name].copy()
        intrusion_frame_data = self.data[self.data["feature_name"] == intrusion_frame_name].copy()

        # -- get variables for intrusion frame interpolation
        gxxgz = kwargs.get("gxxgz", 0)
        gxxgy = kwargs.get("gxxgy", 0)
        gyxgz = kwargs.get("gyxgz", 0)

        interpolatortype = kwargs.get("interpolatortype", "PLI")
        # buffer = kwargs.get("buffer", 0.1)
        nelements = kwargs.get("nelements", 1e2)

        weights = [gxxgz, gxxgy, gyxgz]

        intrusion_frame_builder = IntrusionFrameBuilder(
            interpolatortype=interpolatortype,
            bounding_box=self.bounding_box.with_buffer(kwargs.get("buffer", 0.1)),
            nelements=kwargs.get("nelements", 1e2),
            name=intrusion_frame_name,
            model=self,
            **kwargs,
        )

        self._add_faults(intrusion_frame_builder)
        # intrusion_frame_builder.post_intrusion_faults = faults  # LG unused?

        # -- create intrusion frame using intrusion structures (steps and marginal faults) and flow/inflation measurements
        if len(intrusion_frame_parameters) == 0:
            logger.error("Please specify parameters to build intrusion frame")
        intrusion_frame_builder.set_intrusion_frame_parameters(
            intrusion_data, intrusion_frame_parameters
        )
        intrusion_frame_builder.create_constraints_for_c0()

        intrusion_frame_builder.set_intrusion_frame_data(intrusion_frame_data)

        ## -- create intrusion frame
        intrusion_frame_builder.setup(
            nelements=nelements,
            w2=weights[0],
            w1=weights[1],
            gxygz=weights[2],
        )

        intrusion_frame = intrusion_frame_builder.frame

        # -- create intrusion builder to compute distance thresholds along the frame coordinates
        intrusion_builder = IntrusionBuilder(
            intrusion_frame,
            model=self,
            # interpolator=interpolator,
            name=f"{intrusion_name}_feature",
            lateral_extent_model=intrusion_lateral_extent_model,
            vertical_extent_model=intrusion_vertical_extent_model,
            **kwargs,
        )
        intrusion_builder.set_data_for_extent_calculation(intrusion_data)

        intrusion_builder.build_arguments = {
            "geometric_scaling_parameters": geometric_scaling_parameters,
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
            if f.type == FeatureType.FAULT:
                feature_builder.add_fault(f)

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
        Looks through the feature list and adds any the domain_fault to the features
        that already exist in the stack until an unconformity is reached. domain faults
        to the feature. The domain fault masks everything where the fault scalar field
         is < 0 as being active when added to feature.

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

        if feature.type == FeatureType.FAULT:
            return
        for f in reversed(self.features):
            if f.type == FeatureType.UNCONFORMITY and f.name != feature.name:
                logger.info(f"Adding {f.name} as unconformity to {feature.name}")
                feature.add_region(f)
            if f.type == FeatureType.ONLAPUNCONFORMITY and f.name != feature.name:
                feature.add_region(f)
                break

    def add_unconformity(self, feature: GeologicalFeature, value: float) -> UnconformityFeature:
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
            logger.warning("Cannot add unconformtiy, base feature is None")
            return
        # look backwards through features and add the unconformity as a region until
        # we get to an unconformity
        uc_feature = UnconformityFeature(feature, value)
        feature.add_region(uc_feature.inverse())
        for f in reversed(self.features):
            if f.type == FeatureType.UNCONFORMITY:
                logger.debug(f"Reached unconformity {f.name}")
                break
            logger.debug(f"Adding {uc_feature.name} as unconformity to {f.name}")
            if f.type == FeatureType.FAULT or f.type == FeatureType.INACTIVEFAULT:
                continue
            if f == feature:
                continue
            else:
                f.add_region(uc_feature)
        # now add the unconformity to the feature list
        self._add_feature(uc_feature)
        return uc_feature

    def add_onlap_unconformity(self, feature: GeologicalFeature, value: float) -> GeologicalFeature:
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
        feature.regions = []
        uc_feature = UnconformityFeature(feature, value, False, onlap=True)
        feature.add_region(uc_feature.inverse())
        for f in reversed(self.features):
            if f.type == FeatureType.UNCONFORMITY:
                # f.add_region(uc_feature)
                continue
            if f.type == FeatureType.FAULT:
                continue
            if f != feature:
                f.add_region(uc_feature)
        self._add_feature(uc_feature.inverse())

        return uc_feature

    def create_and_add_domain_fault(
        self, fault_surface_data, nelements=10000, interpolatortype="FDI", **kwargs
    ):
        """
        Parameters
        ----------
        fault_surface_data : string
            name of the domain fault data in the data frame

        Returns
        -------
        domain_Fault : GeologicalFeature
            the created domain fault

        Notes
        -----
        * :meth:`LoopStructural.GeologicalModel.get_interpolator`

        """
        domain_fault_feature_builder = GeologicalFeatureBuilder(
            bounding_box=self.bounding_box,
            interpolatortype=interpolatortype,
            nelements=nelements,
            name=fault_surface_data,
            model=self,
            **kwargs,
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
        domain_fault.type = FeatureType.DOMAINFAULT
        self._add_feature(domain_fault)
        self._add_domain_fault_below(domain_fault)

        domain_fault_uc = UnconformityFeature(domain_fault, 0)
        # iterate over existing features and add the unconformity as a region
        # so the feature is only evaluated where the unconformity is positive
        return domain_fault_uc

    def create_and_add_fault(
        self,
        fault_surface_data,
        displacement,
        interpolatortype="FDI",
        tol=None,
        fault_slip_vector=None,
        fault_normal_vector=None,
        fault_center=None,
        major_axis=None,
        minor_axis=None,
        intermediate_axis=None,
        faultfunction="BaseFault",
        faults=[],
        force_mesh_geometry: bool = False,
        points: bool = False,
        fault_buffer=0.2,
        fault_trace_anisotropy=0.0,
        fault_dip=90,
        fault_dip_anisotropy=0.0,
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

        Notes
        -----
        * :meth:`LoopStructural.GeologicalModel.get_interpolator`
        * :class:`LoopStructural.modelling.features.builders.FaultBuilder`
        * :meth:`LoopStructural.modelling.features.builders.FaultBuilder.setup`
        """
        if "fault_extent" in kwargs and major_axis is None:
            major_axis = kwargs["fault_extent"]
        if "fault_influence" in kwargs and minor_axis is None:
            minor_axis = kwargs["fault_influence"]
        if "fault_vectical_radius" in kwargs and intermediate_axis is None:
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
        if fault_slip_vector is not None:
            fault_slip_vector = np.array(fault_slip_vector, dtype="float")
        if fault_center is not None:
            fault_center = np.array(fault_center, dtype="float")

        for k, v in kwargs.items():
            logger.info(f"{k}: {v}")

        if tol is None:
            tol = self.tol
            # divide the tolerance by half of the minor axis, as this is the equivalent of the distance
            # of the unit vector
            # if minor_axis:
            # tol *= 0.1*minor_axis

        if displacement == 0:
            logger.warning(f"{fault_surface_data} displacement is 0")

        if "data_region" in kwargs:
            kwargs.pop("data_region")
            logger.error("kwarg data_region currently not supported, disabling")
        displacement_scaled = displacement / self.scale_factor
        fault_frame_builder = FaultBuilder(
            interpolatortype,
            bounding_box=self.bounding_box,
            nelements=kwargs.pop("nelements", 1e4),
            name=fault_surface_data,
            model=self,
            **kwargs,
        )
        fault_frame_data = self.data.loc[self.data["feature_name"] == fault_surface_data].copy()
        self._add_faults(fault_frame_builder, features=faults)
        # add data
        fault_frame_data = self.data.loc[self.data["feature_name"] == fault_surface_data].copy()

        if fault_center is not None and ~np.isnan(fault_center).any():
            fault_center = self.scale(fault_center, inplace=False)
        if minor_axis:
            minor_axis = minor_axis / self.scale_factor
        if major_axis:
            major_axis = major_axis / self.scale_factor
        if intermediate_axis:
            intermediate_axis = intermediate_axis / self.scale_factor
        fault_frame_builder.create_data_from_geometry(
            fault_frame_data=fault_frame_data,
            fault_center=fault_center,
            fault_normal_vector=fault_normal_vector,
            fault_slip_vector=fault_slip_vector,
            minor_axis=minor_axis,
            major_axis=major_axis,
            intermediate_axis=intermediate_axis,
            points=points,
            force_mesh_geometry=force_mesh_geometry,
            fault_buffer=fault_buffer,
            fault_trace_anisotropy=fault_trace_anisotropy,
            fault_dip=fault_dip,
            fault_dip_anisotropy=fault_dip_anisotropy,
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
                fault.add_region(f)
                break
        if displacement == 0:
            fault.type = FeatureType.INACTIVEFAULT
        self._add_feature(fault)

        return fault

    # TODO move rescale to bounding box/transformer
    def rescale(self, points: np.ndarray, inplace: bool = False) -> np.ndarray:
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
        if not inplace:
            points = points.copy()
        points *= self.scale_factor
        points += self.origin
        return points

    # TODO move scale to bounding box/transformer
    def scale(self, points: np.ndarray, inplace: bool = False) -> np.ndarray:
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
        points : np.a::rray((N,3),dtype=double)

        """
        points = np.array(points).astype(float)
        if not inplace:
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
        return self.bounding_box.regular_grid(nsteps=nsteps, shuffle=shuffle, order=order)

    def evaluate_model(self, xyz: np.ndarray, scale: bool = True) -> np.ndarray:
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
        # set strat id to -1 to identify which areas of the model aren't covered
        strat_id[:] = -1
        if self.stratigraphic_column is None:
            logger.warning("No stratigraphic column defined")
            return strat_id
        for group in reversed(self.stratigraphic_column.keys()):
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

    def evaluate_model_gradient(self, points: np.ndarray, scale: bool = True) -> np.ndarray:
        """Evaluate the gradient of the stratigraphic column at each location

        Parameters
        ----------
        points : np.ndarray
            location to evaluate
        scale : bool, optional
            whether to scale the points into model domain, by default True

        Returns
        -------
        np.ndarray
            N,3 array of gradient vectors
        """
        xyz = np.array(points)
        if scale:
            xyz = self.scale(xyz, inplace=False)
        grad = np.zeros(xyz.shape)
        for group in reversed(self.stratigraphic_column.keys()):
            if group == "faults":
                continue
            feature_id = self.feature_name_index.get(group, -1)
            if feature_id >= 0:
                feature = self.features[feature_id]
                gradt = feature.evaluate_gradient(xyz)
                grad[~np.isnan(gradt).any(axis=1)] = gradt[~np.isnan(gradt).any(axis=1)]
            if feature_id == -1:
                logger.error(f"Model does not contain {group}")

        return grad

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
        return vals * -self.scale_factor  # convert from restoration magnutude to displacement

    def get_feature_by_name(self, feature_name) -> GeologicalFeature:
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
            raise ValueError(f"{feature_name} does not exist!")

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
        if verbose:
            print(
                f"Updating geological model. There are: \n {nfeatures} \
                    geological features that need to be interpolated\n"
            )

        if progressbar:
            try:
                from tqdm.auto import tqdm

                # Load tqdm with size counter instead of file counter
                with tqdm(total=nfeatures) as pbar:
                    for f in self.features:
                        pbar.set_description(f"Interpolating {f.name}")
                        f.builder.up_to_date(callback=pbar.update)
                return
            except ImportError:
                logger.warning("Failed to import tqdm, disabling progress bar")

        for f in self.features:
            f.builder.up_to_date()

    def stratigraphic_ids(self):
        """Return a list of all stratigraphic ids in the model

        Returns
        -------
        ids : list
            list of unique stratigraphic ids, featurename, unit name and min and max scalar values
        """
        ids = []
        if self.stratigraphic_column is None:
            logger.warning('No stratigraphic column defined')
            return ids
        for group in self.stratigraphic_column.keys():
            if group == "faults":
                continue
            for name, series in self.stratigraphic_column[group].items():
                ids.append([series["id"], group, name, series['min'], series['max']])
        return ids

    def get_fault_surfaces(self, faults: List[str] = []):
        surfaces = []
        if len(faults) == 0:
            faults = self.fault_names()

        for f in faults:
            surfaces.extend(self.get_feature_by_name(f).surfaces([0], self.bounding_box))
        return surfaces

    def get_stratigraphic_surfaces(self, units: List[str] = [], bottoms: bool = True):
        ## TODO change the stratigraphic column to its own class and have methods to get the relevant surfaces
        surfaces = []
        units = []
        if self.stratigraphic_column is None:
            return []
        for group in self.stratigraphic_column.keys():
            if group == "faults":
                continue
            for series in self.stratigraphic_column[group].values():
                series['feature_name'] = group
                units.append(series)
        unit_table = pd.DataFrame(units)
        for u in unit_table['feature_name'].unique():

            values = unit_table.loc[unit_table['feature_name'] == u, 'min' if bottoms else 'max']
            if 'name' not in unit_table.columns:
                unit_table['name'] = unit_table['feature_name']

            names = unit_table[unit_table['feature_name'] == u]['name']
            values = values.loc[~np.logical_or(values == np.inf, values == -np.inf)]
            surfaces.extend(
                self.get_feature_by_name(u).surfaces(
                    values.to_list(), self.bounding_box, name=names.loc[values.index].to_list()
                )
            )

        return surfaces

    def get_block_model(self, name='block model'):
        grid = self.bounding_box.structured_grid(name=name)

        grid.cell_properties['stratigraphy'] = self.evaluate_model(
            self.rescale(self.bounding_box.cell_centers())
        )
        return grid, self.stratigraphic_ids()

    def save(
        self,
        filename: str,
        block_model: bool = True,
        stratigraphic_surfaces=True,
        fault_surfaces=True,
        stratigraphic_data=True,
        fault_data=True,
    ):
        path = pathlib.Path(filename)
        extension = path.suffix
        parent = path.parent
        name = path.stem
        stratigraphic_surfaces = self.get_stratigraphic_surfaces()
        if fault_surfaces:
            for s in self.get_fault_surfaces():
                ## geoh5 can save everything into the same file
                if extension == ".geoh5" or extension == '.omf':
                    s.save(filename)
                else:
                    s.save(f'{parent}/{name}_{s.name}{extension}')
        if stratigraphic_surfaces:
            for s in self.get_stratigraphic_surfaces():
                if extension == ".geoh5" or extension == '.omf':
                    s.save(filename)
                else:
                    s.save(f'{parent}/{name}_{s.name}{extension}')
        if block_model:
            grid, ids = self.get_block_model()
            if extension == ".geoh5" or extension == '.omf':
                grid.save(filename)
            else:
                grid.save(f'{parent}/{name}_block_model{extension}')
        if stratigraphic_data:
            if self.stratigraphic_column is not None:
                for group in self.stratigraphic_column.keys():
                    if group == "faults":
                        continue
                    for data in self.__getitem__(group).get_data():
                        if extension == ".geoh5" or extension == '.omf':
                            data.save(filename)
                        else:
                            data.save(f'{parent}/{name}_{group}_data{extension}')
        if fault_data:
            for f in self.fault_names():
                for d in self.__getitem__(f).get_data():
                    if extension == ".geoh5" or extension == '.omf':

                        d.save(filename)
                    else:
                        d.save(f'{parent}/{name}_{group}{extension}')
