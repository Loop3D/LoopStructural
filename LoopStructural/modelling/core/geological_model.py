"""
Main entry point for creating a geological model
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from LoopStructural import LoopStructuralConfig

from ...geometry import BoundingBox
from ...modelling.features import (
    BaseFeature,
    FeatureType,
    GeologicalFeature,
    StructuralFrame,
    UnconformityFeature,
)
from ...modelling.features.fault import FaultSegment
from ...modelling.features.fold import FoldFrame
from ...utils import LoopValueError, getLogger, public_api, strikedip2vector, timed_stage
from ...utils.helper import (
    all_heading,
    gradient_vec_names,
)
from ..features._feature_converters import (
    add_fold_to_feature as _add_fold_to_feature,
)
from ..features._feature_converters import (
    convert_feature_to_structural_frame as _convert_feature_to_structural_frame,
)
from ._feature_registry import FeatureBuilderRegistry
from ._model_exporter import ModelExporter
from ._model_feature_factory import ModelFeatureFactory
from ._model_relationships import FeatureRelationshipManager
from ._model_serializer import ModelSerializer
from .stratigraphic_column import StratigraphicColumn

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


    """

    def __init__(self, *args):
        """
        Parameters
        ----------
        bounding_box : BoundingBox
            the bounding box of the model
        origin : np.array(3,dtype=doubles)
            the origin of the model
        maximum : np.array(3,dtype=doubles)
            the maximum of the model

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
        args = list(args)
        if len(args) == 0:
            raise ValueError("Must provide either bounding_box or origin and maximum")
        if len(args) == 1:
            bounding_box = args[0]
            if not isinstance(bounding_box, BoundingBox):
                raise ValueError("Must provide a bounding box")
            # A pre-built BoundingBox already carries its own local transform
            # (defaulted to a zero local_origin/identity rotation in
            # BoundingBox.__init__, or explicitly configured by the caller via
            # set_local_transform/from_dict/with_buffer), so we deliberately
            # do not call set_local_transform again here -- doing so would
            # override any anchoring the caller already set up.
            self.bounding_box = bounding_box
        if len(args) == 2:
            origin = np.array(args[0])
            maximum = np.array(args[1])
            if not isinstance(origin, np.ndarray) or not isinstance(maximum, np.ndarray):
                raise ValueError("Must provide origin and maximum as numpy arrays")
            self.bounding_box = BoundingBox(
                dimensions=3,
                origin=origin,
                maximum=maximum,
            )
            # Anchor the interpolation frame near zero for numerical
            # conditioning without leaking the shift into the public
            # origin/maximum, which now stay in world coordinates.
            self.bounding_box.set_local_transform(local_origin=origin)
        logger.info("Initialising geological model")
        self.features = []
        self.feature_name_index = {}
        self._data = pd.DataFrame()  # None

        self._stratigraphic_column = StratigraphicColumn()

        self.tol = 1e-10 * np.max(self.bounding_box.maximum - self.bounding_box.origin)
        self._dtm = None

    @public_api(tier="stable")
    def to_dict(self):
        """
        Convert the geological model to a json string

        Returns
        -------
        json : str
            json string of the geological model
        """
        return ModelSerializer.to_dict(self)

    @public_api(tier="provisional")
    def to_recipe_dict(self, data_reference=None):
        """Return a YAML/JSON-friendly recipe for rebuilding the model.

        This captures the construction inputs needed for stage 3a: bounding
        box, stratigraphic column, and either inline model data or a file
        reference to it.
        """
        return ModelSerializer.to_recipe_dict(self, data_reference=data_reference)

    @classmethod
    @public_api(tier="provisional")
    def from_recipe_dict(cls, recipe):
        """Rebuild a geological model from a recipe dictionary."""
        return ModelSerializer.from_recipe_dict(cls, recipe)

    @public_api(tier="provisional")
    def to_recipe_json(self, data_reference=None, indent=2):
        """Return a JSON-formatted string of the recipe.

        Parameters
        ----------
        data_reference : str, optional
            Path to an external CSV file to reference instead of embedding
            data inline in the JSON. If None, data is embedded.
        indent : int, optional
            JSON indentation level. Default is 2.

        Returns
        -------
        str
            JSON-formatted recipe string.
        """
        return ModelSerializer.to_recipe_json(self, data_reference=data_reference, indent=indent)

    @classmethod
    @public_api(tier="provisional")
    def from_recipe_json(cls, json_str):
        """Rebuild a geological model from a JSON-formatted recipe string.

        Parameters
        ----------
        json_str : str
            JSON-formatted recipe string.

        Returns
        -------
        GeologicalModel
            The reconstructed geological model.

        Raises
        ------
        TypeError
            If json_str is not a string or does not parse as valid JSON.
        """
        return ModelSerializer.from_recipe_json(cls, json_str)

    @public_api(tier="provisional")
    def save_recipe(self, filename, data_reference=None):
        """Save the recipe to a JSON file.

        Parameters
        ----------
        filename : str or Path
            Path to the output JSON file.
        data_reference : str, optional
            Path to an external CSV file to reference instead of embedding
            data inline in the JSON. If None, data is embedded.
        """
        ModelSerializer.save_recipe(self, filename, data_reference=data_reference)

    @classmethod
    @public_api(tier="provisional")
    def load_recipe(cls, filename):
        """Load a geological model from a recipe JSON file.

        Parameters
        ----------
        filename : str or Path
            Path to the recipe JSON file.

        Returns
        -------
        GeologicalModel
            The reconstructed geological model.
        """
        return ModelSerializer.load_recipe(cls, filename)

    def __str__(self):
        return f"GeologicalModel with {len(self.features)} features"

    def _ipython_key_completions_(self):
        return self.feature_name_index.keys()

    def prepare_data(self, data: pd.DataFrame, include_feature_name: bool = True) -> pd.DataFrame:
        data = data.copy()
        # Data is kept in world coordinates end-to-end; the interpolator
        # projects into its local frame when constraints are set.

        if "type" in data:
            logger.warning("'type' is deprecated replace with 'feature_name' \n")
            data.rename(columns={"type": "feature_name"}, inplace=True)
        if "feature_name" not in data and include_feature_name:
            logger.error("Data does not contain 'feature_name' column")
            raise ValueError("Cannot load data")
        for h in all_heading():
            if h not in data:
                data[h] = np.nan
                if h == "w":
                    data[h] = 1.0
                if h == "coord":
                    data[h] = 0
                if h == "polarity":
                    data[h] = 1.0
        # LS wants polarity as -1 or 1, change 0  to -1
        data.loc[data["polarity"] == 0, "polarity"] = -1.0
        data.loc[np.isnan(data["w"]), "w"] = 1.0
        if "strike" in data and "dip" in data:
            logger.info("Converting strike and dip to vectors")
            mask = np.all(~np.isnan(data.loc[:, ["strike", "dip"]]), axis=1)
            data.loc[mask, gradient_vec_names()] = (
                strikedip2vector(data.loc[mask, "strike"], data.loc[mask, "dip"])
                * data.loc[mask, "polarity"].to_numpy()[:, None]
            )
            data.drop(["strike", "dip"], axis=1, inplace=True)
        data[['X', 'Y', 'Z', 'val', 'nx', 'ny', 'nz', 'gx', 'gy', 'gz', 'tx', 'ty', 'tz']] = data[
            ['X', 'Y', 'Z', 'val', 'nx', 'ny', 'nz', 'gx', 'gy', 'gz', 'tx', 'ty', 'tz']
        ].astype(float)
        return data

    @classmethod
    @public_api(tier="stable")
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
                if "angle" in properties and float(properties["angle"]) < 30 and (
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
        for s in processor.stratigraphic_column:
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
    @public_api(tier="stable")
    def from_file(cls, file, allow_pickle: bool = True):
        """Load a geological model from file

        .. warning::
            Model files are loaded using `dill` (an extension of `pickle`).
            Unpickling data is **not safe** against maliciously constructed
            data: loading a file from an untrusted or unauthenticated source
            can execute arbitrary code on your machine. Only call
            ``from_file`` on files you created yourself or that come from a
            source you fully trust. If you need to load model definitions
            from an untrusted source, use the JSON/dictionary-based
            ``GeologicalModel.from_recipe_dict``/``to_recipe_dict`` recipe
            format instead, or pass ``allow_pickle=False`` here to make sure
            pickle-based loading is refused outright.

        Parameters
        ----------
        file : string
            path to the file
        allow_pickle : bool, optional
            whether to allow loading the file using `dill`/`pickle`, by
            default True. Set to False to refuse pickle-based deserialisation
            (e.g. when the file may come from an untrusted source) -- in that
            case a :class:`LoopValueError` is raised instead of attempting to
            unpickle the file. Use ``GeologicalModel.from_recipe_dict`` for a
            safe, JSON-based alternative.

        Returns
        -------
        GeologicalModel
            the geological model object
        """
        return ModelSerializer.from_file(cls, file, allow_pickle=allow_pickle)

    def __getitem__(self, feature_name):
        """Accessor for feature in features using feature_name_index

        Parameters
        ----------
        feature_name : string
            name of the feature to return
        """
        return self.get_feature_by_name(feature_name)

    def __setitem__(self, feature_name, feature):
        """Set a feature in the model using feature_name_index

        Parameters
        ----------
        feature_name : string
            name of the feature to set
        feature : GeologicalFeature
            the geological feature to set
        """
        if not issubclass(type(feature), BaseFeature):
            raise TypeError("feature must be a GeologicalFeature")
        if feature.name != feature_name:
            raise ValueError("feature name does not match key")
        self._add_feature(feature)

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
            raise TypeError("DTM must be a callable function")
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

    @public_api(tier="stable")
    def feature_names(self):
        return self.feature_name_index.keys()

    @public_api(tier="stable")
    def fault_names(self):
        """Get name of all faults in the model

        Returns
        -------
        list
            list of the names of the faults in the model
        """
        return [f.name for f in self.faults]

    @public_api(tier="stable")
    def to_file(self, file):
        """Save a model to a pickle file requires dill

        Parameters
        ----------
        file : string
            path to file location
        """
        ModelSerializer.to_file(self, file)

    def _add_feature(self, feature, index: int | None = None):
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
            if index is not None:
                if index < 0 or index > len(self.features):
                    raise IndexError(f"Index {index} out of bounds for features list")
                self.features.insert(index, feature)
                self.feature_name_index[feature.name] = index
                logger.info(f"Adding {feature.name} to model at location {index}")
                for feature_index, feature_in_list in enumerate(self.features):
                    self.feature_name_index[feature_in_list.name] = feature_index
            else:
                self.features.append(feature)
                self.feature_name_index[feature.name] = len(self.features) - 1
                logger.info(f"Adding {feature.name} to model at location {len(self.features)}")
        FeatureRelationshipManager.add_domain_fault_above(self, feature)
        if feature.type == FeatureType.INTERPOLATED:
            FeatureRelationshipManager.add_unconformity_above(self, feature)
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
         feature_name, X, Y, Z, nx, ny, nz, val, strike, dip, dip_dir, plunge,
         plunge_dir, azimuth

        Returns
        -------
        Note
        ----
        feature_name can be any unique identifier for the feature the data point
        'eg' 'S0', 'S2', 'F1_axis'
        it is then used by the create functions to get the correct data
        """
        if data is None:
            return
        if not issubclass(type(data), pd.DataFrame):
            logger.warning("Data is not a pandas data frame, trying to read data frame " "from csv")
            try:
                data = pd.read_csv(data)
            except (OSError, ValueError, pd.errors.ParserError) as e:
                logger.error(f"Could not load pandas data frame from data: {e}")
                raise LoopValueError("Cannot load data") from e
        logger.info(f"Adding data to GeologicalModel with {len(data)} data points")
        self._data = data.copy()
        # self._data[['X','Y','Z']] = self.bounding_box.project(self._data[['X','Y','Z']].to_numpy())

    def set_model_data(self, data):
        logger.warning("deprecated method. Model data can now be set using the data attribute")
        self.data = data.copy()

    @property
    def stratigraphic_column(self):
        """Get the stratigraphic column of the model

        Returns
        -------
        StratigraphicColumn
            the stratigraphic column of the model
        """
        return self._stratigraphic_column

    @stratigraphic_column.setter
    def stratigraphic_column(self, stratigraphic_column: StratigraphicColumn | dict):
        """Set the stratigraphic column of the model

        Parameters
        ----------
        stratigraphic_column : StratigraphicColumn
            the stratigraphic column to set
        """
        if isinstance(stratigraphic_column, dict):
            self.set_stratigraphic_column(stratigraphic_column)
            return
        elif not isinstance(stratigraphic_column, StratigraphicColumn):
            raise TypeError("stratigraphic_column must be a StratigraphicColumn object")
        self._stratigraphic_column = stratigraphic_column

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
        warnings.warn(
            "set_stratigraphic_column is deprecated, use model.stratigraphic_column.add_units instead",
            DeprecationWarning,
            stacklevel=2,
        )
        self.stratigraphic_column.clear(basement=False)
        # if the colour for a unit hasn't been specified we can just sample from
        # a colour map e.g. tab20
        logger.info("Adding stratigraphic column to model")
        for i, g in enumerate(stratigraphic_column.keys()):
            if g == 'faults':
                logger.info('Not adding faults to stratigraphic column')
                continue
            for u in stratigraphic_column[g]:
                thickness = 0
                if "min" in stratigraphic_column[g][u] and "max" in stratigraphic_column[g][u]:
                    min_val = stratigraphic_column[g][u]["min"]
                    max_val = stratigraphic_column[g][u].get("max", None)
                    thickness = max_val - min_val if max_val is not None else None
                logger.info(
                    f"""
                               model.stratigraphic_column.add_unit({u},
                               colour={stratigraphic_column[g][u].get("colour", None)},
                                 thickness={thickness})"""
                )
                self.stratigraphic_column.add_unit(
                    u,
                    colour=stratigraphic_column[g][u].get("colour", None),
                    thickness=thickness,
                )

            self.stratigraphic_column.add_unconformity(
                name=f"{g}unconformity",
            )
            self.stratigraphic_column.group_mapping[f'Group_{i}'] = g

    @public_api(tier="provisional")
    def create_and_add_feature(self, feature_type: str, name: str, **params):
        """Create a feature of the given type and add it to the model.

        Generic dispatch entry point backed by
        :class:`LoopStructural.modelling.core._feature_registry.FeatureBuilderRegistry`
        (see ``API.md``). The 7 built-in feature types (``foliation``,
        ``fold_frame``, ``folded_foliation``, ``folded_fold_frame``,
        ``intrusion``, ``domain_fault``, ``fault``) are registered against
        the same builder logic the ``create_and_add_*`` convenience
        methods use - this is the extension point for new feature types
        without modifying ``GeologicalModel``.

        Parameters
        ----------
        feature_type : str
            registered feature type, see
            :meth:`FeatureBuilderRegistry.registered_types`
        name : str
            corresponding to the feature_name in the data
        **params
            forwarded to the registered factory

        Returns
        -------
        feature : BaseFeature
            the created geological feature, or None if it could not be built
        """
        return FeatureBuilderRegistry.create(feature_type, self, name, **params)

    @public_api(tier="stable")
    def create_and_add_foliation(
        self,
        series_surface_name: str,
        *,
        index: int | None = None,
        data: pd.DataFrame | None = None,
        interpolatortype: str = "FDI",
        nelements: int = LoopStructuralConfig.nelements,
        tol=None,
        faults=None,
        **kwargs,
    ):
        """Create a foliation feature and add it to the model.

        See :meth:`~._model_feature_factory.ModelFeatureFactory.build_foliation`
        for parameter documentation. Thin
        wrapper around :meth:`create_and_add_feature` (see ``API.md``);
        kept as a stable, unchanged entry point.
        """
        return self.create_and_add_feature(
            "foliation",
            series_surface_name,
            index=index,
            data=data,
            interpolatortype=interpolatortype,
            nelements=nelements,
            tol=tol,
            faults=faults,
            **kwargs,
        )

    @public_api(tier="stable")
    def create_and_add_fold_frame(
        self,
        fold_frame_name: str,
        *,
        index: int | None = None,
        data=None,
        interpolatortype="FDI",
        nelements=LoopStructuralConfig.nelements,
        tol=None,
        buffer=0.1,
        **kwargs,
    ):
        """Create a fold frame and add it to the model.

        See :meth:`~._model_feature_factory.ModelFeatureFactory.build_fold_frame`
        for parameter documentation. Thin
        wrapper around :meth:`create_and_add_feature` (see ``API.md``);
        kept as a stable, unchanged entry point.
        """
        return self.create_and_add_feature(
            "fold_frame",
            fold_frame_name,
            index=index,
            data=data,
            interpolatortype=interpolatortype,
            nelements=nelements,
            tol=tol,
            buffer=buffer,
            **kwargs,
        )

    @public_api(tier="stable")
    def create_and_add_folded_foliation(
        self,
        foliation_name,
        *,
        index: int | None = None,
        data=None,
        interpolatortype="DFI",
        nelements=LoopStructuralConfig.nelements,
        buffer=0.1,
        fold_frame=None,
        svario=True,
        tol=None,
        invert_fold_norm=False,
        **kwargs,
    ):
        """Create a folded foliation and add it to the model.

        See :meth:`~._model_feature_factory.ModelFeatureFactory.build_folded_foliation`
        for parameter documentation.
        Thin wrapper around :meth:`create_and_add_feature` (see
        ``API.md``); kept as a stable, unchanged entry point.
        """
        return self.create_and_add_feature(
            "folded_foliation",
            foliation_name,
            index=index,
            data=data,
            interpolatortype=interpolatortype,
            nelements=nelements,
            buffer=buffer,
            fold_frame=fold_frame,
            svario=svario,
            tol=tol,
            invert_fold_norm=invert_fold_norm,
            **kwargs,
        )

    @public_api(tier="stable")
    def create_and_add_folded_fold_frame(
        self,
        fold_frame_name: str,
        *,
        index: int | None = None,
        data: pd.DataFrame | None = None,
        interpolatortype="FDI",
        nelements=LoopStructuralConfig.nelements,
        fold_frame=None,
        tol=None,
        **kwargs,
    ):
        """Create a folded fold frame and add it to the model.

        See :meth:`~._model_feature_factory.ModelFeatureFactory.build_folded_fold_frame`
        for parameter documentation.
        Thin wrapper around :meth:`create_and_add_feature` (see
        ``API.md``); kept as a stable, unchanged entry point.
        """
        return self.create_and_add_feature(
            "folded_fold_frame",
            fold_frame_name,
            index=index,
            data=data,
            interpolatortype=interpolatortype,
            nelements=nelements,
            fold_frame=fold_frame,
            tol=tol,
            **kwargs,
        )

    @public_api(tier="stable")
    def create_and_add_intrusion(
        self,
        intrusion_name,
        intrusion_frame_name,
        *,
        intrusion_frame_parameters=None,
        intrusion_lateral_extent_model=None,
        intrusion_vertical_extent_model=None,
        geometric_scaling_parameters=None,
        **kwargs,
    ):
        """Create an intrusion and add it to the model.

        See :meth:`~._model_feature_factory.ModelFeatureFactory.build_intrusion`
        for parameter documentation. Thin
        wrapper around :meth:`create_and_add_feature` (see ``API.md``);
        kept as a stable, unchanged entry point.
        """
        if intrusion_frame_parameters is None:
            intrusion_frame_parameters = {}
        if geometric_scaling_parameters is None:
            geometric_scaling_parameters = {}
        return self.create_and_add_feature(
            "intrusion",
            intrusion_name,
            intrusion_frame_name=intrusion_frame_name,
            intrusion_frame_parameters=intrusion_frame_parameters,
            intrusion_lateral_extent_model=intrusion_lateral_extent_model,
            intrusion_vertical_extent_model=intrusion_vertical_extent_model,
            geometric_scaling_parameters=geometric_scaling_parameters,
            **kwargs,
        )

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
        FeatureRelationshipManager.add_faults(self, feature_builder, features=features)

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
        FeatureRelationshipManager.add_domain_fault_below(self, domain_fault)

    @public_api(tier="stable")
    def add_unconformity(
        self, feature: GeologicalFeature, value: float, index: int | None = None
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
        return FeatureRelationshipManager.add_unconformity(self, feature, value, index=index)

    @public_api(tier="stable")
    def add_onlap_unconformity(
        self, feature: GeologicalFeature, value: float, index: int | None = None
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
        return FeatureRelationshipManager.add_onlap_unconformity(self, feature, value, index=index)

    @public_api(tier="provisional")
    def add_fold_to_feature(
        self, feature_name: str, fold_frame: FoldFrame, **kwargs
    ) -> GeologicalFeature:
        """Add a fold to an already-built feature, replacing it in the model.

        Promoted (``API.md``) from the previously-private
        ``LoopStructural.modelling.features._feature_converters.add_fold_to_feature``,
        which the QGIS plugin imports directly today; that private module
        is unchanged, this is a public, tested entry point for the same
        behaviour.

        Parameters
        ----------
        feature_name : str
            name of an existing feature already in the model
        fold_frame : FoldFrame
            the fold frame to fold the feature around
        **kwargs
            forwarded to :class:`FoldEvent` / ``FoldedFeatureBuilder.from_feature_builder``

        Returns
        -------
        feature : GeologicalFeature
            the folded feature, replacing the original in the model
        """
        feature = self.get_feature_by_name(feature_name)
        folded_feature = _add_fold_to_feature(feature, fold_frame, **kwargs)
        self._add_feature(folded_feature)
        return folded_feature

    @public_api(tier="provisional")
    def convert_feature_to_structural_frame(self, feature_name: str, **kwargs) -> StructuralFrame:
        """Convert an already-built feature into a structural frame, replacing it in the model.

        Promoted (``API.md``) from the previously-private
        ``LoopStructural.modelling.features._feature_converters.convert_feature_to_structural_frame``.

        Parameters
        ----------
        feature_name : str
            name of an existing feature already in the model
        **kwargs
            forwarded to ``StructuralFrameBuilder.from_feature_builder``

        Returns
        -------
        frame : StructuralFrame
            the structural frame, replacing the original feature in the model
        """
        feature = self.get_feature_by_name(feature_name)
        frame = _convert_feature_to_structural_frame(feature, **kwargs)
        self._add_feature(frame)
        return frame

    @public_api(tier="stable")
    def create_and_add_domain_fault(
        self,
        fault_surface_data,
        *,
        nelements=LoopStructuralConfig.nelements,
        interpolatortype="FDI",
        index: int | None = None,
        **kwargs,
    ):
        """Create a domain fault and add it to the model.

        See :meth:`~._model_feature_factory.ModelFeatureFactory.build_domain_fault`
        for parameter documentation. Thin
        wrapper around :meth:`create_and_add_feature` (see ``API.md``);
        kept as a stable, unchanged entry point.
        """
        return self.create_and_add_feature(
            "domain_fault",
            fault_surface_data,
            nelements=nelements,
            interpolatortype=interpolatortype,
            index=index,
            **kwargs,
        )

    @public_api(tier="stable")
    def create_and_add_fault(
        self,
        fault_name: str,
        displacement: float,
        *,
        index: int | None = None,
        data: pd.DataFrame | None = None,
        interpolatortype="FDI",
        tol=None,
        fault_slip_vector=None,
        fault_normal_vector=None,
        fault_center=None,
        major_axis=None,
        minor_axis=None,
        intermediate_axis=None,
        faultfunction="BaseFault",
        faults=None,
        force_mesh_geometry: bool = False,
        points: bool = False,
        fault_buffer=0.2,
        fault_trace_anisotropy=0.0,
        fault_dip=90,
        fault_dip_anisotropy=0.0,
        fault_pitch=None,
        **kwargs,
    ):
        """Create a fault and add it to the model.

        See :meth:`~._model_feature_factory.ModelFeatureFactory.build_fault`
        for parameter documentation. Thin
        wrapper around :meth:`create_and_add_feature` (see ``API.md``);
        kept as a stable, unchanged entry point.
        """
        if faults is None:
            faults = []
        return self.create_and_add_feature(
            "fault",
            fault_name,
            displacement=displacement,
            index=index,
            data=data,
            interpolatortype=interpolatortype,
            tol=tol,
            fault_slip_vector=fault_slip_vector,
            fault_normal_vector=fault_normal_vector,
            fault_center=fault_center,
            major_axis=major_axis,
            minor_axis=minor_axis,
            intermediate_axis=intermediate_axis,
            faultfunction=faultfunction,
            faults=faults,
            force_mesh_geometry=force_mesh_geometry,
            points=points,
            fault_buffer=fault_buffer,
            fault_trace_anisotropy=fault_trace_anisotropy,
            fault_dip=fault_dip,
            fault_dip_anisotropy=fault_dip_anisotropy,
            fault_pitch=fault_pitch,
            **kwargs,
        )

    # TODO move rescale to bounding box/transformer
    @public_api(tier="stable")
    def rescale(self, points: np.ndarray, *, inplace: bool = False) -> np.ndarray:
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

        return self.bounding_box.reproject(points, inplace=inplace)

    # TODO move scale to bounding box/transformer
    @public_api(tier="stable")
    def scale(self, points: np.ndarray, *, inplace: bool = False) -> np.ndarray:
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
        return self.bounding_box.project(np.array(points).astype(float), inplace=inplace)

    @public_api(tier="stable")
    def regular_grid(self, *, nsteps=None, shuffle=True, rescale=False, order="C"):
        """
        Return a regular grid within the model bounding box

        Parameters
        ----------
        nsteps : tuple
            number of cells in x,y,z

        Returns
        -------
        xyz : np.array((N,3),dtype=float)
            locations of points in regular grid, in world coordinates
        """
        return self.bounding_box.regular_grid(
            nsteps=nsteps, shuffle=shuffle, order=order, local=False
        )

    @public_api(tier="stable")
    def evaluate_model(self, xyz: np.ndarray, *, scale: bool = True) -> np.ndarray:
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
        # `scale` is retained for API-signature compatibility only: features
        # now project world -> local coordinates internally (via the
        # interpolator's bounding_box), so xyz is always treated as world
        # coordinates here.
        xyz = np.array(xyz)
        strat_id = np.zeros(xyz.shape[0], dtype=int)
        # set strat id to -1 to identify which areas of the model aren't covered
        strat_id[:] = -1
        if self.stratigraphic_column is None:
            logger.warning("No stratigraphic column defined")
            return strat_id

        s_id = 0
        for g in reversed(self.stratigraphic_column.get_groups()):
            feature_id = self.feature_name_index.get(g.name, -1)
            if feature_id >= 0:
                vals = self.features[feature_id].evaluate_value(xyz)
                for u in g.units:
                    strat_id[np.logical_and(vals < u.max(), vals > u.min())] = s_id
                    s_id += 1
            if feature_id == -1:
                logger.error(f"Model does not contain {g.name}")

        return strat_id

    @public_api(tier="stable")
    def evaluate_model_gradient(self, points: np.ndarray, *, scale: bool = True) -> np.ndarray:
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
        # `scale` is retained for API-signature compatibility only -- see
        # evaluate_model.
        xyz = np.array(points)
        grad = np.zeros(xyz.shape)
        for g in reversed(self.stratigraphic_column.get_groups()):
            feature_id = self.feature_name_index.get(g.name, -1)
            if feature_id >= 0:
                gradt = self.features[feature_id].evaluate_gradient(xyz)
                grad[~np.isnan(gradt).any(axis=1)] = gradt[~np.isnan(gradt).any(axis=1)]
            if feature_id == -1:
                logger.error(f"Model does not contain {g.name}")

        return grad

    @public_api(tier="stable")
    def evaluate_fault_displacements(self, points, scale=True):
        """Evaluate the fault displacement magnitude at each location


        Parameters
        ----------
        points : np.array((N,3),dtype=float)
            locations
        scale : bool
            whether to rescale the points before evaluating model

        Returns
        -------
        fault_displacement : np.array(N,dtype=float)
            the fault displacement magnitude
        """
        # `scale` is retained for API-signature compatibility only -- see
        # evaluate_model.
        vals = np.zeros(points.shape[0])
        for f in self.features:
            if f.type == FeatureType.FAULT:
                disp = f.displacementfeature.evaluate_value(points)
                vals[~np.isnan(disp)] += disp[~np.isnan(disp)]
        return vals  # convert from restoration magnutude to displacement

    @public_api(tier="stable")
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

    @public_api(tier="stable")
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
        # `scale` is retained for API-signature compatibility only -- see
        # evaluate_model.
        feature = self.get_feature_by_name(feature_name)
        if feature:
            return feature.evaluate_value(xyz)
        else:
            return np.zeros(xyz.shape[0])

    @public_api(tier="stable")
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
        # `scale` is retained for API-signature compatibility only -- see
        # evaluate_model.
        feature = self.get_feature_by_name(feature_name)
        if feature:
            return feature.evaluate_gradient(xyz)
        else:
            return np.zeros(xyz.shape[0])

    @public_api(tier="stable")
    def update(self, verbose=False, progressbar=True):
        total_dof = 0
        nfeatures = 0
        for f in self.features:
            if f.type == FeatureType.FAULT:
                nfeatures += 3
                total_dof += f[0].interpolator.dof * 3
                continue
            if isinstance(f, StructuralFrame):
                nfeatures += 3
                total_dof += f[0].interpolator.dof * 3
                continue
            if f.type == FeatureType.INTERPOLATED:
                nfeatures += 1
                total_dof += f.interpolator.dof
                continue
        if verbose:
            logger.info(
                f"Updating geological model. There are: \n {nfeatures} \
                    geological features that need to be interpolated\n"
            )

        with timed_stage(logger, "update", nfeatures=nfeatures, total_dof=total_dof):
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

    @public_api(tier="stable")
    def stratigraphic_ids(self):
        """Return a list of all stratigraphic ids in the model

        Returns
        -------
        ids : list
            list of unique stratigraphic ids, featurename, unit name and min and max scalar values
        """
        return self.stratigraphic_column.get_stratigraphic_ids()

    @public_api(tier="stable")
    def get_fault_surfaces(self, faults: list[str] | None = None):
        return ModelExporter.get_fault_surfaces(self, faults=faults)

    @public_api(tier="stable")
    def get_stratigraphic_surfaces(self, units: list[str] | None = None, bottoms: bool = True):
        return ModelExporter.get_stratigraphic_surfaces(self, units=units, bottoms=bottoms)

    @public_api(tier="stable")
    def get_block_model(self, name='block model'):
        return ModelExporter.get_block_model(self, name=name)

    @public_api(tier="stable")
    def save(
        self,
        filename: str,
        block_model: bool = True,
        stratigraphic_surfaces=True,
        fault_surfaces=True,
        stratigraphic_data=True,
        fault_data=True,
    ):
        ModelExporter.save(
            self,
            filename,
            block_model=block_model,
            stratigraphic_surfaces=stratigraphic_surfaces,
            fault_surfaces=fault_surfaces,
            stratigraphic_data=stratigraphic_data,
            fault_data=fault_data,
        )


# Wire the built-in feature types up to GeologicalModel.create_and_add_feature
# (see FeatureBuilderRegistry / API.md). Each factory calls the corresponding
# ModelFeatureFactory.build_* staticmethod; new feature types register here
# without modifying GeologicalModel's source.
FeatureBuilderRegistry.register(
    "foliation", lambda model, name, **params: ModelFeatureFactory.build_foliation(model, name, **params)
)
FeatureBuilderRegistry.register(
    "fold_frame",
    lambda model, name, **params: ModelFeatureFactory.build_fold_frame(model, name, **params),
)
FeatureBuilderRegistry.register(
    "folded_foliation",
    lambda model, name, **params: ModelFeatureFactory.build_folded_foliation(model, name, **params),
)
FeatureBuilderRegistry.register(
    "folded_fold_frame",
    lambda model, name, **params: ModelFeatureFactory.build_folded_fold_frame(model, name, **params),
)
FeatureBuilderRegistry.register(
    "intrusion",
    lambda model, name, **params: ModelFeatureFactory.build_intrusion(
        model, name, params.pop("intrusion_frame_name"), **params
    ),
)
FeatureBuilderRegistry.register(
    "domain_fault",
    lambda model, name, **params: ModelFeatureFactory.build_domain_fault(model, name, **params),
)
FeatureBuilderRegistry.register(
    "fault",
    lambda model, name, **params: ModelFeatureFactory.build_fault(
        model, name, params.pop("displacement"), **params
    ),
)
