from abc import ABCMeta, abstractmethod
from LoopStructural.modelling.features import FeatureType
from LoopStructural.utils import getLogger
from LoopStructural.utils.typing import NumericInput


import numpy as np

logger = getLogger(__name__)


class BaseFeature(metaclass=ABCMeta):
    """
    Base class for geological features.
    """

    def __init__(self, name: str, model=None, faults: list = [], regions: list = [], builder=None):
        """Base geological feature, this is a virtual class and should not be
        used directly. Inheret from this to implement a new type of geological
        feature or use one of the exisitng implementations

        Parameters
        ----------
        name :
            Name of the geological feature to add
        model : GeologicalModel, optional
            the model the feature is associated with, by default None
        faults : list, optional
            any faults that fault this feature, by default []
        regions : list, optional
            any regions that affect this feature, by default []
        builder : GeologicalFeatureBuilder, optional
            the builder of the feature, by default None
        """
        self.name = name
        self.type = FeatureType.BASE
        self.regions = regions
        self._faults = []
        if faults:
            self.faults = faults
        self._model = model
        self.builder = builder
        self.faults_enabled = True
        self._min = None
        self._max = None

    @property
    def faults(self):
        return self._faults

    @faults.setter
    def faults(self, faults: list):
        _faults = []
        try:
            for f in faults:
                if not issubclass(type(f), BaseFeature):
                    raise TypeError("Faults must be a list of BaseFeature")
                _faults.append(f)
        except TypeError:
            logger.error(
                f'Faults must be a list of BaseFeature \n Trying to set using {type(faults)}'
            )
            raise TypeError("Faults must be a list of BaseFeature")

        self._faults = _faults

    def to_json(self):
        """
        Returns a json representation of the geological feature

        Returns
        -------
        json : dict
            json representation of the geological feature
        """
        json = {}
        json["name"] = self.name
        json["regions"] = [r.to_json() for r in self.regions]
        json["faults"] = [f.name for f in self.faults]
        json["type"] = self.type
        return json

    def __str__(self):
        _str = "-----------------------------------------------------\n"
        _str += f"{self.name} {self.type} \n"
        _str += "-----------------------------------------------------\n"
        _str += f"\t{len(self.regions)} regions\n"
        for r in self.regions:
            _str += f"\t \t{r.__str__}\n"
        _str += f"\t{len(self.faults)} faults.\n"
        _str += f"\tFault enabled {self.faults_enabled}\n"

        for f in self.faults:
            _str += f"\t \t{f.__str__}\n"
        return _str

    def __repr__(self):
        return self.__str__()

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        from LoopStructural import GeologicalModel

        # causes circular import, could delay import?
        if type(model) == GeologicalModel:
            self._model = model
        elif not model:
            self._model = None
            logger.error("Model not set")
        else:
            raise TypeError("Model must be a GeologicalModel")

    def toggle_faults(self):
        """
        Turn the fault off for a feature
        This function is only really used for debugging or creating methods
        explanation figures

        Returns
        -------

        """
        logger.warning(f"Toggling faults for feature {self.name}")
        self.faults_enabled = not self.faults_enabled

    def add_region(self, region):
        """
        Adds a region where the geological feature is active to the model.

        Parameters
        ----------
        region : boolean function(x,y,z)
                returns true if inside region, false if outside
                can be passed as a lambda function e.g.
                lambda pos : feature.evaluate_value(pos) > 0

        Returns
        -------

        """
        self.regions.append(region)

    def __call__(self, xyz):
        """Calls evaluate_value method

        Parameters
        ----------
        xyz : np.ndarray
            location to evaluate feature

        Returns
        -------
        np.ndarray
            the value of the feature at the locations
        """
        return self.evaluate_value(xyz)

    @abstractmethod
    def evaluate_value(self, pos):
        """
        Evaluate the feature at a given position.
        """
        raise NotImplementedError

    def evaluate_normalised_value(self, pos: NumericInput):
        """Evaluate the feature value scaling between 0 and 1

        Parameters
        ----------
        pos : NumericInput
            An array or arraylike object with locations
        """
        value = self.evaluate_value(pos)
        return (value - self.min()) / (self.max() - self.min())

    def _calculate_mask(self, evaluation_points: np.ndarray) -> np.ndarray:
        """Calculate the mask for which evaluation points need to be calculated

        Parameters
        ----------
        evaluation_points : np.ndarray
            location to be evaluated, Nx3 array

        Returns
        -------
        np.ndarray
            bool mask Nx1 ndarray
        """
        mask = np.zeros(evaluation_points.shape[0]).astype(bool)

        mask[:] = True
        # check regions
        for r in self.regions:
            # try:
            mask = np.logical_and(mask, r(evaluation_points))
        return mask

    def _apply_faults(self, evaluation_points: np.ndarray, reverse: bool = False) -> np.ndarray:
        """Calculate the restored location of the points given any faults if faults are enabled

        Parameters
        ----------
        evaluation_points : np.ndarray
            location to be evaluated, Nx3 array

        Returns
        -------
        np.ndarray
            faulted value Nx1 ndarray
        """

        if self.faults_enabled:
            # check faults
            for f in self.faults:
                evaluation_points = f.apply_to_points(evaluation_points, reverse=reverse)
        return evaluation_points

    @abstractmethod
    def evaluate_gradient(self, pos):
        """
        Evaluate the gradient of the feature at a given position.
        """

        raise NotImplementedError

    def min(self):
        """Calculate the min value of the geological feature
        in the model

        Returns
        -------
        minimum, float
            min value of the feature evaluated on a regular grid in the model domain
        """
        if self.model is None:
            return 0

        return np.nanmin(self.evaluate_value(self.model.regular_grid((10, 10, 10))))

    def max(self):
        """Calculate the maximum value of the geological feature
        in the model

        Returns
        -------
        maximum, float
            max value of the feature evaluated on a regular grid in the model domain
        """
        if self.model is None:
            return 0
        return np.nanmax(self.evaluate_value(self.model.regular_grid((10, 10, 10))))

    def __tojson__(self):
        regions = [r.name for r in self.regions]
        faults = [f.name for f in self.faults]
        return {
            "name": self.name,
            "type": self.type,
            "regions": regions,
            "faults": faults,
        }
