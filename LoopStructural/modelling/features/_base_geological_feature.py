from LoopStructural.modelling.features import FeatureType
from LoopStructural.utils import getLogger

# from LoopStructural import GeologicalModel
import numpy as np

logger = getLogger(__name__)


class BaseFeature:
    """
    Base class for geological features.
    """

    def __init__(self, name, model, faults, regions, builder):
        self.name = name
        self.type = FeatureType.BASE
        self.regions = regions
        self.faults = faults
        self._model = model
        self.builder = builder
        self.faults_enabled = True

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
        # causes circular import, could delay import?
        # if type(model) == GeologicalModel:
        self._model = model

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
        return self.evaluate_value(xyz)

    def evaluate_value(self, pos):
        """
        Evaluate the feature at a given position.
        """
        raise NotImplementedError

    def evaluate_gradient(self, pos):
        """
        Evaluate the gradient of the feature at a given position.
        """
        raise NotImplementedError

    def min(self):
        if self.model is None:
            return 0
        return np.nanmin(self.evaluate_value(self.model.regular_grid((10, 10, 10))))

    def max(self):
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
