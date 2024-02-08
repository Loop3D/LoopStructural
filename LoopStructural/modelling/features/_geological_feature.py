"""
Geological features
"""

from ...modelling.features import BaseFeature
from ...utils import getLogger
from ...modelling.features import FeatureType
from ...interpolators import GeologicalInterpolator
import numpy as np

from ...utils import LoopValueError, rotate

logger = getLogger(__name__)


class GeologicalFeature(BaseFeature):
    """
    Geological feature is class that is used to represent a geometrical element in a geological
    model. For example foliations, fault planes, fold rotation angles etc.

    Attributes
    ----------
    name : string
        should be a unique name for the geological feature
    support : a ScalarField
        holds the property values for the feature and links to the
        support geometry
    data : list
        list containing geological data
    region : list
        list of boolean functions defining whether the feature is
        active
    faults : list
        list of FaultSegments that affect this feature
    """

    def __init__(
        self,
        name: str,
        interpolator: GeologicalInterpolator,
        builder=None,
        regions: list = [],
        faults: list = [],
        model=None,
    ):
        """Default constructor for geological feature

        Parameters
        ----------
        name: str
        interpolator : GeologicalInterpolator
        builder : GeologicalFeatureBuilder
        region : list
        faults : list
        model : GeologicalModel


        """
        BaseFeature.__init__(self, name, model, faults, regions, builder)
        self.name = name
        self.interpolator = interpolator
        self.builder = builder
        self.type = FeatureType.INTERPOLATED

    def to_json(self):
        """
        Returns a json representation of the geological feature

        Returns
        -------
        json : dict
            json representation of the geological feature
        """
        json = super().to_json()
        print(self.name, json)
        json["interpolator"] = self.interpolator.to_json()
        return json

    def is_valid(self):
        return self.interpolator.valid

    def __getitem__(self, key):
        return self._attributes[key]

    def __setitem__(self, key, item):
        self._attributes[key] = item

    def set_model(self, model):
        self.model = model

    def evaluate_value(self, evaluation_points: np.ndarray) -> np.ndarray:
        """
        Evaluate the scalar field value of the geological feature at the locations
        specified

        Parameters
        ----------
        evaluation_points : np.ndarray
            location to evaluate the scalar value

        Returns
        -------
        values : numpy array
            numpy array containing evaluated values

        """
        if evaluation_points.shape[1] != 3:
            raise LoopValueError("Need Nx3 array of xyz points to evaluate value")
        # TODO need to add a generic type checker for all methods
        # if evaluation_points is not a numpy array try and convert
        # otherwise error
        evaluation_points = np.asarray(evaluation_points)
        self.builder.up_to_date()
        # check if the points are within the display region
        v = np.zeros(evaluation_points.shape[0])
        v[:] = np.nan
        mask = self._calculate_mask(evaluation_points)
        evaluation_points, axis, angle = self._apply_faults(evaluation_points)
        if mask.dtype not in [int, bool]:
            logger.error(f"Unable to evaluate value for {self.name}")
        else:
            v[mask] = self.interpolator.evaluate_value(evaluation_points[mask, :])
        return v

    def evaluate_gradient(self, pos: np.ndarray) -> np.ndarray:
        """

        Parameters
        ----------
        locations : numpy array
            location where the gradient is being evaluated

        Returns
        -------


        """
        if pos.shape[1] != 3:
            raise LoopValueError("Need Nx3 array of xyz points to evaluate gradient")
        self.builder.up_to_date()
        v = np.zeros(pos.shape)
        v[:] = np.nan
        mask = self._calculate_mask(pos)
        original_pos = pos.copy()
        pos, axis, angle = self._apply_faults(pos)
        if mask.dtype not in [int, bool]:
            logger.error(f"Unable to evaluate gradient for {self.name}")
        else:
            v[mask, :] = self.interpolator.evaluate_gradient(pos[mask, :])
            v_norm = np.linalg.norm(v[mask, :], axis=1)
            v[mask, :] /= v_norm[:, None]
            if axis is not None:
                for i in reversed(range(axis.shape[1])):
                    v[mask, :] = rotate(v[mask, :], axis[:, i, :], angle[:, i])
            v[mask, :] *= v_norm[:, None]
        return v

    def evaluate_gradient_misfit(self):
        """

        Returns
        -------
        misfit : np.array(N,dtype=double)
            dot product between interpolated gradient and constraints
        """
        self.builder.up_to_date()
        grad = self.interpolator.get_gradient_constraints()
        norm = self.interpolator.get_norm_constraints()

        dot = []
        if grad.shape[0] > 0:
            grad /= np.linalg.norm(grad, axis=1)[:, None]
            model_grad = self.evaluate_gradient(grad[:, :3])
            dot.append(np.einsum("ij,ij->i", model_grad, grad[:, :3:6]).tolist())

        if norm.shape[0] > 0:
            norm /= np.linalg.norm(norm, axis=1)[:, None]
            model_norm = self.evaluate_gradient(norm[:, :3])
            dot.append(np.einsum("ij,ij->i", model_norm, norm[:, :3:6]))

        return np.array(dot)

    def evaluate_value_misfit(self):
        """

        Returns
        -------
        misfit : np.array(N,dtype=double)
            difference between interpolated scalar field and value constraints
        """
        self.builder.up_to_date()

        locations = self.interpolator.get_value_constraints()
        diff = np.abs(locations[:, 3] - self.evaluate_value(locations[:, :3]))
        diff /= self.max() - self.min()
        return diff

    def copy(self, name=None):
        if not name:
            name = f"{self.name}_copy"
        feature = GeologicalFeature(
            name=name,
            faults=self.faults,
            regions=[],  # feature.regions.copy(),  # don't want to share regionsbetween unconformity and # feature.regions,
            builder=self.builder,
            model=self.model,
            interpolator=self.interpolator,
        )
        return feature
