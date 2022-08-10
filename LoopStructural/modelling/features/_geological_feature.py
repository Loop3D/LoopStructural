"""
Geological features
"""
from LoopStructural.modelling.features import BaseFeature
from LoopStructural.utils import getLogger
from LoopStructural.modelling.features import FeatureType
import numpy as np

from LoopStructural.utils import getLogger, LoopValueError

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
        name,
        interpolator,
        builder=None,
        regions=[],
        faults=[],
        model=None,
    ):
        """Default constructor for geological feature

        Parameters
        ----------
        name: string
        interpolator : GeologicalInterpolator
        builder : GeologicalFeatureBuilder
        data :
        region :
        type :
        faults : list


        """
        BaseFeature.__init__(self, name, model, faults, regions, builder)
        self.name = name
        self.interpolator = interpolator
        self.builder = builder
        self.type = FeatureType.INTERPOLATED

    def is_valid(self):
        return self.interpolator.valid

    def __getitem__(self, key):
        return self._attributes[key]

    def __setitem__(self, key, item):
        self._attributes[key] = item

    def set_model(self, model):
        self.model = model

    def evaluate_value(self, evaluation_points):
        """
        Evaluate the scalar field value of the geological feature at the locations
        specified

        Parameters
        ----------
        evaluation_points : numpy array
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
        self.builder.up_to_date()
        # check if the points are within the display region
        v = np.zeros(evaluation_points.shape[0])
        v[:] = np.nan
        mask = np.zeros(evaluation_points.shape[0]).astype(bool)
        mask[:] = True
        # check regions
        for r in self.regions:
            # try:
            mask = np.logical_and(mask, r(evaluation_points))
            # except:
            #     logger.error("nan slicing")
        # apply faulting after working out which regions are visible
        if self.faults_enabled:
            for f in self.faults:
                evaluation_points = f.apply_to_points(evaluation_points)
        if mask.dtype not in [int, bool]:
            logger.error(f"Unable to evaluate value for {self.name}")
        else:
            v[mask] = self.interpolator.evaluate_value(evaluation_points[mask, :])
        return v

    def evaluate_gradient(self, evaluation_points):
        """

        Parameters
        ----------
        locations : numpy array
            location where the gradient is being evaluated

        Returns
        -------

        """
        if evaluation_points.shape[1] != 3:
            raise LoopValueError("Need Nx3 array of xyz points to evaluate gradient")
        self.builder.up_to_date()
        v = np.zeros(evaluation_points.shape)
        v[:] = np.nan
        mask = np.zeros(evaluation_points.shape[0]).astype(bool)
        mask[:] = True
        # check regions
        for r in self.regions:
            try:
                mask = np.logical_and(mask, r(evaluation_points))
            except:
                logger.error("nan slicing caught")
        # apply faulting after working out which regions are visible
        if self.faults_enabled:
            for f in self.faults:
                evaluation_points = f.apply_to_points(evaluation_points)
        if mask.dtype not in [int, bool]:
            logger.error(f"Unable to evaluate gradient for {self.name}")
        else:
            v[mask, :] = self.interpolator.evaluate_gradient(evaluation_points[mask, :])

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

    def mean(self):
        """
        Calculate average of the support values

        Returns
        -------
        mean : double
            average value of the feature evaluated on a (10,10,10) grid in model area

        """
        if self.model is None:
            return 0
        return np.mean(self.evaluate_value(self.model.regular_grid((10, 10, 10))))

    def min(self):
        """

        Returns
        -------
        min : double
            min value of the feature evaluated on a (10,10,10) grid in model area
        """
        if self.model is None:
            return 0
        return np.nanmin(self.evaluate_value(self.model.regular_grid((10, 10, 10))))

    def max(self):
        """
        Calculate average of the support values

        Returns
        -------
        max : double
            max value of the feature evaluated on a (10,10,10) grid in model area
        """
        if self.model is None:
            return 0
        return np.nanmax(self.evaluate_value(self.model.regular_grid((10, 10, 10))))
