"""Geological features for LoopStructural modeling.

This module contains classes for representing geometrical elements in geological
models such as foliations, fault planes, and fold rotation angles.
"""

from LoopStructural.utils.maths import regular_tetraherdron_for_points, gradient_from_tetrahedron
from ...modelling.features import BaseFeature
from ...utils import getLogger
from ...modelling.features import FeatureType
import numpy as np
from typing import Optional, List, Union
from ...datatypes import ValuePoints, VectorPoints

from ...utils import LoopValueError

logger = getLogger(__name__)


class GeologicalFeature(BaseFeature):
    """A geological feature representing a geometrical element in a geological model.
    
    This class provides the foundation for representing various geological structures
    such as foliations, fault planes, fold rotation angles, and other geometrical
    elements within a geological model.

    Parameters
    ----------
    name : str
        Unique name for the geological feature
    builder : GeologicalFeatureBuilder
        Builder object used to construct the feature
    regions : list, optional
        List of boolean functions defining where the feature is active, by default []
    faults : list, optional
        List of FaultSegments that affect this feature, by default []
    interpolator : GeologicalInterpolator, optional
        Interpolator for the feature. If None, uses builder's interpolator, by default None
    model : GeologicalModel, optional
        The geological model containing this feature, by default None

    Attributes
    ----------
    name : str
        Unique name for the geological feature
    builder : GeologicalFeatureBuilder
        Builder object used to construct the feature
    interpolator : GeologicalInterpolator
        Interpolator used to compute field values
    type : FeatureType
        Type of the feature (set to INTERPOLATED)
    """

    def __init__(
        self,
        name: str,
        builder,
        regions: list = [],
        faults: list = [],
        interpolator=None,
        model=None,
    ):
        """Initialize the geological feature.

        Parameters
        ----------
        name : str
            Unique name for the geological feature
        builder : GeologicalFeatureBuilder
            Builder object used to construct the feature
        regions : list, optional
            List of boolean functions defining where the feature is active, by default []
        faults : list, optional
            List of FaultSegments that affect this feature, by default []
        interpolator : GeologicalInterpolator, optional
            Interpolator for the feature. If None, uses builder's interpolator, by default None
        model : GeologicalModel, optional
            The geological model containing this feature, by default None
        """
        BaseFeature.__init__(self, name, model, faults, regions, builder)
        self.name = name
        self.builder = builder
        self.interpolator = self.builder.interpolator if self.builder is not None else interpolator
        self.type = FeatureType.INTERPOLATED

    def to_json(self):
        """Return a JSON representation of the geological feature.

        Returns
        -------
        dict
            Dictionary containing the feature's state suitable for JSON serialization,
            including interpolator configuration
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

    def evaluate_value(self, pos: np.ndarray, ignore_regions=False, fillnan=None) -> np.ndarray:
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
        if pos.shape[1] != 3:
            raise LoopValueError("Need Nx3 array of xyz points to evaluate value")
        # TODO need to add a generic type checker for all methods
        # if evaluation_points is not a numpy array try and convert
        # otherwise error
        evaluation_points = np.asarray(pos)
        # if there is a builder lets make sure that the feature is up to date
        if self.builder is not None:
            self.builder.up_to_date()
        # check if the points are within the display region
        v = np.zeros(evaluation_points.shape[0])
        v[:] = np.nan
        mask = self._calculate_mask(pos, ignore_regions=ignore_regions)
        evaluation_points = self._apply_faults(evaluation_points)
        if mask.dtype not in [int, bool]:
            logger.error(f"Unable to evaluate value for {self.name}")
        else:
            v[mask] = self.interpolator.evaluate_value(evaluation_points[mask, :])
        if fillnan == 'nearest':
            import scipy.spatial as spatial

            nanmask = np.isnan(v)
            tree = spatial.cKDTree(evaluation_points[~nanmask, :])
            _d, i = tree.query(evaluation_points[nanmask, :])
            v[nanmask] = v[~nanmask][i]
        return v

    def evaluate_gradient(
        self, pos: np.ndarray, ignore_regions=False, element_scale_parameter=None
    ) -> np.ndarray:
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
        logger.info(f'Calculating gradient for {self.name}')
        if element_scale_parameter is None:
            if self.model is not None:
                element_scale_parameter = np.min(self.model.bounding_box.step_vector) / 10
            else:
                element_scale_parameter = 1
        else:
            try:
                element_scale_parameter = float(element_scale_parameter)
            except ValueError:
                logger.error("element_scale_parameter must be a float")
                element_scale_parameter = 1

        self.builder.up_to_date()

        v = np.zeros(pos.shape)
        v[:] = np.nan
        mask = self._calculate_mask(pos, ignore_regions=ignore_regions)
        # evaluate the faults on the nodes of the faulted feature support
        # then evaluate the gradient at these points
        if len(self.faults) > 0:
            # generate a regular tetrahedron for each point
            # we will then move these points by the fault and then recalculate the gradient.
            # this should work...
            resolved = False
            tetrahedron = regular_tetraherdron_for_points(pos, element_scale_parameter)

            while resolved:
                for f in self.faults:
                    v = (
                        f[0]
                        .evaluate_value(tetrahedron.reshape(-1, 3), fillnan='nearest')
                        .reshape(tetrahedron.shape[0], 4)
                    )
                    flag = np.logical_or(np.all(v > 0, axis=1), np.all(v < 0, axis=1))
                    if np.any(~flag):
                        logger.warning(
                            f"Points are too close to fault {f[0].name}. Refining the tetrahedron"
                        )
                        element_scale_parameter *= 0.5
                        tetrahedron = regular_tetraherdron_for_points(pos, element_scale_parameter)

                resolved = True

            tetrahedron_faulted = self._apply_faults(np.array(tetrahedron.reshape(-1, 3))).reshape(
                tetrahedron.shape
            )

            values = self.interpolator.evaluate_value(tetrahedron_faulted.reshape(-1, 3)).reshape(
                (-1, 4)
            )
            v[mask, :] = gradient_from_tetrahedron(tetrahedron[mask, :, :], values[mask])

            return v
        pos = self._apply_faults(pos)
        if mask.dtype not in [int, bool]:
            logger.error(f"Unable to evaluate gradient for {self.name}")
        else:
            v[mask, :] = self.interpolator.evaluate_gradient(pos[mask, :])
        logger.info(f'Gradient calculated for {self.name}')
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
        )
        return feature

    def get_data(self, value_map: Optional[dict] = None) -> List[Union[ValuePoints, VectorPoints]]:
        """Return the data associated with this geological feature

        Parameters
        ----------
        value_map : Optional[dict], optional
            A dictionary to map scalar values to another property, by default None

        Returns
        -------
        List[Union[ValuePoints, VectorPoints]]
            A container of either ValuePoints or VectorPoints
        """

        if self.builder is None:
            return []
        value_constraints = self.builder.get_value_constraints()
        gradient_constraints = self.builder.get_gradient_constraints()
        norm_constraints = self.builder.get_norm_constraints()
        # inequality_pair_constraints = self.builder.get_inequality_pair_constraints()
        # inequality_constraints = self.builder.get_inequality_constraints()
        data = []
        if gradient_constraints.shape[0] > 0:

            data.append(
                VectorPoints(
                    locations=self.model.rescale(gradient_constraints[:, :3]),
                    vectors=gradient_constraints[:, 3:6],
                    name=f'{self.name}+_gradient',
                )
            )
        if norm_constraints.shape[0] > 0:
            data.append(
                VectorPoints(
                    locations=self.model.rescale(norm_constraints[:, :3]),
                    vectors=norm_constraints[:, 3:6],
                    name=f'{self.name}_norm',
                )
            )
        if value_constraints.shape[0] > 0:
            if value_map is not None:
                for name, v in value_map.items():
                    data.append(
                        ValuePoints(
                            locations=self.model.rescale(
                                value_constraints[value_constraints == v, :3]
                            ),
                            values=value_constraints[value_constraints == v, 3],
                            name=f"{name}_value",
                        )
                    )
            else:
                data.append(
                    ValuePoints(
                        locations=self.model.rescale(value_constraints[:, :3]),
                        values=value_constraints[:, 3],
                        name=f"{self.name}_value",
                    )
                )
        # if inequality_constraints.shape[0] > 0:

        #     data.append(
        #         ValuePoints(
        #             locations=self.model.rescale(
        #                 inequality_constraints[:, :3]
        #             ),
        #             values=value_constraints[:, 3],
        #             name=f"{name}_inequality",
        #             properties = {'l':inequality_constraints[:,3],'u':inequality_constraints[:,4]}

        #         )
        #     )
        return data
