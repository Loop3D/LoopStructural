"""Geological features for LoopStructural modeling.

This module contains classes for representing geometrical elements in geological
models such as foliations, fault planes, and fold rotation angles.
"""

from __future__ import annotations

import dataclasses
import warnings

import numpy as np
from loop_interpolation import DiscreteInterpolator
from loop_interpolation._diagnostics import (
    ConstraintDiagnosticsReport,
    ResidualSummary,
    build_exclusion_records,
    direction_cosines,
    summarise_residuals,
)

from LoopStructural.utils.maths import gradient_from_tetrahedron, regular_tetraherdron_for_points

from ...geometry import ValuePoints, VectorPoints
from ...utils import LoopValueError, getLogger
from . import FeatureType
from ._base_geological_feature import BaseFeature

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
        regions: list | None = None,
        faults: list | None = None,
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
        if regions is None:
            regions = []
        if faults is None:
            faults = []
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
        logger.debug("%s %s", self.name, json)
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
        pos : np.ndarray
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
            from scipy import spatial

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

        if self.builder is not None:
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

            while not resolved:
                resolved = True
                for f in self.faults:
                    side = (
                        f[0]
                        .evaluate_value(tetrahedron.reshape(-1, 3), fillnan='nearest')
                        .reshape(tetrahedron.shape[0], 4)
                    )
                    flag = np.logical_or(np.all(side > 0, axis=1), np.all(side < 0, axis=1))
                    if np.any(~flag):
                        logger.warning(
                            f"Points are too close to fault {f[0].name}. Refining the tetrahedron"
                        )
                        element_scale_parameter *= 0.5
                        tetrahedron = regular_tetraherdron_for_points(pos, element_scale_parameter)
                        resolved = False

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

    @staticmethod
    def _first_stale(builder, label: str) -> str | None:
        """Find the first stale builder or interpolator below a builder.

        The method only reads flags. It follows the rule of
        ``BaseBuilder.up_to_date``. A structural frame builder or a fault
        builder has no flag of its own, so the method reads its component
        builders. It also reads the faults of each builder.

        Parameters
        ----------
        builder : object or None
            Builder to check.
        label : str
            Name of the builder in the text of the result.

        Returns
        -------
        str or None
            Text that names the first stale object, or None if none is stale.
        """
        if builder is None:
            return None
        for fault in getattr(builder, "faults", None) or []:
            fault_builder = getattr(fault, "builder", fault)
            fault_name = getattr(fault, "name", "a fault")
            reason = GeologicalFeature._first_stale(fault_builder, f"the fault {fault_name}")
            if reason is not None:
                return reason
        components = getattr(builder, "builders", None)
        if components is not None:
            for index, component in enumerate(components):
                name = getattr(component, "name", index)
                reason = GeologicalFeature._first_stale(component, f"{label} ({name})")
                if reason is not None:
                    return reason
            return None
        if getattr(builder, "_up_to_date", True) is False:
            return f"The builder of {label} is stale. Build the model first."
        interpolator = getattr(builder, "_interpolator", None)
        if getattr(interpolator, "up_to_date", True) is False:
            return f"The interpolator of {label} is stale. Build the model first."
        return None

    def _stale_reason(self) -> str | None:
        """Find the first stale object of the model.

        The method only reads flags. It does not build the model.

        Returns
        -------
        str or None
            Text that names the first stale object and says to build the model
            first, or None if no object is stale.
        """
        reason = self._first_stale(self.builder, self.name)
        if reason is not None:
            return reason
        for fault in getattr(self, "faults", None) or []:
            reason = self._first_stale(
                getattr(fault, "builder", None), f"the fault {getattr(fault, 'name', 'a fault')}"
            )
            if reason is not None:
                return reason
        if getattr(self.interpolator, "up_to_date", True) is False:
            return f"The interpolator of {self.name} is stale. Build the model first."
        return None

    def _constraint_row_counts(self) -> list:
        """Count the stored rows of the gradient, normal and value constraints.

        Returns
        -------
        list of int
            Row counts of the gradient, normal and value arrays.
        """
        return [
            int(np.shape(self.interpolator.get_gradient_constraints())[0]),
            int(np.shape(self.interpolator.get_norm_constraints())[0]),
            int(np.shape(self.interpolator.get_value_constraints())[0]),
        ]

    def _gradient_cosines_fresh(self) -> list:
        """Compute the gradient and normal cosines without a stale check.

        The method copies the stored arrays. It does not change them.
        The stored positions are already restored by the faults of the feature.
        The method evaluates the interpolator, not the feature, so no fault
        acts a second time.

        Returns
        -------
        list of numpy.ndarray
            One array of cosines for each non-empty family (gradient first,
            then normal).
        """
        bounding_box = self.interpolator.bounding_box
        cosines = []
        for array in (
            self.interpolator.get_gradient_constraints(),
            self.interpolator.get_norm_constraints(),
        ):
            if np.shape(array)[0] == 0:
                continue
            rows = np.array(array, dtype=float, copy=True)
            positions = rows[:, :3]
            data = rows[:, 3:6]
            if bounding_box is not None:
                positions = bounding_box.reproject(positions)
                data = self.interpolator.local_to_world_vectors(data)
            model = self.interpolator.evaluate_gradient(positions)
            cosines.append(direction_cosines(model, data))
        return cosines

    def _value_residuals_fresh(self) -> np.ndarray:
        """Compute the value residuals without a stale check.

        The stored positions are already restored by the faults of the feature.
        The method evaluates the interpolator, not the feature, so no fault
        acts a second time.

        Returns
        -------
        numpy.ndarray
            Absolute difference in field units, shape (N,).
        """
        rows = np.array(self.interpolator.get_value_constraints(), dtype=float, copy=True)
        positions = rows[:, :3]
        if self.interpolator.bounding_box is not None:
            positions = self.interpolator.bounding_box.reproject(positions)
        return np.abs(rows[:, 3] - self.interpolator.evaluate_value(positions))

    def _physical_residual_summaries(self) -> dict[str, ResidualSummary]:
        """Summarise the physical residuals of the value, gradient and normal data.

        The method does not check if the model is stale. The caller must check.

        Returns
        -------
        dict of str to ResidualSummary
            One summary for each non-empty family. The value summary uses field
            units. The gradient and normal summaries use degrees.
        """
        summaries: dict[str, ResidualSummary] = {}
        counts = self._constraint_row_counts()
        if counts[2] > 0:
            summaries["value"] = summarise_residuals(
                "value", "physical", "field", self._value_residuals_fresh()
            )
        # The cosine list has one array for each non-empty family, in this order.
        names = [name for name, count in zip(("gradient", "normal"), counts[:2]) if count > 0]
        for name, cosine in zip(names, self._gradient_cosines_fresh()):
            angles = np.degrees(np.arccos(np.clip(cosine, -1, 1)))
            summaries[name] = summarise_residuals(name, "physical", "degree", angles)
        return summaries

    def verification_report(self) -> ConstraintDiagnosticsReport:
        """Build a report of the constraints, residuals and exclusions.

        The method returns a new object. It does not store the report, build the
        model or change the state of the model.

        Returns
        -------
        ConstraintDiagnosticsReport
            The report. The residual blocks are empty if the model is not
            solved. On a stale model the method warns, leaves the residual
            blocks empty and still gives the exclusions.
        """
        base = self.interpolator._build_constraint_diagnostics_report()
        physical: dict[str, ResidualSummary] = {}
        assembled: dict[str, ResidualSummary] = {}
        reason = self._stale_reason()
        if reason is not None:
            warnings.warn(reason, UserWarning, stacklevel=2)
        elif getattr(self.interpolator, "c", None) is not None:
            physical = self._physical_residual_summaries()
            if hasattr(self.interpolator, "calculate_residual_for_constraints"):
                residuals = self.interpolator.calculate_residual_for_constraints()
                for name, residual in residuals.items():
                    assembled[name] = summarise_residuals(name, "assembled", "assembled", residual)
        excluded = getattr(self.builder, "excluded_rows", {}) if self.builder is not None else {}
        rows_per_point = (
            {"gradient": 2} if isinstance(self.interpolator, DiscreteInterpolator) else None
        )
        return dataclasses.replace(
            base,
            physical_residuals=physical,
            assembled_residuals=assembled,
            exclusions=build_exclusion_records(excluded, base.families, rows_per_point),
        )

    def evaluate_gradient_cosines(self) -> np.ndarray:
        """Compute the cosine between the model gradient and the gradient data.

        The cosine uses world coordinates (see COMPAT.md). The gradient rows come
        first and the normal rows come second. The method does not change the
        model or the stored constraints.

        Returns
        -------
        numpy.ndarray
            Cosines, shape (N,). The array is empty if there are no rows. On a
            stale model the method warns and returns NaN for each row.
        """
        reason = self._stale_reason()
        if reason is not None:
            warnings.warn(reason, UserWarning, stacklevel=2)
            n = sum(self._constraint_row_counts()[:2])
            return np.full(n, np.nan)
        cosines = self._gradient_cosines_fresh()
        if len(cosines) == 0:
            return np.zeros(0)
        return np.concatenate(cosines)

    def evaluate_value_residuals(self) -> np.ndarray:
        """Compute the residual between the model field and the value data.

        The method does not change the model or the stored constraints.

        Returns
        -------
        numpy.ndarray
            Absolute difference in field units, shape (N,), in the order of the
            value rows. On a stale model the method warns and returns NaN for
            each row.
        """
        reason = self._stale_reason()
        if reason is not None:
            warnings.warn(reason, UserWarning, stacklevel=2)
            return np.full(self._constraint_row_counts()[2], np.nan)
        return self._value_residuals_fresh()

    @staticmethod
    def _misfit_shape(counts: list) -> tuple:
        """Give the shape of the old gradient misfit array for the row counts."""
        if len(counts) == 0:
            return (0,)
        if len(counts) == 2 and counts[0] == counts[1]:
            return (2, counts[0])
        return (1, sum(counts))

    def evaluate_gradient_misfit(self) -> np.ndarray:
        """Compute the cosine misfit of the gradient and normal constraints.

        The value is a cosine in world coordinates (see COMPAT.md).

        Returns
        -------
        numpy.ndarray
            Shape (1, N) for one family, (2, N) for two families of equal
            length, (1, N1 + N2) for two families of different length and (0,)
            with no rows. On a stale model the method warns and returns NaN in
            the same shape.
        """
        reason = self._stale_reason()
        if reason is not None:
            warnings.warn(reason, UserWarning, stacklevel=2)
            counts = [n for n in self._constraint_row_counts()[:2] if n > 0]
            return np.full(self._misfit_shape(counts), np.nan)
        cosines = self._gradient_cosines_fresh()
        if len(cosines) == 0:
            return np.zeros(0)
        if len(cosines) == 2 and len(cosines[0]) == len(cosines[1]):
            return np.vstack(cosines)
        return np.concatenate(cosines)[None, :]

    def evaluate_value_misfit(self) -> np.ndarray:
        """Compute the value misfit.

        Returns
        -------
        numpy.ndarray
            Absolute difference between the model field and the value data,
            shape (N,). The difference is divided by the range of the field
            only if the range is finite and larger than zero. On a stale model
            the method warns and returns NaN for each row.
        """
        reason = self._stale_reason()
        if reason is not None:
            warnings.warn(reason, UserWarning, stacklevel=2)
            return np.full(self._constraint_row_counts()[2], np.nan)
        diff = self._value_residuals_fresh()
        value_range = self.max() - self.min()
        if np.isfinite(value_range) and value_range > 0:
            diff = diff / value_range
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

    def get_data(self, value_map: dict | None = None) -> list[ValuePoints | VectorPoints]:
        """Return the data associated with this geological feature

        The locations and the vectors are in world coordinates, as the builder stores them.

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
                    locations=gradient_constraints[:, :3],
                    vectors=gradient_constraints[:, 3:6],
                    name=f'{self.name}+_gradient',
                )
            )
        if norm_constraints.shape[0] > 0:
            data.append(
                VectorPoints(
                    locations=norm_constraints[:, :3],
                    vectors=norm_constraints[:, 3:6],
                    name=f'{self.name}_norm',
                )
            )
        if value_constraints.shape[0] > 0:
            if value_map is not None:
                for name, v in value_map.items():
                    data.append(
                        ValuePoints(
                            locations=value_constraints[value_constraints == v, :3],
                            values=value_constraints[value_constraints == v, 3],
                            name=f"{name}_value",
                        )
                    )
            else:
                data.append(
                    ValuePoints(
                        locations=value_constraints[:, :3],
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
