"""Base geological interpolator for LoopStructural.

This module contains the abstract base class for all geological interpolators
used in LoopStructural geological modelling framework.
"""

from abc import ABCMeta, abstractmethod
from ._interpolatortype import InterpolatorType
import numpy as np
import json

from typing import Dict, Optional, Union
from loop_common.interfaces.representation import BaseRepresentation
from loop_common.logging import get_logger as getLogger
from ._diagnostics import (
    ConstraintDiagnosticsReport,
    ConstraintFamilyDiagnostics,
    RegionCoverageDiagnostics,
)
from .constraints import (
    ValueConstraint,
    GradientConstraint,
    InterfaceConstraint,
    InequalityConstraint,
    InequalityPair,
)
from ._validation import (
    check_unsupported_combinations,
    ValidationError,
)

logger = getLogger(__name__)


class LoopTypeError(TypeError):
    pass


class GeologicalInterpolator(BaseRepresentation, metaclass=ABCMeta):
    """Abstract base class for geological interpolators.

    This class defines the interface for all geological interpolators in
    LoopStructural, providing methods for setting constraints and evaluating
    the interpolated scalar field.

    Attributes
    ----------
    data : dict
        Dictionary containing numpy arrays for gradient, value, normal, and tangent data
    n_g : int
        Number of gradient constraints
    n_i : int
        Number of interface/value constraints
    n_n : int
        Number of normal constraints
    n_t : int
        Number of tangent constraints
    type : InterpolatorType
        The type of interpolator
    up_to_date : bool
        Whether the interpolator needs to be rebuilt
    constraints : list
        List of applied constraints
    valid : bool
        Whether the interpolator is in a valid state
    dimensions : int
        Number of spatial dimensions (default 3)
    support : object
        The support structure used by the interpolator
    """

    @abstractmethod
    def __init__(self, data=None, up_to_date=False):
        """Initialize the geological interpolator.

        This method sets up the basic data structures and parameters required
        for geological interpolation.

        Parameters
        ----------
        data : dict, optional
            Dictionary containing constraint data arrays, by default {}
        up_to_date : bool, optional
            Whether the interpolator is already built and up to date, by default False

        Notes
        -----
        This is an abstract method that must be implemented by subclasses.
        All subclasses should call this parent constructor to ensure proper
        initialization of the base data structures.
        """
        self._data = {}
        self.data = data  # None
        self.clean()  # init data structure

        self.n_g = 0
        self.n_i = 0
        self.n_n = 0
        self.n_t = 0

        self.type = InterpolatorType.BASE
        self.up_to_date = up_to_date
        self.constraints = []
        self.__str = "Base Geological Interpolator"
        self.valid = True
        self.dimensions = 3  # default to 3d
        self.support = None
        self.latest_diagnostics_report: Optional[ConstraintDiagnosticsReport] = None

    @abstractmethod
    def set_nelements(self, nelements: int) -> int:
        """Set the number of elements for the interpolation support.

        Parameters
        ----------
        nelements : int
            Target number of elements

        Returns
        -------
        int
            Actual number of elements set

        Notes
        -----
        This is an abstract method that must be implemented by subclasses.
        The actual number of elements may differ from the requested number
        depending on the interpolator's constraints.
        """
        pass

    @property
    @abstractmethod
    def n_elements(self) -> int:
        """Get the number of elements in the interpolation support.

        Returns
        -------
        int
            Number of elements

        Notes
        -----
        This is an abstract property that must be implemented by subclasses.
        """
        pass

    @property
    def data(self):
        """Get the constraint data dictionary.

        Returns
        -------
        dict
            Dictionary containing constraint data arrays
        """
        return self._data

    @data.setter
    def data(self, data):
        """Set the constraint data dictionary.

        Parameters
        ----------
        data : dict or None
            Dictionary containing constraint data arrays. If None, an empty dict is used.
        """
        if data is None:
            data = {}
        for k, v in data.items():
            self._data[k] = np.array(v)

    def __str__(self):
        """Return string representation of the interpolator.

        Returns
        -------
        str
            String describing the interpolator type and constraint counts
        """
        name = f"{self.type} \n"
        name += f"{self.n_g} gradient points\n"
        name += f"{self.n_i} interface points\n"
        name += f"{self.n_n} normal points\n"
        name += f"{self.n_t} tangent points\n"
        name += f"{self.n_g + self.n_i + self.n_n + self.n_t} total points\n"
        return name

    def check_array(self, array: np.ndarray):
        """Validate and convert input to numpy array.

        Parameters
        ----------
        array : array_like
            Input array to validate and convert

        Returns
        -------
        np.ndarray
            Validated numpy array

        Raises
        ------
        LoopTypeError
            If the array cannot be converted to a numpy array
        """
        try:
            return np.array(array)
        except Exception as e:
            raise LoopTypeError(str(e))

    def _coerce_value_constraint(
        self, points: Union[np.ndarray, ValueConstraint]
    ) -> ValueConstraint:
        if isinstance(points, ValueConstraint):
            return points
        return ValueConstraint.from_array(points, dimensions=self.dimensions)

    def _coerce_gradient_constraint(
        self, points: Union[np.ndarray, GradientConstraint], is_normal: bool = False
    ) -> GradientConstraint:
        if isinstance(points, GradientConstraint):
            return points
        return GradientConstraint.from_array(
            points, dimensions=self.dimensions, is_normal=is_normal
        )

    def _coerce_interface_constraint(
        self, points: Union[np.ndarray, InterfaceConstraint]
    ) -> InterfaceConstraint:
        if isinstance(points, InterfaceConstraint):
            return points
        return InterfaceConstraint.from_array(points, dimensions=self.dimensions)

    def _coerce_inequality_constraint(
        self, points: Union[np.ndarray, InequalityConstraint]
    ) -> InequalityConstraint:
        if isinstance(points, InequalityConstraint):
            return points
        return InequalityConstraint.from_array(points, dimensions=self.dimensions)

    def _coerce_inequality_pair_constraint(
        self, points: Union[np.ndarray, InequalityPair]
    ) -> InequalityPair:
        if isinstance(points, InequalityPair):
            return points
        return InequalityPair.from_array(points, dimensions=self.dimensions)

    @abstractmethod
    def set_region(self, **kwargs):
        """Set the interpolation region.

        Parameters
        ----------
        **kwargs : dict
            Region parameters specific to the interpolator implementation

        Notes
        -----
        This is an abstract method that must be implemented by subclasses.
        The specific parameters depend on the interpolator type.
        """
        pass

    def set_value_constraints(self, points: Union[np.ndarray, ValueConstraint]):
        """Set value constraints for the interpolation.

        Parameters
        ----------
        points : np.ndarray
            Array containing the value constraints with shape (n_points, 4-5).
            Columns should be [X, Y, Z, value, weight]. If weight is not provided,
            a weight of 1.0 is assumed for all points.

        Raises
        ------
        ValidationError
            If points array fails shape, dtype, finiteness, or logical checks

        Notes
        -----
        Value constraints specify known scalar field values at specific locations.
        These are typically used for interface points or measured data values.
        All values must be finite (not NaN or inf), and array must be convertible
        to float64.
        """
        try:
            check_unsupported_combinations(self.data, "value")
            points = self._coerce_value_constraint(points).to_array()
            self.data["value"] = points.copy()
            self.n_i = points.shape[0]
            self.up_to_date = False
        except ValidationError as e:
            raise ValidationError(f"Failed to set value constraints: {e}") from e

    def set_gradient_constraints(self, points: Union[np.ndarray, GradientConstraint]):
        """Set gradient constraints for the interpolation.

        Parameters
        ----------
        points : np.ndarray
            Array containing gradient constraints with shape (n_points, 7-8).
            Columns should be [X, Y, Z, gx, gy, gz, weight]. If weight is not
            provided, a weight of 1.0 is assumed for all points.

        Raises
        ------
        ValidationError
            If points array fails shape, dtype, finiteness, or logical checks.
            Also raised if gradient vectors have zero magnitude.

        Notes
        -----
        Gradient constraints specify the direction and magnitude of the scalar
        field gradient at specific locations. These are typically derived from
        structural measurements like bedding or foliation orientations.
        All values must be finite (not NaN or inf), and gradient vectors must
        have non-zero magnitude.
        """
        try:
            check_unsupported_combinations(self.data, "gradient")
            points = self._coerce_gradient_constraint(points).to_array()
            self.n_g = points.shape[0]
            self.data["gradient"] = points.copy()
            self.up_to_date = False
        except ValidationError as e:
            raise ValidationError(f"Failed to set gradient constraints: {e}") from e

    def set_normal_constraints(self, points: Union[np.ndarray, GradientConstraint]):
        """Set normal constraints for the interpolation.

        Parameters
        ----------
        points : np.ndarray
            Array containing normal constraints with shape (n_points, 7-8).
            Columns should be [X, Y, Z, nx, ny, nz, weight]. If weight is not
            provided, a weight of 1.0 is assumed for all points.

        Raises
        ------
        ValidationError
            If points array fails shape, dtype, finiteness, or logical checks.
            Also raised if normal vectors have zero magnitude.

        Notes
        -----
        Normal constraints specify surface normal directions at specific locations.
        All values must be finite (not NaN or inf), and normal vectors must
        have non-zero magnitude.
        If no weights are provided, w = 1 is assigned to each normal constraint.
        """
        try:
            check_unsupported_combinations(self.data, "normal")
            points = self._coerce_gradient_constraint(points, is_normal=True).to_array()
            self.n_n = points.shape[0]
            self.data["normal"] = points.copy()
            self.up_to_date = False
        except ValidationError as e:
            raise ValidationError(f"Failed to set normal constraints: {e}") from e

    def set_tangent_constraints(self, points: Union[np.ndarray, GradientConstraint]):
        """Set tangent constraints for the interpolation.

        Parameters
        ----------
        points : np.ndarray
            Array containing tangent constraints with shape (n_points, 7-8).
            Columns should be [X, Y, Z, tx, ty, tz, weight]. If weight is not
            provided, a weight of 1.0 is assumed for all points.

        Raises
        ------
        ValidationError
            If points array fails shape, dtype, finiteness, or logical checks.
            Also raised if tangent vectors have zero magnitude.

        Notes
        -----
        Tangent constraints specify tangent directions at specific locations.
        All values must be finite (not NaN or inf), and tangent vectors must
        have non-zero magnitude. If no weights are provided, w = 1 is assigned
        to each tangent constraint.
        """
        try:
            check_unsupported_combinations(self.data, "tangent")
            points = self._coerce_gradient_constraint(points).to_array()
            self.n_t = points.shape[0]
            self.data["tangent"] = points.copy()
            self.up_to_date = False
        except ValidationError as e:
            raise ValidationError(f"Failed to set tangent constraints: {e}") from e

    def set_interface_constraints(self, points: Union[np.ndarray, InterfaceConstraint]):
        """Set interface constraints for the interpolation.

        Parameters
        ----------
        points : np.ndarray
            Array containing interface constraints with shape (n_points, 4-5).
            Columns should be [X, Y, Z, interface_id, weight]. If weight is not
            provided, a weight of 1.0 is assumed for all points.

        Raises
        ------
        ValidationError
            If points array fails shape, dtype, or finiteness checks.

        Notes
        -----
        Interface constraints mark surface boundaries between different units.
        All values must be finite (not NaN or inf).
        """
        try:
            check_unsupported_combinations(self.data, "interface")
            points = self._coerce_interface_constraint(points).to_array()
            self.data["interface"] = points.copy()
            self.up_to_date = False
        except ValidationError as e:
            raise ValidationError(f"Failed to set interface constraints: {e}") from e

    def set_value_inequality_constraints(self, points: Union[np.ndarray, InequalityConstraint]):
        """Set inequality value constraints for the interpolation.

        Parameters
        ----------
        points : np.ndarray
            Array containing inequality constraints with shape (n_points, 5-6).
            Columns should be [X, Y, Z, lower_bound, upper_bound, weight].
            If weight is not provided, a weight of 1.0 is assumed for all points.

        Raises
        ------
        ValidationError
            If points array fails shape, dtype, finiteness, or logical checks.
            Also raised if lower_bound >= upper_bound for any constraint.

        Notes
        -----
        Inequality constraints specify bounds on scalar field values at points.
        For each constraint, lower_bound must be strictly less than upper_bound.
        All values must be finite (not NaN or inf).
        """
        try:
            check_unsupported_combinations(self.data, "inequality")
            points = self._coerce_inequality_constraint(points).to_array()
            self.data["inequality"] = points.copy()
            self.up_to_date = False
        except ValidationError as e:
            raise ValidationError(f"Failed to set inequality value constraints: {e}") from e

    def set_inequality_pairs_constraints(self, points: Union[np.ndarray, InequalityPair]):
        """Set inequality pairs constraints for the interpolation.

        Parameters
        ----------
        points : np.ndarray
            Array containing inequality pairs constraints with shape (n_points, 4-5).
            Columns should be [X, Y, Z, rock_id, weight]. If weight is not
            provided, a weight of 1.0 is assumed for all points.

        Raises
        ------
        ValidationError
            If points array fails shape, dtype, or finiteness checks.

        Notes
        -----
        Inequality pairs constraints enforce ordering relationships between pairs
        of points. All values must be finite (not NaN or inf).
        """
        try:
            check_unsupported_combinations(self.data, "inequality_pairs")
            points = self._coerce_inequality_pair_constraint(points).to_array()
            self.data["inequality_pairs"] = points.copy()
            self.up_to_date = False
        except ValidationError as e:
            raise ValidationError(f"Failed to set inequality pairs constraints: {e}") from e

    def get_value_constraints(self):
        """

        Returns
        -------
        numpy array
        """
        return self.data["value"]

    def get_gradient_constraints(self):
        """

        Returns
        -------
        numpy array
        """
        return self.data["gradient"]

    def get_tangent_constraints(self):
        """

        Returns
        -------
        numpy array
        """

        return self.data["tangent"]

    def get_norm_constraints(self):
        """

        Returns
        -------
        numpy array
        """
        return self.data["normal"]

    def get_data_locations(self):
        """Get the location of all data points

        Returns
        -------
        numpy array
            Nx3 - X,Y,Z location of all data points
        """
        return np.vstack([d[:, :3] for d in self.data.values()])

    def get_interface_constraints(self):
        """Get the location of interface constraints

        Returns
        -------
        numpy array
            Nx4 - X,Y,Z,id location of all interface constraints
        """
        return self.data["interface"]

    def get_inequality_value_constraints(self):
        return self.data["inequality"]

    def get_inequality_pairs_constraints(self):
        return self.data["inequality_pairs"]

    def _outside_model_points_from_data(self) -> Dict[str, int]:
        if self.support is None or not hasattr(self.support, "inside"):
            return {}

        mapping = {
            "value": "value",
            "gradient": "gradient",
            "normal": "normal",
            "tangent": "tangent",
            "interface": "interface",
            "inequality": "inequality_value",
            "inequality_pairs": "inequality_pairs",
        }

        outside = {}
        for data_key, family_name in mapping.items():
            points = self.data.get(data_key)
            if points is None or points.shape[0] == 0:
                outside[family_name] = 0
                continue
            xyz = np.asarray(points[:, : self.dimensions], dtype=float)
            try:
                inside = self.support.inside(xyz)
                outside[family_name] = int((~inside).sum())
            except (ValueError, TypeError) as e:
                logger.warning(
                    f"Failed to compute outside-model coverage for constraint family '{family_name}': {e}. "
                    f"Using zero outside points."
                )
                outside[family_name] = 0
        return outside

    def _weight_stats_from_points(self, points: np.ndarray):
        if points.shape[1] <= self.dimensions + 1:
            return None, None, None
        weights = np.asarray(points[:, -1], dtype=float)
        if weights.size == 0:
            return None, None, None
        return (
            float(np.mean(weights)),
            float(np.min(weights)),
            float(np.max(weights)),
        )

    def _build_constraint_diagnostics_report(self) -> ConstraintDiagnosticsReport:
        mapping = {
            "value": "value",
            "gradient": "gradient",
            "normal": "normal",
            "tangent": "tangent",
            "interface": "interface",
            "inequality": "inequality_value",
            "inequality_pairs": "inequality_pairs",
        }

        outside = self._outside_model_points_from_data()
        families = {}
        for data_key, family_name in mapping.items():
            points = self.data.get(data_key)
            if points is None:
                points = np.zeros((0, self.dimensions + 1), dtype=float)
            row_count = int(points.shape[0])
            mean_w, min_w, max_w = self._weight_stats_from_points(points)
            families[family_name] = ConstraintFamilyDiagnostics(
                name=family_name,
                active=row_count > 0,
                row_count=row_count,
                dropped_rows=outside.get(family_name, 0),
                effective_weight_mean=mean_w,
                effective_weight_min=min_w,
                effective_weight_max=max_w,
                source_point_count=row_count,
                outside_model_point_count=outside.get(family_name, 0),
            )

        region_coverage = None
        if self.support is not None and hasattr(self.support, "n_nodes"):
            total_nodes = int(self.support.n_nodes)
            region_mask = np.ones(total_nodes, dtype=bool)
            if hasattr(self, "region"):
                try:
                    region_mask = np.asarray(self.region, dtype=bool)
                except (ValueError, TypeError) as e:
                    logger.warning(
                        f"Failed to interpret region mask (expected array-like bool): {e}. "
                        f"Using full domain (all nodes active)."
                    )
                    region_mask = np.ones(total_nodes, dtype=bool)
            active_nodes = int(np.sum(region_mask))
            region_coverage = RegionCoverageDiagnostics(
                total_support_nodes=total_nodes,
                active_region_nodes=active_nodes,
                inactive_region_nodes=total_nodes - active_nodes,
                active_fraction=0.0 if total_nodes == 0 else float(active_nodes / total_nodes),
            )

        return ConstraintDiagnosticsReport(
            interpolator_type=self.type.name,
            families=families,
            region_coverage=region_coverage,
            outside_model_points=outside,
        )

    def get_constraint_diagnostics_report(
        self, refresh: bool = False
    ) -> ConstraintDiagnosticsReport:
        if self.latest_diagnostics_report is None or refresh:
            self.latest_diagnostics_report = self._build_constraint_diagnostics_report()
        return self.latest_diagnostics_report

    # @abstractmethod
    def setup(self, **kwargs) -> ConstraintDiagnosticsReport:
        """Run setup and return a diagnostics report."""
        report = self.setup_interpolator(**kwargs)
        if isinstance(report, ConstraintDiagnosticsReport):
            self.latest_diagnostics_report = report
            return report
        return self.get_constraint_diagnostics_report(refresh=True)

    @abstractmethod
    def setup_interpolator(self, **kwargs):
        """
        Runs all of the required setting up stuff
        """
        raise NotImplementedError("setup_interpolator must be implemented by subclasses")

    @abstractmethod
    def solve_system(self, solver, solver_kwargs: Optional[dict] = None) -> bool:
        """
        Solves the interpolation equations
        """
        pass

    @abstractmethod
    def update(self) -> bool:
        return False

    @abstractmethod
    def evaluate_value(self, locations: np.ndarray):
        raise NotImplementedError("evaluate_value not implemented")

    @abstractmethod
    def evaluate_gradient(self, locations: np.ndarray):
        raise NotImplementedError("evaluate_gradient not implemented")

    def surfaces(self, value):
        raise NotImplementedError("Surface extraction not implemented for this representation")

    @abstractmethod
    def reset(self):
        pass

    @abstractmethod
    def add_value_constraints(self, w: float = 1.0):
        pass

    @abstractmethod
    def add_gradient_constraints(self, w: float = 1.0):
        pass

    @abstractmethod
    def add_norm_constraints(self, w: float = 1.0):
        pass

    @abstractmethod
    def add_tangent_constraints(self, w: float = 1.0):
        pass

    @abstractmethod
    def add_interface_constraints(self, w: float = 1.0):
        pass

    @abstractmethod
    def add_value_inequality_constraints(self, w: float = 1.0):
        pass

    @abstractmethod
    def add_inequality_pairs_constraints(
        self,
        w: float = 1.0,
        upper_bound=np.finfo(float).eps,
        lower_bound=-np.inf,
        pairs: Optional[list] = None,
    ):
        pass

    def to_dict(self):
        def _to_json_safe(value):
            if isinstance(value, np.ndarray):
                return value.tolist()
            if isinstance(value, (np.floating, np.integer)):
                return value.item()
            if isinstance(value, dict):
                return {k: _to_json_safe(v) for k, v in value.items()}
            if isinstance(value, (list, tuple)):
                return [_to_json_safe(v) for v in value]
            return value

        payload = {
            "type": self.type.value,
            "data": _to_json_safe({k: np.asarray(v) for k, v in self.data.items()}),
            "up_to_date": self.up_to_date,
            "valid": self.valid,
        }
        if self.support is not None and hasattr(self.support, "to_dict"):
            support_dict = _to_json_safe(self.support.to_dict())
            if (
                isinstance(support_dict, dict)
                and "nsteps" in support_dict
                and hasattr(self.support, "nsteps_cells")
            ):
                support_dict["nsteps"] = _to_json_safe(np.asarray(self.support.nsteps_cells))
            if isinstance(support_dict, dict) and "type" not in support_dict:
                support_type = getattr(self.support, "type", None)
                if support_type is not None:
                    support_dict["type"] = getattr(support_type, "numerator", support_type)
            payload["support"] = support_dict
        return payload

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    @classmethod
    def from_json(cls, json_str: str) -> "GeologicalInterpolator":
        return cls.from_dict(json.loads(json_str))

    def to_yaml(self, file_path: Optional[str] = None) -> None | str:
        try:
            import yaml
        except ImportError as exc:
            raise ImportError("PyYAML is required for YAML export: pip install pyyaml") from exc
        if file_path is None:
            return yaml.safe_dump(self.to_dict(), sort_keys=False, allow_unicode=True)
        with open(file_path, "w") as f:
            yaml.safe_dump(self.to_dict(), f, sort_keys=False, allow_unicode=True)

    @classmethod
    def from_yaml(cls, yaml_str: str) -> "GeologicalInterpolator":
        try:
            import yaml
        except ImportError as exc:
            raise ImportError("PyYAML is required for YAML import: pip install pyyaml") from exc
        payload = yaml.safe_load(yaml_str)
        return cls.from_dict(payload)

    @classmethod
    def from_dict(cls, data):
        from ._interpolator_factory import InterpolatorFactory

        return InterpolatorFactory.from_dict(data.copy())

    def clean(self):
        """
        Removes all of the data from an interpolator

        Returns
        -------

        """
        self.data = {
            "gradient": np.zeros((0, 7)),
            "value": np.zeros((0, 5)),
            "normal": np.zeros((0, 7)),
            "tangent": np.zeros((0, 7)),
            "interface": np.zeros((0, 5)),
            "inequality": np.zeros((0, 6)),
            "inequality_pairs": np.zeros((0, 4)),
        }
        self.up_to_date = False
        self.n_g = 0
        self.n_i = 0
        self.n_n = 0
        self.n_t = 0

    def debug(self):
        """Helper function for debugging when the interpolator isn't working"""
        error_string = ""
        error_code = 0
        if (
            self.type > InterpolatorType.BASE_DISCRETE
            and self.type < InterpolatorType.BASE_DATA_SUPPORTED
        ):

            def mask(xyz):
                return self.support.inside(xyz)

        else:

            def mask(xyz):
                return np.ones(xyz.shape[0], dtype=bool)

        if (
            len(
                np.unique(
                    self.get_value_constraints()[mask(self.get_value_constraints()[:, :3]), 3]
                )
            )
            == 1
        ):
            error_code += 1
            error_string += "There is only one unique value in the model interpolation support \n"
            error_string += "Try increasing the model bounding box \n"
        if len(self.get_norm_constraints()[mask(self.get_norm_constraints()[:, :3]), :]) == 0:
            error_code += 1
            error_string += "There are no norm constraints in the model interpolation support \n"
            error_string += "Try increasing the model bounding box or adding more data\n"
        if error_code > 1:
            print(error_string)
