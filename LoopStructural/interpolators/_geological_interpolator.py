"""Base geological interpolator for LoopStructural.

This module contains the abstract base class for all geological interpolators
used in LoopStructural geological modelling framework.
"""

from abc import ABCMeta, abstractmethod
from LoopStructural.utils.exceptions import LoopTypeError
from ..interpolators import InterpolatorType
import numpy as np

from typing import Optional
from ..utils import getLogger

logger = getLogger(__name__)


class GeologicalInterpolator(metaclass=ABCMeta):
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
    def __init__(self, data={}, up_to_date=False):
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

    def to_json(self):
        """Return a JSON representation of the geological interpolator.

        Returns
        -------
        dict
            Dictionary containing the interpolator's state and configuration
            suitable for JSON serialization

        Notes
        -----
        This method packages the essential state of the interpolator including
        its type, constraints, data, and build status for serialization.
        """
        json = {}
        json["type"] = self.type
        # json["name"] = self.propertyname
        json["constraints"] = self.constraints
        json["data"] = self.data
        json["type"] = self.type
        # json["dof"] = self.dof
        json["up_to_date"] = self.up_to_date
        return json

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

    def set_value_constraints(self, points: np.ndarray):
        """Set value constraints for the interpolation.

        Parameters
        ----------
        points : np.ndarray
            Array containing the value constraints with shape (n_points, 4-5).
            Columns should be [X, Y, Z, value, weight]. If weight is not provided,
            a weight of 1.0 is assumed for all points.

        Raises
        ------
        ValueError
            If points array doesn't have the minimum required columns

        Notes
        -----
        Value constraints specify known scalar field values at specific locations.
        These are typically used for interface points or measured data values.
        """
        points = self.check_array(points)
        if points.shape[1] == self.dimensions + 1:
            points = np.hstack([points, np.ones((points.shape[0], 1))])
        if points.shape[1] < self.dimensions + 2:
            raise ValueError("Value points must at least have X,Y,Z,val,w")
        self.data["value"] = points
        self.n_i = points.shape[0]
        self.up_to_date = False

    def set_gradient_constraints(self, points: np.ndarray):
        """Set gradient constraints for the interpolation.

        Parameters
        ----------
        points : np.ndarray
            Array containing gradient constraints with shape (n_points, 7-8).
            Columns should be [X, Y, Z, gx, gy, gz, weight]. If weight is not 
            provided, a weight of 1.0 is assumed for all points.

        Raises
        ------
        ValueError
            If points array doesn't have the minimum required columns

        Notes
        -----
        Gradient constraints specify the direction and magnitude of the scalar
        field gradient at specific locations. These are typically derived from
        structural measurements like bedding or foliation orientations.
        """
        if points.shape[1] == self.dimensions * 2:
            points = np.hstack([points, np.ones((points.shape[0], 1))])
        if points.shape[1] < self.dimensions * 2 + 1:
            raise ValueError("Gradient constraints must at least have X,Y,Z,gx,gy,gz")
        self.n_g = points.shape[0]
        self.data["gradient"] = points
        self.up_to_date = False

    def set_normal_constraints(self, points: np.ndarray):
        """

        Parameters
        ----------
        points : np.ndarray
            array containing the value constraints usually 7-8 columns.
            X,Y,Z,nx,ny,nz,(weight, default : 1 for each row)

        Returns
        -------

        Notes
        -------
        If no weights are provided, w = 1 is assigned to each normal constraint.

        """
        if points.shape[1] == self.dimensions * 2:
            points = np.hstack([points, np.ones((points.shape[0], 1))])
            logger.info("No weight provided for normal constraints, all weights are set to 1")
        if points.shape[1] < self.dimensions * 2 + 1:
            raise ValueError("Normal constraints must at least have X,Y,Z,nx,ny,nz")
        self.n_n = points.shape[0]
        self.data["normal"] = points
        self.up_to_date = False

    def set_tangent_constraints(self, points: np.ndarray):
        """

        Parameters
        ----------
        points : np.ndarray
            array containing the value constraints usually 7-8 columns.
            X,Y,Z,nx,ny,nz,weight

        Returns
        -------

        """
        if points.shape[1] == self.dimensions * 2:
            points = np.hstack([points, np.ones((points.shape[0], 1))])
        if points.shape[1] < self.dimensions * 2 + 1:
            raise ValueError("Tangent constraints must at least have X,Y,Z,tx,ty,tz")
        self.data["tangent"] = points
        self.up_to_date = False

    def set_interface_constraints(self, points: np.ndarray):
        self.data["interface"] = points
        self.up_to_date = False

    def set_value_inequality_constraints(self, points: np.ndarray):
        if points.shape[1] < self.dimensions + 2:
            raise ValueError("Inequality constraints must at least have X,Y,Z,lower,upper")
        self.data["inequality"] = points
        self.up_to_date = False

    def set_inequality_pairs_constraints(self, points: np.ndarray):
        if points.shape[1] < self.dimensions + 1:
            raise ValueError("Inequality pairs constraints must at least have X,Y,Z,rock_id")

        self.data["inequality_pairs"] = points
        self.up_to_date = False

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

    # @abstractmethod
    def setup(self, **kwargs):
        """
        Runs all of the required setting up stuff
        """
        self.setup_interpolator(**kwargs)

    @abstractmethod
    def setup_interpolator(self, **kwargs):
        """
        Runs all of the required setting up stuff
        """
        self.setup_interpolator(**kwargs)

    @abstractmethod
    def solve_system(self, solver, solver_kwargs: dict = {}) -> bool:
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
        return {
            "type": self.type,
            "data": self.data,
            "up_to_date": self.up_to_date,
            "valid": self.valid,
        }

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
