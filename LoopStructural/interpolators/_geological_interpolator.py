"""
Base geological interpolator
"""

from abc import ABCMeta, abstractmethod
from LoopStructural.utils.exceptions import LoopTypeError
from ..interpolators import InterpolatorType
import numpy as np

from typing import Optional
from ..utils import getLogger

logger = getLogger(__name__)


class GeologicalInterpolator(metaclass=ABCMeta):
    """
    Attributes
    ----------
    data : dict
        a dictionary with np.arrays for gradient, value, normal, tangent data
    """

    @abstractmethod
    def __init__(self, data={}, up_to_date=False):
        """
        This class is the base class for a geological interpolator and contains
        all of the main interface functions. Any class that is inheriting from
        this should be callable by using any of these functions. This will
        enable interpolators to be interchanged.
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

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        if data is None:
            data = {}
        for k, v in data.items():
            self._data[k] = np.array(v)

    def __str__(self):
        name = f"{self.type} \n"
        name += f"{self.n_g} gradient points\n"
        name += f"{self.n_i} interface points\n"
        name += f"{self.n_n} normal points\n"
        name += f"{self.n_t} tangent points\n"
        name += f"{self.n_g + self.n_i + self.n_n + self.n_t} total points\n"
        return name

    def check_array(self, array: np.ndarray):
        try:
            return np.array(array)
        except Exception as e:
            raise LoopTypeError(str(e))

    def to_json(self):
        """
        Returns a json representation of the geological interpolator

        Returns
        -------
        json : dict
            json representation of the geological interpolator
        """
        json = {}
        json["type"] = self.type
        # json["name"] = self.propertyname
        json["constraints"] = self.constraints
        json["data"] = self.data
        json["type"] = self.type
        # json["dof"] = self.nx
        json["up_to_date"] = self.up_to_date
        return json

    @abstractmethod
    def set_region(self, **kwargs):
        pass

    def set_value_constraints(self, points: np.ndarray):
        """

        Parameters
        ----------
        points : np.ndarray
            array containing the value constraints usually 4-5 columns.
            X,Y,Z,val,weight

        Returns
        -------

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
        """

        Parameters
        ----------
        points : np.ndarray
            array containing the value constraints usually 7-8 columns.
            X,Y,Z,gx,gy,gz,weight

        Returns
        -------

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
            X,Y,Z,nx,ny,nz,weight

        Returns
        -------

        """
        if points.shape[1] == self.dimensions * 2:
            points = np.hstack([points, np.ones((points.shape[0], 1))])
        if points.shape[1] < self.dimensions * 2 + 1:
            raise ValueError("Nonrmal constraints must at least have X,Y,Z,nx,ny,nz")
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
