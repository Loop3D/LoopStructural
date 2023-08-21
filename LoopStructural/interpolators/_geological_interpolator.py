"""
Base geological interpolator
"""
from ..interpolators import InterpolatorType
import numpy as np

from ..utils import getLogger

logger = getLogger(__name__)


class GeologicalInterpolator:
    """
    Attributes
    ----------
    data : dict
        a dictionary with np.arrays for gradient, value, normal, tangent data
    """

    def __init__(self):
        """
        This class is the base class for a geological interpolator and contains
        all of the main interface functions. Any class that is inheriting from
        this should be callable by using any of these functions. This will
        enable interpolators to be interchanged.
        """
        self.data = {}  # None
        self.clean()  # init data structure

        self.n_g = 0
        self.n_i = 0
        self.n_n = 0
        self.n_t = 0

        self.type = InterpolatorType.BASE
        self.up_to_date = False
        self.constraints = []
        self.propertyname = "defaultproperty"
        self.__str = "Base Geological Interpolator"
        self.valid = False

    def __str__(self):
        name = f"{self.type} \n"
        name += f"{self.n_g} gradient points\n"
        name += f"{self.n_i} interface points\n"
        name += f"{self.n_n} normal points\n"
        name += f"{self.n_t} tangent points\n"
        name += f"{self.n_g + self.n_i + self.n_n + self.n_t} total points\n"
        return name

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
        json["name"] = self.propertyname
        json["constraints"] = self.constraints
        json["data"] = self.data
        json["type"] = self.type
        json["dof"] = self.nx
        json["up_to_date"] = self.up_to_date
        return json

    def set_region(self, **kwargs):
        pass

    def set_property_name(self, name):
        """
        Set the name of the interpolated property
        Parameters
        ----------
        name : string
            name of the property to be saved on a mesh

        Returns
        -------

        """
        self.propertyname = name

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
        if points.shape[1] < 4:
            raise ValueError("Value points must at least have X,Y,Z,val")
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
        if points.shape[1] < 7:
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
        if points.shape[1] < 7:
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
        if points.shape[1] < 7:
            raise ValueError("Tangent constraints must at least have X,Y,Z,tx,ty,tz")
        self.data["tangent"] = points
        self.up_to_date = False

    def set_interface_constraints(self, points: np.ndarray):
        self.data["interface"] = points
        self.up_to_date = False

    def set_inequality_constraints(self, points: np.ndarray):
        self.data["inequality"] = points
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

    def get_inequality_constraints(self):
        return self.data["inequality"]

    def setup_interpolator(self, **kwargs):
        """
        Runs all of the required setting up stuff
        """
        self._setup_interpolator(**kwargs)

    def solve_system(self, **kwargs):
        """
        Solves the interpolation equations
        """
        self._solve(**kwargs)
        self.up_to_date = True

    def update(self):
        return False

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
                    self.get_value_constraints()[
                        mask(self.get_value_constraints()[:, :3]), 3
                    ]
                )
            )
            == 1
        ):
            error_code += 1
            error_string += (
                "There is only one unique value in the model interpolation support \n"
            )
            error_string += "Try increasing the model bounding box \n"
        if (
            len(
                self.get_norm_constraints()[mask(self.get_norm_constraints()[:, :3]), :]
            )
            == 0
        ):
            error_code += 1
            error_string += (
                "There are no norm constraints in the model interpolation support \n"
            )
            error_string += (
                "Try increasing the model bounding box or adding more data\n"
            )
        if error_code > 1:
            print(error_string)
