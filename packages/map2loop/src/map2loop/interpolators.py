from abc import ABC, abstractmethod
from typing import Any, Union
import beartype
import numpy
from numpy import ndarray
from scipy.interpolate import Rbf, LinearNDInterpolator
from sklearn.cluster import DBSCAN
import pandas


from .utils import strike_dip_vector, generate_grid

from .logging import getLogger
logger = getLogger(__name__)  

class Interpolator(ABC):
    """
    Base Class of Interpolator used to force structure of Interpolator

    Args:
        ABC (ABC): Derived from Abstract Base Class
    """

    def __init__(self):
        """
        Initialiser of for Interpolator
        """
        self.interpolator_label = "InterpolatorBaseClass"

    def type(self):
        """
        Getter for subclass type label

        Returns:
            str: Name of subclass
        """
        return self.interpolator_label

    @beartype.beartype
    @abstractmethod
    def setup_interpolation(self, structure_data: pandas.DataFrame):
        """
        abstract method to setup interpolation (abstract method)

        Args:
            structure_data (pandas.DataFrame): sampled structural data
        """
        pass

    @beartype.beartype
    @abstractmethod
    def setup_grid(self, bounding_box: dict):
        """
        abstract method to setup an XY grid (abstract method)

        Args:
            bounding_box (dict): a dictionary containing the bounding box of the map data.
            The bounding box dictionary should comply with the following format: {
                "minx": value,
                "maxx": value,
                "miny": value,
                "maxy": value,
            }
        """
        pass

    @beartype.beartype
    @abstractmethod
    def interpolate(self, ni: Union[list, numpy.ndarray]):
        """
        Interpolator method

        Args:

            ni (int): number of points


        Returns:
            float: interpolated value
        """

        pass

    @beartype.beartype
    @abstractmethod
    def __call__(
        self, bounding_box: dict, structure_data: pandas.DataFrame, interpolator: Any = None
    ) -> Any:
        """
        Execute interpolate method (abstract method)

        Args:
            bounding_box (dict): a dictionary containing the bounding box of the map data
            structure_data (pandas.DataFrame): sampled structural data
            interpolator: type of interpolator to use by default SciPy Rbf interpolator

        Returns:
            list: sorted list of unit names

        """
        pass


class NormalVectorInterpolator(Interpolator):
    """
    This class is a subclass of the Interpolator abstract base class. It implements the normal vector interpolation
    method for a given set of data points.

    Attributes:
        dataframe (pandas.DataFrame): A DataFrame that stores the processed data points for interpolation.
        x (numpy.ndarray): A numpy array that stores the x-coordinates of the data points.
        y (numpy.ndarray): A numpy array that stores the y-coordinates of the data points.
        xi (numpy.ndarray): A numpy array that stores the x-coordinates of the grid points for interpolation.
        yi (numpy.ndarray): A numpy array that stores the y-coordinates of the grid points for interpolation.
        interpolator_label (str): A string that stores the label of the interpolator. For this class, it is
        "NormalVectorInterpolator".

    Methods:
        type(): Returns the label of the interpolator.
        setup_interpolation(structure_data: pandas.DataFrame): Sets up the interpolation by preparing the data points for interpolation.
        setup_grid(bounding_box: dict): Sets up the grid for interpolation.
        interpolator(ni: Any) -> numpy.ndarray: Performs the interpolation for a given set of values.
        interpolate(bounding_box: dict, structure_data: pandas.DataFrame, interpolator: Any) -> numpy.ndarray: Executes the interpolation method.
    """

    def __init__(self):
        """
        Initialiser of for NormalVectorInterpolator class
        """
        self.dataframe = None
        self.x = None
        self.y = None
        self.xi = None
        self.yi = None
        self.interpolator_label = "NormalVectorInterpolator"

    def type(self):
        """
        Getter for subclass type label

        Returns:
            str: Name of subclass
        """
        return self.interpolator_label

    @beartype.beartype
    def setup_interpolation(self, structure_data: pandas.DataFrame):
        """
        Setup the interpolation method (abstract method)

        Args:
            structure_data (pandas.DataFrame): sampled structural data
        """
        # the following code is a slightly modified version from LoopStructural's InputDataProcessor
        if (
            "nx" not in structure_data.columns
            or "ny" not in structure_data.columns
            or "nz" not in structure_data.columns
        ):
            if "strike" not in structure_data.columns and "azimuth" in structure_data.columns:
                structure_data["strike"] = structure_data["azimuth"] - 90
            if "strike" not in structure_data.columns and "dipdir" in structure_data.columns:
                structure_data["strike"] = structure_data["dipdir"] - 90

            if "strike" not in structure_data.columns and "dipDir" in structure_data.columns:
                structure_data["strike"] = structure_data["dipDir"] - 90
            if "strike" not in structure_data.columns and "DIPDIR" in structure_data.columns:
                structure_data["strike"] = structure_data["DIPDIR"] - 90

            if "strike" in structure_data.columns and "dip" in structure_data.columns:
                structure_data["nx"] = numpy.nan
                structure_data["ny"] = numpy.nan
                structure_data["nz"] = numpy.nan
                structure_data[["nx", "ny", "nz"]] = strike_dip_vector(
                    structure_data["strike"], structure_data["dip"]
                )
            if "strike" in structure_data.columns and "DIP" in structure_data.columns:
                structure_data["nx"] = numpy.nan
                structure_data["ny"] = numpy.nan
                structure_data["nz"] = numpy.nan
                structure_data[["nx", "ny", "nz"]] = strike_dip_vector(
                    structure_data["strike"], structure_data["DIP"]
                )

            if (
                "nx" not in structure_data.columns
                or "ny" not in structure_data.columns
                or "nz" not in structure_data.columns
            ):
                raise ValueError(
                    "Contact orientation data must contain either strike/dipdir, dip, or nx, ny, nz"
                )

            if (
                "X" not in structure_data.columns
                or "Y" not in structure_data.columns
                or "Z" not in structure_data.columns
            ):
                if "X" not in structure_data.columns and "easting" in structure_data.columns:
                    structure_data["X"] = structure_data["easting"]
                if "X" not in structure_data.columns and "easting" not in structure_data.columns:
                    structure_data["X"] = structure_data["geometry"].apply(lambda geom: geom.x)

                if "Y" not in structure_data.columns and "northing" in structure_data.columns:
                    structure_data["Y"] = structure_data["northing"]
                if "Y" not in structure_data.columns and "northing" not in structure_data.columns:
                    structure_data["Y"] = structure_data["geometry"].apply(lambda geom: geom.y)
                if "Z" not in structure_data.columns and "altitude" in structure_data.columns:
                    structure_data["Z"] = structure_data["altitude"]
                if "Z" not in structure_data.columns and "altitude" not in structure_data.columns:
                    structure_data["Z"] = 0

                if (
                    "X" not in structure_data.columns
                    or "Y" not in structure_data.columns
                    or "Z" not in structure_data.columns
                ):
                    raise ValueError(
                        "Contact orientation data must contain either X, Y, Z or easting, northing, altitude"
                    )

        self.dataframe = structure_data
        self.x = structure_data["X"].to_numpy()
        self.y = structure_data["Y"].to_numpy()

    @beartype.beartype
    def setup_grid(self, bounding_box: dict):
        """
        Setup the grid for interpolation

        Args:
            bounding_box (dict): a dictionary containing the bounding box of the map data.
            The bounding box dictionary should comply with the following format: {
                "minx": value,
                "maxx": value,
                "miny": value,
                "maxy": value,
            }
        """
        self.xi, self.yi = generate_grid(bounding_box)

    @beartype.beartype
    def interpolate(self, ni: Union[ndarray, list], interpolator: Any = Rbf) -> numpy.ndarray:
        # TODO: 1. add code to use LoopStructural interpolators
        """
        Inverse Distance Weighting interpolation method

        Args:
            ni (int): value to interpolate
            interpolator: type of interpolator to use by default SciPy Rbf interpolator

        Returns:
            Rbf: radial basis function object
        """
        if interpolator is Rbf:
            rbf = Rbf(self.x, self.y, ni, function="linear")
            return rbf(self.xi, self.yi)

        if interpolator is LinearNDInterpolator:
            lnd_interpolator = LinearNDInterpolator(list(zip(self.x, self.y)), ni)
            return lnd_interpolator(self.xi, self.yi)

    @beartype.beartype
    def __call__(
        self, bounding_box: dict, structure_data: pandas.DataFrame, interpolator: Any = Rbf
    ) -> numpy.ndarray:
        """
        Execute interpolation method

        Args:
            bounding_box (dict): a dictionary containing the bounding box of the map data
            structure_data (pandas.DataFrame): sampled structural data
            interpolator: type of interpolator to use by default SciPy Rbf interpolator

        Returns:
            list: sorted list of unit names

        """
        self.setup_interpolation(structure_data)
        self.setup_grid(bounding_box)
        # get normal vector components
        nx = self.dataframe["nx"].to_numpy()
        ny = self.dataframe["ny"].to_numpy()
        nz = self.dataframe["nz"].to_numpy()

        # interpolate each component of the normal vector nx, ny, nz
        nx_interp = self.interpolate(nx)
        ny_interp = self.interpolate(ny)
        nz_interp = self.interpolate(nz)

        vecs = numpy.array([nx_interp, ny_interp, nz_interp]).T
        # normalize the vectors
        vecs /= numpy.linalg.norm(vecs, axis=1)[:, None]

        return vecs


class DipDipDirectionInterpolator(Interpolator):
    """


    Args:
        Interpolator(ABC): Derived from Abstract Base Class
    """

    def __init__(self, data_type=None):
        """
        Initialiser of for IDWInterpolator
        """
        if data_type is None:
            self.data_type = ["dip", "dipdir"]
        else:
            self.data_type = data_type
        self.x = None
        self.y = None
        self.xi = None
        self.yi = None
        self.dip = None
        self.dipdir = None
        self.cell_size = None
        self.interpolator_label = "DipDipDirectionInterpolator"

    def type(self):
        """
        Getter for subclass type label

        Returns:
            str: Name of subclass
        """
        return self.interpolator_label

    @beartype.beartype
    def setup_interpolation(self, structure_data: pandas.DataFrame):
        """
        Setup the interpolation method

        Args:
            structure_data (pandas.DataFrame): sampled structural data
        """
        # Check for collocated points and remove them
        structure_data = structure_data.drop_duplicates(subset=['X', 'Y']).copy()
        # Check for collocated point clusters and average them
        coords = structure_data[["X", "Y"]].values
        db = DBSCAN(eps=self.cell_size, min_samples=1).fit(coords)
        structure_data["cluster"] = db.labels_
        # Check if there are any clusters with more than one point (indicating collocated points)
        collocated_clusters = structure_data['cluster'].value_counts()
        collocated_clusters = collocated_clusters[collocated_clusters > 1]

        if not collocated_clusters.empty:
            # Log a warning if collocated points are detected
            logger.warning(
                f"Detected {len(collocated_clusters)} collocated point clusters. Aggregating these points.\n " 
            )

        # Aggregate data for collocated points by taking the mean of X, Y, DIP, and DIPDIR within each cluster
        aggregated_data = (
            structure_data.groupby("cluster")
            .agg({"X": "mean", "Y": "mean", "DIP": "mean", "DIPDIR": "mean"})
            .reset_index(drop=True)
        )

        # setup variables for interpolation
        self.x = aggregated_data["X"].to_numpy()
        self.y = aggregated_data["Y"].to_numpy()
        if "dip" in self.data_type:
            self.dip = aggregated_data["DIP"].to_numpy()
        if "dipdir" in self.data_type:
            self.dipdir = aggregated_data["DIPDIR"].to_numpy()

    @beartype.beartype
    def setup_grid(self, bounding_box: dict):
        """
        Setup the grid for interpolation

        Args:
            bounding_box (dict): a dictionary containing the bounding box of the map data.
            The bounding box dictionary should comply with the following format: {
                "minx": value,
                "maxx": value,
                "miny": value,
                "maxy": value,
            }
        """
        self.xi, self.yi, self.cell_size = generate_grid(bounding_box)

    @beartype.beartype
    def interpolate(self, ni: Union[ndarray, list], interpolator: Any = Rbf):
        # TODO: 1. add code to use LoopStructural interpolators
        """
        Inverse Distance Weighting interpolation method

        Args:
            ni (int): value to interpolate
            interpolator: type of interpolator to use by default SciPy Rbf interpolator

        Returns:
            Rbf: radial basis function object
        """
        if interpolator is Rbf:
            rbf = Rbf(self.x, self.y, ni, function="linear")
            return rbf(self.xi, self.yi)

        if interpolator is LinearNDInterpolator:
            lnd_interpolator = LinearNDInterpolator(list(zip(self.x, self.y)), ni)
            return lnd_interpolator(self.xi, self.yi)

    @beartype.beartype
    def __call__(
        self, bounding_box: dict, structure_data: pandas.DataFrame, interpolator: Any = Rbf
    ):
        """
        Execute interpolation method (abstract method)

        Args:
            bounding_box (dict): a dictionary containing the bounding box of the map data
            structure_data (pandas.DataFrame): sampled structural data
            interpolator (Union[Rbf, LinearNDInterpolator]): type of interpolator to use by default SciPy Rbf interpolator

        Returns:
            numpy.ndarray: interpolated dip and dip direction values

        """
        self.setup_grid(bounding_box)
        self.setup_interpolation(structure_data)

        # interpolate dip and dip direction
        if self.dip is not None and self.dipdir is not None:
            interpolated_dip = self.interpolate(self.dip, interpolator)
            interpolated_dipdir = self.interpolate(self.dipdir, interpolator)
            interpolated = numpy.array([interpolated_dip, interpolated_dipdir]).T
            return interpolated

        if self.dip is not None and self.dipdir is None:
            return self.interpolate(self.dip, interpolator)

        if self.dipdir is not None and self.dip is None:
            return self.interpolate(self.dipdir, interpolator)
