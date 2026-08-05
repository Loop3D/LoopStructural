from abc import ABC, abstractmethod
import beartype
import pandas
import geopandas
from .mapdata import MapData
import numpy as np

from .logging import getLogger

logger = getLogger(__name__)


class FaultOrientation(ABC):
    """
    Base Class of Fault Orientation assigner to force structure of FaultOrientation

    Args:
        ABC (ABC): Derived from Abstract Base Class
    """

    def __init__(self):
        """
        Initialiser of for orientation assigner
        """
        self.label = "FaultOrientationBaseClass"

    def type(self):
        """
        Getter for subclass type label

        Returns:
            str: Name of subclass
        """
        return self.label

    @beartype.beartype
    @abstractmethod
    def calculate(
        self,
        fault_trace: geopandas.GeoDataFrame,
        fault_orientations: pandas.DataFrame,
        map_data: MapData,
    ) -> pandas.DataFrame:
        """
        Execute fault orientation assigned method (abstract method)

        Args:
            faults (pandas.DataFrame): the data frame of the faults to add throw values to
            fault_orientation (pandas.DataFrame): data frame with fault orientations to assign to faults
            map_data (map2loop.MapData): a catchall so that access to all map data is available

        Returns:
            pandas.DataFrame: fault orientations assigned to a fault label
        """
        pass


class FaultOrientationNearest(FaultOrientation):
    """
    FaultOrientation class which estimates fault orientation based on nearest orientation
    """

    def __init__(self):
        """
        Initialiser for nearest version of the fault orientation assigner
        """
        self.label = "FaultOrientationNearest"

    @beartype.beartype
    def calculate(
        self,
        fault_trace: geopandas.GeoDataFrame,
        fault_orientations: pandas.DataFrame,
        map_data: MapData,
    ) -> pandas.DataFrame:
        """
        Assigns the nearest fault orientation to a fault

        Args:
            fault_trace (geopandas.GeoDataFrame): the data frame of the fault traces
            fault_orientation (pandas.DataFrame): data frame with fault orientations to assign to faults
            map_data (map2loop.MapData): a catchall so that access to all map data is available

        Returns:
            pandas.DataFrame: fault orientations assigned to a fault label
        """
        logger.info("Assigning fault orientations to fault traces from nearest orientation")
        orientations = fault_orientations.copy()
        logger.info(f'There are {len(orientations)} fault orientations to assign')

        orientations["ID"] = -1

        for i in orientations.index:
            p = orientations.loc[i, :].geometry
            orientations.loc[i, "ID"] = fault_trace.loc[
                fault_trace.index[np.argmin(fault_trace.distance(p))], "ID"
            ]
            orientations.loc[i, "X"] = p.x
            orientations.loc[i, "Y"] = p.y

        return orientations.drop(columns="geometry")
