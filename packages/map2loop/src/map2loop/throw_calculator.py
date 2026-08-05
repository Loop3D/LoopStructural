from abc import ABC, abstractmethod
import beartype
import pandas
import geopandas
from .mapdata import MapData


class ThrowCalculator(ABC):
    """
    Base Class of Throw Calculator used to force structure of ThrowCalculator

    Args:
        ABC (ABC): Derived from Abstract Base Class
    """

    def __init__(self):
        """
        Initialiser of for Sorter
        """
        self.throw_calculator_label = "ThrowCalculatorBaseClass"

    def type(self):
        """
        Getter for subclass type label

        Returns:
            str: Name of subclass
        """
        return self.throw_calculator_label

    @beartype.beartype
    @abstractmethod
    def compute(
        self,
        faults: pandas.DataFrame,
        stratigraphic_order: list,
        basal_contacts: geopandas.GeoDataFrame,
        map_data: MapData,
    ) -> pandas.DataFrame:
        """
        Execute throw calculator method (abstract method)

        Args:
            faults (pandas.DataFrame): the data frame of the faults to add throw values to
            stratigraphic_order (list): a list of unit names sorted from youngest to oldest
            basal_contacts (geopandas.GeoDataFrame): basal contact geo data with locations and unit names of the contacts (columns must contain ["ID","basal_unit","type","geometry"])
            map_data (map2loop.MapData): a catchall so that access to all map data is available

        Returns:
            pandas.DataFrame: fault data frame with throw values (avgDisplacement and avgDownthrowDir) filled in
        """
        pass


class ThrowCalculatorAlpha(ThrowCalculator):
    """
    ThrowCalculator class which estimates fault throw values based on units, basal_contacts and stratigraphic order
    """

    def __init__(self):
        """
        Initialiser for alpha version of the throw calculator
        """
        self.throw_calculator_label = "ThrowCalculatorAlpha"

    @beartype.beartype
    def compute(
        self,
        faults: pandas.DataFrame,
        stratigraphic_order: list,
        basal_contacts: pandas.DataFrame,
        map_data: MapData,
    ) -> pandas.DataFrame:
        """
        Execute throw calculator method takes fault data, basal_contacts and stratigraphic order and attempts to estimate fault throw.

        Args:
            faults (pandas.DataFrame): the data frame of the faults to add throw values to
            stratigraphic_order (list): a list of unit names sorted from youngest to oldest
            basal_contacts (geopandas.GeoDataFrame): basal contact geo data with locations and unit names of the contacts (columns must contain ["ID","basal_unit","type","geometry"])
            map_data (map2loop.MapData): a catchall so that access to all map data is available

        Returns:
            pandas.DataFrame: fault data frame with throw values (avgDisplacement and avgDownthrowDir) filled in
        """
        # For each fault take the geometric join of all contact lines and that fault line

        # For each contact join take the length of that join as an approximation of the minimum throw of the fault

        # Take all the min throw approximations and set the largest one as the avgDisplacement

        # If a fault has no contact lines the maximum throw should be less than the thickness of the containing
        # unit (if we exclude map height changes and fault angle)

        # Set any remaining displacement values to default value
        faults["avgDisplacement"] = faults.apply(
            lambda row: 100 if row["avgDisplacement"] == -1 else row["avgDisplacement"], axis=1
        )
        return faults
