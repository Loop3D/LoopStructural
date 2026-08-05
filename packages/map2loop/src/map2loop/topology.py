from typing import Optional
# internal imports
from .m2l_enums import VerboseLevel
from .contact_extractor import ContactExtractor

# external imports
import geopandas
import pandas
import numpy

from .logging import getLogger

logger = getLogger(__name__)


class Topology:
    """
    Topology calculator for units and faults relationships

    Attributes
    ----------
    sorted_units: None or list
        the stratigraphic column
    fault_fault_relationships: None or pandas.DataFrame
        data frame of fault to fault relationships with columns ["Fault1", "Fault2", "Type", "Angle"]
    unit_fault_relationships: None or pandas.DataFrame
        data frame of unit fault relationships with columns ["Unit", "Fault"]
    unit_unit_relationships: None or pandas.DataFrame
        data frame of unit unit relationships with columns ["Index1", "UnitName1", "Index2", "UnitName2"]
    map_data: MapData
        A pointer to the map data structure in project
    verbose_level: m2l_enum.VerboseLevel
        A selection that defines how much console logging is output
    """

    def __init__(
        self, geology_data: geopandas.GeoDataFrame, 
        fault_data: Optional[geopandas.GeoDataFrame] = None,  
        verbose_level: VerboseLevel = VerboseLevel.NONE,
        id_field: str = 'ID',
    ):
        """
        The initialiser for the map2model wrapper

        Args:
            map_data (MapData):
                The project map data structure to reference
            verbose_level (VerboseLevel, optional):
                How much console output is sent. Defaults to VerboseLevel.ALL.
        """
        self.sorted_units = None
        self._fault_fault_relationships = None
        self._unit_fault_relationships = None
        self._unit_unit_relationships = None
        self.geology_data = geology_data
        self.fault_data = fault_data
        self.verbose_level = verbose_level
        self.buffer_radius = 500
        self.fault_id_field = id_field
    @property
    def fault_fault_relationships(self):
        if self._fault_fault_relationships is None:
            self._calculate_fault_fault_relationships()
            
        return self._fault_fault_relationships

    @property
    def unit_fault_relationships(self):
        if self._unit_fault_relationships is None:
            self._calculate_fault_unit_relationships()
            
        return self._unit_fault_relationships

    @property
    def unit_unit_relationships(self):
        if self._unit_unit_relationships is None:
            self._calculate_unit_unit_relationships()
            
        return self._unit_unit_relationships

    def reset(self):
        """
        Reset before the next run
        """
        logger.info("Resetting topology calculator")
        self.sorted_units = None
        self._fault_fault_relationships = None
        self._unit_fault_relationships = None
        self._unit_unit_relationships = None

    def get_sorted_units(self):
        """
        Getter for the map2model sorted units

        Returns:
            list: The map2model stratigraphic column estimate
        """
        raise NotImplementedError("This method is not implemented")
        

    def get_fault_fault_relationships(self):
        """
        Getter for the fault fault relationships

        Returns:
            pandas.DataFrame: The fault fault relationships
        """

        return self.fault_fault_relationships

    def get_unit_fault_relationships(self):
        """
        Getter for the unit fault relationships

        Returns:
            pandas.DataFrame: The unit fault relationships
        """

        return self.unit_fault_relationships

    def get_unit_unit_relationships(self):
        """
        Getter for the unit unit relationships

        Returns:
            pandas.DataFrame: The unit unit relationships
        """

        return self.unit_unit_relationships

    def _calculate_fault_fault_relationships(self):

        faults = self.fault_data.copy()
        # reset index so that we can index the adjacency matrix with the index
        faults.reset_index(inplace=True)
        buffers = faults.buffer(self.buffer_radius)
        # create the adjacency matrix
        intersection = geopandas.sjoin(
            geopandas.GeoDataFrame(geometry=buffers), geopandas.GeoDataFrame(geometry=faults["geometry"])
        )
        intersection["index_left"] = intersection.index
        intersection.reset_index(inplace=True)

        adjacency_matrix = numpy.zeros((faults.shape[0], faults.shape[0]), dtype=bool)
        adjacency_matrix[
            intersection.loc[:, "index_left"], intersection.loc[:, "index_right"]
        ] = True
        f1, f2 = numpy.where(numpy.tril(adjacency_matrix, k=-1))
        df = pandas.DataFrame(
            {'Fault1': faults.loc[f1, 'ID'].to_list(), 'Fault2': faults.loc[f2, 'ID'].to_list()}
        )
        df['Angle'] = 60  # make it big to prevent LS from making splays
        df['Type'] = 'T'
        self._fault_fault_relationships = df

    def _calculate_fault_unit_relationships(self):
        """Calculate unit/fault relationships using geopandas sjoin.
        This will return
        """
        units = self.geology_data["UNITNAME"].unique()
        faults = self.fault_data.copy().reset_index().drop(columns=['index'])
        adjacency_matrix = numpy.zeros((len(units), faults.shape[0]), dtype=bool)
        for i, u in enumerate(units):
            unit = self.geology_data[self.geology_data["UNITNAME"] == u]
            intersection = geopandas.sjoin(
                geopandas.GeoDataFrame(geometry=faults["geometry"]),
                geopandas.GeoDataFrame(geometry=unit["geometry"]),
            )
            intersection["index_left"] = intersection.index
            intersection.reset_index(inplace=True)
            adjacency_matrix[i, intersection.loc[:, "index_left"]] = True
        u, f = numpy.where(adjacency_matrix)
        df = pandas.DataFrame({"Unit": units[u].tolist(), "Fault": faults.loc[f, "ID"].to_list()})
        self._unit_fault_relationships = df

    def _calculate_unit_unit_relationships(self):
        extractor = ContactExtractor(
            self.geology_data,
            self.fault_data,
        )
        contacts = extractor.extract_all_contacts()
        self._unit_unit_relationships = contacts.copy().drop(
            columns=['length', 'geometry']
        )
        return self._unit_unit_relationships

    def run(self, verbose_level: VerboseLevel = None):
        """
        The main execute function that prepares, runs and parse the output of the map2model process

        Args:
            verbose_level (VerboseLevel, optional):
                How much console output is sent. Defaults to None (which uses the wrapper attribute).
        """

        self.get_fault_fault_relationships()
        self.get_unit_fault_relationships()
        self.get_unit_unit_relationships()
        
