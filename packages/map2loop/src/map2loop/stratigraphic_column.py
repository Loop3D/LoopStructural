import pandas
import numpy
import geopandas


class StratigraphicColumn:
    """
    A class containing all the fault and fold summaries and relationships

    Attributes
    ----------
    column: list
        List of stratigraphic units in time order
    groups: list
        List of stratigraphic groups in time order
    stratigraphicUnitColumns: numpy.dtype
        Column names and types for stratigraphic unit summary
    stratigraphicUnits: pandas.DataFrame
        The stratigraphic units
    lithologyUnitColumns: numpy.dtype
        Column names and types for lithology layer summary
    lithologyUnits: pandas.DataFrame
        The lithology units

    """

    def __init__(self):
        """
        The initialiser for the stratigraphic units. All attributes are defaulted
        """
        self.column = []
        self.groups = []

        # Create empty dataframes for units
        self.stratigraphicUnitColumns = numpy.dtype(
            [
                ("layerId", int),
                ("name", str),
                ("minAge", float),
                ("maxAge", float),
                ("group", str),
                ("supergroup", str),
                ("colour", str),
            ]
        )
        self.stratigraphicUnits = pandas.DataFrame(
            numpy.empty(0, dtype=self.stratigraphicUnitColumns)
        )
        self.stratigraphicUnits = self.stratigraphicUnits.set_index("name")

        self.lithologyUnitColumns = numpy.dtype(
            [
                ("layerId", int),
                ("name", str),
                ("minAge", float),
                ("maxAge", float),
                ("group", str),
                ("colour", str),
            ]
        )
        self.lithologyUnits = pandas.DataFrame(numpy.empty(0, dtype=self.lithologyUnitColumns))
        self.lithologyUnits = self.lithologyUnits.set_index("name")

    def findStratigraphicUnit(self, id):
        """
        Find the unit in the units list based on its layerId or name

        Args:
            id (int or str):
                The layerId or name to look for

        Returns:
            pandas.DataFrame: The sliced data frame containing the requested unit
        """
        if issubclass(type(id), int):
            return self.stratigraphicUnits[self.stratigraphicUnits["layerId"] == id]
        elif issubclass(type(id), str):
            return self.stratigraphicUnits[self.stratigraphicUnits["name"] == id]
        else:
            print("ERROR: Unknown identifier type used to find stratigraphic unit")

    def findLithologyUnit(self, id):
        """
        Find the lithology unit in the units list based on its layerId or name

        Args:
            id (int or str):
                The layerId or name to look for

        Returns:
            pandas.DataFrame: The sliced data frame containing the requested unit
        """
        if issubclass(type(id), int):
            return self.lithologyUnits[self.lithologyUnits["layerId"] == id]
        elif issubclass(type(id), str):
            return self.lithologyUnits[self.lithologyUnits["name"] == id]
        else:
            print("ERROR: Unknown identifier type used to find lithology unit")

    def addStratigraphicUnit(self, unit):
        """
        Add stratigraphic unit to the units list

        Args:
            fault (pandas.DataFrame or dict):
                The unit information to add
        """
        if issubclass(type(unit), pandas.DataFrame) or issubclass(type(unit), dict):
            if "name" in unit.keys():
                if unit["name"] in self.stratigraphicUnits.index:
                    print("Replacing stratigraphic unit", unit["name"])
                self.stratigraphicUnits.loc[unit["name"]] = unit
            else:
                print("No name field in stratigraphic unit", unit)
        else:
            print("Cannot add unit to dataframe with type", type(unit))

    def addLithologyUnit(self, unit):
        """
        Add lithology unit to the units list

        Args:
            fault (pandas.DataFrame or dict):
                The unit information to add
        """
        if issubclass(type(unit), pandas.DataFrame) or issubclass(type(unit), dict):
            if "name" in unit.keys():
                if unit["name"] in self.lithologyUnits.index:
                    print("Replacing lithology unit", unit["name"])
                self.lithologyUnits.loc[unit["name"]] = unit
            else:
                print("No name field in lithology unit", unit)
        else:
            print("Cannot add unit to dataframe with type", type(unit))

    def populate(self, geology_map_data: geopandas.GeoDataFrame):
        """
        Parse the geodataframe data into the stratigraphic units list

        Args:
            geology_map_data (geopandas.GeoDataFrame):
                The geodataframe with the unit data
        """
        if geology_map_data.shape[0] == 0:
            return
        geology_data = geology_map_data.copy()
        geology_data = geology_data.drop_duplicates(subset=["UNITNAME"])
        geology_data = geology_data.reset_index(drop=True)
        # geology_data = geology_data.dropna(subset=["UNITNAME"])

        self.stratigraphicUnits = pandas.DataFrame(
            numpy.empty(geology_data.shape[0], dtype=self.stratigraphicUnitColumns)
        )
        self.stratigraphicUnits["layerId"] = numpy.arange(geology_data.shape[0])
        self.stratigraphicUnits["name"] = geology_data["UNITNAME"]
        self.stratigraphicUnits["minAge"] = geology_data["MIN_AGE"]
        self.stratigraphicUnits["maxAge"] = geology_data["MAX_AGE"]
        self.stratigraphicUnits["group"] = geology_data["GROUP"]
        self.stratigraphicUnits["supergroup"] = geology_data["SUPERGROUP"]
        self.stratigraphicUnits["colour"] = "#000000"
        # self.stratigraphicUnits["indexInGroup"] = -1

        self.groups = list(self.stratigraphicUnits['group'].unique())

    def set_stratigraphic_unit_parameter_by_name(self, name: str, parameter: str, value):
        """
        Set a specific parameter on a specific stratigraphic unit

        Args:
            name (str): The name of the stratigraphic unit
            parameter (str): The colmn name of the parameters
            value (str or int or float): The value to set
        """
        self.stratigraphicUnits.iloc[self.stratigraphicUnits["name"] == name][parameter] = value

    def sort_from_relationship_list(self, relationshipList: list):
        """
        Sort the stratigraphic column based on the list of name

        Args:
            relationshipList (list):
                The order of the units by name
        """
        sorter = dict(zip(relationshipList, range(len(relationshipList))))
        self.stratigraphicUnits["stratigraphic_Order"] = self.stratigraphicUnits["name"].map(sorter)
        self.stratigraphicUnits.sort_values(["stratigraphic_Order"], inplace=True)
