import pandas
import numpy
import beartype
import geopandas
import math

from .utils import calculate_minimum_fault_length


from .logging import getLogger

logger = getLogger(__name__)


class DeformationHistory:
    """
    A class containing all the fault and fold summaries and relationships

    Attributes
    ----------
    history: list
        The time ordered list of deformation events
    faultColumns: numpy.dtype
        Column names and types for fault summary
    foldColumns: numpy.dtype
        Column names and types for fold summary
    faults: pandas.DataFrame
        The fault summary
    folds: pandas.DataFrame
        The fold summary

    """

    def __init__(self, project):
        """
        The initialiser for the deformation history. All attributes are defaulted
        """
        self.history = []
        self.fault_fault_relationships = []
        self.project = project
        
        # Create empty fault and fold dataframes
        self.faultColumns = numpy.dtype(
            [
                ("eventId", int),
                ("name", str),
                ("minAge", float),
                ("maxAge", float),
                ("group", str),
                ("supergroup", str),
                ("avgDisplacement", float),
                ("avgDownthrowDir", float),
                ("influenceDistance", float),
                ("verticalRadius", float),
                ("horizontalRadius", float),
                ("colour", str),
                ("centreX", float),
                ("centreY", float),
                ("centreZ", float),
                ("avgSlipDirX", float),
                ("avgSlipDirY", float),
                ("avgSlipDirZ", float),
                ("avgNormalX", float),
                ("avgNormalY", float),
                ("avgNormalZ", float),
                ("length", float),
            ]
        )
        self.faults = pandas.DataFrame(numpy.empty(0, dtype=self.faultColumns))
        # self.faults = self.faults.set_index("name")
        
        self.foldColumns = numpy.dtype(
            [
                ("eventId", int),
                ("name", str),
                ("minAge", float),
                ("maxAge", float),
                ("periodic", bool),
                ("wavelength", float),
                ("amplitude", float),
                ("asymmetry", bool),
                ("asymmetryShift", float),
                ("secondaryWavelength", float),
                ("secondaryAmplitude", float),
            ]
        )
        self.folds = pandas.DataFrame(numpy.empty(0, dtype=self.foldColumns))
        # self.folds = self.folds.set_index("name")
        
        
        
    def findfault(self, id):
        """
        Find the fault in the summary based on its eventId

        Args:
            id (int or str):
                The eventId or name to look for

        Returns:
            pandas.DataFrame: The sliced data frame containing the requested fault
        """
        if issubclass(type(id), int):
            logger.info(f"Finding fault with eventId {id}")
            return self.faults[self.faults["eventId"] == id]
        elif issubclass(type(id), str):
            logger.info(f"Finding fault with name {id}")
            return self.faults[self.faults["name"] == id]
        else:
            logger.error("ERROR: Unknown identifier type used to find fault")

    def findfold(self, id):
        """
        Find the fold in the summary based on its eventId

        Args:
            id (int or str):
                The eventId or name to look for

        Returns:
            pandas.DataFrame: The sliced data frame containing the requested fold
        """
        if issubclass(type(id), int):
            logger.info(f"Finding fold with eventId {id}")
            return self.folds[self.folds["foldId"] == id]
        elif issubclass(type(id), str):
            logger.info(f"Finding fold with name {id}")
            return self.folds[self.folds["name"] == id]
        else:
            logger.error("ERROR: Unknown identifier type used to find fold")

    def addFault(self, fault):
        """
        Add fault to the fault summary

        Args:
            fault (pandas.DataFrame or dict):
                The fault information to add
        """
        if issubclass(type(fault), pandas.DataFrame) or issubclass(type(fault), dict):
            if "name" in fault.keys():
                if fault["name"] in self.faults.index:
                    logger.warning("Replacing fault", fault["name"])
                self.faults[fault["name"]] = fault
                logger.info("Adding fault", fault["name"])
            else:
                logger.error("No name field in fault", fault)
        else:
            logger.error("Cannot add fault to dataframe with type", type(fault))

    def removeFaultByName(self, name: str):
        """
        Remove the fault from the summary by name

        Args:
            name (str):
                The name of the fault(s) to remove
        """
        logger.info(f"Removing fault with name {name}")
        self.faults = self.faults.drop[self.faults.index[self.faults["name"] != name]]

    def removeFaultByEventId(self, eventId: int):
        """
        Remove the fault from the summary by eventId

        Args:
            eventId (int):
                The eventId of the fault to remove
        """
        logger.info(f"Removing fault with eventId {eventId}")
        self.faults = self.faults[self.faults["eventId"] != eventId].copy()

    def addFold(self, fold):
        """
        Add fold to the fold summary

        Args:
            fold (pandas.DataFrame or dict):
                The fold information to add
        """
        if issubclass(type(fold), pandas.DataFrame) or issubclass(type(fold), dict):
            if "name" in fold.keys():
                if fold["name"] in self.folds.index:
                    logger.warning("Replacing fold", fold["name"])
                logger.info("Adding fold", fold["name"])
                self.folds[fold["name"]] = fold
            else:
                logger.error("No name field in fold", fold)
        else:
            logger.error("Cannot add fold to dataframe with type", type(fold))

    @beartype.beartype
    def populate(self, faults_map_data: geopandas.GeoDataFrame):
        """
        Populate the fault (and fold) summaries from a geodataframe

        Args:
            faults_map_data (geopandas.GeoDataFrame):
                The parsed data frame from the map
        """
        logger.info("Populating fault/fold summary")
        if faults_map_data.shape[0] == 0:
            return
        faults_data = faults_map_data.copy()
        faults_data = faults_data.dissolve(by="NAME", as_index=False)
        faults_data = faults_data.reset_index(drop=True)

        self.stratigraphicUnits = pandas.DataFrame(
            numpy.empty(faults_data.shape[0], dtype=self.faultColumns)
        )
        self.faults["eventId"] = faults_data["ID"]
        self.faults["name"] = faults_data["NAME"]
        self.faults["minAge"] = -1.0
        self.faults["maxAge"] = -1.0
        self.faults["group"] = ""
        self.faults["supergroup"] = ""
        self.faults["avgDisplacement"] = -1.0
        self.faults["avgDownthrowDir"] = numpy.nan
        self.faults["influenceDistance"] = numpy.nan
        self.faults["verticalRadius"] = numpy.nan
        self.faults["horizontalRadius"] = numpy.nan
        self.faults["colour"] = "#000000"
        self.faults["centreX"] = numpy.nan
        self.faults["centreY"] = numpy.nan
        self.faults["centreZ"] = numpy.nan
        self.faults["avgSlipDirX"] = numpy.nan
        self.faults["avgSlipDirY"] = numpy.nan
        self.faults["avgSlipDirZ"] = numpy.nan
        self.faults["avgNormalX"] = numpy.nan
        self.faults["avgNormalY"] = numpy.nan
        self.faults["avgNormalZ"] = numpy.nan
        self.faults["length"] = faults_data.geometry.length
        for index, fault in self.faults.iterrows():
            bounds = faults_map_data[faults_map_data["ID"] == fault["eventId"]].geometry.bounds
            xdist = float(bounds.maxx.iloc[0] - bounds.minx.iloc[0])
            ydist = float(bounds.maxy.iloc[0] - bounds.miny.iloc[0])
            length = math.sqrt(xdist * xdist + ydist * ydist)
            self.faults.at[index, "verticalRadius"] = length
            self.faults.at[index, "horizontalRadius"] = length / 2.0
            self.faults.at[index, "influenceDistance"] = length / 4.0

    @beartype.beartype
    def summarise_data(self, fault_observations: pandas.DataFrame):
        """
        Use fault observations data to add summary data for each fault

        Args:
            fault_observations (pandas.DataFrame):
                The fault observations data
        """
        logger.info("Summarising fault data")
        id_list = self.faults["eventId"].unique()
        for id in id_list:
            observations = fault_observations[fault_observations["ID"] == id]
            if len(observations) < 2:
                self.removeFaultByEventId(id)

        # id_list = self.faults["eventId"].unique()
        for index, fault in self.faults.iterrows():
            observations = fault_observations[fault_observations["ID"] == fault["eventId"]]
            # calculate centre point
            self.faults.at[index, "centreX"] = numpy.mean(observations["X"])
            self.faults.at[index, "centreY"] = numpy.mean(observations["Y"])
            self.faults.at[index, "centreZ"] = numpy.mean(observations["Z"])
    
    
    def get_faults_for_export(self):
        """
        Get the faults for export (removes any fault that is shorter than the cutoff)

        Returns:
            pandas.DataFrame: The filtered fault summary
        """
        # if no minimum fault length is set, calculate it
        if self.project.get_minimum_fault_length() < 0:
            self.project.set_minimum_fault_length( calculate_minimum_fault_length(
                bbox=self.project.bounding_box, area_percentage=0.05
            ))
        return self.faults[self.faults["length"] >= self.project.get_minimum_fault_length()].copy()

    @beartype.beartype
    def get_fault_relationships_with_ids(self, fault_fault_relationships: pandas.DataFrame):
        """
        Ammend the fault relationships DataFrame with the fault eventIds

        Args:
            fault_fault_relationships (pandas.DataFrame): The fault_fault_relationships

        Returns:
            pandas.DataFrame: The fault_relationships with the correct eventIds
        """
        logger.info("Getting fault relationships with eventIds")

        faultIds = self.get_faults_for_export()[["eventId", "name"]].copy()
        rel = fault_fault_relationships.copy()
        rel = rel.merge(faultIds, left_on="Fault1", right_on="eventId")
        rel.rename(columns={"eventId": "eventId1"}, inplace=True)
        rel.drop(columns=["name"], inplace=True)
        rel = rel.merge(faultIds, left_on="Fault2", right_on="eventId")
        rel.rename(columns={"eventId": "eventId2"}, inplace=True)
        rel.drop(columns=["name"], inplace=True)
        return rel
