from abc import ABC, abstractmethod
import copy
import beartype
import pandas
import numpy as np
import math
from typing import Callable, Dict, Union, Optional, List
from map2loop.topology import Topology
import geopandas
from osgeo import gdal
from map2loop.utils import value_from_raster
from .logging import getLogger
import networkx as nx

logger = getLogger(__name__)


class Sorter(ABC):
    """
    Base Class of Sorter used to force structure of Sorter

    Args:
        ABC (ABC): Derived from Abstract Base Class
    """

    def __init__(self):
        """
        Initialiser for Sorter

        Args:
            unit_relationships (pandas.DataFrame): the relationships between units (columns must contain ["Index1", "Unitname1", "Index2", "Unitname2"])
            contacts (pandas.DataFrame): unit contacts with length of the contacts in metres
            geology_data (geopandas.GeoDataFrame): the geology data
            structure_data (geopandas.GeoDataFrame): the structure data
            dtm_data (gdal.Dataset): the dtm data
        """
        self.sorter_label = "SorterBaseClass"

    def type(self):
        """
        Getter for subclass type label

        Returns:
            str: Name of subclass
        """
        return self.sorter_label

    @beartype.beartype
    @abstractmethod
    def sort(self, units: pandas.DataFrame) -> list:
        """
        Execute sorter method (abstract method)

        Args:
            units (pandas.DataFrame): the data frame to sort (columns must contain ["layerId", "name", "minAge", "maxAge", "group"])

        Returns:
            list: sorted list of unit names
        """
        pass

    def __call__(self, **kwargs):
        return self.sort(**kwargs)


class SorterUseNetworkX(Sorter):
    """
    Sorter class which returns a sorted list of units based on the unit relationships using a topological graph sorting algorithm
    """

    required_arguments: List[str] = ['geology_data', 'unit_name_column']

    def __init__(
        self,
        *,
        unit_name_column: Optional[str] = 'name',
        unit_relationships: Optional[pandas.DataFrame] = None,
        geology_data: Optional[geopandas.GeoDataFrame] = None,
    ):
        """
        Initialiser for networkx graph sorter

        Args:
            unit_relationships (pandas.DataFrame): the relationships between units
        """
        super().__init__()
        self.sorter_label = "SorterUseNetworkX"
        self.unit_name_column = unit_name_column
        if geology_data is not None:
            self.set_geology_data(geology_data)
        elif unit_relationships is not None:
            self.unit_relationships = unit_relationships
        else:
            self.unit_relationships = None

    def set_geology_data(self, geology_data: geopandas.GeoDataFrame):
        """
        Set geology data and calculate topology and unit relationships

        Args:
            geology_data (geopandas.GeoDataFrame): the geology data
        """
        self._calculate_topology(geology_data)

    def _calculate_topology(self, geology_data: geopandas.GeoDataFrame):
        if geology_data is None:
            raise ValueError("geology_data is required")

        if isinstance(geology_data, geopandas.GeoDataFrame) is False:
            raise TypeError("geology_data must be a geopandas.GeoDataFrame")

        if 'UNITNAME' not in geology_data.columns:
            raise ValueError("geology_data must contain 'UNITNAME' column")

        self.topology = Topology(geology_data=geology_data)
        self.unit_relationships = self.topology.get_unit_unit_relationships()

    @beartype.beartype
    def sort(self, units: pandas.DataFrame) -> list:
        """
        Execute sorter method takes unit data and returns the sorted unit names based on this algorithm.

        Args:
            units (pandas.DataFrame): the data frame to sort

        Returns:
            list: the sorted unit names
        """
        import networkx as nx

        if self.unit_relationships is None:
            raise ValueError("SorterUseNetworkX requires 'unit_relationships' argument")
        graph = nx.DiGraph()
        name_to_index = {}
        for row in units.iterrows():
            graph.add_node(int(row[1]["layerId"]), name=row[1]["name"])
            name_to_index[row[1]["name"]] = int(row[1]["layerId"])
        for row in self.unit_relationships.iterrows():
            graph.add_edge(name_to_index[row[1]["UNITNAME_1"]], name_to_index[row[1]["UNITNAME_2"]])

        cycles = list(nx.simple_cycles(graph))
        for i in range(0, len(cycles)):
            if graph.has_edge(cycles[i][0], cycles[i][1]):
                graph.remove_edge(cycles[i][0], cycles[i][1])
                logger.warning(
                    " SorterUseNetworkX: Cycle found and contact edge removed:",
                    units["name"][cycles[i][0]],
                    units["name"][cycles[i][1]],
                )

        indexes = list(nx.topological_sort(graph))
        order = [units["name"][i] for i in list(indexes)]
        logger.info("Stratigraphic order calculated using networkx topological sort")
        logger.info(','.join(order))
        return order


class SorterUseHint(SorterUseNetworkX):
    required_arguments: List[str] = ['unit_relationships']

    def __init__(self, *, geology_data: Optional[geopandas.GeoDataFrame] = None):
        logger.warning("SorterUseHint is deprecated in v3.2. Using SorterUseNetworkX instead")
        super().__init__(geology_data=geology_data)


class SorterAgeBased(Sorter):
    """
    Sorter class which returns a sorted list of units based on the min and max ages of the units
    """

    required_arguments = ['min_age_column', 'max_age_column', 'unit_name_column']

    def __init__(
        self,
        *,
        unit_name_column: Optional[str] = 'name',
        min_age_column: Optional[str] = 'minAge',
        max_age_column: Optional[str] = 'maxAge',
    ):
        """
        Initialiser for age based sorter
        """
        super().__init__()
        self.unit_name_column = unit_name_column
        self.min_age_column = min_age_column
        self.max_age_column = max_age_column
        self.sorter_label = "SorterAgeBased"

    def sort(self, units: pandas.DataFrame) -> list:
        """
        Execute sorter method takes unit data and returns the sorted unit names based on this algorithm.

        Args:
            units (pandas.DataFrame): the data frame to sort

        Returns:
            list: the sorted unit names
        """
        logger.info("Calling age based sorter")
        sorted_units = units.copy()

        if self.min_age_column in units.columns and self.max_age_column in units.columns:
            # print(sorted_units["minAge"], sorted_units["maxAge"])
            sorted_units["meanAge"] = sorted_units.apply(
                lambda row: (row[self.min_age_column] + row[self.max_age_column]) / 2.0, axis=1
            )
        else:
            logger.error(
                f"Columns {self.min_age_column} and {self.max_age_column} must be present in units DataFrame"
            )
            logger.error(f"Available columns are: {units.columns.tolist()}")
            raise ValueError(
                f"Columns {self.min_age_column} and {self.max_age_column} must be present in units DataFrame"
            )
        if "group" in units.columns:
            sorted_units = sorted_units.sort_values(by=["group", "meanAge"])
        else:
            sorted_units = sorted_units.sort_values(by=["meanAge"])
        logger.info("Stratigraphic order calculated using age based sorting")
        for _i, row in sorted_units.iterrows():
            logger.info(
                f"{row[self.unit_name_column]} - {row[self.min_age_column]} - {row[self.max_age_column]}"
            )

        return list(sorted_units[self.unit_name_column])


class SorterAlpha(Sorter):
    """
    Sorter class which returns a sorted list of units based on the adjacency of units
    prioritising the units with lower number of contacting units
    """

    required_arguments = ['contacts', 'unit_name_column', 'unitname1_column', 'unitname2_column']

    def __init__(
        self,
        *,
        contacts: Optional[geopandas.GeoDataFrame] = None,
        unit_name_column: Optional[str] = 'name',
        unitname1_column: Optional[str] = 'UNITNAME_1',
        unitname2_column: Optional[str] = 'UNITNAME_2',
    ):
        """
        Initialiser for adjacency based sorter

        Args:
            contacts (geopandas.GeoDataFrame): unit contacts with length of the contacts in metres
        """
        super().__init__()
        self.contacts = contacts
        self.unit_name_column = unit_name_column
        self.sorter_label = "SorterAlpha"
        self.unitname1_column = unitname1_column
        self.unitname2_column = unitname2_column
        if (
            self.unitname1_column not in contacts.columns
            or self.unitname2_column not in contacts.columns
            or 'length' not in contacts.columns
        ):
            raise ValueError(
                f"contacts GeoDataFrame must contain '{self.unitname1_column}', '{self.unitname2_column}' and 'length' columns"
            )

    def sort(self, units: pandas.DataFrame) -> list:
        """
        Execute sorter method takes unit data and returns the sorted unit names based on this algorithm.

        Args:
            units (pandas.DataFrame): the data frame to sort

        Returns:
            list: the sorted unit names
        """
        if self.contacts is None:
            raise ValueError(
                "contacts must be set (not None) before calling sort() in SorterAlpha."
            )
        if len(self.contacts) == 0:
            raise ValueError("contacts GeoDataFrame is empty in SorterAlpha.")
        if 'length' not in self.contacts.columns:
            self.contacts['length'] = self.contacts.geometry.length
        self.contacts['length'] = self.contacts['length'].astype(float)
        sorted_contacts = self.contacts.sort_values(by="length", ascending=False)[
            [self.unitname1_column, self.unitname2_column, "length"]
        ]
        unit_names = list(units[self.unit_name_column].unique())
        graph = nx.Graph()
        for unit in unit_names:
            graph.add_node(unit, name=unit)
        max_weight = max(list(sorted_contacts["length"])) + 1
        for _, row in sorted_contacts.iterrows():
            graph.add_edge(
                row[self.unitname1_column],
                row[self.unitname2_column],
                weight=int(max_weight - row["length"]),
            )

        cnode = None
        new_graph = nx.DiGraph()
        while graph.number_of_nodes() > 0:
            if cnode is None:
                df = pandas.DataFrame(columns=["unit", "num_neighbours"])
                df["unit"] = list(graph.nodes)
                df["num_neighbours"] = df.apply(
                    lambda row: len(list(graph.neighbors(row["unit"]))), axis=1
                )
                df.sort_values(by=["num_neighbours"], inplace=True)
                df.reset_index(inplace=True, drop=True)
                cnode = df["unit"][0]
                new_graph.add_node(cnode)
            neighbour_edge_count = {}
            for neighbour in list(graph.neighbors(cnode)):
                neighbour_edge_count[neighbour] = len(list(graph.neighbors(neighbour)))
            if len(neighbour_edge_count) == 0:
                graph.remove_node(cnode)
                cnode = None
            else:
                node_with_min_edges = min(neighbour_edge_count, key=neighbour_edge_count.get)
                if neighbour_edge_count[node_with_min_edges] < 2:
                    new_graph.add_node(node_with_min_edges)
                    new_graph.add_edge(cnode, node_with_min_edges)
                    graph.remove_node(node_with_min_edges)
                else:
                    new_graph.add_node(node_with_min_edges)
                    new_graph.add_edge(cnode, node_with_min_edges)
                    graph.remove_node(cnode)
                    cnode = node_with_min_edges
        order = list(reversed(list(nx.topological_sort(new_graph))))
        logger.info("Stratigraphic order calculated using adjacency based sorting")
        logger.info(','.join(order))
        return order


class SorterMaximiseContacts(Sorter):
    """
    Sorter class which returns a sorted list of units based on the adjacency of units
    prioritising the maximum length of each contact
    """

    required_arguments = ['contacts', 'unit_name_column', 'unitname1_column', 'unitname2_column']

    def __init__(
        self,
        *,
        contacts: Optional[geopandas.GeoDataFrame] = None,
        unit_name_column: str = 'name',
        unitname1_column: str = 'UNITNAME_1',
        unitname2_column: str = 'UNITNAME_2',
    ):
        """
        Initialiser for adjacency based sorter

        Args:
            contacts (pandas.DataFrame): unit contacts with length of the contacts in metres
        """
        super().__init__()
        self.sorter_label = "SorterMaximiseContacts"
        # variables for visualising/interrogating the sorter
        self.graph = None
        self.route = None
        self.directed_graph = None
        self.contacts = contacts
        self.unit_name_column = unit_name_column
        self.unitname1_column = unitname1_column
        self.unitname2_column = unitname2_column
        if (
            self.unitname1_column not in contacts.columns
            or self.unitname2_column not in contacts.columns
            or 'length' not in contacts.columns
        ):
            raise ValueError(
                f"contacts GeoDataFrame must contain '{self.unitname1_column}', '{self.unitname2_column}' and 'length' columns"
            )

    def sort(self, units: pandas.DataFrame) -> list:
        """
        Execute sorter method takes unit data and returns the sorted unit names based on this algorithm.

        Args:
            units (pandas.DataFrame): the data frame to sort

        Returns:
            list: the sorted unit names
        """
        import networkx as nx
        import networkx.algorithms.approximation as nx_app

        if self.contacts is None:
            raise ValueError("SorterMaximiseContacts requires 'contacts' argument")
        if len(self.contacts) == 0:
            raise ValueError("contacts GeoDataFrame is empty in SorterMaximiseContacts.")
        if "length" not in self.contacts.columns:
            self.contacts['length'] = self.contacts.geometry.length
        self.contacts['length'] = self.contacts['length'].astype(float)
        sorted_contacts = self.contacts.sort_values(by="length", ascending=False)
        self.graph = nx.Graph()
        unit_names = list(units[self.unit_name_column].unique())
        for unit in unit_names:
            ## some units may not have any contacts e.g. if they are intrusives or sills. If we leave this then the
            ## sorter crashes
            if (
                unit not in sorted_contacts[self.unitname1_column].values
                or unit not in sorted_contacts[self.unitname2_column].values
            ):
                continue
            self.graph.add_node(unit, name=unit)

        max_weight = max(list(sorted_contacts["length"])) + 1
        sorted_contacts['length'] /= max_weight
        for _, row in sorted_contacts.iterrows():
            self.graph.add_edge(
                row[self.unitname1_column], row[self.unitname2_column], weight=(1 - row["length"])
            )

        self.route = nx_app.traveling_salesman_problem(self.graph)
        edge_list = list(nx.utils.pairwise(self.route))
        self.directed_graph = nx.DiGraph()
        self.directed_graph.add_node(edge_list[0][0])
        for edge in edge_list:
            if edge[1] not in self.directed_graph.nodes():
                self.directed_graph.add_node(edge[1])
                self.directed_graph.add_edge(edge[0], edge[1])

        # we need to reverse the order of the graph to get the correct order
        order = list(
            reversed(
                list(
                    nx.dfs_preorder_nodes(
                        self.directed_graph, source=list(self.directed_graph.nodes())[0]
                    )
                )
            )
        )
        logger.info("Stratigraphic order calculated using adjacency based sorting")
        logger.info(','.join(order))
        return order


class SorterObservationProjections(Sorter):
    """
    Sorter class which returns a sorted list of units based on the adjacency of units
    using the direction of observations to predict which unit is adjacent to the current one
    """

    required_arguments = [
        'contacts',
        'geology_data',
        'structure_data',
        'dtm_data',
        'unit_name_column',
        'unitname1_column',
        'unitname2_column',
    ]

    def __init__(
        self,
        *,
        unitname1_column: Optional[str] = 'UNITNAME_1',
        unitname2_column: Optional[str] = 'UNITNAME_2',
        unit_name_column: Optional[str] = 'name',
        contacts: Optional[geopandas.GeoDataFrame] = None,
        geology_data: Optional[geopandas.GeoDataFrame] = None,
        structure_data: Optional[geopandas.GeoDataFrame] = None,
        dtm_data: Optional[gdal.Dataset] = None,
        length: Union[float, int] = 1000,
    ):
        """
        Initialiser for adjacency based sorter

        Args:
            contacts (pandas.DataFrame): unit contacts with length of the contacts in metres
            geology_data (geopandas.GeoDataFrame): the geology data
            structure_data (geopandas.GeoDataFrame): the structure data
            dtm_data (gdal.Dataset): the dtm data
            length (int): the length of the projection in metres
        """
        super().__init__()
        self.contacts = contacts
        self.geology_data = geology_data
        self.structure_data = structure_data
        self.dtm_data = dtm_data
        self.unit_name_column = unit_name_column
        self.sorter_label = "SorterObservationProjections"
        self.length = length
        self.lines = []
        self.unit1name_column = unitname1_column
        self.unit2name_column = unitname2_column

    def sort(self, units: pandas.DataFrame) -> list:
        """
        Execute sorter method takes unit data and returns the sorted unit names based on this algorithm.

        Args:
            units (pandas.DataFrame): the data frame to sort

        Returns:
            list: the sorted unit names
        """
        import networkx as nx
        import networkx.algorithms.approximation as nx_app
        from shapely.geometry import LineString, Point

        if self.contacts is None:
            raise ValueError("SorterObservationProjections requires 'contacts' argument")
        if self.geology_data is None:
            raise ValueError("SorterObservationProjections requires 'geology_data' argument")
        geol = self.geology_data.copy()
        if "INTRUSIVE" in geol.columns:
            geol = geol.drop(geol.index[geol["INTRUSIVE"]])
        if "SILL" in geol.columns:
            geol = geol.drop(geol.index[geol["SILL"]])
        if self.structure_data is None:
            raise ValueError("structure_data is required for sorting but is None.")
        orientations = self.structure_data.copy()
        if self.dtm_data is None:
            raise ValueError("DTM data (self.dtm_data) is not set. Cannot proceed with sorting.")
        inv_geotransform = gdal.InvGeoTransform(self.dtm_data.GetGeoTransform())
        dtm_array = np.array(self.dtm_data.GetRasterBand(1).ReadAsArray().T)

        # Create a map of maps to store younger/older observations
        ordered_unit_observations = []
        for _, row in orientations.iterrows():
            # get containing unit
            containing_unit = geol[geol.contains(row.geometry)]
            if len(containing_unit) > 1:
                logger.info(f"Orientation {row.ID} is within multiple units")
                logger.info(f"Check geology map around coordinates {row.geometry}")

            if len(containing_unit) < 1:
                logger.info(f"Orientation {row.ID} is not in a unit")
                logger.info(f"Check geology map around coordinates {row.geometry}")
            else:
                first_unit_name = containing_unit.iloc[0]["UNITNAME"]
                # Get units that a projected line passes through
                length = self.length
                dipDirRadians = row.DIPDIR * math.pi / 180.0
                dipRadians = row.DIP * math.pi / 180.0
                start = row.geometry
                end = Point(
                    start.x + math.sin(dipDirRadians) * length,
                    start.y + math.cos(dipDirRadians) * length,
                )
                line = LineString([start, end])
                self.lines.append(line)
                inter = geol[line.intersects(geol.geometry)]

                if len(inter) > 1:
                    intersect = line.intersection(inter.geometry.boundary)
                    # # Remove containing unit
                    intersect = intersect.drop(containing_unit.index)

                    # sort by distance from start point
                    sub = geol.loc[intersect.index].copy()
                    sub["distance"] = geol.distance(start)
                    sub = sub.sort_values(by="distance")

                    # Get first unit it hits and the point of intersection
                    second_unit_name = sub.iloc[0].UNITNAME

                    if intersect.loc[sub.index[0]].geom_type == "MultiPoint":
                        second_intersect_point = intersect.loc[sub.index[0]].geoms[0]
                    elif intersect.loc[sub.index[0]].geom_type == "Point":
                        second_intersect_point = intersect.loc[sub.index[0]]
                    else:
                        continue

                    # Get heights for intersection point and start of ray
                    height = value_from_raster(inv_geotransform, dtm_array, start.x, start.y)
                    first_intersect_point = Point(start.x, start.y, height)
                    height = value_from_raster(
                        inv_geotransform,
                        dtm_array,
                        second_intersect_point.x,
                        second_intersect_point.y,
                    )
                    second_intersect_point = Point(second_intersect_point.x, start.y, height)

                    # Check vertical difference between points and compare to projected dip angle
                    horizontal_dist = (
                        first_intersect_point.x - first_intersect_point.x,
                        second_intersect_point.y - first_intersect_point.y,
                    )
                    horizontal_dist = math.sqrt(horizontal_dist[0] ** 2 + horizontal_dist[1] ** 2)
                    projected_height = first_intersect_point.z + horizontal_dist * math.cos(
                        dipRadians
                    )

                    if second_intersect_point.z < projected_height:
                        ordered_unit_observations += [(first_unit_name, second_unit_name)]
                    else:
                        ordered_unit_observations += [(second_unit_name, first_unit_name)]
        self.ordered_unit_observations = ordered_unit_observations
        # Create a matrix of older versus younger frequency from observations
        unit_names = geol.UNITNAME.unique()
        df = pandas.DataFrame(0, index=unit_names, columns=unit_names)
        for younger, older in ordered_unit_observations:
            df.loc[younger, older] += 1
        print(df, df.max())
        max_value = max(df.max())

        # Using the older/younger matrix create a directed graph
        g = nx.DiGraph()
        remaining_units = unit_names
        for unit1 in unit_names:
            g.add_node(unit1)
        for unit1 in unit_names:
            remaining_units = remaining_units[1:]
            for unit2 in remaining_units:
                if unit1 != unit2:
                    weight = df.loc[unit1, unit2] - df.loc[unit2, unit1]
                    if weight < 0:
                        g.add_edge(unit1, unit2, weight=max_value + weight)
                    elif weight > 0:
                        g.add_edge(unit2, unit1, weight=max_value - weight)
                    if df.loc[unit1, unit2] > 0 and df.loc[unit2, unit1] > 0 and weight == 0:
                        # if both units have the same weight add a bidirectional edge
                        pass
                        print('')
                        g.add_edge(unit2, unit1, weight=max_value)
                        g.add_edge(unit1, unit2, weight=max_value)
        self.G = g
        # Link in unlinked units from contacts with max weight
        g_undirected = g.to_undirected()
        for unit in unit_names:
            if len(list(g_undirected.neighbors(unit))) < 1:
                mask1 = self.contacts[self.unit1name_column] == unit
                mask2 = self.contacts[self.unit2name_column] == unit
                for _, row in self.contacts[mask1 | mask2].iterrows():
                    if unit == row[self.unit1name_column]:
                        g.add_edge(row[self.unit2name_column], unit, weight=max_value * 10)
                    else:
                        g.add_edge(row[self.unit1name_column], unit, weight=max_value * 10)

        # Run travelling salesman using the observation evidence as weighting
        route = nx_app.traveling_salesman_problem(g.to_undirected())
        self.route = route
        edge_list = list(nx.utils.pairwise(route))
        self.edge_list = edge_list
        dd = nx.DiGraph()
        dd.add_node(edge_list[0][0])
        for edge in edge_list:
            if edge[1] not in dd.nodes():
                dd.add_node(edge[1])
                dd.add_edge(edge[0], edge[1])
        self.directed = dd
        logger.info("Stratigraphic order calculated using observation based sorting")
        order = list(nx.dfs_preorder_nodes(dd, source=list(dd.nodes())[0]))
        logger.info(','.join(order))
        return order


def _is_missing(value) -> bool:
    """
    Check if a group or supergroup value is empty

    Args:
        value: the value from the group or supergroup column

    Returns:
        bool: True if the value is None, NaN, an empty string or "None"/"nan"
    """
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return str(value).strip() in ("", "None", "nan")


def relabel_contacts(
    contacts: pandas.DataFrame,
    labels: Dict[str, str],
    unitname1_column: str = 'UNITNAME_1',
    unitname2_column: str = 'UNITNAME_2',
) -> pandas.DataFrame:
    """
    Change the unit names of the contacts to labels (for example the group of each unit)

    A contact with a unit that is not in labels is removed. A contact between two units
    with the same label is removed. The lengths of the contacts between the same two
    labels are added together.

    Args:
        contacts (pandas.DataFrame): the contacts, with a 'length' column or a geometry
        labels (dict): the label of each unit name
        unitname1_column (str): the name of the column with the first unit name
        unitname2_column (str): the name of the column with the second unit name

    Returns:
        pandas.DataFrame: the contacts between the labels, with the columns
        [unitname1_column, unitname2_column, 'length']
    """
    columns = [unitname1_column, unitname2_column, 'length']
    if contacts is None or len(contacts) == 0:
        return pandas.DataFrame(columns=columns)
    if 'length' in contacts.columns:
        lengths = contacts['length'].astype(float)
    elif isinstance(contacts, geopandas.GeoDataFrame):
        lengths = contacts.geometry.length
    else:
        lengths = pandas.Series(1.0, index=contacts.index)
    label1 = contacts[unitname1_column].map(labels)
    label2 = contacts[unitname2_column].map(labels)
    keep = label1.notna() & label2.notna() & (label1 != label2)
    pairs = pandas.DataFrame(
        {unitname1_column: label1[keep], unitname2_column: label2[keep], 'length': lengths[keep]}
    )
    if len(pairs) == 0:
        return pandas.DataFrame(columns=columns)
    # (a, b) and (b, a) are the same contact
    swap = pairs[unitname1_column] > pairs[unitname2_column]
    pairs.loc[swap, [unitname1_column, unitname2_column]] = pairs.loc[
        swap, [unitname2_column, unitname1_column]
    ].values
    return pairs.groupby([unitname1_column, unitname2_column], as_index=False)['length'].sum()


def relabel_unit_relationships(
    unit_relationships: pandas.DataFrame, labels: Dict[str, str]
) -> pandas.DataFrame:
    """
    Change the unit names of the unit relationships to labels (for example the group of each unit)

    A relationship with a unit that is not in labels is removed. A relationship between
    two units with the same label is removed. The direction of each relationship is kept.

    Args:
        unit_relationships (pandas.DataFrame): the relationships, with the columns
            'UNITNAME_1' and 'UNITNAME_2'
        labels (dict): the label of each unit name

    Returns:
        pandas.DataFrame: the relationships between the labels
    """
    columns = ['UNITNAME_1', 'UNITNAME_2']
    if unit_relationships is None or len(unit_relationships) == 0:
        return pandas.DataFrame(columns=columns)
    label1 = unit_relationships['UNITNAME_1'].map(labels)
    label2 = unit_relationships['UNITNAME_2'].map(labels)
    keep = label1.notna() & label2.notna() & (label1 != label2)
    relationships = pandas.DataFrame({'UNITNAME_1': label1[keep], 'UNITNAME_2': label2[keep]})
    return relationships.drop_duplicates().reset_index(drop=True)


class SorterHierarchical(Sorter):
    """
    Sorter class which keeps the units of each supergroup and of each group together

    This sorter uses a different sorter (for example SorterAlpha) at each level:
    1. It sorts the supergroups.
    2. It sorts the groups in each supergroup.
    3. It sorts the units in each group.

    To sort the supergroups (or the groups), the data of the sorter (contacts, unit
    relationships and geology) is changed so that each supergroup (or group) is one
    unit. At each step, the sorter uses only the data of the units in that step.
    Thus each supergroup has its own stratigraphic order, and the contacts between
    two supergroups only set the order of the two supergroups.

    A group with no supergroup is one item at the supergroup level, the same as a
    supergroup. A unit with no group (and no supergroup) is one item at that level.
    """

    required_arguments: List[str] = ['sorter']

    def __init__(
        self,
        *,
        sorter: Sorter,
        group_column: Optional[str] = 'group',
        supergroup_column: Optional[str] = 'supergroup',
        postprocess: Optional[Callable[[list, Dict[str, str]], list]] = None,
    ):
        """
        Initialiser for hierarchical sorter

        Args:
            sorter (Sorter): the sorter to use at each level
            group_column (str, optional): the column of the units with the group of each
                unit. Set to None to not use groups. Defaults to 'group'.
            supergroup_column (str, optional): the column of the units with the
                supergroup of each unit. Set to None to not use supergroups.
                Defaults to 'supergroup'.
            postprocess (callable, optional): a function that is applied to the result
                of each step, postprocess(order, labels) -> order. labels is the label of
                each unit name in that step (the group, the supergroup or the unit name).
                For example, use it to repair the closed route of a travelling salesman
                sorter. Defaults to None.
        """
        super().__init__()
        if isinstance(sorter, SorterHierarchical):
            raise TypeError("sorter must not be a SorterHierarchical")
        self.sorter = sorter
        self.group_column = group_column
        self.supergroup_column = supergroup_column
        self.postprocess = postprocess
        self.unit_name_column = getattr(sorter, 'unit_name_column', None) or 'name'
        self.sorter_label = f"SorterHierarchical({sorter.sorter_label})"

    def sort(self, units: pandas.DataFrame) -> list:
        """
        Execute sorter method takes unit data and returns the sorted unit names based on this algorithm.

        Args:
            units (pandas.DataFrame): the data frame to sort

        Returns:
            list: the sorted unit names
        """
        if self.unit_name_column not in units.columns:
            raise ValueError(f"Column {self.unit_name_column} must be present in units DataFrame")
        units = units.drop_duplicates(subset=[self.unit_name_column]).reset_index(drop=True)
        levels = [
            column
            for column in (self.supergroup_column, self.group_column)
            if column and column in units.columns
        ]
        if not levels:
            logger.warning(
                f"{self.sorter_label}: no group or supergroup column in the units, "
                "so the units are sorted with no hierarchy"
            )
        order = self._sort_level(units, levels)
        logger.info("Stratigraphic order calculated using hierarchical sorting")
        logger.info(','.join(order))
        return order

    def _sort_level(self, units: pandas.DataFrame, levels: List[str]) -> list:
        """
        Sort the units with the first level in levels, then each part with the next levels

        Args:
            units (pandas.DataFrame): the units to sort
            levels (list): the group columns to use, from the highest level

        Returns:
            list: the sorted unit names
        """
        names = list(units[self.unit_name_column])
        if not levels:
            return self._sort_step(self._units_table(units), {name: name for name in names})
        column, next_levels = levels[0], levels[1:]

        # Find the label of each unit at this level. A unit with no value in this
        # column uses the value of the next level (or its name), so that a group
        # with no supergroup is one item at the supergroup level.
        level_values = {
            str(value).strip() for value in units[column] if not _is_missing(value)
        }
        labels = {}
        parts = {}
        for _, row in units.iterrows():
            name = row[self.unit_name_column]
            label_column, label = None, name
            for level in levels:
                if not _is_missing(row[level]):
                    label_column, label = level, str(row[level]).strip()
                    break
            if label_column != column and label in level_values:
                label = f"{label} ({label_column or 'unit'})"
            labels[name] = label
            parts.setdefault(label, []).append(name)
        if len(parts) == 1:
            return self._sort_level(units, next_levels)

        label_order = self._sort_step(self._label_table(units, labels), labels)
        order = []
        for label in label_order:
            part = units[units[self.unit_name_column].isin(parts[label])]
            order += self._sort_level(part, next_levels)
        return order

    def _units_table(self, units: pandas.DataFrame) -> pandas.DataFrame:
        """
        Make the units table for the sorter, for the units in one step

        Args:
            units (pandas.DataFrame): the units

        Returns:
            pandas.DataFrame: the units, with 'layerId' equal to the index
        """
        hierarchy_columns = [
            column for column in (self.group_column, self.supergroup_column) if column
        ]
        table = units.drop(columns=[c for c in hierarchy_columns if c in units.columns])
        table = table.reset_index(drop=True)
        # SorterUseNetworkX reads units["name"][layerId]
        table['layerId'] = table.index
        if 'name' not in table.columns:
            table['name'] = table[self.unit_name_column]
        return table

    def _label_table(self, units: pandas.DataFrame, labels: Dict[str, str]) -> pandas.DataFrame:
        """
        Make a units table for the sorter where each label is one unit

        The minimum age of a label is the minimum of the minimum ages of its units and
        the maximum age is the maximum of the maximum ages.

        Args:
            units (pandas.DataFrame): the units
            labels (dict): the label of each unit name

        Returns:
            pandas.DataFrame: one row for each label
        """
        unit_labels = units[self.unit_name_column].map(labels)
        table = pandas.DataFrame({self.unit_name_column: list(dict.fromkeys(unit_labels))})
        for attribute, default, aggregate in (
            ('min_age_column', 'minAge', 'min'),
            ('max_age_column', 'maxAge', 'max'),
        ):
            column = getattr(self.sorter, attribute, None) or default
            if column in units.columns:
                ages = pandas.to_numeric(units[column], errors='coerce')
                ages = ages.groupby(unit_labels).agg(aggregate)
                table[column] = table[self.unit_name_column].map(ages)
        return self._units_table(table)

    def _relabelled_sorter(self, labels: Dict[str, str]) -> Sorter:
        """
        Make a copy of the sorter that uses only the units in labels, with their labels as the unit names

        Args:
            labels (dict): the label of each unit name

        Returns:
            Sorter: the copy of the sorter
        """
        sorter = copy.copy(self.sorter)
        contacts = getattr(sorter, 'contacts', None)
        if contacts is not None:
            unitname1_column = (
                getattr(sorter, 'unitname1_column', None)
                or getattr(sorter, 'unit1name_column', None)
                or 'UNITNAME_1'
            )
            unitname2_column = (
                getattr(sorter, 'unitname2_column', None)
                or getattr(sorter, 'unit2name_column', None)
                or 'UNITNAME_2'
            )
            sorter.contacts = relabel_contacts(
                contacts, labels, unitname1_column, unitname2_column
            )
        unit_relationships = getattr(sorter, 'unit_relationships', None)
        if unit_relationships is not None:
            sorter.unit_relationships = relabel_unit_relationships(unit_relationships, labels)
        geology_data = getattr(sorter, 'geology_data', None)
        if geology_data is not None and 'UNITNAME' in geology_data.columns:
            geology_data = geology_data[geology_data['UNITNAME'].isin(list(labels))].copy()
            geology_data['UNITNAME'] = geology_data['UNITNAME'].map(labels)
            sorter.geology_data = geology_data.reset_index(drop=True)
        # copy.copy shares lists, so do not add to the list of the original sorter
        if isinstance(getattr(sorter, 'lines', None), list):
            sorter.lines = []
        return sorter

    def _sort_step(self, table: pandas.DataFrame, labels: Dict[str, str]) -> list:
        """
        Sort the labels in the table with a relabelled copy of the sorter

        If the sorter fails, the labels are kept in the order of the table. If the sorter
        does not return some labels, they are added at the end.

        Args:
            table (pandas.DataFrame): one row for each label
            labels (dict): the label of each unit name

        Returns:
            list: the sorted labels
        """
        expected = list(table[self.unit_name_column])
        if len(expected) < 2:
            return expected
        try:
            order = self._relabelled_sorter(labels).sort(table)
        except Exception as e:
            logger.warning(
                f"{self.sorter_label}: could not sort {expected} ({e}). "
                "The order of these is not changed."
            )
            return expected
        expected_set = set(expected)
        order = [label for label in dict.fromkeys(order) if label in expected_set]
        missing = [label for label in expected if label not in order]
        if missing:
            logger.warning(
                f"{self.sorter_label}: the sorter did not give a position for {missing}. "
                "They are added at the end."
            )
            order += missing
        if self.postprocess is not None:
            order = list(self.postprocess(order, labels))
        return order
