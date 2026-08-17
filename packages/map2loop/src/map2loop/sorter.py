from abc import ABC, abstractmethod
import beartype
import pandas
import numpy as np
import math
from typing import Union, Optional, List
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
