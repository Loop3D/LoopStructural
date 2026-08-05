# internal imports
from .utils import (
    create_points,
    rebuild_sampled_basal_contacts,
    calculate_endpoints,
    multiline_to_line,
    find_segment_strike_from_pt,
    set_z_values_from_raster_df,
    value_from_raster,
    segment_measure_range,
    clean_line_geometry,
    nearest_orientation_to_line,
    iter_line_segments
)
from .interpolators import DipDipDirectionInterpolator

from .logging import getLogger
logger = getLogger(__name__)  

# external imports
from abc import ABC, abstractmethod
from typing import Optional, List
from collections import defaultdict
import scipy.interpolate
import beartype
import numpy
import scipy
from scipy.spatial import cKDTree
import pandas
import geopandas
import shapely
from shapely.geometry import LineString
import math
from osgeo import gdal
from shapely.errors import UnsupportedGEOSVersionError

class ThicknessCalculator(ABC):
    """
    Base Class of Thickness Calculator used to force structure of ThicknessCalculator

    Args:
        ABC (ABC): Derived from Abstract Base Class
    """

    def __init__(
        self, 
        dtm_data: Optional[gdal.Dataset] = None, 
        bounding_box: Optional[dict] = None, 
        max_line_length: Optional[float] = None,
        is_strike: Optional[bool] = False,
        ):
        """
        Initialiser of for ThicknessCalculator
        """
        self.thickness_calculator_label = "ThicknessCalculatorBaseClass"
        self.max_line_length = max_line_length
        self.dtm_data = dtm_data
        self.bounding_box = bounding_box
        self.is_strike = is_strike

    def type(self):
        """
        Getter for subclass type label

        Returns:
            str: Name of subclass
        """
        return self.thickness_calculator_label

    @beartype.beartype
    @abstractmethod
    def compute(
        self,
        units: pandas.DataFrame,
        stratigraphic_order: list,
        basal_contacts: geopandas.GeoDataFrame,
        structure_data: pandas.DataFrame,
        geology_data: geopandas.GeoDataFrame,
        sampled_contacts: pandas.DataFrame,
    ) -> pandas.DataFrame:
        """
        Execute thickness calculator method (abstract method)

        Args:
            units (pandas.DataFrame): the data frame of units to add thicknesses to
            stratigraphic_order (list): a list of unit names sorted from youngest to oldest
            basal_contacts (geopandas.GeoDataFrame): basal contact geo data with locations and unit names of the contacts (columns must contain ["ID","basal_unit","type","geometry"])
            structure_data (pandas.DataFrame): sampled structural data
            map_data (map2loop.MapData): a catchall so that access to all map data is available

        Returns:
            pandas.DataFrame: units dataframe with added thickness column for calculated thickness values

        Note:
        -----
        Geological units' thicknesses are returned in meters. The thickness is calculated based on the stratigraphic order of the units.
        If the thickness is not calculated for a given unit, the assigned thickness is -1.
        For the bottom and top units of the stratigraphic sequence, the assigned thickness is -1.
        """

        if len(stratigraphic_order) < 3:
            logger.warning(
                f"Cannot make any thickness calculations with only {len(stratigraphic_order)} units"
            )
            return units

    def _check_thickness_percentage_calculations(self, thicknesses: pandas.DataFrame):
        units_with_no_thickness = len(thicknesses[thicknesses['ThicknessMean'] == -1])
        total_units = len(thicknesses)

        if total_units > 0 and (units_with_no_thickness / total_units) >= 0.75:
            logger.warning(
                f"More than {int(0.75 * 100)}% of units ({units_with_no_thickness}/{total_units}) "
                f"have a calculated thickness of -1. This may indicate that {self.thickness_calculator_label} "
                f"is not suitable for this dataset."
            )
    def __call__(self, **kwargs):
        return self.compute(**kwargs)

class ThicknessCalculatorAlpha(ThicknessCalculator):
    """
    ThicknessCalculator class which estimates unit thickness based on units, basal_contacts and stratigraphic order
    """

    def __init__(self):
        """
        Initialiser for alpha version of the thickness calculator
        """
        self.thickness_calculator_label = "ThicknessCalculatorAlpha"

    @beartype.beartype
    def compute(
        self,
        units: pandas.DataFrame,
        stratigraphic_order: list,
        basal_contacts: geopandas.GeoDataFrame,
        structure_data: pandas.DataFrame,
        geology_data: geopandas.GeoDataFrame,
        sampled_contacts: pandas.DataFrame,
    ) -> pandas.DataFrame:
        """
        Execute thickness calculator method takes unit data, basal_contacts and stratigraphic order and attempts to estimate unit thickness.
        Note: Thicknesses of the top and bottom units are not possible with this data and so are assigned the average of all other calculated unit thicknesses.

        Args:
            units (pandas.DataFrame): the data frame of units to add thicknesses to
            stratigraphic_order (list): a list of unit names sorted from youngest to oldest
            basal_contacts (geopandas.GeoDataFrame): basal contact geo data with locations and unit names of the contacts (columns must contain ["ID","basal_unit","type","geometry"])
            structure_data (pandas.DataFrame): sampled structural data
            map_data (map2loop.MapData): a catchall so that access to all map data is available

        Returns:
            pandas.DataFrame: units dataframe with added thickness column for calculated thickness values
        """
        # TODO: If we have orientation data near basal contact points we can estimate
        # the actual distance between contacts rather than just using the horizontal distance
        no_distance = -1.0
        basal_contacts = basal_contacts[basal_contacts["type"] == "BASAL"]
        thicknesses = units.copy()

        # Set default values
        thicknesses["ThicknessMean"] = no_distance
        thicknesses["ThicknessMedian"] = no_distance
        thicknesses["ThicknessStdDev"] = no_distance

        basal_unit_list = basal_contacts["basal_unit"].to_list()
        sampled_basal_contacts = rebuild_sampled_basal_contacts(
            basal_contacts=basal_contacts, sampled_contacts=sampled_contacts
        )

        if len(stratigraphic_order) < 3:
            logger.warning(
                f"ThicknessCalculatorAlpha: Cannot make any thickness calculations with only {len(stratigraphic_order)} units"
            )
            return thicknesses

        for i in range(1, len(stratigraphic_order) - 1):
            # Compare basal contacts of adjacent units
            if (
                stratigraphic_order[i] in basal_unit_list
                and stratigraphic_order[i + 1] in basal_unit_list
            ):
                contact1 = sampled_basal_contacts[
                    sampled_basal_contacts["basal_unit"] == stratigraphic_order[i]
                ]["geometry"].to_list()[0]
                contact2 = sampled_basal_contacts[
                    sampled_basal_contacts["basal_unit"] == stratigraphic_order[i + 1]
                ]["geometry"].to_list()[0]
                if contact1 is not None and contact2 is not None:
                    distance = contact1.distance(contact2)
                else:
                    logger.warning(
                        f"ThicknessCalculatorAlpha: Cannot calculate thickness between {stratigraphic_order[i]} and {stratigraphic_order[i + 1]} \n"
                    )
                    distance = no_distance
            else:
                logger.warning(
                    f"ThicknessCalculatorAlpha: Cannot calculate thickness between {stratigraphic_order[i]} and {stratigraphic_order[i + 1]} \n"
                )

                distance = no_distance

            # Maximum thickness is the horizontal distance between the minimum of these distances
            # Find row in unit_dataframe corresponding to unit and replace thickness value if it is -1 or larger than distance
            idx = thicknesses.index[thicknesses["name"] == stratigraphic_order[i]].tolist()[0]
            if thicknesses.loc[idx, "ThicknessMean"] == -1:
                val = distance
            else:
                val = min(distance, thicknesses.at[idx, "ThicknessMean"])
            thicknesses.loc[idx, "ThicknessMean"] = val

        self._check_thickness_percentage_calculations(thicknesses)
        
        return thicknesses


class InterpolatedStructure(ThicknessCalculator):
    """
    This class is a subclass of the ThicknessCalculator abstract base class. It implements the thickness calculation
    method for a given set of interpolated data points.

    Attributes:
        thickness_calculator_label (str): A string that stores the label of the thickness calculator.
        For this class, it is "InterpolatedStructure".

    Methods:
        compute(units: pandas.DataFrame, stratigraphic_order: list, basal_contacts: pandas.DataFrame, map_data: MapData)
        -> pandas.DataFrame: Calculates a thickness map for the overall map area.
    """

    def __init__(
        self, 
        dtm_data: Optional[gdal.Dataset] = None, 
        bounding_box: Optional[dict] = None, 
        max_line_length: Optional[float] = None,
        is_strike: Optional[bool] = False
        ):
        """
        Initialiser for interpolated structure version of the thickness calculator
        """
        super().__init__(dtm_data, bounding_box, max_line_length, is_strike)
        self.thickness_calculator_label = "InterpolatedStructure"
        self.lines = None

    @beartype.beartype
    def compute(
        self,
        units: pandas.DataFrame,
        stratigraphic_order: list,
        basal_contacts: geopandas.GeoDataFrame,
        structure_data: pandas.DataFrame,
        geology_data: geopandas.GeoDataFrame,
        sampled_contacts: pandas.DataFrame,
    ) -> pandas.DataFrame:
        """
        Execute thickness calculator method takes unit data, basal_contacts, stratigraphic order, orientation data and
        DTM data to estimate unit thickness.

        The method works by iterating over the stratigraphic order of units. For each unit, it finds the basal contact
        points and the top contact line. It then calculates the shortest line between these points. For each basal
        contact point, it finds the interpolated points that are within 10% of the length of the shortest line. It then
        calculates the true thickness of the unit at these points using the formula t = L . sin dip, where L is the
        length of the line orthogonal to both contacts and dip is the dip of the interpolated points.
        The method then calculates the median thickness and standard deviation for the unit.

        Args:
            units (pandas.DataFrame): the data frame of units to add thicknesses to
            stratigraphic_order (list): a list of unit names sorted from youngest to oldest
            basal_contacts (geopandas.GeoDataFrame): basal contact geo data with locations and unit names of
            the contacts (columns must contain ["ID","basal_unit","type","geometry"])
            structure_data (pandas.DataFrame): sampled structural data
            map_data (map2loop.MapData): a catchall so that access to all map data is available

        Returns:
            pandas.DataFrame: units dataframe with added thickness columns:
            "ThicknessMedian" is the median thickness of the unit,
            "ThicknessStdDev" is the standard deviation of the thickness of the unit
        """

        basal_contacts = basal_contacts[basal_contacts["type"] == "BASAL"].copy()

        thicknesses = units.copy()
        # Set default value
        # thicknesses["ThicknessMedian"] is the median thickness of the unit
        thicknesses["ThicknessMedian"] = -1.0
        thicknesses['ThicknessMean'] = -1.0
        # thicknesses["ThicknessStdDev"] is the standard deviation of the thickness of the unit
        thicknesses["ThicknessStdDev"] = -1.0
        thicknesses['ThicknessStdDev'] = thicknesses['ThicknessStdDev'].astype('float64')
        basal_unit_list = basal_contacts["basal_unit"].to_list()
        # increase buffer around basal contacts to ensure that the points are included as intersections
        basal_contacts["geometry"] = basal_contacts["geometry"].buffer(0.01)
        # get the sampled contacts
        contacts = geopandas.GeoDataFrame(sampled_contacts)
        # build points from x and y coordinates
        geometry2 = geopandas.points_from_xy(contacts['X'], contacts['Y'])
        contacts.set_geometry(geometry2, inplace=True)

        # set the crs of the contacts to the crs of the units
        contacts = contacts.set_crs(crs=basal_contacts.crs)
        if self.dtm_data is not None:
            # get the elevation Z of the contacts
            contacts = set_z_values_from_raster_df(self.dtm_data, contacts)
            # update the geometry of the contact points to include the Z value
            contacts["geometry"] = contacts.apply(
                lambda row: shapely.geometry.Point(row.geometry.x, row.geometry.y, row["Z"]), axis=1
            )
        else:
            contacts["geometry"] = contacts.apply(
                lambda row: shapely.geometry.Point(row.geometry.x, row.geometry.y,), axis=1
            )
        # spatial join the contact points with the basal contacts to get the unit for each contact point
        contacts = contacts.sjoin(basal_contacts, how="inner", predicate="intersects")
        # keep only necessary columns
        if 'Z' not in contacts.columns:
            contacts = contacts[["X", "Y", "geometry", "basal_unit"]].copy()
        if 'Z' in contacts.columns:
            contacts = contacts[["X", "Y", "Z", "geometry", "basal_unit"]].copy()
        # Interpolate the dip of the contacts
        interpolator = DipDipDirectionInterpolator(data_type="dip")
        # Interpolate the dip of the contacts
        dip = interpolator(self.bounding_box, structure_data, interpolator=scipy.interpolate.Rbf)
        # create a GeoDataFrame of the interpolated orientations
        interpolated_orientations = geopandas.GeoDataFrame()
        # add the dip and dip direction to the GeoDataFrame
        interpolated_orientations["dip"] = dip
        # interpolated_orientations["dip_direction"] = dip_direction
        # get the x and y coordinates of the interpolated points
        interpolated_orientations["X"] = interpolator.xi
        interpolated_orientations["Y"] = interpolator.yi
        xy = numpy.array([interpolator.xi, interpolator.yi]).T
        # create Point objects from the x and y coordinates
        interpolated_orientations.set_geometry(create_points(xy), inplace=True)
        # set the crs of the interpolated orientations to the crs of the units
        interpolated_orientations = interpolated_orientations.set_crs(crs=basal_contacts.crs)
        if self.dtm_data is not None:
            # get the elevation Z of the interpolated points
            interpolated_orientations = set_z_values_from_raster_df(self.dtm_data, interpolated_orientations)
            # update the geometry of the interpolated points to include the Z value
            interpolated_orientations["geometry"] = interpolated_orientations.apply(
                lambda row: shapely.geometry.Point(row.geometry.x, row.geometry.y, row["Z"]), axis=1
            )
        # for each interpolated point, assign name of unit using spatial join
        units = geology_data.copy()
        interpolated_orientations = interpolated_orientations.sjoin(
            units, how="inner", predicate="within"
        )
        interpolated_orientations = interpolated_orientations[
            ["geometry", "dip", "UNITNAME"]
        ].copy()

        _lines = []
        _dips = []
        _location_tracking = []

        for i in range(0, len(stratigraphic_order) - 1):
            if (
                stratigraphic_order[i] in basal_unit_list
                and stratigraphic_order[i + 1] in basal_unit_list
            ):
                basal_contact = contacts.loc[
                    contacts["basal_unit"] == stratigraphic_order[i]
                ].copy()
                top_contact = basal_contacts.loc[
                    basal_contacts["basal_unit"] == stratigraphic_order[i + 1]
                ].copy()
                top_contact_geometry = [
                    shapely.geometry.shape(geom.__geo_interface__) for geom in top_contact.geometry
                ]
                if basal_contact is not None and top_contact is not None:
                    interp_points = interpolated_orientations.loc[
                        interpolated_orientations["UNITNAME"] == stratigraphic_order[i], "geometry"
                    ].copy()
                    dip = interpolated_orientations.loc[
                        interpolated_orientations["UNITNAME"] == stratigraphic_order[i], "dip"
                    ].to_numpy()

                    _thickness = []

                    for _, row in basal_contact.iterrows():
                        # find the shortest line between the basal contact points and top contact points
                        short_line = shapely.shortest_line(row.geometry, top_contact_geometry)
                        # check if the short line is
                        if self.max_line_length is not None and short_line[0].length > self.max_line_length:
                            continue
                        if self.dtm_data is not None:
                            inv_geotransform = gdal.InvGeoTransform(self.dtm_data.GetGeoTransform())
                            data_array = numpy.array(self.dtm_data.GetRasterBand(1).ReadAsArray().T)

                        # extract the end points of the shortest line
                        p1 = numpy.zeros(3)
                        p1[0] = numpy.asarray(short_line[0].coords[0][0])
                        p1[1] = numpy.asarray(short_line[0].coords[0][1])
                        if self.dtm_data is not None:
                            # get the elevation Z of the end point p1
                            p1[2] = value_from_raster(inv_geotransform, data_array, p1[0], p1[1])
                        # create array to store xyz coordinates of the end point p2
                        p2 = numpy.zeros(3)
                        p2[0] = numpy.asarray(short_line[0].coords[-1][0])
                        p2[1] = numpy.asarray(short_line[0].coords[-1][1])
                        if self.dtm_data is not None:
                            # get the elevation Z of the end point p2
                            p2[2] = value_from_raster(inv_geotransform, data_array, p2[0], p2[1])
                        # calculate the length of the shortest line
                        line_length = scipy.spatial.distance.euclidean(p1, p2)
                        # find the indices of the points that are within 5% of the length of the shortest line
                        try:
                            # GEOS 3.10.0+
                            indices = shapely.dwithin(short_line[0], interp_points, line_length * 0.25)
                        except UnsupportedGEOSVersionError:
                            indices= numpy.array([shapely.distance(short_line[0],point)<= (line_length * 0.25) for point in interp_points])
                        # get the dip of the points that are within
                        _dip = numpy.deg2rad(dip[indices])
                        if len(_dip) > 0:
                            _lines.extend([short_line[0]]*len(_dip))
                            _dips.extend(_dip)
                        # calculate the true thickness t = L * sin(dip)
                        thickness = line_length * numpy.sin(_dip)

                        # add location tracking
                        location_tracking = pandas.DataFrame(
                            {
                                "p1_x": [p1[0]], "p1_y": [p1[1]], "p1_z": [p1[2]],
                                "p2_x": [p2[0]], "p2_y": [p2[1]], "p2_z": [p2[2]],
                                "thickness": [thickness],
                                "unit": [stratigraphic_order[i]]
                            }
                        )
                        _location_tracking.append(location_tracking)

                        # Average thickness along the shortest line
                        if all(numpy.isnan(thickness)):
                            pass
                        else:
                            _thickness.append(numpy.nanmean(thickness))

                    # calculate the median thickness and standard deviation for the unit
                    _thickness = numpy.asarray(_thickness, dtype=numpy.float64)

                    median = numpy.nanmedian(_thickness)
                    mean = numpy.nanmean(_thickness)
                    std_dev = numpy.nanstd(_thickness, dtype=numpy.float64)

                    idx = thicknesses.index[
                        thicknesses["name"] == stratigraphic_order[i + 1]
                    ].tolist()[0]
                    thicknesses.loc[idx, "ThicknessMean"] = mean
                    thicknesses.loc[idx, "ThicknessMedian"] = median
                    thicknesses.loc[idx, "ThicknessStdDev"] = std_dev

            else:
                logger.warning(
                    f"Thickness Calculator InterpolatedStructure: Cannot calculate thickness between {stratigraphic_order[i]} and {stratigraphic_order[i + 1]}\n"
                )

        # Combine all location_tracking DataFrames into a single DataFrame
        if _location_tracking and len(_location_tracking) > 0:
            combined_location_tracking = pandas.concat(_location_tracking, ignore_index=True)
            combined_location_tracking['geometry'] = combined_location_tracking.apply(
                lambda row: LineString(
                    [
                        (row['p1_x'], row['p1_y'], row['p1_z']),
                        (row['p2_x'], row['p2_y'], row['p2_z']),
                    ]
                ),
                axis=1,
            )
        else:
            combined_location_tracking = pandas.DataFrame()
        # Save the combined DataFrame as an attribute of the class

        # Convert to GeoDataFrame and set CRS to match basal_contacts
        if not combined_location_tracking.empty:
            self.location_tracking = geopandas.GeoDataFrame(combined_location_tracking, geometry='geometry', crs=basal_contacts.crs)
        else:
            self.location_tracking = geopandas.GeoDataFrame(columns=['p1_x', 'p1_y', 'p1_z', 'p2_x', 'p2_y', 'p2_z', 'thickness', 'unit', 'geometry'], crs=basal_contacts.crs)
        # Create GeoDataFrame for lines
   
        self.lines = geopandas.GeoDataFrame(geometry=_lines, crs=basal_contacts.crs)
        self.lines['dip'] = _dips

        # Check thickness calculation
        self._check_thickness_percentage_calculations(thicknesses)

        return thicknesses

class StructuralPoint(ThicknessCalculator):
    '''
    This class is a subclass of the ThicknessCalculator abstract base class. It implements the thickness calculation using a deterministic workflow based on stratigraphic measurements.

    Attributes:
        thickness_calculator_label (str): A string that stores the label of the thickness calculator.
        For this class, it is "StrucuturalPoint".

    Methods:
        compute(units: pandas.DataFrame, stratigraphic_order: list, basal_contacts: pandas.DataFrame, map_data: MapData)
        -> pandas.DataFrame: Calculates the thickness in meters for each unit in the stratigraphic column.

    '''

    def __init__(
        self, 
        dtm_data: Optional[gdal.Dataset] = None, 
        bounding_box: Optional[dict] = None, 
        max_line_length: Optional[float] = None,
        is_strike: Optional[bool] = False
        ):
        super().__init__(dtm_data, bounding_box, max_line_length, is_strike)
        self.thickness_calculator_label = "StructuralPoint"
        self.strike_allowance = 30
        self.lines = None

    @beartype.beartype
    def compute(
        self,
        units: pandas.DataFrame,
        stratigraphic_order: list,
        basal_contacts: geopandas.GeoDataFrame,
        structure_data: pandas.DataFrame,
        geology_data: geopandas.GeoDataFrame,
        sampled_contacts: pandas.DataFrame,
    ) -> pandas.DataFrame:
        """
        Method overview:
        - define perpendicular line, with strike perpendicular to the stratigraphic measurement's strike
        - find intersection points between the perpendicular line and the geological contacts
        - Perform the following checks:
            1) is there more than one intersection?
            2) are the intersections between two different lithologies, and if so, grab the neighboring lithologies only.
            3) is the distance between the two intersections less than half the map dimensions? (avoids incorrect intersections to be picked)
            4) is the stratigraphic measurement strike within 30 degrees of the strike of the geological contacts?

        - once the intersections pass the checks, calculate the thickness of the unit at the intersection points, using the general formula L*sin(dip)

        Attributes:
        -----------
            units (pandas.DataFrame): the data frame of units to add thicknesses to
            stratigraphic_order (list): a list of unit names sorted from youngest to oldest
            basal_contacts (geopandas.GeoDataFrame): basal contact geo data with locations and unit names of
            the contacts (columns must contain ["ID","basal_unit","type","geometry"])
            structure_data (pandas.DataFrame): sampled structural data
            map_data (map2loop.MapData): a catchall so that access to all map data is available

        Returns:
        --------
            pandas.DataFrame: units dataframe with added thickness columns:
            "ThicknessMedian" is the median thickness of the unit,
            "ThicknessStdDev" is the standard deviation of the thickness of the unit

        Note:
        -----
            This method is highly dependent on the existence of stratigraphic measurements that follow the strike of the geological contacts.
            If an unit does not contain a stratigraphic measurement, the thickness will not be calculated. Interpolated Structure may be used for such situations and future versions of map2loop will attempt to solve for this.
            If the thickness is not calculated for a given unit, the assigned thickness will be -1.
            For the bottom and top units of the stratigraphic sequence, the assigned thickness will also be -1.
        """
        # input sampled data
        sampled_structures = structure_data
        basal_contacts = basal_contacts.copy()

        # grab geology polygons and calculate bounding boxes for each lithology
        geology = geology_data.copy()
        geology[['minx', 'miny', 'maxx', 'maxy']] = geology.bounds

        # create a GeoDataFrame of the sampled structures
        sampled_structures = geopandas.GeoDataFrame(
            sampled_structures,
            geometry=geopandas.points_from_xy(sampled_structures.X, sampled_structures.Y),
            crs=basal_contacts.crs,
        )
        if 'UNITNAME' in sampled_structures.columns:
            sampled_structures = sampled_structures.drop(columns=['UNITNAME'])
        # add unitname to the sampled structures
        sampled_structures['unit_name'] = geopandas.sjoin(
            sampled_structures, geology, how='inner', predicate='within'
        )['UNITNAME']

        # remove nans from sampled structures
        # this happens when there are strati measurements within intrusions. If intrusions are removed from the geology map, unit_name will then return a nan
        logger.info(
            f"skipping row(s) {sampled_structures[sampled_structures['unit_name'].isnull()].index.to_numpy()} in sampled structures dataset, as they do not spatially coincide with a valid geology polygon \n"
        )
        sampled_structures = sampled_structures.dropna(subset=['unit_name'])

        # rebuild basal contacts lines based on sampled dataset
        sampled_basal_contacts = rebuild_sampled_basal_contacts(
            basal_contacts, sampled_contacts
        )

        # calculate map dimensions
        map_dx = geology.total_bounds[2] - geology.total_bounds[0]
        map_dy = geology.total_bounds[3] - geology.total_bounds[1]

        # create empty lists to store thicknesses and lithologies
        thicknesses = []
        lis = []
        _lines = []
        _dip = []

        # loop over each sampled structural measurement
        for s in range(0, len(sampled_structures)):

            # make a shapely point from the measurement
            measurement = sampled_structures.iloc[s]
            measurement_pt = shapely.geometry.Point(measurement.X, measurement.Y)

            # find unit and strike
            litho_in = measurement['unit_name']
            if self.is_strike:
                strike = measurement['DIPDIR']
            else:
                strike = (measurement['DIPDIR'] - 90) % 360

            # find bounding box of the lithology
            bbox_poly = geology[geology['UNITNAME'] == litho_in][['minx', 'miny', 'maxx', 'maxy']]

            # check if litho_in is in geology
            # for a special case when the litho_in is not in the geology
            if len(geology[geology['UNITNAME'] == litho_in]) == 0:
                logger.info(
                    f"There are structural measurements in unit - {litho_in} - that are not in the geology shapefile. Skipping this structural measurement"
                )
                continue
            else:
                # make a subset of the geology polygon & find neighbour units
                GEO_SUB = geology[geology['UNITNAME'] == litho_in]['geometry'].values[0]

            neighbor_list = list(
                basal_contacts[GEO_SUB.intersects(basal_contacts.geometry)]['basal_unit']
            )

            # draw orthogonal line to the strike (default value 10Km), and clip it by the bounding box of the lithology
            if self.max_line_length is None:
                self.max_line_length = 10000
            B = calculate_endpoints(measurement_pt, float(strike), float(self.max_line_length), bbox_poly)
            b = geopandas.GeoDataFrame({'geometry': [B]}).set_crs(basal_contacts.crs)

            # find all intersections
            all_intersections = sampled_basal_contacts.overlay(
                b, how='intersection', keep_geom_type=False
            )
            all_intersections = all_intersections[
                all_intersections['geometry'].geom_type == 'Point'
            ]

            # clip intersections by the neighbouring geology polygons
            final_intersections = all_intersections[
                all_intersections['basal_unit'].isin(neighbor_list)
            ]

            # sometimes the intersections will return as MultiPoint, so we need to convert them to nearest point
            if 'MultiPoint' in final_intersections['geometry'].geom_type.values:
                multi = final_intersections[
                    final_intersections['geometry'].geom_type == 'MultiPoint'
                ].index
                for m in multi:
                    nearest_ = shapely.ops.nearest_points(
                        final_intersections.loc[m, :].geometry, measurement_pt
                    )[0]
                    final_intersections.at[m, 'geometry'] = nearest_
                    final_intersections.at[m, 'geometry'] = nearest_

            # check to see if there's less than 2 intersections
            if len(final_intersections) < 2:
                continue

            # check to see if the intersections cross two lithologies"
            if len(final_intersections['basal_unit'].unique()) == 1:
                continue

            # declare the two intersection points
            int_pt1 = final_intersections.iloc[0].geometry
            int_pt2 = final_intersections.iloc[1].geometry

            # if the intersections are too far apart, skip
            if (
                math.sqrt(((int_pt1.x - int_pt2.x) ** 2) + ((int_pt1.y - int_pt2.y) ** 2))
                > map_dx / 2
                or math.sqrt(((int_pt1.x - int_pt2.x) ** 2) + ((int_pt1.y - int_pt2.y) ** 2))
                > map_dy / 2
            ):
                continue

            # find the segments that the intersections belong to
            seg1 = sampled_basal_contacts[
                sampled_basal_contacts['basal_unit'] == final_intersections.iloc[0]['basal_unit']
            ].geometry.iloc[0]
            seg2 = sampled_basal_contacts[
                sampled_basal_contacts['basal_unit'] == final_intersections.iloc[1]['basal_unit']
            ].geometry.iloc[0]

            # simplify the geometries to LineString
            if seg1.geom_type == 'MultiLineString':
                seg1 = multiline_to_line(seg1)
            if seg2.geom_type == 'MultiLineString':
                seg2 = multiline_to_line(seg2)

            # find the strike of the segments
            strike1 = find_segment_strike_from_pt(seg1, int_pt1, measurement)
            strike2 = find_segment_strike_from_pt(seg2, int_pt2, measurement)

            # check to see if the strike of the stratigraphic measurement is within the strike allowance of the strike of the geological contact
            b_s = strike - self.strike_allowance, strike + self.strike_allowance
            if not (b_s[0] < strike1 < b_s[1] and b_s[0] < strike2 < b_s[1]):
                continue

            # build the debug info
            line = shapely.geometry.LineString([int_pt1, int_pt2])
            _lines.append(line)
            _dip.append(measurement['DIP'])  

            # find the lenght of the segment
            L = math.sqrt(((int_pt1.x - int_pt2.x) ** 2) + ((int_pt1.y - int_pt2.y) ** 2))

            # if length is higher than max_line_length, skip
            if self.max_line_length is not None and L > self.max_line_length:
                continue

            # calculate thickness
            thickness = L * math.sin(math.radians(measurement['DIP']))

            thicknesses.append(thickness)
            lis.append(litho_in)

        # create the debug gdf
        self.lines = geopandas.GeoDataFrame(geometry=_lines, crs=basal_contacts.crs)
        self.lines["DIP"] = _dip

        # create a DataFrame of the thicknesses median and standard deviation by lithology
        result = pandas.DataFrame({'unit': lis, 'thickness': thicknesses})
        result = result.groupby('unit')['thickness'].agg(['median', 'mean', 'std']).reset_index()
        result.rename(columns={'thickness': 'ThicknessMedian'}, inplace=True)

        output_units = units.copy()
        # remove the old thickness column
        output_units['ThicknessMedian'] = numpy.full(len(output_units), numpy.nan)
        output_units['ThicknessMean'] = numpy.full(len(output_units), numpy.nan)
        output_units['ThicknessStdDev'] = numpy.full(len(output_units), numpy.nan)

        # find which units have no thickness calculated
        names_not_in_result = units[~units['name'].isin(result['unit'])]['name'].to_list()
        # assign the thicknesses to the each unit
        for _, unit in result.iterrows():
            idx = units.index[units['name'] == unit['unit']].tolist()[0]
            output_units.loc[idx, 'ThicknessMedian'] = unit['median']
            output_units.loc[idx, 'ThicknessMean'] = unit['mean']
            output_units.loc[idx, 'ThicknessStdDev'] = unit['std']

        output_units["ThicknessMean"] = output_units["ThicknessMean"].fillna(-1)
        output_units["ThicknessMedian"] = output_units["ThicknessMedian"].fillna(-1)
        output_units["ThicknessStdDev"] = output_units["ThicknessStdDev"].fillna(-1)

        # handle the units that have no thickness
        for unit in names_not_in_result:
            # if no thickness has been calculated for the unit
            if (
                # not a top//bottom unit
                (output_units[output_units['name'] == unit]['ThicknessMedian'].all() == 0)
                and (unit != stratigraphic_order[-1])
                and (unit != stratigraphic_order[0])
            ):
                idx = stratigraphic_order.index(unit)
                # throw warning to user
                logger.warning(
                    'Thickness Calculator StructuralPoint: Cannot calculate thickness between',
                    unit,
                    "and ",
                    stratigraphic_order[idx + 1],
                    "\n",
                )
                # assign -1 as thickness
                output_units.loc[output_units["name"] == unit, "ThicknessMedian"] = -1
                output_units.loc[output_units["name"] == unit, "ThicknessMean"] = -1
                output_units.loc[output_units["name"] == unit, "ThicknessStdDev"] = -1

            # if top//bottom unit assign -1
            if unit in [stratigraphic_order[-1], stratigraphic_order[0]]:
                output_units.loc[output_units["name"] == unit, "ThicknessMedian"] = -1
                output_units.loc[output_units["name"] == unit, "ThicknessMean"] = -1
                output_units.loc[output_units["name"] == unit, "ThicknessStdDev"] = -1

        self._check_thickness_percentage_calculations(output_units)

        return output_units


class AlongSection(ThicknessCalculator):
    """Thickness calculator that estimates true thicknesses along supplied section lines."""

    def __init__(
        self,
        sections: geopandas.GeoDataFrame,
        dtm_data: Optional[gdal.Dataset] = None,
        bounding_box: Optional[dict] = None,
        max_line_length: Optional[float] = None,
        is_strike: Optional[bool] = False,
    ):
        super().__init__(dtm_data, bounding_box, max_line_length, is_strike)  # initialise base calculator bits
        self.thickness_calculator_label = "AlongSection"  # label used externally to identify this strategy
        self.sections = (
            sections.copy()  # keep a copy so editing outside does not mutate our internal state
            if sections is not None
            else geopandas.GeoDataFrame({"geometry": []}, geometry="geometry")  # ensure predictable structure when no sections supplied
        )
        self.section_thickness_records: geopandas.GeoDataFrame = geopandas.GeoDataFrame()  # populated as GeoDataFrame during compute
        self.section_intersection_points: dict = {}

    def type(self):
        """Return the calculator label."""
        return self.thickness_calculator_label

    @beartype.beartype
    def compute(
        self,
        units: pandas.DataFrame,
        stratigraphic_order: list,
        basal_contacts: geopandas.GeoDataFrame,
        structure_data: pandas.DataFrame,
        geology_data: geopandas.GeoDataFrame,
        sampled_contacts: pandas.DataFrame,
    ) -> pandas.DataFrame:
        """Estimate unit thicknesses along cross sections.

        Args:
            units (pandas.DataFrame): Table of stratigraphic units that will be annotated with
                thickness statistics. Must expose a ``name`` column used to map aggregated
                thickness values back to each unit.
            stratigraphic_order (list): Accepted for API compatibility but not used by this
                calculator. Retained so the method signature matches the other calculators.
            basal_contacts (geopandas.GeoDataFrame): Unused placeholder maintained for parity
                with the other calculators.
            structure_data (pandas.DataFrame): Orientation measurements containing X/Y columns
                and one dip column (``DIP``, ``dip`` or ``Dip``). A KD-tree is built from these
                points so each split section segment can grab the nearest dip when converting its
                length to true thickness.
            geology_data (geopandas.GeoDataFrame): Mandatory polygon dataset describing map
                units. Each section is intersected with these polygons to generate unit-aware
                line segments whose lengths underpin the thickness calculation.
            sampled_contacts (pandas.DataFrame): Part of the public interface but not used in
                this implementation.

        Returns:
            pandas.DataFrame: Copy of ``units`` with three new columns – ``ThicknessMean``,
            ``ThicknessMedian`` and ``ThicknessStdDev`` – populated from the accumulated segment
            thickness samples. Units that never receive a measurement retain the sentinel value
            ``-1`` in each column.

        Workflow:
            1. Validate the presence of section geometries, geology polygons, and a recognizable
               unit-name column. Reproject sections to match the geology CRS when required.
            2. Pre-build helpers: a spatial index over geology polygons for fast section/polygon
               queries, and (if possible) a KD-tree of orientation points so the closest dip to a
               split segment can be queried in O(log n).
            3. For every section, reduce complex/multi-part geometries to a single `LineString`,
               intersect it with candidate geology polygons, and collect the resulting split line
               segments along with their length and measure positions along the parent section.
            4. Walk the ordered segments and keep those bounded by two distinct neighbouring units
               (i.e., the segment sits between different units on either side). Fetch the nearest
               dip (fallback to 90° once if no structures exist), convert the segment length to
               true thickness using ``length * sin(dip)``, and store the sample plus provenance in
               ``self.section_thickness_records``.
            5. Aggregate all collected samples per unit (mean/median/std) and return the enriched
               ``units`` table. The helper ``self.section_intersection_points`` captures the raw
               split segments grouped by section to aid downstream inspection.

        Notes:
            - ``stratigraphic_order``, ``basal_contacts`` and ``sampled_contacts`` are retained for
              interface compatibility but do not influence this algorithm.
            - When no nearby orientation measurements exist, the method emits a single warning and
              assumes a vertical dip (90°) for affected segments to avoid silently dropping units.
            - ``self.section_thickness_records`` is materialised as a ``GeoDataFrame`` so the split
              segments (geometry column) can be visualised directly in notebooks or GIS clients.
        """

        if self.sections is None or self.sections.empty:
            logger.warning("AlongSection: No sections provided; skipping thickness calculation.")
            return units

        if geology_data is None or geology_data.empty:
            logger.warning(
                "AlongSection: Geology polygons are required to split sections; skipping thickness calculation."
            )
            return units
        unit_column = "UNITNAME"
        sections = self.sections.copy()
        geology = geology_data.copy()

        if geology.crs is not None:
            if sections.crs is None:
                sections = sections.set_crs(geology.crs)
            elif sections.crs != geology.crs:
                sections = sections.to_crs(geology.crs)
        elif sections.crs is not None:
            geology = geology.set_crs(sections.crs)

        try:
            geology_sindex = geology.sindex
        except Exception as e:
            logger.error("Failed to create spatial index for geology data: %s", e, exc_info=True)
            geology_sindex = None

        sections_bounds = sections.total_bounds
        geology_bounds = geology.total_bounds
        if (
            sections_bounds[2] < geology_bounds[0]
            or sections_bounds[0] > geology_bounds[2]
            or sections_bounds[3] < geology_bounds[1]
            or sections_bounds[1] > geology_bounds[3]
        ):
            logger.warning("AlongSection: Sections do not overlap geology extent; skipping thickness calculation.")
            return units

        if geology_sindex is None:
            overlaps = False
            for geom in sections.geometry:
                if geom is None or geom.is_empty:
                    continue
                if geology.intersects(geom).any():
                    overlaps = True
                    break
            if not overlaps:
                logger.warning("AlongSection: Sections do not overlap geology; skipping thickness calculation.")
                return units
        else:
            overlaps = False
            for geom in sections.geometry:
                if geom is None or geom.is_empty:
                    continue
                candidate_idx = list(geology_sindex.intersection(geom.bounds))
                if not candidate_idx:
                    continue
                if geology.iloc[candidate_idx].intersects(geom).any():
                    overlaps = True
                    break
            if not overlaps:
                logger.warning("AlongSection: Sections do not overlap geology; skipping thickness calculation.")
                return units

        dip_column = next((col for col in ("DIP", "dip", "Dip") if col in structure_data.columns), None)
        orientation_tree = None
        orientation_dips = None
        orientation_coords = None
        if dip_column is not None and {"X", "Y"}.issubset(structure_data.columns):
            orient_df = structure_data.dropna(subset=["X", "Y", dip_column]).copy()
            if not orient_df.empty:
                orient_df["_dip_value"] = pandas.to_numeric(orient_df[dip_column], errors="coerce")
                orient_df = orient_df.dropna(subset=["_dip_value"])
                if not orient_df.empty:
                    try:
                        orientation_coords = orient_df[["X", "Y"]].astype(float).to_numpy()
                    except ValueError:
                        logger.debug(
                            "Failed to convert orientation coordinates to float for %d rows. Data quality issue likely. Example rows: %s",
                            len(orient_df),
                            orient_df[["X", "Y"]].head(3).to_dict(orient="records")
                        )
                        orientation_coords = numpy.empty((0, 2))
                    if orientation_coords.size:
                        orientation_dips = orient_df["_dip_value"].astype(float).to_numpy()
                        try:
                            orientation_tree = cKDTree(orientation_coords)
                        except (ValueError, TypeError) as e:
                            logger.error(f"Failed to construct cKDTree for orientation data: {e}")
                            orientation_tree = None
        default_dip_warning_emitted = False

        units_lookup = dict(zip(units["name"], units.index))
        thickness_by_unit = {name: [] for name in units_lookup.keys()}
        thickness_records: List[dict] = []
        split_segments_by_section: dict = defaultdict(list)

        for section_idx, section_row in sections.iterrows():
            line = clean_line_geometry(section_row.geometry)
            if line is None or line.length == 0:
                continue
            section_id = section_row.get("ID", section_idx)

            candidate_idx = (
                list(geology_sindex.intersection(line.bounds))
                if geology_sindex is not None
                else list(range(len(geology)))
            )
            if not candidate_idx:
                continue

            split_segments = []
            for idx in candidate_idx:
                poly = geology.iloc[idx]
                polygon_geom = poly.geometry
                if polygon_geom is None or polygon_geom.is_empty:
                    continue
                intersection = line.intersection(polygon_geom)
                if intersection.is_empty:
                    continue
                for segment in iter_line_segments(intersection):
                    seg_length = segment.length
                    if seg_length <= 0:
                        continue
                    start_measure, end_measure = segment_measure_range(line, segment)
                    segment_record = {
                        "section_id": section_id,
                        "geometry": segment,
                        "unit": poly[unit_column],
                        "length": seg_length,
                        "start_measure": start_measure,
                        "end_measure": end_measure,
                    }
                    split_segments.append(segment_record)
                    split_segments_by_section[section_id].append(segment_record)

            if len(split_segments) < 1:
                continue

            split_segments.sort(key=lambda item: (item["start_measure"], item["end_measure"]))

            for idx_segment, segment in enumerate(split_segments):
                prev_unit = split_segments[idx_segment - 1]["unit"] if idx_segment > 0 else None
                next_unit = (
                    split_segments[idx_segment + 1]["unit"]
                    if idx_segment + 1 < len(split_segments)
                    else None
                )

                if prev_unit is None or next_unit is None:
                    continue
                if prev_unit == segment["unit"] or next_unit == segment["unit"]:
                    continue
                if prev_unit == next_unit:
                    continue

                dip_value, orientation_idx, orientation_distance = nearest_orientation_to_line(
                    orientation_tree, 
                    orientation_dips,
                    orientation_coords,
                    segment["geometry"]
                )
                if numpy.isnan(dip_value):
                    dip_value = 90.0
                    orientation_idx = None
                    orientation_distance = None
                    if not default_dip_warning_emitted:
                        logger.warning(
                            "AlongSection: Missing structure measurements near some sections; assuming vertical dip (90°) for those segments."
                        )
                        default_dip_warning_emitted = True

                thickness = segment["length"] * math.sin(math.radians(dip_value))
                if thickness <= 0:
                    continue

                unit_name = segment["unit"]
                if unit_name not in thickness_by_unit:
                    thickness_by_unit[unit_name] = []
                thickness_by_unit[unit_name].append(thickness)

                thickness_records.append(
                    {
                        "section_id": section_id,
                        "unit_id": unit_name,
                        "thickness": thickness,
                        "segment_length": segment["length"],
                        "dip_used_deg": dip_value,
                        "prev_unit": prev_unit,
                        "next_unit": next_unit,
                        "orientation_index": orientation_idx,
                        "orientation_distance": orientation_distance,
                        "geometry": segment["geometry"],
                    }
                )

        output_units = units.copy()
        output_units["ThicknessMean"] = -1.0
        output_units["ThicknessMedian"] = -1.0
        output_units["ThicknessStdDev"] = -1.0

        for unit_name, values in thickness_by_unit.items():
            if not values:
                continue
            idx = units_lookup.get(unit_name)
            if idx is None:
                continue
            arr = numpy.asarray(values, dtype=numpy.float64)
            output_units.at[idx, "ThicknessMean"] = float(numpy.nanmean(arr))
            output_units.at[idx, "ThicknessMedian"] = float(numpy.nanmedian(arr))
            output_units.at[idx, "ThicknessStdDev"] = float(numpy.nanstd(arr))

        if thickness_records:
            self.section_thickness_records = geopandas.GeoDataFrame(
                thickness_records,
                geometry="geometry",
                crs=sections.crs,
            )
        else:
            self.section_thickness_records = geopandas.GeoDataFrame(
                columns=[
                    "section_id",
                    "unit_id",
                    "thickness",
                    "segment_length",
                    "dip_used_deg",
                    "prev_unit",
                    "next_unit",
                    "orientation_index",
                    "orientation_distance",
                    "geometry",
                ],
                geometry="geometry",
                crs=sections.crs,
            )
        self.section_intersection_points = dict(split_segments_by_section)
        self._check_thickness_percentage_calculations(output_units)

        return output_units
