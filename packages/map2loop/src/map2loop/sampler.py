# internal imports
from .utils import set_z_values_from_raster_df

# external imports
from abc import ABC, abstractmethod
import beartype
import geopandas
import pandas
import shapely
import numpy
from typing import Optional
from osgeo import gdal


class Sampler(ABC):
    """
    Base Class of Sampler used to force structure of Sampler

    Args:
        ABC (ABC): Derived from Abstract Base Class
    """

    def __init__(self, dtm_data=None, geology_data=None):
        """
        Initialiser of for Sampler
        """
        self.sampler_label = "SamplerBaseClass"
        self.dtm_data = dtm_data
        self.geology_data = geology_data

    def type(self):
        """
        Getter for subclass type label

        Returns:
            str: Name of subclass
        """
        return self.sampler_label

    @beartype.beartype
    @abstractmethod
    def sample(
        self, spatial_data: geopandas.GeoDataFrame
    ) -> pandas.DataFrame:
        """
        Execute sampling method (abstract method)

        Args:
            spatial_data (geopandas.GeoDataFrame): data frame to sample

        Returns:
            pandas.DataFrame: data frame containing samples
        """
        pass
    def __call__(self, **kwargs):
        return self.sample(**kwargs)

class SamplerDecimator(Sampler):
    """
    Decimator sampler class which decimates the geo data frame based on the decimation value
    ie. decimation = 10 means take every tenth point
    Note: This only works on data frames with lists of points with columns "X" and "Y"
    """

    @beartype.beartype
    def __init__(self, decimation: int = 1, dtm_data: Optional[gdal.Dataset] = None, geology_data: Optional[geopandas.GeoDataFrame] = None):
        """
        Initialiser for decimator sampler

        Args:
            decimation (int, optional): stride of the points to sample. Defaults to 1.
            dtm_data (Optional[gdal.Dataset], optional): digital terrain map data. Defaults to None.
            geology_data (Optional[geopandas.GeoDataFrame], optional): geology data. Defaults to None.
        """
        super().__init__(dtm_data, geology_data)
        self.sampler_label = "SamplerDecimator"
        decimation = max(decimation, 1)
        self.decimation = decimation


    @beartype.beartype
    def sample(
        self, spatial_data: geopandas.GeoDataFrame
    ) -> pandas.DataFrame:
        """
        Execute sample method takes full point data, samples the data and returns the decimated points

        Args:
            spatial_data (geopandas.GeoDataFrame): the data frame to sample

        Returns:
            pandas.DataFrame: the sampled data points
        """
        data = spatial_data.copy()
        data["X"] = data.geometry.x
        data["Y"] = data.geometry.y
        if self.dtm_data is not None:
            result = set_z_values_from_raster_df(self.dtm_data, data)
            if result is not None:
                data["Z"] = result["Z"]
            else:
                data["Z"] = None
        else:
            data["Z"] = None
        if self.geology_data is not None:
            data["layerID"] = geopandas.sjoin(
                data, self.geology_data, how='left'
            )['index_right']
        else:
            data["layerID"] = None
        data.reset_index(drop=True, inplace=True)

        return pandas.DataFrame(data[:: self.decimation].drop(columns="geometry"))


class SamplerSpacing(Sampler):
    """
    Spacing based sampler which decimates the geo data frame based on the distance between points along a line or
    in the case of a polygon along the boundary of that polygon
    ie. spacing = 500 means take a sample every 500 metres
    Note: This only works on data frames that contain MultiPolgon, Polygon, MultiLineString and LineString geometry
    """

    @beartype.beartype
    def __init__(self, spacing: float = 50.0, dtm_data: Optional[gdal.Dataset] = None, geology_data: Optional[geopandas.GeoDataFrame] = None):
        """
        Initialiser for spacing sampler

        Args:
            spacing (float, optional): The distance between samples. Defaults to 50.0.
            dtm_data (Optional[gdal.Dataset], optional): digital terrain map data. Defaults to None.
            geology_data (Optional[geopandas.GeoDataFrame], optional): geology data. Defaults to None.
        """
        super().__init__(dtm_data, geology_data)
        self.sampler_label = "SamplerSpacing"
        spacing = max(spacing, 1.0)
        self.spacing = spacing

    def __call__(self, **kwargs):
        return self.sample(**kwargs)
    @beartype.beartype
    def sample(
        self, spatial_data: geopandas.GeoDataFrame
    ) -> pandas.DataFrame:
        """
        Execute sample method takes full point data, samples the data and returns the sampled points

        Args:
            spatial_data (geopandas.GeoDataFrame): the data frame to sample (must contain column ["ID"])

        Returns:
            pandas.DataFrame: the sampled data points
        """
        schema = {"ID": str, "X": float, "Y": float, "featureId": str}
        df = pandas.DataFrame(columns=schema.keys()).astype(schema)
        for _, row in spatial_data.iterrows():
            if type(row.geometry) is shapely.geometry.multipolygon.MultiPolygon:
                targets = row.geometry.boundary.geoms
            elif type(row.geometry) is shapely.geometry.polygon.Polygon:
                targets = [row.geometry.boundary]
            elif type(row.geometry) is shapely.geometry.multilinestring.MultiLineString:
                targets = row.geometry.geoms
            elif type(row.geometry) is shapely.geometry.linestring.LineString:
                targets = [row.geometry]
            else:
                targets = []

            # For the main cases Polygon and LineString the list 'targets' has one element
            for a, target in enumerate(targets):
                df2 = pandas.DataFrame(columns=schema.keys()).astype(schema)
                distances = numpy.arange(0, target.length, self.spacing)[:-1]
                points = [target.interpolate(distance) for distance in distances]
                df2["X"] = [point.x for point in points]
                df2["Y"] = [point.y for point in points]

                # # account for holes//rings in polygons
                df2["featureId"] = str(a)
                # 1. check if line is "closed"
                if target.is_ring:
                    target_polygon = shapely.geometry.Polygon(target)
                    if target_polygon.exterior.is_ccw:  # if counterclockwise --> hole
                        for j, target2 in enumerate(targets):
                            # skip if line or point
                            if len(target2.coords) >= 2:
                                continue
                            # which poly is the hole in? assign featureId of the same poly
                            t2_polygon = shapely.geometry.Polygon(target2)
                            if target.within(t2_polygon):  #
                                df2['featureId'] = str(j)

                df2["ID"] = row["ID"] if "ID" in spatial_data.columns else 0
                df = df2 if len(df) == 0 else pandas.concat([df, df2])

        df.reset_index(drop=True, inplace=True)
        return df
