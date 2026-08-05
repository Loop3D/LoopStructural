import geopandas
import map2loop
import pathlib
from osgeo import gdal
gdal.UseExceptions()
def map2loop_dir(folder)-> pathlib.Path:
    path = pathlib.Path(map2loop.__file__).parent
    path = path / "_datasets"/"geodata_files"/f'{folder}'
    return path
def load_hamersley_geology():
    """
    Loads Hamersley geology data from a shapefile

    Args:
        path (str):
            The path to the shapefile

    Returns:
        geopandas.GeoDataFrame: The geology data
    """
    
    path = map2loop_dir('hamersley') / "geology.geojson"
    return geopandas.read_file(str(path))


def load_hamersley_structure():
    """
    Loads Hamersley structure data from a shapefile

    Args:
        path (str):
            The path to the shapefile

    Returns:
        geopandas.GeoDataFrame: The structure data
    """

    path = map2loop_dir('hamersley') / "structure.geojson"
    return geopandas.read_file(str(path))


def load_hamersley_dtm():
    """
    Load DTM data from a raster file

    Returns:
        gdal.Dataset: The DTM data
    """
    path = map2loop_dir('hamersley') / "dtm_rp.tif"
    return gdal.Open(str(path))
