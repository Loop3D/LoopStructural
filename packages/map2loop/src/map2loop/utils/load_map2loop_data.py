import geopandas
import pathlib
from osgeo import gdal
gdal.UseExceptions()

# Use the path of this file to locate the datasets directory
def map2loop_dir(folder) -> pathlib.Path:
    path = pathlib.Path(__file__).parent.parent.parent / 'map2loop' / '_datasets' / 'geodata_files' / f'{folder}'
    return path

def load_hamersley_geology():
    path = map2loop_dir('hamersley') / "geology.geojson"
    return geopandas.read_file(str(path))

def load_hamersley_structure():
    path = map2loop_dir('hamersley') / "structure.geojson"
    return geopandas.read_file(str(path))

def load_hamersley_dtm():
    path = map2loop_dir('hamersley') / "dtm_rp.tif"
    return gdal.Open(str(path))
