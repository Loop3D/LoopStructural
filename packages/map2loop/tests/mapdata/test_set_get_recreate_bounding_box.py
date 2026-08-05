import pytest
from map2loop.mapdata import MapData
import geopandas
import shapely


@pytest.fixture
def md():
    return MapData()


def test_set_bounding_box_with_tuple(md):
    bounding_box = (0, 10, 0, 10)
    md.set_bounding_box(bounding_box)

    assert md.bounding_box == {
        "minx": 0,
        "maxx": 10,
        "miny": 0,
        "maxy": 10,
        "top": 0,
        "base": 2000,
    }, "MapData.set_bounding_box() not working as expected"


def test_set_bounding_box_with_dict(md):
    bounding_box = {"minx": 0, "maxx": 10, "miny": 0, "maxy": 10, "top": 5, "base": 15}
    md.set_bounding_box(bounding_box)

    assert md.bounding_box == bounding_box, "MapData.set_bounding_box() not working as expected"


def test_bounding_box_keys(md):
    bounding_box = (0, 10, 0, 10)
    md.set_bounding_box(bounding_box)

    for key in ["minx", "maxx", "miny", "maxy", "top", "base"]:
        assert key in md.bounding_box, f"MapData.bounding_box missing key: {key}"


def test_bounding_box_polygon(md):
    bounding_box = (0, 10, 0, 10)
    md.set_bounding_box(bounding_box)

    minx, miny, maxx, maxy = 0, 0, 10, 10
    lat_point_list = [miny, miny, maxy, maxy, miny]
    lon_point_list = [minx, maxx, maxx, minx, minx]
    expected_polygon = geopandas.GeoDataFrame(
        index=[0],
        crs=md.working_projection,
        geometry=[shapely.geometry.Polygon(zip(lon_point_list, lat_point_list))],
    )

    assert md.bounding_box_polygon.equals(
        expected_polygon
    ), "MapData.bounding_box_polygon not returning the correct GeoDataFrame"


def test_get_bounding_box_as_dict(md):
    bounding_box = {"minx": 0, "maxx": 10, "miny": 0, "maxy": 10, "top": 5, "base": 15}
    md.set_bounding_box(bounding_box)
    result = md.get_bounding_box()

    assert (
        result == bounding_box
    ), "MapData.get_bounding_box() not returning the correct bounding box"


def test_get_bounding_box_as_polygon(md):
    bounding_box = (0, 10, 0, 10)
    md.set_bounding_box(bounding_box)
    result = md.get_bounding_box(polygon=True)

    assert isinstance(
        result, geopandas.GeoDataFrame
    ), "MapData.get_bounding_box(polygon=True) not returning a GeoDataFrame"
    assert result.equals(
        md.bounding_box_polygon
    ), "MapData.get_bounding_box(polygon=True) not returning the correct GeoDataFrame"


def test_recreate_bounding_box_str(md):
    bounding_box = (0, 10, 0, 10)
    md.set_working_projection("EPSG:4326")
    md.set_bounding_box(bounding_box)
    md.recreate_bounding_box_str()

    expected_str = "0,0,10,10,EPSG:4326"
    assert (
        md.bounding_box_str == expected_str
    ), "MapData.recreate_bounding_box_str() not working as expected"


def test_set_bounding_box_with_missing_keys(md):
    bounding_box = {
        "minx": 0,
        "maxx": 10,
        "miny": 0
        # Missing "maxy", "top", "base"
    }
    with pytest.raises(KeyError):
        md.set_bounding_box(
            bounding_box
        ), "MapData.set_bounding_box accepting wrong argument, but should raise KeyError"
