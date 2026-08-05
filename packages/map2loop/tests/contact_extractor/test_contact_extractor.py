import sys
sys.path.append('/usr/lib/python3/dist-packages')
from map2loop.contact_extractor import ContactExtractor
import geopandas as gpd
from shapely.geometry import Polygon

def simple_geology():
    poly1 = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    poly2 = Polygon([(1, 0), (2, 0), (2, 1), (1, 1)])
    return gpd.GeoDataFrame(
        {
            "UNITNAME": ["A", "B"],
            "INTRUSIVE": [False, False],
            "SILL": [False, False],
            "geometry": [poly1, poly2],
        },
        geometry="geometry",
        crs="EPSG:28350",
    )

def test_extract_all_contacts():
    geology = simple_geology()
    extractor = ContactExtractor(geology, None)
    contacts = extractor.extract_all_contacts()
    assert len(contacts) == 1
    assert {contacts.loc[0, "UNITNAME_1"], contacts.loc[0, "UNITNAME_2"]} == {"A", "B"}

def test_extract_basal_contacts():
    geology = simple_geology()
    extractor = ContactExtractor(geology, None)
    extractor.extract_all_contacts()
    basal = extractor.extract_basal_contacts(["A", "B"], save_contacts=True)
    assert len(basal) == 1
    assert basal.loc[0, "basal_unit"] == "A"
    assert basal.loc[0, "type"] == "BASAL"
    assert extractor.basal_contacts is not None
    assert extractor.all_basal_contacts is not None
