import pathlib
from map2loop.project import Project
from map2loop.m2l_enums import Datatype
import map2loop
from unittest.mock import patch


# Sample test function for lithology and fault ignore codes
def test_set_get_ignore_codes():
    # Set up a sample bounding box and other required data
    bbox_3d = {
        "minx": 515687.31005864,
        "miny": 7493446.76593407,
        "maxx": 562666.860106543,
        "maxy": 7521273.57407786,
        "base": -3200,
        "top": 3000,
    }

    # Set up the config dictionary
    config_dictionary = {
        "structure": {"dipdir_column": "azimuth2", "dip_column": "dip"},
        "geology": {"unitname_column": "unitname", "alt_unitname_column": "code"},
        "fault": {'structtype_column': 'feature', 'fault_text': 'Fault'},
    }
    with patch.object(Project, 'validate_required_inputs', return_value=None):
        project = Project(
            working_projection='EPSG:28350',
            bounding_box=bbox_3d,
            geology_filename=str(
                pathlib.Path(map2loop.__file__).parent
                / pathlib.Path('_datasets/geodata_files/hamersley/geology.geojson')
            ),
            fault_filename=str(
                pathlib.Path(map2loop.__file__).parent
                / pathlib.Path('_datasets/geodata_files/hamersley/faults.geojson')
            ),
            dtm_filename=str(
                pathlib.Path(map2loop.__file__).parent
                / pathlib.Path('_datasets/geodata_files/hamersley/dtm_rp.tif')
            ),
            config_dictionary=config_dictionary,
            structure_filename="", 
        )

    # Define test ignore codes for lithology and faults
    lithology_codes = ["cover", "Fortescue_Group", "A_FO_od"]
    fault_codes = ['Fault_9', "NotAFault"]

    # Test Lithology ignore codes
    project.set_ignore_lithology_codes(lithology_codes)
    assert (
        project.map_data.get_ignore_lithology_codes() == lithology_codes
    ), "Lithology ignore codes mismatch"

    # Test Fault ignore codes
    project.set_ignore_fault_codes(fault_codes)
    assert project.map_data.get_ignore_fault_codes() == fault_codes, "Fault ignore codes mismatch"

    project.map_data.parse_fault_map()
    # check if the skipped fault has been removed now:
    assert not project.map_data.get_map_data(Datatype.FAULT)['NAME'].str.contains('Fault_9').any()

    project.map_data.parse_geology_map()
    # check if the skipped lithology has been removed now:
    assert (
        not project.map_data.get_map_data(Datatype.GEOLOGY)['UNITNAME']
        .str.contains('Fortescue_Group')
        .any()
    )
    assert (
        not project.map_data.get_map_data(Datatype.GEOLOGY)['UNITNAME'].str.contains('cover').any()
    )
