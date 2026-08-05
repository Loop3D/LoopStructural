import pytest
import pathlib
from map2loop.project import Project
import map2loop

# ------------------------------------------------------------------------------
# Common fixtures or helper data (bounding box, minimal filenames, etc.)
# ------------------------------------------------------------------------------

@pytest.fixture
def minimal_bounding_box():
    return {
        "minx": 515687.31005864,
        "miny": 7493446.76593407,
        "maxx": 562666.860106543,
        "maxy": 7521273.57407786,
        "base": -3200,
        "top": 3000,
    }

@pytest.fixture
def geology_file():
    return str(
        pathlib.Path(map2loop.__file__).parent
        / pathlib.Path('_datasets/geodata_files/hamersley/geology.geojson')
    )

@pytest.fixture
def structure_file():
    return str(
        pathlib.Path(map2loop.__file__).parent
        / pathlib.Path('_datasets/geodata_files/hamersley/structure.geojson')
    )

@pytest.fixture
def dtm_file():
    return str(
        pathlib.Path(map2loop.__file__).parent
        / pathlib.Path('_datasets/geodata_files/hamersley/dtm_rp.tif')
    )

@pytest.fixture
def valid_config_dictionary():
    """
    A valid config dictionary that meets the 'structure' and 'geology' requirements
    """
    return {
        "structure": {
            "dipdir_column": "azimuth2",
            "dip_column": "dip"
        },
        "geology": {
            "unitname_column": "unitname",
            "alt_unitname_column": "code",
        }
    }



# 1) config_filename and config_dictionary both present should raise ValueError
def test_config_filename_and_dictionary_raises_error(
    minimal_bounding_box, geology_file, dtm_file, structure_file, valid_config_dictionary
):

    with pytest.raises(ValueError, match="Both 'config_filename' and 'config_dictionary' were provided"):
        Project(
            bounding_box=minimal_bounding_box,
            working_projection="EPSG:28350",
            geology_filename=geology_file,
            dtm_filename=dtm_file,
            structure_filename=structure_file,
            config_filename="dummy_config.json",
            config_dictionary=valid_config_dictionary,
        )
        
# 2) No config_filename or config_dictionary should raise ValueError
def test_no_config_provided_raises_error(
    minimal_bounding_box, geology_file, dtm_file, structure_file
):

    with pytest.raises(ValueError, match="A config file is required to run map2loop"):
        Project(
            bounding_box=minimal_bounding_box,
            working_projection="EPSG:28350",
            geology_filename=geology_file,
            dtm_filename=dtm_file,
            structure_filename=structure_file,
        )

# 3) Passing an unexpected argument should raise TypeError
def test_unexpected_argument_raises_error(
    minimal_bounding_box, geology_file, dtm_file, structure_file, valid_config_dictionary
):
   
    with pytest.raises(TypeError, match="unexpected keyword argument 'config_file'"):
        Project(
            bounding_box=minimal_bounding_box,
            working_projection="EPSG:28350",
            geology_filename=geology_file,
            dtm_filename=dtm_file,
            structure_filename=structure_file,
            config_dictionary=valid_config_dictionary,
            config_file="wrong_kwarg.json", 
        )

# 4) Dictionary missing a required key should raise ValueError

def test_dictionary_missing_required_key_raises_error(
    minimal_bounding_box, geology_file, dtm_file, structure_file
):

    invalid_dictionary = {
        "structure": {"dipdir_column": "azimuth2", "dip_column": "dip"},
        "geology": {"unitname_column": "unitname"}  # alt_unitname_column missing
    }

    with pytest.raises(ValueError, match="Missing required key 'alt_unitname_column' for 'geology'"):
        Project(
            bounding_box=minimal_bounding_box,
            working_projection="EPSG:28350",
            geology_filename=geology_file,
            dtm_filename=dtm_file,
            structure_filename=structure_file,
            config_dictionary=invalid_dictionary,
        )

# 5) All good => The Project should be created without errors
def test_good_config_runs_successfully(
    minimal_bounding_box, geology_file, dtm_file, structure_file, valid_config_dictionary
):
    project = None
    try:
        project = Project(
            bounding_box=minimal_bounding_box,
            working_projection="EPSG:28350",
            geology_filename=geology_file,
            dtm_filename=dtm_file,
            structure_filename=structure_file,
            config_dictionary=valid_config_dictionary,
        )
    except Exception as e:
        pytest.fail(f"Project initialization raised an unexpected exception: {e}")

    assert project is not None, "Project was not created."
    assert project.map_data.config.structure_config["dipdir_column"] == "azimuth2"
    assert project.map_data.config.structure_config["dip_column"] == "dip"
    assert project.map_data.config.geology_config["unitname_column"] == "unitname"
    assert project.map_data.config.geology_config["alt_unitname_column"] == "code"