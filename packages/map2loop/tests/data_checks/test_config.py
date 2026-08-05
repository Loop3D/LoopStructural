import pytest
from map2loop.data_checks import validate_config_dictionary


@pytest.fixture
def valid_config():
    return {
        "structure": {
            "orientation_type": "dip direction",
            "dipdir_column": "azimuth",
            "dip_column": "inclinatn",
            "description_column": "DESCRIPTION",
            "bedding_text": "bed",
            "overturned_column": "no_col",
            "overturned_text": "blah",
            "objectid_column": "geographic",
            "desciption_column": "sub_type"
        },
        "geology": {
            "unitname_column": "formatted_",
            "alt_unitname_column": "abbreviate",
            "group_column": "no_col",
            "supergroup_column": "interpreta",
            "description_column": "text_descr",
            "minage_column": "no_col",
            "maxage_column": "no_col",
            "rocktype_column": "rank",
            "alt_rocktype_column": "type",
            "sill_text": "sill",
            "intrusive_text": "intrusion",
            "volcanic_text": "volc",
            "objectid_column": "ID",
            "ignore_lithology_codes": ["cover"]
        },
        "fault": {
            "structtype_column": "featuretyp",
            "fault_text": "s",
            "dip_null_value": "0",
            "dipdir_flag": "num",
            "dipdir_column": "no_col",
            "dip_column": "no_col",
            "orientation_type": "dip direction",
            "dipestimate_column": "no_col",
            "dipestimate_text": "no_col",
            "name_column": "no_col",
            "objectid_column": "geographic",
            "minimum_fault_length": 100.0, 
            "ignore_fault_codes": []
        },
        "fold": {
            "structtype_column": "featuretyp",
            "fold_text": "fold",
            "description_column": "no_col",
            "synform_text": "syn",
            "foldname_column": "NAME",
            "objectid_column": "geographic"
        }
    }


def test_valid_config_no_errors(valid_config):
    # Should not raise any error
    validate_config_dictionary(valid_config)


def test_missing_required_section(valid_config):

    config_missing_structure = dict(valid_config)
    del config_missing_structure["structure"]  # remove required section

    with pytest.raises(ValueError) as exc_info:
        validate_config_dictionary(config_missing_structure)
    assert "Missing required section 'structure'" in str(exc_info.value)


def test_missing_required_key(valid_config):
    
    config_missing_dip = dict(valid_config)
    
    del config_missing_dip["structure"]["dip_column"] # remove required key

    with pytest.raises(ValueError) as exc_info:
        validate_config_dictionary(config_missing_dip)
    assert "Missing required key 'dip_column' for 'structure'" in str(exc_info.value)


def test_unrecognized_section(valid_config):

    config_extra_section = dict(valid_config)
    config_extra_section["random_section"] = {"random_key": "random_value"}

    with pytest.raises(ValueError) as exc_info:
        validate_config_dictionary(config_extra_section)
    assert "Unrecognized section 'random_section'" in str(exc_info.value)


def test_unrecognized_key_in_section(valid_config):
    
    config_extra_key = dict(valid_config)
    config_extra_key["structure"]["random_key"] = "random_value"

    with pytest.raises(ValueError) as exc_info:
        validate_config_dictionary(config_extra_key)
    assert "Key 'random_key' is not an allowed key in the 'structure' section." in str(exc_info.value)


def test_legacy_key_detected(valid_config):

    config_with_legacy = dict(valid_config)
    config_with_legacy["structure"]["otype"] = "legacy_value"  # 'otype' --> legacy key
    with pytest.raises(ValueError) as exc_info:
        validate_config_dictionary(config_with_legacy)
    assert "Legacy key found in config - 'otype'" in str(exc_info.value)


def test_minimum_fault_length_wrong_type(valid_config):

    config_wrong_mfl = dict(valid_config)
    config_wrong_mfl["fault"]["minimum_fault_length"] = "one_hundred"  # invalid type

    with pytest.raises(ValueError) as exc_info:
        validate_config_dictionary(config_wrong_mfl)
    assert "minimum_fault_length must be a number" in str(exc_info.value)


def test_minimum_fault_length_missing(valid_config):
    """
    Remove minimum_fault_length entirely. That should be fine (None -> no check).
    """
    config_no_mfl = dict(valid_config)
    del config_no_mfl["fault"]["minimum_fault_length"]

    # Should not raise any error, as it's optional
    validate_config_dictionary(config_no_mfl)

