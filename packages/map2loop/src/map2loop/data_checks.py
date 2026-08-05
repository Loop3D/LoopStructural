# internal imports
from .m2l_enums import Datatype

# external imports
import beartype as beartype
from beartype.typing import Tuple, Optional, List, Dict, Type, Union
import geopandas
import shapely
import pandas

from .logging import getLogger
logger = getLogger(__name__)  

@beartype.beartype
def check_geology_fields_validity(mapdata) -> tuple[bool, str]:
    #TODO (AR) - add check for gaps in geology data (inspo here: https://medium.com/@achm.firmansyah/an-approach-for-checking-overlaps-and-gaps-in-polygons-using-geopandas-ebd6606e7f70 )
    """
    Validate the columns in GEOLOGY geodataframe

    Several checks to ensure that the geology data:
    - Is loaded and valid.
    - Contains required columns with appropriate types and no missing or blank values.
    - Has optional columns with valid types, if present.
    - Does not contain duplicate in IDs.
    - Ensures the geometry column has valid geometries.

    Returns:
        Tuple[bool, str]: A tuple indicating success (False) or failure (True)
    """
    # Check if geology data is loaded and valid
    if (
        mapdata.raw_data[Datatype.GEOLOGY] is None
        or type(mapdata.raw_data[Datatype.GEOLOGY]) is not geopandas.GeoDataFrame
    ):
        logger.error("GEOLOGY data is not loaded or is not a valid GeoDataFrame")
        return (True, "GEOLOGY data is not loaded or is not a valid GeoDataFrame")
    
    geology_data = mapdata.raw_data[Datatype.GEOLOGY]
    config = mapdata.config.geology_config
    
    # 2. Validate geometry
    failed, message = validate_geometry(
        geodata=geology_data,
        expected_geom_types=[shapely.geometry.Polygon, shapely.geometry.MultiPolygon],
        datatype_name="GEOLOGY"
    )
    if failed:
        return (failed, message)
    
    
    # check required columns in geology
    required_columns = ["unitname_column", "alt_unitname_column"]
    
    failed, message = validate_required_columns(
        geodata=geology_data,
        config=config,
        required_columns=required_columns,
        expected_type=str,
        check_blank=True,
        datatype_name="GEOLOGY"
    )
    if failed:
        return (failed, message)
    
    # check optional columns
    optional_string_columns = [
        "group_column", "supergroup_column", "description_column",
        "rocktype_column", "alt_rocktype_column",
    ]
    
    string_warnings = validate_optional_columns(
        geodata=geology_data,
        config=config,
        optional_columns=optional_string_columns,
        expected_type=str,
        check_blank=True,  
        datatype_name="GEOLOGY"
    )
    ### only emit warnings for optional columns
    for warning in string_warnings:
        logger.warning(warning)
    
    # 5. Validate Optional Numeric Columns
    optional_numeric_columns = ["minage_column", "maxage_column"]
    numeric_warnings = validate_optional_columns(
        geodata=geology_data,
        config=config,
        optional_columns=optional_numeric_columns,
        expected_type=(int, float),
        check_blank=False,
        datatype_name="GEOLOGY"
    )
    
    ### only emit warnings for optional columns
    for warning in numeric_warnings:
        logger.warning(warning)
    
    # # 4. check ID column
    if "objectid_column" in config:
        failed, message = validate_id_column(
            geodata=geology_data,
            config=config,
            id_config_key="objectid_column", 
            geodata_name="GEOLOGY")
        
        if failed:
            return (failed, message)

    logger.info("Geology fields validation passed.")
    return (False, "")


@beartype.beartype
def check_structure_fields_validity(mapdata) -> Tuple[bool, str]:
    """
    Validate the structure data for required and optional fields.

    Performs the following checks:
    - Ensures the structure map is loaded, valid, and contains at least two structures.
    - Validates the geometry column
    - Checks required numeric columns (`dip_column`, `dipdir_column`) for existence, dtype, range, and null values.
    - Checks optional string columns (`description_column`, `overturned_column`) for type and null/empty values.
    - Validates the optional numeric `objectid_column` for type, null values, and duplicates.

    Returns:
        Tuple[bool, str]: A tuple where the first value indicates if validation failed (True = failed),
                        and the second value provides a message describing the issue.
    """
    
    # Check type and size of loaded structure map
    if (
        mapdata.raw_data[Datatype.STRUCTURE] is None
        or type(mapdata.raw_data[Datatype.STRUCTURE]) is not geopandas.GeoDataFrame
    ):
        logger.warning("Structure map is not loaded or valid")
        return (True, "Structure map is not loaded or valid")

    if len(mapdata.raw_data[Datatype.STRUCTURE]) < 2:
        logger.warning(
            "Datatype STRUCTURE: map does with not enough orientations to complete calculations (need at least 2), projection may be inconsistent"
        )
    
    structure_data = mapdata.raw_data[Datatype.STRUCTURE]
    config = mapdata.config.structure_config

    # 2. Validate geometry
    failed, message = validate_geometry(
        geodata=structure_data,
        expected_geom_types=[shapely.Point, shapely.MultiPoint],
        datatype_name="STRUCTURE"
    )
    if failed:
        return (failed, message)
    
    
    # check required columns in structure (numeric dips & dip dir)
    required_columns = ["dipdir_column", "dip_column"]
    failed, message = validate_required_columns(
        geodata=structure_data,
        config=config,
        required_columns=required_columns,
        expected_type=(int, float),
        check_blank=False,
        datatype_name="STRUCTURE"
    )
    if failed:
        return (failed, message)

    # 4. Validate Dip and Dip Direction value ranges
    dip_columns = ["dip_column", "dipdir_column"]
    dip_validation_failed, dip_message = validate_dip_columns(
        geodata=structure_data,
        config=config,
        dip_columns=dip_columns,
        datatype_name="STRUCTURE",
        allow_nulls=False  # Dip and dipdir cannot have nulls in structure data
    )
    if dip_validation_failed:
        logger.warning(dip_message)
    
    # check optional columns
    optional_string_columns = ["description_column", "overturned_column"]
    string_warnings = validate_optional_columns(
        geodata=structure_data,
        config=config,
        optional_columns=optional_string_columns,
        expected_type=str,
        check_blank=True,  
        datatype_name="STRUCTURE"
    )
    
    ## only emit warnings for optional columns
    for warning in string_warnings:
        logger.warning(warning)

    # check ID column 
    if "objectid_column" in config:
        failed, id_message = validate_id_column(
            geodata=structure_data,
            config=config,
            id_config_key="objectid_column", 
            geodata_name="STRUCTURE")
        
        if failed:
            return (failed, id_message)
        
    return (False, "")

@beartype.beartype
def check_fault_fields_validity(mapdata) -> Tuple[bool, str]:
    
    # Check type of loaded fault map
    if (
        mapdata.raw_data[Datatype.FAULT] is None
        or type(mapdata.raw_data[Datatype.FAULT]) is not geopandas.GeoDataFrame
    ):
        logger.warning("Fault map is not loaded or valid")
        return (True, "Fault map is not loaded or valid")
    
    fault_data = mapdata.raw_data[Datatype.FAULT]
    config = mapdata.config.fault_config
    
    # 2. Validate geometry
    failed, message = validate_geometry(
        geodata=fault_data,
        expected_geom_types=[shapely.LineString, shapely.MultiLineString],
        datatype_name="FAULT"
    )
    if failed:
        return (failed, message)
    
    # # Check "structtype_column" if it exists
    text_keys = {
        "fault_text": "fault_text"
    }
    structtype_validation_failed, structtype_message = validate_structtype_column(
        geodata=fault_data,
        config=config,
        datatype_name="FAULT",
        required=True,  # Assuming structtype_column is required in FAULT
        text_keys=text_keys
    )
    if structtype_validation_failed:
        return (structtype_validation_failed, structtype_message)
    
    
    
    #checks on name column
    name_column = config.get("name_column")
    if name_column not in fault_data.columns:
        logger.warning(
            f"Datatype FAULT: Column '{name_column}' (config key 'name_column') is missing from the fault data."
            "Please ensure it is present, or remove that key from the config."
        )
    
    if name_column and name_column in fault_data.columns:
        # Check if the column contains non-string values
        if not fault_data[name_column].apply(lambda x: isinstance(x, str) or pandas.isnull(x)).all():
            logger.error(
                f"Datatype FAULT: Column '{name_column}' (config key 'name_column') contains non-string values. Ensure all values are valid strings."
            )
            return (True, f"Datatype FAULT: Column '{name_column}' (config key 'name_column') contains non-string values.")
        
        # Check for NaN values
        if fault_data[name_column].isnull().any():
            logger.warning(
                f"Datatype FAULT: Column '{name_column}' (config key 'name_column') contains NaN or empty values. This may affect processing."
            )
        
        # Check for duplicate values
        if fault_data[name_column].duplicated().any():
            logger.warning(
                f"Datatype FAULT: Column '{name_column}' contains duplicate values. This may affect processing."
            )

    # # dips & strikes
    dip_columns = ["dip_column", "dipdir_column"]
    dip_validation_failed, dip_message = validate_dip_columns(
        geodata=fault_data,
        config=config,
        dip_columns=dip_columns,
        datatype_name="FAULT",
        allow_nulls=True  # Dip fields can be empty
    )
    if dip_validation_failed:
        logger.warning(dip_message)
    
    # dip estimates
    dip_estimate_column = config.get("dip_estimate_column")
    valid_directions = [
        "north_east", "south_east", "south_west", "north_west",
        "north", "east", "south", "west"
    ]

    if dip_estimate_column:  
        if dip_estimate_column in fault_data.columns:
            # Ensure all values are in the set of valid directions or are NaN
            invalid_values = fault_data[dip_estimate_column][
                ~fault_data[dip_estimate_column].apply(lambda x: x in valid_directions or pandas.isnull(x))
            ]

            if not invalid_values.empty:
                logger.error(
                    f"Datatype FAULT: Column '{dip_estimate_column}' contains invalid values not in the set of allowed dip estimates: {valid_directions}."
                )
                return (
                    True,
                    f"Datatype FAULT: Column '{dip_estimate_column}' contains invalid values. Allowed values: {valid_directions}.",
                )

            # Warn if there are NaN or empty values
            if fault_data[dip_estimate_column].isnull().any():
                logger.warning(
                    f"Datatype FAULT: Column '{dip_estimate_column}' contains NaN or empty values. This may affect processing."
                )
        else:
            logger.error(
                f"Datatype FAULT: Column '{dip_estimate_column}' is missing from the fault data. Please ensure the column name is correct or remove that key from the config."
            )
            return (True, f"Datatype FAULT: Column '{dip_estimate_column}' is missing from the fault data.")

    
    # # 4. check ID column
    if "objectid_column" in config:
        id_validation_failed, id_message = validate_id_column(
            geodata=fault_data,
            config=config,
            id_config_key="objectid_column", 
            geodata_name="FAULT")
        
        if id_validation_failed:
            return (id_validation_failed, id_message)
    
    return (False, "")

@beartype.beartype
def check_fold_fields_validity(mapdata) -> Tuple[bool, str]:
    # Check type of loaded fold map
    if (
        mapdata.raw_data[Datatype.FOLD] is None
        or type(mapdata.raw_data[Datatype.FOLD]) is not geopandas.GeoDataFrame
    ):
        logger.warning("Fold map is not loaded or valid")
        return (True, "Fold map is not loaded or valid")

    folds = mapdata.raw_data[Datatype.FOLD]
    config = mapdata.config.fold_config

    # Debugging: Print column names in the fold_data
    logger.debug(f"Fold data columns: {folds.columns.tolist()}")

    # 2. Validate geometry
    failed, message = validate_geometry(
        geodata=folds,
        expected_geom_types=[shapely.LineString, shapely.MultiLineString],
        datatype_name="FOLD"
    )
    if failed:
        return (failed, message)
    
    ## check structtype column if it exists
    text_keys = {
        "fold_text": "fold_text",
        "synform_text": "synform_text"
    }
    structtype_validation_failed, structtype_message = validate_structtype_column(
        geodata=folds,
        config=config,
        datatype_name="FOLD",
        required=True,  # Assuming structtype_column is required in FOLD
        text_keys=text_keys
    )
    if structtype_validation_failed:
        return (structtype_validation_failed, structtype_message)
                        
    # check description column
    description_column = config.get("description_column", None)
    if description_column:
        # Ensure the column exists in the data
        if description_column not in folds.columns:
            logger.warning(
                f"Datatype FOLD: Column '{description_column}' (config key: 'description_column') is missing from the fold data. Consider removing that key from the config."
            )
        else:
            # Check if all entries in the column are strings
            if not folds[description_column].apply(lambda x: isinstance(x, str) or pandas.isnull(x)).all():
                logger.error(
                    f"Datatype FOLD: Column '{description_column}' (config key: 'description_column') contains non-string values. Please ensure all values in this column are strings."
                )
                return (True, f"Datatype FOLD: Column '{description_column}' (config key: 'description_column') contains non-string values.")

            # Warn about empty or null cells
            if folds[description_column].isnull().any() or folds[description_column].str.strip().eq("").any():
                logger.warning(
                    f"Datatype FOLD: Column '{description_column}' contains NaN, empty, or blank values. Processing might not work as expected."
                )


    # # 4. check ID column
    if "objectid_column" in config:
        id_validation_failed, id_message = validate_id_column(
            geodata=folds,
            config=config,
            id_config_key="objectid_column", 
            geodata_name="FOLD")
        
        if id_validation_failed:
            return (id_validation_failed, id_message)

    return (False, "")


@beartype.beartype
def validate_config_dictionary(config_dict: dict) -> None:
    
    # 1)  check mandatory keys for "structure" and "geology"
    required_keys = {
        "structure": {"dipdir_column", "dip_column"},
        "geology": {"unitname_column", "alt_unitname_column"},
    }

    # Loop over "structure" and "geology"
    for section, keys in required_keys.items():

        # 1) Check that "section" exists
        if section not in config_dict:
            logger.error(f"Missing required section '{section}' in config dictionary.")
            raise ValueError(f"Missing required section '{section}' in config dictionary.")

        # 2) Check that each required key is in config_dict[section]
        for key in keys:
            if key not in config_dict[section]:
                logger.error(f"Missing required key '{key}' for '{section}' section of the config dictionary.")
                raise ValueError(f"Missing required key '{key}' for '{section}' section of the config dictionary.")
    
    # 2) check for legacy keys first:
    legacy_keys = {
        "otype", "dd", "d", "sf", "bedding", "bo", "btype", "gi", "c", "u",
        "g", "g2", "ds", "min", "max", "r1", "r2", "sill", "intrusive", "volcanic",
        "f", "fdipnull", "fdipdip_flag", "fdipdir", "fdip", "fdipest",
        "fdipest_vals", "n", "ff", "t", "syn"
    }

    def check_keys(d: dict, parent_key=""):
        for key, value in d.items():
            if key in legacy_keys:
                logger.error(
                    f"Legacy key found in config - '{key}' at '{parent_key}'. Please use the new config format. Use map2loop.utils.update_from_legacy_file to convert between the formats if needed"
                )
                raise ValueError(
                    f"Legacy key found in config - '{key}' at '{parent_key}'. Please use the new config format. Use map2loop.utils.update_from_legacy_file to convert between the formats if needed"
                )
            if isinstance(value, dict):
                check_keys(value, parent_key=f"{parent_key}{key}.")
    
    check_keys(config_dict)

    # 3) check if all keys are valid:
    allowed_keys_by_section = {
        "structure": {
            "orientation_type", "dipdir_column", "dip_column",
            "description_column", "bedding_text", "overturned_column", "overturned_text",
            "objectid_column", "desciption_column", "interp_source_column"
        },
        "geology": {
            "unitname_column", "alt_unitname_column", "group_column",
            "supergroup_column", "description_column", "minage_column",
            "maxage_column", "rocktype_column",  "alt_rocktype_column",
            "sill_text", "intrusive_text",  "volcanic_text",   "objectid_column", "ignore_lithology_codes",
        },
        "fault": {
            "structtype_column",  "fault_text",  "dip_null_value",
            "dipdir_flag", "dipdir_column",  "dip_column",  "orientation_type",
            "dipestimate_column",  "dipestimate_text","displacement_column", 
            "displacement_text",  "name_column","objectid_column", "featureid_column", "minimum_fault_length",
            "fault_length_column", "fault_length_text", "ignore_fault_codes",
        },
        "fold": {
            "structtype_column", "fold_text", "axial_plane_dipdir_column", 
            "axial_plane_dip_column","tightness_column", "description_column",
            "synform_text", "foldname_column","objectid_column",
        },
    }
    
    for section_name, section_dict in config_dict.items():
        # check section
        if section_name not in allowed_keys_by_section:
            logger.error(f"Unrecognized section '{section_name}' in config dictionary.")
            raise ValueError(f"Unrecognized section '{section_name}' in config dictionary.")

        # check keys
        allowed_keys = allowed_keys_by_section[section_name]
        for key in section_dict.keys():
            if key not in allowed_keys:
                logger.error(f"Key '{key}' is not an allowed key in the '{section_name}' section.")
                raise ValueError(f"Key '{key}' is not an allowed key in the '{section_name}' section.")
    
    # 4) check if minimum fault length is a number
    mfl = config_dict.get("fault", {}).get("minimum_fault_length", None)
    if mfl is not None and not isinstance(mfl, (int, float)):
        logger.error("minimum_fault_length must be a number.")
        raise ValueError(f"minimum_fault_length must be a number, instead got: {type(mfl)}")

@beartype.beartype
def validate_geometry(
    geodata: geopandas.GeoDataFrame,
    expected_geom_types: List[type],
    datatype_name: str
) -> Tuple[bool, str]:
    geodata.geometry = geodata.geometry.make_valid()
    # 1. Check if all geometries are valid
    if not geodata.geometry.is_valid.all():
        logger.error(f"Invalid geometries found in datatype {datatype_name}. Please fix them before proceeding.")
        return True, f"Invalid geometries found in datatype {datatype_name}"
    
    # 2. Check if all geometries are of the expected types
    if not geodata.geometry.apply(lambda geom: isinstance(geom, tuple(expected_geom_types))).all():
        invalid_types = geodata[~geodata.geometry.apply(lambda geom: isinstance(geom, tuple(expected_geom_types)))]
        invalid_indices = invalid_types.index.tolist()
        expected_types_names = ', '.join([geom_type.__name__ for geom_type in expected_geom_types])
        logger.error(
            f"Datatype {datatype_name}: Invalid geometry types found. Expected types: {expected_types_names}. "
            f"Rows with invalid types: {invalid_indices}"
        )
        return True, (
            f"Invalid geometry types found in datatype {datatype_name}. "
            f"All geometries must be {expected_types_names}."
        )
    
    # If all checks pass
    logger.debug(f"Geometry validation passed for datatype {datatype_name}")
    return False, ""


@beartype.beartype
def validate_id_column(
    geodata: geopandas.GeoDataFrame,
    config: dict,
    id_config_key: str, 
    geodata_name: str
) -> Tuple[bool, str]:

    # Retrieve the ID column name from the configuration
    id_column = config.get(id_config_key)
    
    if not id_column:
        error_msg = f"Configuration key '{id_config_key}' is missing."
        logger.error(error_msg)
        return (True, error_msg)
    
    if id_column in geodata.columns:
        geodata[id_column] = pandas.to_numeric(geodata[id_column], errors='coerce')
    
        # Check for non-numeric values (which are now NaN after coercion)
        if geodata[id_column].isnull().any():
            error_msg = (
                f"Datatype {geodata_name}: Column '{id_column}' "
                f"(config key: '{id_config_key}') contains non-numeric or NaN values. "
                "Please rectify the values, or remove this key from the config dictionary to let map2loop assign IDs."
            )
            logger.error(error_msg)
            return (True, error_msg)
        
        if not (geodata[id_column] == geodata[id_column].astype(int)).all():
            error_msg = (
                f"Datatype {geodata_name}: Column '{id_column}' "
                f"(config key: '{id_config_key}') contains non-integer values."
            )
            logger.error(error_msg)
            return (True, error_msg)

        if geodata[id_column].duplicated().any():
            error_msg = (
                f"Datatype {geodata_name}: Column '{id_column}' "
                f"(config key: '{id_config_key}') contains duplicate values."
            )
            logger.error(error_msg)
            return (True, error_msg)
    
    
    elif id_column not in geodata.columns:
        msg = (
            f"Datatype {geodata_name}: Column '{id_column}' "
            f"(config key: '{id_config_key}') is missing from the data. "
            "Map2loop will automatically generate IDs."
        )
        logger.warning(msg)

    return (False, "")

@beartype.beartype
def validate_required_columns(
    geodata: geopandas.GeoDataFrame,
    config: dict,
    required_columns: List[str],
    expected_type: Union[Type, Tuple[Type, ...]],
    check_blank: bool = False,
    datatype_name: str = "UNKNOWN"
) -> Tuple[bool, str]:

    for config_key in required_columns:
        column_name = config.get(config_key)
        
        if not column_name:
            error_msg = (
                f"Configuration key '{config_key}' is missing for datatype '{datatype_name}'."
            )
            logger.error(error_msg)
            return (True, error_msg)
        
        if column_name not in geodata.columns:
            error_msg = (
                f"Datatype {datatype_name.upper()}: Required column with config key '{config_key}' "
                f"(column: '{column_name}')  is missing from the data."
            )
            logger.error(error_msg)
            return (True, error_msg)
        
        # Check data type
        if not geodata[column_name].apply(lambda x: isinstance(x, expected_type)).all():
            error_msg = (
                f"Datatype {datatype_name.upper()}: Column '{config_key}' (column: '{column_name}') "
                f"must contain only {expected_type if isinstance(expected_type, type) else 'numeric'} values."
            )
            logger.error(error_msg)
            return (True, error_msg)
        
        # Check for null values
        if geodata[column_name].isnull().any():
            error_msg = (
                f"Datatype {datatype_name.upper()}: Column '{config_key}' (column: '{column_name}')  "
                f"contains null values. Please ensure all values are present."
            )
            logger.error(error_msg)
            return (True, error_msg)
        
        # Optionally check for blank strings
        if check_blank and issubclass(expected_type, str):
            if geodata[column_name].str.strip().eq("").any():
                error_msg = (
                    f"Datatype {datatype_name.upper()}: Column '{config_key}' (column: '{column_name}') "
                    f"contains blank (empty) values. Please ensure all values are populated."
                )
                logger.error(error_msg)
                return (True, error_msg)
    
    # If all required columns pass validation
    logger.info(f"Datatype {datatype_name.upper()}: All required columns validated successfully.")
    return (False, "")


def validate_optional_columns(
    geodata: geopandas.GeoDataFrame,
    config: Dict[str, str],
    optional_columns: List[str],
    expected_type: Union[Type, Tuple[Type, ...]],
    check_blank: bool = False,
    datatype_name: str = "UNKNOWN"
) -> List[str]:

    warnings = []

    for config_key in optional_columns:
        column_name = config.get(config_key)

        if not column_name:
            warning_msg = (
                f"Configuration key '{config_key}' is missing for datatype '{datatype_name}'. "
                f"Optional column validation for this key is skipped."
            )
            logger.warning(warning_msg)
            warnings.append(warning_msg)
            continue  

        if column_name in geodata.columns:
            # Type Check
            if not geodata[column_name].apply(lambda x: isinstance(x, expected_type) or pandas.isnull(x)).all():
                warning_msg = (
                    f"Datatype {datatype_name.upper()}: Optional column '{column_name}' "
                    f"(config key: '{config_key}') contains values that are not of type "
                    f"{expected_type if isinstance(expected_type, type) else expected_type}. "
                    "Map2loop processing might not work as expected."
                )
                logger.warning(warning_msg)
                warnings.append(warning_msg)

            # Blank String Check (if applicable)
            if check_blank and issubclass(expected_type, str):
                if geodata[column_name].str.strip().eq("").any():
                    warning_msg = (
                        f"Datatype {datatype_name.upper()}: Optional column '{column_name}' "
                        f"(config key: '{config_key}') contains blank (empty) string values. "
                        "Map2loop processing might not work as expected."
                    )
                    logger.warning(warning_msg)
                    warnings.append(warning_msg)

            # Null Value Check
            if geodata[column_name].isnull().any():
                warning_msg = (
                    f"Datatype {datatype_name.upper()}: Optional column '{column_name}' "
                    f"(config key: '{config_key}') contains NaN or null values. "
                    "Map2loop processing might not work as expected."
                )
                logger.warning(warning_msg)
                warnings.append(warning_msg)

        else:
            info_msg = (
                f"Datatype {datatype_name.upper()}: Optional column '{column_name}' "
                f"(config key: '{config_key}') is missing from the data. "
            )
            ###### this might be taking it a bit too far

            logger.info(info_msg)

    return warnings


@beartype.beartype
def validate_dip_columns(
    geodata: geopandas.GeoDataFrame,
    config: Dict[str, str],
    dip_columns: List[str],
    datatype_name: str = "UNKNOWN",
    allow_nulls: bool = False
) -> Tuple[bool, str]:

    validation_failed = False
    messages = []
    
    # Define fixed ranges
    fixed_ranges = {
        "dip_column": (0, 90),
        "dipdir_column": (0, 360)
    }
    
    for key in dip_columns:
        column_name = config.get(key)
        if not column_name and datatype_name == "STRUCTURE": # only mandatory for structure, not faults!
            warning_msg = (
                f"Configuration key '{key}' is missing for datatype '{datatype_name}'. "
                f"Dip column validation for this key is skipped."
            )
            logger.warning(warning_msg)
            messages.append(warning_msg)
            validation_failed = True
            continue

        if column_name in geodata.columns:
            # Coerce to numeric
            geodata[column_name] = pandas.to_numeric(geodata[column_name], errors='coerce')

            # Check for non-numeric or NaN values
            if geodata[column_name].isnull().any():
                if not allow_nulls:
                    warning_msg = (
                        f"Datatype {datatype_name.upper()}: Column '{column_name}' "
                        f"(config key: '{key}') contains non-numeric or NaN values."
                    )
                    logger.warning(warning_msg)
                    messages.append(warning_msg)
                    validation_failed = True

            # Check if all values are numeric
            if not geodata[column_name].apply(lambda x: isinstance(x, (int, float)) or pandas.isnull(x)).all():
                warning_msg = (
                    f"Datatype {datatype_name.upper()}: Column '{column_name}' "
                    f"(config key: '{key}') must contain only numeric values."
                )
                logger.warning(warning_msg)
                messages.append(warning_msg)
                validation_failed = True

            # Range validation
            min_val, max_val = fixed_ranges.get(key, (None, None))
            if min_val is not None and max_val is not None:
                invalid_values = ~geodata[column_name].between(min_val, max_val, inclusive='both')
                if invalid_values.any():
                    warning_msg = (
                        f"Datatype {datatype_name.upper()}: Column '{column_name}' "
                        f"(config key: '{key}') contains values outside the range [{min_val}, {max_val}]. "
                        "Is this intentional?"
                    )
                    logger.warning(warning_msg)
                    messages.append(warning_msg)

    summary_message = "\n".join(messages)
    return (validation_failed, summary_message)


@beartype.beartype
def validate_structtype_column(
    geodata: geopandas.GeoDataFrame,
    config: Dict[str, str],
    datatype_name: str,
    required: bool = True,
    text_keys: Optional[Dict[str, str]] = None
) -> Tuple[bool, str]:

    structtype_key = "structtype_column"
    structtype_column = config.get(structtype_key)

    if not structtype_column:
        if required:
            error_msg = (
                f"Configuration key '{structtype_key}' is missing for datatype '{datatype_name}'. "
                f"Validation for 'structtype_column' is skipped."
            )
            logger.warning(error_msg)
            return (True, error_msg)
        else:
            warning_msg = (
                f"Configuration key '{structtype_key}' is missing for datatype '{datatype_name}'. "
                f"Optional 'structtype_column' validation is skipped."
            )
            logger.warning(warning_msg)
            return (False, "")
    
    if structtype_column not in geodata.columns:
        if required:
            error_msg = (
                f"Datatype {datatype_name.upper()}: '{structtype_column}' (config key: '{structtype_key}') "
                f"is missing from the data. Consider removing that key from the config."
            )
            logger.error(error_msg)
            return (True, error_msg)
        else:
            warning_msg = (
                f"Datatype {datatype_name.upper()}: '{structtype_column}' (config key: '{structtype_key}') "
                f"is missing from the data. Consider removing that key from the config."
            )
            logger.warning(warning_msg)
            return (False, "")
    
    # Check if all entries are strings or nulls
    if not geodata[structtype_column].apply(lambda x: isinstance(x, str) or pandas.isnull(x)).all():
        error_msg = (
            f"Datatype {datatype_name.upper()}: Column '{structtype_column}' "
            f"(config key: '{structtype_key}') contains non-string values. "
            "Please ensure all values in this column are strings."
        )
        logger.error(error_msg)
        return (True, error_msg)
    
    # Warn about empty or null cells
    if geodata[structtype_column].isnull().any() or geodata[structtype_column].str.strip().eq("").any():
        warning_msg = (
            f"Datatype {datatype_name.upper()}: Column '{structtype_column}' contains NaN, empty, or blank values. "
            "Processing might not work as expected."
        )
        logger.warning(warning_msg)
    
    # Check for specific text keys
    if text_keys:
        for text_key, config_key in text_keys.items():
            text_value = config.get(config_key, None)
            text_pattern = '|'.join(text_value) if text_value else None
            if text_value:
                if not isinstance(text_value, str):
                    error_msg = (
                        f"Datatype {datatype_name.upper()}: '{config_key}' must be a string. "
                        "Please ensure it is defined correctly in the config."
                    )
                    logger.error(error_msg)
                    return (True, error_msg)
                
                if not geodata[structtype_column].str.contains(text_pattern, case=False, regex=True, na=False).any():
                    if text_key == "synform_text":
                        warning_msg = (
                            f"Datatype {datatype_name.upper()}: The '{text_key}' value '{text_value}' is not found in column '{structtype_column}'. "
                            "This may impact processing."
                        )
                        logger.warning(warning_msg)
                    else:
                        error_msg = (
                            f"Datatype {datatype_name.upper()}: The '{text_key}' value '{text_value}' is not found in column '{structtype_column}'. "
                            "Project might end up with no faults."
                        )
                        logger.error(error_msg)
                        return (True, error_msg)
    
    return (False, "")
