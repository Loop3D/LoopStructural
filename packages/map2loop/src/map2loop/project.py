# internal imports
from map2loop.fault_orientation import FaultOrientationNearest
from .utils import hex_to_rgb, set_z_values_from_raster_df
from .m2l_enums import VerboseLevel, ErrorState, Datatype
from .mapdata import MapData
from .contact_extractor import ContactExtractor
from .sampler import Sampler, SamplerDecimator, SamplerSpacing
from .thickness_calculator import InterpolatedStructure, ThicknessCalculator
from .throw_calculator import ThrowCalculator, ThrowCalculatorAlpha
from .fault_orientation import FaultOrientation
from .sorter import Sorter, SorterAgeBased, SorterAlpha, SorterUseNetworkX, SorterUseHint
from .stratigraphic_column import StratigraphicColumn
from .deformation_history import DeformationHistory
from .topology import Topology
from .data_checks import validate_config_dictionary

# external imports
import LoopProjectFile as LPF
from osgeo import gdal
gdal.UseExceptions()
import geopandas
import beartype
from beartype.typing import Union, List, Dict, Any
import pathlib
import numpy
import pandas
import os
import re

from .logging import getLogger

logger = getLogger(__name__)


class Project(object):
    """
    The main entry point into using map2loop

    Attributes
    -----------
    verbose_level: m2l_enums.VerboseLevel
        A selection that defines how much console logging is output
    samplers: Sampler
        A list of samplers used to extract point samples from polygons or line segments. Indexed by m2l_enum
    sorter: Sorter
        The sorting algorithm to use for calculating the stratigraphic column
    loop_filename: str
        The name of the loop project file used in this project
    map_data: MapData
        The structure that holds all map and dtm data
    map2model: Topology
        A wrapper around the map2model module that extracts unit and fault adjacency
    stratigraphic_column: StratigraphicColumn
        The structure that holds the unit information and ordering
    deformation_history: DeformationHistory
        The structura that holds the fault and fold information and interactions
    """

    @beartype.beartype
    def __init__(
        self,
        verbose_level: VerboseLevel = VerboseLevel.ALL,
        tmp_path: str = "m2l_data_tmp",
        working_projection=None,
        bounding_box=None,
        use_australian_state_data: str = "",
        geology_filename: str = "",
        structure_filename: str = "",
        fault_filename: str = "",
        fault_orientation_filename: str = "",
        fold_filename: str = "",
        dtm_filename: str = "",
        config_filename: Union[pathlib.Path, str] = "",
        config_dictionary: dict = {},
        clut_filename: Union[pathlib.Path, str] = "",
        save_pre_checked_map_data: bool = False,
        loop_project_filename: str = "",
        overwrite_loopprojectfile: bool = False,
    ):  
        """
        The initialiser for the map2loop project

        Args:
            verbose_level (VerboseLevel, optional):
                How much console output is sent. Defaults to VerboseLevel.ALL.
            tmp_path (str, optional):
                The directory for storing temporary files. Defaults to "./m2l_data_tmp".
            working_projection (str or int, optional):
                The EPSG projection to use. Defaults to None.
            bounding_box (list or dict, optional):
                The boundary extents of the project (m). Defaults to None.
            use_australian_state_data (str, optional):
                Whether to use default state data and which state to use. Defaults to "".
            geology_filename (str, optional):
                The filename of the geology shapefile. Defaults to "".
            structure_filename (str, optional):
                The filename of the structure shapefile. Defaults to "".
            fault_filename (str, optional):
                The filename of the fault shapefile. Defaults to "".
            fault_orientation_filename (str, optional):
                The filename of the fault orientation shapefile. Defaults to "".
            fold_filename (str, optional):
                The filename fo the fold shapefile. Defaults to "".
            dtm_filename (str, optional):
                The filename of the digital terrain map to use. Defaults to "".
            config_filename (str, optional):
                The filename of the configuration json file to use (if not using config_dictionary). Defaults to "".
            config_dictionary (dict, optional):
                A dictionary version of the configuration file. Defaults to {}.
            clut_filename (str, optional):
                The filename of the colour look up table to use. Defaults to "".
            save_pre_checked_map_data (bool, optional):
                A flag to save all map data to file before use. Defaults to False.
            loop_project_filename (str, optional):
                The filename of the loop project file. Defaults to "".

        Raises:
            TypeError: Type of working_projection not a str or int
            ValueError: bounding_box not length 4 or 6
            TypeError: Type of bounding_box not a dict or tuple
            ValueError: use_australian_state_data not in state list ['WA', 'SA', 'QLD', 'NSW', 'TAS', 'VIC', 'ACT', 'NT']
        """
        
        # make sure all the needed arguments are provided
        if not use_australian_state_data: # this check has to skip if using Loop server data
            self.validate_required_inputs(
                bounding_box=bounding_box,
                working_projection=working_projection,
                geology_filename=geology_filename,
                structure_filename=structure_filename,
                dtm_filename=dtm_filename,
                config_dictionary=config_dictionary,
                config_filename=config_filename,
            )
        
        self._error_state = ErrorState.NONE
        self._error_state_msg = ""
        self.verbose_level = verbose_level
        self.samplers = [SamplerDecimator()] * len(Datatype)
        self.set_default_samplers()
        self.bounding_box = bounding_box
        self.contact_extractor = None
        self.sorter = SorterUseHint
        self.throw_calculator = ThrowCalculatorAlpha()
        self.fault_orientation = FaultOrientationNearest()
        self.map_data = MapData(verbose_level=verbose_level)
        self.stratigraphic_column = StratigraphicColumn()
        self.deformation_history = DeformationHistory(project=self)
        self.loop_filename = loop_project_filename
        self.overwrite_lpf = overwrite_loopprojectfile
        self.active_thickness = None
        
        
        # initialise the dataframes to store data in
        self.fault_orientations = pandas.DataFrame(
            columns=["ID", "DIPDIR", "DIP", "X", "Y", "Z", "featureId"]
        )
        self.fault_samples = pandas.DataFrame(columns=["ID", "X", "Y", "Z", "featureId"])
        self.fold_samples = pandas.DataFrame(columns=["ID", "X", "Y", "Z", "featureId"])
        self.geology_samples = pandas.DataFrame(columns=["ID", "X", "Y", "Z", "featureId"])

        # Sanity check on working projection parameter
        if issubclass(type(working_projection), str) or issubclass(type(working_projection), int):
            self.map_data.set_working_projection(working_projection)
        elif type(working_projection) is None:
            if verbose_level != VerboseLevel.NONE:
                logger.warning(
                    "No working projection set, will attempt to use the projection of the geology map"
                )
        else:
            logger.error(f"Invalid type for working_projection {type(working_projection)}")
            raise TypeError(f"Invalid type for working_projection {type(working_projection)}")

        # Sanity check bounding box
        if issubclass(type(bounding_box), dict) or issubclass(type(bounding_box), tuple):
            if len(bounding_box) == 4 or len(bounding_box) == 6:
                self.map_data.set_bounding_box(bounding_box)
            else:
                logger.error(
                    f"Length of bounding_box {len(bounding_box)} is neither 4 (map boundary) nor 6 (volumetric boundary)"
                )
                raise ValueError(
                    f"Length of bounding_box {len(bounding_box)} is neither 4 (map boundary) nor 6 (volumetric boundary)"
                )
        else:
            logger.error(f"Invalid type for bounding_box {type(bounding_box)}")
            raise TypeError(f"Invalid type for bounding_box {type(bounding_box)}")

        # Assign filenames
        if use_australian_state_data != "":
            # Sanity check on state string
            if use_australian_state_data in {"WA", "SA", "QLD", "NSW", "TAS", "VIC", "ACT", "NT"}:
                self.map_data.set_filenames_from_australian_state(use_australian_state_data)
            else:
                logger.error(
                    f"Australian state {use_australian_state_data} not in state url database"
                )
                raise ValueError(
                    f"Australian state {use_australian_state_data} not in state url database"
                )
        # set the data filenames
        if geology_filename != "":
            self.map_data.set_filename(Datatype.GEOLOGY, geology_filename)
        if structure_filename != "":
            self.map_data.set_filename(Datatype.STRUCTURE, structure_filename)
        if fault_filename != "":
            self.map_data.set_filename(Datatype.FAULT, fault_filename)
        if fold_filename != "":
            self.map_data.set_filename(Datatype.FOLD, fold_filename)
        if dtm_filename != "":
            self.map_data.set_filename(Datatype.DTM, dtm_filename)
        if fault_orientation_filename != "":
            self.map_data.set_filename(Datatype.FAULT_ORIENTATION, fault_orientation_filename)
       
        if config_filename != "":
            self.map_data.set_config_filename(config_filename)

        if config_dictionary != {}:
            validate_config_dictionary(config_dictionary)
            self.map_data.config.update_from_dictionary(config_dictionary)
            # print(self.map_data.config)
            # self.map_data.config.validate_config_dictionary(config_dictionary)
            
        if clut_filename != "":
            self.map_data.set_colour_filename(clut_filename)
            
        
        # Load all data (both shape and raster)
        self.map_data.load_all_map_data()

        # If flag to save out data is check do so
        tmp_path = pathlib.Path(tmp_path)

        if save_pre_checked_map_data:
            # check if the path exists, and if not, create
            if not tmp_path.exists():
                tmp_path.mkdir()
            self.map_data.save_all_map_data(tmp_path)

        # Populate the stratigraphic column and deformation history from map data
        self.stratigraphic_column.populate(self.map_data.get_map_data(Datatype.GEOLOGY))
        self.deformation_history.populate(self.map_data.get_map_data(Datatype.FAULT))
        self.topology = Topology(
            self.map_data.get_map_data(Datatype.GEOLOGY), 
            self.map_data.get_map_data(Datatype.FAULT)
            )
        self.thickness_calculator = [InterpolatedStructure(
            dtm_data=self.map_data.get_map_data(Datatype.DTM),
            bounding_box=self.bounding_box,
            is_strike=False
        )]


    @beartype.beartype
    def validate_required_inputs(
        self,
        bounding_box: Dict[str, Union[float, int]],
        working_projection: str,
        geology_filename: str,
        structure_filename: str,
        dtm_filename: str,
        config_filename: str = None,
        config_dictionary: Dict[str, Any] = {},
    ) -> None:

        required_inputs = {
            "bounding_box": bounding_box,
            "working_projection": working_projection, # this may be removed when fix is added for https://github.com/Loop3D/map2loop/issues/103
            "geology_filename": geology_filename,
            "structure_filename": structure_filename,
            "dtm_filename": dtm_filename,
        }

        # Check for missing required inputs in project
        missing_inputs = [key for key, value in required_inputs.items() if not value]

        if missing_inputs:
            missing_list = ", ".join(missing_inputs)
            logger.error(
                f"Project construction is missing required inputs: {missing_list}. "
                "Please add them to the Project()."
            )
            raise ValueError(
                f"Project construction is missing required inputs: {missing_list}. "
                "Please add them to the Project()."
            )

        # Either config_filename or config_dictionary must be provided (but not both or neither)
        if not config_filename and not config_dictionary:
            logger.error(
                "A config file is required to run map2loop - use either 'config_filename' or 'config_dictionary' to initialise the project."
            )
            raise ValueError(
                 "A config file is required to run map2loop - use either 'config_filename' or 'config_dictionary' to initialise the project."
            )
        if config_filename and config_dictionary:
            logger.error(
                "Both 'config_filename' and 'config_dictionary' were provided. Please specify only one config."
            )
            raise ValueError(
                "Both 'config_filename' and 'config_dictionary' were provided. Please specify only one config."
            )

            
    
    # Getters and Setters
    @beartype.beartype
    def set_ignore_lithology_codes(self, codes: list):
        """
        Set the lithology unit names to be ignored in the geology shapefile.

        This method sets the lithology codes that should be excluded from the geology shapefile
        and triggers the re-population of the stratigraphic column using the updated data
        from the geological map, ensuring the excluded lithologies are not considered.

        Args:
            codes (list):
                A list of strings representing the lithology unit names to be ignored
                in the geological shapefile.
        """
        self.map_data.set_ignore_lithology_codes(codes)
        # Re-populate the units in the column with the new set of ignored geographical units
        self.stratigraphic_column.populate(self.map_data.get_map_data(Datatype.GEOLOGY))

    @beartype.beartype
    def set_ignore_fault_codes(self, codes: list):
        """
        Set the fault names to be ignored in the fault map.

        This method sets the fault codes to be ignored from the fault map and triggers
        re-parsing of the fault map to exclude the ignored faults during subsequent processing.

        Args:
            codes (list):
                A list of strings representing the fault unit names to be ignored
                in the fault map.
        """
        self.map_data.set_ignore_fault_codes(codes)
        # Re-populate the units in the column with the new set of ignored geographical units
        self.map_data.parse_fault_map()  # remove the ignored faults

    @beartype.beartype
    def set_sorter(self, sorter: Sorter):
        """
        Set the sorter for determining stratigraphic order of units

        Args:
            sorter (Sorter):
                The sorter to use.  Must be of base class Sorter
        """
        logger.info(f"Setting sorter to {sorter.sorter_label}")
        self.sorter = sorter

    def get_sorter(self):
        """
        Get the name of the sorter being used

        Returns:
            str: The name of the sorter used
        """
        return self.sorter.sorter_label

    @beartype.beartype
    def set_fault_orientation(self, fault_orientation: FaultOrientation):
        """
        Set the fault orientation calculator that estimates fault orientation values for all faults

        Args:
            fault_orientation (FaultOrientation):
                The calculator to use. Must be of base class FaultOrientation
        """
        logger.info(f"Setting fault orientation calculator to {fault_orientation.label}")
        self.fault_orientation = fault_orientation

    def get_fault_orientation(self):
        """
        Get the name of the fault orientation calculator being used

        Returns:
            str: The name of the fault orientation calculator used
        """
        return self.fault_orientation.label

    @beartype.beartype
    def set_throw_calculator(self, throw_calculator: ThrowCalculator):
        """
        Set the throw calculator that estimates fault throw values for all faults

        Args:
            throw_calculator (ThrowCalculator):
                The calculator to use. Must be of base class ThrowCalculator
        """
        logger.info(f"Setting throw calculator to {throw_calculator.throw_calculator_label}")
        self.throw_calculator = throw_calculator

    def get_throw_calculator(self):
        """
        Get the name of the throw calculator being used

        Returns:
            str: The name of the throw calculator used
        """
        return self.throw_calculator.throw_calculator_label

    def set_default_samplers(self):
        """
        Initialisation function to set or reset the point samplers
        """
        logger.info("Setting default samplers")
        self.samplers[Datatype.STRUCTURE] = SamplerDecimator(1)
        self.samplers[Datatype.GEOLOGY] = SamplerSpacing(50.0)
        self.samplers[Datatype.FAULT] = SamplerSpacing(50.0)
        self.samplers[Datatype.FOLD] = SamplerSpacing(50.0)
        self.samplers[Datatype.DTM] = SamplerSpacing(50.0)

    @beartype.beartype
    def set_sampler(self, datatype: Datatype, sampler: Sampler):
        """
        Set the point sampler for a specific datatype

        Args:
            datatype (Datatype):
                The datatype to use this sampler on
            sampler (Sampler):
                The sampler to use
        """
        allowed_samplers = {
            Datatype.STRUCTURE: SamplerDecimator,
            Datatype.GEOLOGY: SamplerSpacing,
            Datatype.FAULT: SamplerSpacing,
            Datatype.FOLD: SamplerSpacing,
            Datatype.DTM: SamplerSpacing,
        }

        # Check for wrong sampler
        if datatype in allowed_samplers:
            allowed_sampler_type = allowed_samplers[datatype]
            if not isinstance(sampler, allowed_sampler_type):
                raise ValueError(
                    f"Got wrong argument for this datatype: {type(sampler).__name__}, please use {allowed_sampler_type.__name__} instead"
                )
        ## does the enum print the number or the label?
        logger.info(f"Setting sampler for {datatype} to {sampler.sampler_label}")
        self.samplers[datatype] = sampler

    @beartype.beartype
    def get_sampler(self, datatype: Datatype):
        """
        Get the sampler name being used for a datatype

        Args:
            datatype (Datatype): The datatype of the sampler

        Returns:
            str: The name of the sampler being used on the specified datatype
        """
        return self.samplers[datatype].sampler_label

    @beartype.beartype
    def set_minimum_fault_length(self, length: Union[float, int]):
        """
        Set the cutoff length for faults to ignore

        Args:
            length (float):
                The cutoff length
        """
        logger.info(f"Setting minimum fault length to {length}")
        self.map_data.config.fault_config['minimum_fault_length'] = length
        self.map_data.parse_fault_map()

    @beartype.beartype
    def get_minimum_fault_length(self) -> float:
        """
        Get the cutoff length for faults

        Returns:
            float: The cutoff length
        """
        return float(self.map_data.config.fault_config['minimum_fault_length'])

    @beartype.beartype    
    def set_active_thickness(self, thickness_calculator: ThicknessCalculator):
        """
        Sets the active_thickness attribute based on the provided thickness_calculator.
        Args:
            thickness_calculator (object or str): The thickness calculator object or its label.
                If an object is provided, it should have a 'thickness_calculator_label' attribute.
        Returns:
            None
        Raises:
            ValueError: If the thickness calculator label cannot be determined.
        """
        
        try:
            label = thickness_calculator.thickness_calculator_label
        except AttributeError:
            raise ValueError("The provided thickness calculator object does not have a 'thickness_calculator_label' attribute.")
        self.active_thickness = label
    
    @beartype.beartype
    def get_active_thickness(self) -> str:
        """
        Retrieves the active_thickness attribute.

        Returns:
            str: The label of the active thickness calculator.
        """
        return self.active_thickness
            
    # Processing functions
    def sample_map_data(self):
        """
        Use the samplers to extract points along polylines or unit boundaries
        """
        geology_data = self.map_data.get_map_data(Datatype.GEOLOGY)
        dtm_data = self.map_data.get_map_data(Datatype.DTM)
        
        logger.info(f"Sampling geology map data using {self.samplers[Datatype.GEOLOGY].sampler_label}")
        self.geology_samples = self.samplers[Datatype.GEOLOGY].sample(geology_data)

        logger.info(f"Sampling structure map data using {self.samplers[Datatype.STRUCTURE].sampler_label}")
        self.samplers[Datatype.STRUCTURE].dtm_data = dtm_data
        self.samplers[Datatype.STRUCTURE].geology_data = geology_data
        self.structure_samples = self.samplers[Datatype.STRUCTURE].sample(self.map_data.get_map_data(Datatype.STRUCTURE))

        logger.info(f"Sampling fault map data using {self.samplers[Datatype.FAULT].sampler_label}")
        self.fault_samples = self.samplers[Datatype.FAULT].sample(self.map_data.get_map_data(Datatype.FAULT))

        logger.info(f"Sampling fold map data using {self.samplers[Datatype.FOLD].sampler_label}")
        self.fold_samples = self.samplers[Datatype.FOLD].sample(self.map_data.get_map_data(Datatype.FOLD))

    def extract_geology_contacts(self):
        """
        Use the stratigraphic column, and fault and geology data to extract points along contacts
        """
        # Use stratigraphic column to determine basal contacts
        if self.contact_extractor is None:
            self.contact_extractor = ContactExtractor(
                self.map_data.get_map_data(Datatype.GEOLOGY),
                self.map_data.get_map_data(Datatype.FAULT),
            )
            self.contact_extractor.extract_all_contacts()

        basal_contacts = self.contact_extractor.extract_basal_contacts(self.stratigraphic_column.column)

        # sample the contacts
        self.map_data.sampled_contacts = self.samplers[Datatype.GEOLOGY].sample(basal_contacts)
        dtm_data = self.map_data.get_map_data(Datatype.DTM)
        set_z_values_from_raster_df(dtm_data, self.map_data.sampled_contacts)

    def calculate_stratigraphic_order(self, take_best=False):
        """
        Use unit relationships, unit ages and the sorter to create a stratigraphic column
        """
        if self.contact_extractor is None:
            self.contact_extractor = ContactExtractor(
                self.map_data.get_map_data(Datatype.GEOLOGY),
                self.map_data.get_map_data(Datatype.FAULT),
            )
            self.contact_extractor.extract_all_contacts()
        if take_best:
            sorters = [
                SorterAgeBased(),
                SorterAlpha(
                    contacts=self.contact_extractor.contacts,
                ),
                SorterUseNetworkX(
                    unit_relationships=self.topology.get_unit_unit_relationships(),
                ),
            ]
            logger.info(
                f"Calculating best stratigraphic column from {[sorter.sorter_label for sorter in sorters]}"
            )

            columns = [
                sorter.sort(
                    self.stratigraphic_column.stratigraphicUnits,
                )
                for sorter in sorters
            ]
            basal_contacts = [
                self.contact_extractor.extract_basal_contacts(
                    column, save_contacts=False
                )
                for column in columns
            ]
            basal_lengths = [
                sum(list(contacts[contacts["type"] == "BASAL"]["geometry"].length))
                for contacts in basal_contacts
            ]
            max_length = -1
            column = columns[0]
            best_sorter = sorters[0]
            for i in range(len(sorters)):
                if basal_lengths[i] > max_length:
                    max_length = basal_lengths[i]
                    column = columns[i]
                    best_sorter = sorters[i]
            logger.info(
                f"Best sorter {best_sorter.sorter_label} calculated contact length of {max_length}"
            )
            self.stratigraphic_column.column = column
        else:
            logger.info(f'Calculating stratigraphic column using sorter {self.sorter.sorter_label}')
            # Update sorter with current data based on what it needs
            if hasattr(self.sorter, 'unit_relationships') and self.sorter.unit_relationships is None:
                self.sorter.unit_relationships = self.topology.get_unit_unit_relationships()
            if hasattr(self.sorter, 'contacts') and self.sorter.contacts is None:
                self.sorter.contacts = self.contact_extractor.contacts
            if hasattr(self.sorter, 'geology_data') and self.sorter.geology_data is None:
                self.sorter.geology_data = self.map_data.get_map_data(Datatype.GEOLOGY)
            if hasattr(self.sorter, 'structure_data') and self.sorter.structure_data is None:
                self.sorter.structure_data = self.map_data.get_map_data(Datatype.STRUCTURE)
            if hasattr(self.sorter, 'dtm_data') and self.sorter.dtm_data is None:
                self.sorter.dtm_data = self.map_data.get_map_data(Datatype.DTM)
            if hasattr(self.sorter, 'min_age_column') and self.sorter.min_age_column is None:
                self.sorter.min_age_column = self.stratigraphic_column.get_min_age_column()
            if hasattr(self.sorter, 'max_age_column') and self.sorter.max_age_column is None:
                self.sorter.max_age_column = self.stratigraphic_column.get_max_age_column()
            
            self.stratigraphic_column.column = self.sorter.sort(
                self.stratigraphic_column.stratigraphicUnits,
            )

    @beartype.beartype
    def set_thickness_calculator(
        self, thickness_calculator: Union['ThicknessCalculator', List['ThicknessCalculator']], is_strike: bool = False
    ) -> None:
        """
        Sets the thickness_calculator attribute for the object.

        If a single instance of ThicknessCalculator is passed, it wraps it in a list.
        If a list of ThicknessCalculator instances is passed, it validates that all elements
        are instances of ThicknessCalculator before setting the attribute.

        Args:
            thickness_calculator (ThicknessCalculator or list of ThicknessCalculator):
            An instance or a list of ThicknessCalculator objects.

        Raises:
            TypeError: If the provided thickness_calculator is not an instance of
                    ThicknessCalculator or a list of such instances.
        """
        if isinstance(thickness_calculator, ThicknessCalculator):
            thickness_calculator.is_strike = is_strike
            thickness_calculator = [thickness_calculator]

        # Now check if thickness_calculator is a list of valid instances
        if not isinstance(thickness_calculator, list) or not all(
            isinstance(tc, ThicknessCalculator) for tc in thickness_calculator
        ):
            raise TypeError(
                "All items must be instances of ThicknessCalculator or a single ThicknessCalculator instance."
            )

        # Finally, set the calculators
        for tc in thickness_calculator:
            tc.is_strike = is_strike
        self.thickness_calculator = thickness_calculator

    def get_thickness_calculator(self) -> List[str]:
        """
        Retrieves the thickness_calculator_label from the thickness_calculator attribute.

        This method checks if the thickness_calculator attribute is a list or a single object:
        - If it's a list of ThicknessCalculator objects, it returns a list of their labels.
        - If it's a single ThicknessCalculator object, it returns the label as a single-item list.
        - If neither, it raises a TypeError.

        Returns:
            list: A list of thickness_calculator_label(s) from the ThicknessCalculator object(s).

        Raises:
            TypeError: If thickness_calculator is neither a list of objects nor a single object
                    with a 'thickness_calculator_label' attribute.
        """

        if isinstance(self.thickness_calculator, list):
            # If it's a list, return labels from all items
            return [
                calculator.thickness_calculator_label for calculator in self.thickness_calculator
            ]
        elif hasattr(self.thickness_calculator, 'thickness_calculator_label'):
            # If it's a single object, return the label as a list
            return [self.thickness_calculator.thickness_calculator_label]
        else:
            raise TypeError(
                "self.thickness_calculator must be either a list of objects or a single object with a thickness_calculator_label attribute"
            )

    def calculate_unit_thicknesses(self, basal_contacts):
        """
        Calculates the unit thickness statistics (mean, median, standard deviation) for each stratigraphic unit
        in the stratigraphic column using the provided thickness calculators.

        For each calculator in the `thickness_calculator` list:
        - Computes the thickness statistics using the `compute()` method of each calculator.
        - Repeats the computed results to match the number of rows in the stratigraphic units.
        - Appends these results as new columns to the `stratigraphicUnits` dataframe.

        The new columns added for each calculator will be named in the format:
        - {calculator_label}_mean
        - {calculator_label}_median
        - {calculator_label}_stddev

        Additionally, stores the labels of the calculators in the `thickness_calculator_labels` attribute.

        Returns:
            None

        Raises:
            None
        """

        labels = []

        for calculator in self.thickness_calculator:
            calculator.dtm_data = self.map_data.get_map_data(Datatype.DTM)
            calculator.bounding_box = self.bounding_box
            result = calculator.compute(
                self.stratigraphic_column.stratigraphicUnits,
                self.stratigraphic_column.column,
                basal_contacts,
                self.structure_samples,
                self.map_data.get_map_data(Datatype.GEOLOGY),
                self.map_data.sampled_contacts,
            )[['ThicknessMean', 'ThicknessMedian', 'ThicknessStdDev']].to_numpy()

            label = calculator.thickness_calculator_label
            labels.append(label)

            # Repeat the results for the number of rows in stratigraphicUnits
            num_rows = self.stratigraphic_column.stratigraphicUnits.shape[0]
            repeated_result = numpy.tile(result, (num_rows // result.shape[0], 1))

            # Append the repeated results to the lists
            mean_col_name = f"{label}_mean"
            median_col_name = f"{label}_median"
            stddev_col_name = f"{label}_stddev"

            # Attach the results as new columns to the stratigraphic column dataframe
            self.stratigraphic_column.stratigraphicUnits[mean_col_name] = repeated_result[:, 0]
            self.stratigraphic_column.stratigraphicUnits[median_col_name] = repeated_result[:, 1]
            self.stratigraphic_column.stratigraphicUnits[stddev_col_name] = repeated_result[:, 2]

        self.thickness_calculator_labels = labels
        if self.active_thickness is None:
            self.active_thickness = labels[0]
            
    def calculate_fault_orientations(self):
        if self.map_data.get_map_data(Datatype.FAULT_ORIENTATION) is not None:
            logger.info(f"Calculating fault orientations using {self.fault_orientation.label}")
            self.fault_orientations = self.fault_orientation.calculate(
                self.map_data.get_map_data(Datatype.FAULT),
                self.map_data.get_map_data(Datatype.FAULT_ORIENTATION),
                self.map_data,
            )
            dtm_data = self.map_data.get_map_data(Datatype.DTM)
            set_z_values_from_raster_df(dtm_data, self.fault_orientations)
        else:
            logger.warning(
                "No fault orientation data found, skipping fault orientation calculation"
            )

    def apply_colour_to_units(self):
        """
        Apply the clut file to the units in the stratigraphic column
        """
        self.stratigraphic_column.stratigraphicUnits = self.map_data.colour_units(
            self.stratigraphic_column.stratigraphicUnits
        )

    def sort_stratigraphic_column(self):
        """
        Sort the units in the stratigraphic column data structure to match the column order
        """
        logger.info('Sorting stratigraphic column')
        self.stratigraphic_column.sort_from_relationship_list(self.stratigraphic_column.column)

    def summarise_fault_data(self):
        """
        Use the fault shapefile to make a summary of each fault by name
        """
        dtm_data = self.map_data.get_map_data(Datatype.DTM)
        set_z_values_from_raster_df(dtm_data, self.fault_samples)

        self.deformation_history.summarise_data(self.fault_samples)
        self.deformation_history.faults = self.throw_calculator.compute(
            self.deformation_history.faults,
            self.stratigraphic_column.column,
            self.contact_extractor.basal_contacts,
            self.map_data,
        )
        logger.info(f'There are {self.deformation_history.faults.shape[0]} faults in the dataset')

    def run_all(self, user_defined_stratigraphic_column=None, take_best=False):
        """
        Runs the full map2loop process

        Args:
            user_defined_stratigraphic_column (None or list, optional):
                A user fed list that overrides the stratigraphic column sorter. Defaults to None.
        """
        logger.info('Running all map2loop processes')
        if user_defined_stratigraphic_column is not None:
            logger.info(f'User defined stratigraphic column: {user_defined_stratigraphic_column}')

        # Calculate contacts before stratigraphic column
        self.contact_extractor = ContactExtractor(
            self.map_data.get_map_data(Datatype.GEOLOGY),
            self.map_data.get_map_data(Datatype.FAULT),
        )
        self.map_data.contacts = self.contact_extractor.extract_all_contacts()

        # Calculate the stratigraphic column
        if issubclass(type(user_defined_stratigraphic_column), list):
            self.stratigraphic_column.column = user_defined_stratigraphic_column
            self.topology.run()  # if we use a user defined stratigraphic column, we still need to calculate the results of map2model
        else:
            if user_defined_stratigraphic_column is not None:
                logger.warning(
                    f"user_defined_stratigraphic_column is not of type list and is {type(user_defined_stratigraphic_column)}. Attempting to calculate column"
                )  # why not try casting to a list?
            self.calculate_stratigraphic_order(take_best)
        self.sort_stratigraphic_column()

        # Calculate basal contacts based on stratigraphic column
        self.extract_geology_contacts()
        self.sample_map_data()
        self.calculate_unit_thicknesses(self.contact_extractor.basal_contacts)
        self.calculate_fault_orientations()
        self.summarise_fault_data()
        self.apply_colour_to_units()
        self.save_into_projectfile()

    def save_into_projectfile(self):
        """
        Creates or updates a loop project file with all the data extracted from the map2loop process
        """
        # Open project file
        logger.info('Saving data into loop project file')
        if not self.loop_filename:
            logger.info('No loop project file specified, creating a new one')
            output_dir = pathlib.Path.cwd()  
            output_dir.mkdir(parents=True, exist_ok=True) 
            filename = "new_project.loop3d"
            self.loop_filename = str(output_dir / filename)

        file_exists = os.path.isfile(self.loop_filename)

        if file_exists:
            if self.overwrite_lpf:
                logger.info('Overwriting existing loop project file')
                try:
                    os.remove(self.loop_filename)
                    file_exists = False
                    logger.info(f"\nExisting file '{self.loop_filename}' was successfully deleted.")
                except Exception as e:
                    logger.error(f"\nFailed to delete existing file '{self.loop_filename}': {e}")
                    raise e
            else:
                logger.error(
                    f"\nThere is an existing '{self.loop_filename}' with the same name as specified in project. map2loop process may fail. Set 'overwrite_loopprojectfile' to True to avoid this"
                )
                return

        # Initialize the LoopProjectFile
        if not file_exists:
            LPF.CreateBasic(self.loop_filename)

        version_mismatch = False
        existing_extents = None

        if file_exists:
            file_version = LPF.Get(self.loop_filename, "version", verbose=False)
            if file_version["errorFlag"] is True:
                logger.error(f"Error: {file_version['errorString']}")
                logger.error(
                    f"       Cannot export loop project file as current file of name {self.loop_filename} is not a loop project file"
                )
                raise Exception(
                    f"Cannot export loop project file as current file of name {self.loop_filename} is not a loop project file"
                )

            else:
                version_mismatch = file_version["value"] != LPF.LoopVersion()
                if version_mismatch:
                    logger.warning(
                        f"Mismatched loop project file versions {LPF.LoopVersion()} and {file_version}, old version will be replaced"
                    )
            resp = LPF.Get(self.loop_filename, "extents")
            if not resp["errorFlag"]:
                existing_extents = resp["value"]

        # Save extents
        if existing_extents is None:
            LPF.Set(
                self.loop_filename,
                "extents",
                geodesic=[-180, -179, 0, 1],
                utm=[
                    1,
                    1,
                    self.map_data.bounding_box["minx"],
                    self.map_data.bounding_box["maxx"],
                    self.map_data.bounding_box["miny"],
                    self.map_data.bounding_box["maxy"],
                ],
                depth=[self.map_data.bounding_box["top"], self.map_data.bounding_box["base"]],
                spacing=[1000, 1000, 500],
                preference="utm",
                epsg=self.map_data.get_working_projection(),
            )
        else:
            # TODO: Check loopfile extents match project extents before continuing
            # if mismatch on extents warn the user and create new file
            LPF.Set(self.loop_filename, "extents", **existing_extents)

        # Save unit information
        stratigraphic_data = numpy.zeros(
            len(self.stratigraphic_column.stratigraphicUnits), LPF.stratigraphicLayerType
        )
        stratigraphic_thicknesses = numpy.zeros(
            len(self.stratigraphic_column.stratigraphicUnits), LPF.stratigraphicThicknessType)
        
        stratigraphic_data["layerId"] = self.stratigraphic_column.stratigraphicUnits["layerId"]
        stratigraphic_data["minAge"] = self.stratigraphic_column.stratigraphicUnits["minAge"]
        stratigraphic_data["maxAge"] = self.stratigraphic_column.stratigraphicUnits["maxAge"]
        stratigraphic_data["name"] = self.stratigraphic_column.stratigraphicUnits["name"]
        stratigraphic_data["group"] = self.stratigraphic_column.stratigraphicUnits["group"]
        stratigraphic_data["enabled"] = 1

        stratigraphic_thicknesses['name']= self.stratigraphic_column.stratigraphicUnits["name"]
        
        # store all of the thickness estimates in a separate table
        for i, label in enumerate(self.thickness_calculator_labels):
            stratigraphic_thicknesses[f'thickness{i+1}_mean'] = self.stratigraphic_column.stratigraphicUnits.get(f'{label}_mean',0)
            stratigraphic_thicknesses[f'thickness{i+1}_median'] = self.stratigraphic_column.stratigraphicUnits.get(f'{label}_median',0)
            stratigraphic_thicknesses[f'thickness{i+1}_stddev'] = self.stratigraphic_column.stratigraphicUnits.get(f'{label}_stddev',0)
        
        # store the active thickness calculator as the default thickness
        stratigraphic_data["ThicknessMean"] = self.stratigraphic_column.stratigraphicUnits.get(f'{self.active_thickness}_mean',0)
        stratigraphic_data["ThicknessMedian"] = self.stratigraphic_column.stratigraphicUnits.get(f'{self.active_thickness}_median',0)
        stratigraphic_data["ThicknessStdDev"] = self.stratigraphic_column.stratigraphicUnits.get(f'{self.active_thickness}_stddev',0) 

        # Assign colours to startigraphic data
        stratigraphic_data["colour1Red"] = [
            int(a[1:3], 16) for a in self.stratigraphic_column.stratigraphicUnits["colour"]
        ]
        stratigraphic_data["colour1Green"] = [
            int(a[3:5], 16) for a in self.stratigraphic_column.stratigraphicUnits["colour"]
        ]
        stratigraphic_data["colour1Blue"] = [
            int(a[5:7], 16) for a in self.stratigraphic_column.stratigraphicUnits["colour"]
        ]
        stratigraphic_data["colour2Red"] = [int(a * 0.95) for a in stratigraphic_data["colour1Red"]]

        stratigraphic_data["colour2Green"] = [
            int(a * 0.95) for a in stratigraphic_data["colour1Green"]
        ]
        stratigraphic_data["colour2Blue"] = [
            int(a * 0.95) for a in stratigraphic_data["colour1Blue"]
        ]
        
        n_thick_calcs = len(self.thickness_calculator_labels)
        # get thickness calculator labels, and fill up with None if empty values up to 5 placeholders
        while len(self.thickness_calculator_labels) < 5:
            self.thickness_calculator_labels.append("None")

        headers = 'name;'+';'.join([f'{l}_mean;{l}_median;{l}_stddev' for l in self.thickness_calculator_labels[:5]])
        headers = headers.split(';') # split into list

        # save into LPF
        LPF.Set(
            self.loop_filename,
            "stratigraphicLog",
            data=stratigraphic_data,
            verbose=False,
        )
        LPF.Set(self.loop_filename,
                "stratigraphicThicknesses",
                data=stratigraphic_thicknesses,
                headers=headers,
                ncols=1+3*n_thick_calcs, # index and mean, median, stddev for each thickness calculator
                verbose=False)

        # Save contacts
        contacts_data = numpy.zeros(len(self.map_data.sampled_contacts), LPF.contactObservationType)
        contacts_data["layerId"] = self.map_data.sampled_contacts["ID"]
        contacts_data["easting"] = self.map_data.sampled_contacts["X"]
        contacts_data["northing"] = self.map_data.sampled_contacts["Y"]
        contacts_data["altitude"] = self.map_data.sampled_contacts["Z"]
        contacts_data["featureId"] = self.map_data.sampled_contacts["featureId"]
        LPF.Set(self.loop_filename, "contacts", data=contacts_data)

        # Save fault trace information
        faults_obs_data = numpy.zeros(
            len(self.fault_samples) + len(self.fault_orientations), LPF.faultObservationType
        )
        faults_obs_data["val"] = numpy.nan
        faults_obs_data["eventId"][0 : len(self.fault_samples)] = self.fault_samples["ID"]
        faults_obs_data["easting"][0 : len(self.fault_samples)] = self.fault_samples["X"]
        faults_obs_data["northing"][0 : len(self.fault_samples)] = self.fault_samples["Y"]
        faults_obs_data["altitude"][0 : len(self.fault_samples)] = self.fault_samples["Z"]
        faults_obs_data["featureId"][0 : len(self.fault_samples)] = self.fault_samples["featureId"]
        faults_obs_data["dipDir"][0 : len(self.fault_samples)] = numpy.nan
        faults_obs_data["dip"][0 : len(self.fault_samples)] = numpy.nan
        faults_obs_data["posOnly"][0 : len(self.fault_samples)] = 1
        faults_obs_data[
            "displacement"
        ] = 100  # self.fault_samples["DISPLACEMENT"] #TODO remove note needed

        faults_obs_data["eventId"][len(self.fault_samples) :] = self.fault_orientations["ID"]
        faults_obs_data["easting"][len(self.fault_samples) :] = self.fault_orientations["X"]
        faults_obs_data["northing"][len(self.fault_samples) :] = self.fault_orientations["Y"]
        faults_obs_data["altitude"][len(self.fault_samples) :] = self.fault_orientations["Z"]
        faults_obs_data["featureId"][len(self.fault_samples) :] = self.fault_orientations[
            "featureId"
        ]
        faults_obs_data["dipDir"][len(self.fault_samples) :] = self.fault_orientations["DIPDIR"]
        faults_obs_data["dip"][len(self.fault_samples) :] = self.fault_orientations["DIP"]
        faults_obs_data["posOnly"][len(self.fault_samples) :] = 0
        LPF.Set(self.loop_filename, "faultObservations", data=faults_obs_data)

        faults = self.deformation_history.get_faults_for_export()
        faults_data = numpy.zeros(len(faults), LPF.faultEventType)
        faults_data["eventId"] = faults["eventId"]
        faults_data["name"] = faults["name"]
        faults_data["minAge"] = faults["minAge"]
        faults_data["maxAge"] = faults["maxAge"]
        faults_data["group"] = faults["group"]
        faults_data["supergroup"] = faults["supergroup"]
        faults_data["avgDisplacement"] = faults["avgDisplacement"]
        faults_data["avgDownthrowDir"] = faults["avgDownthrowDir"]
        faults_data["influenceDistance"] = faults["influenceDistance"]
        faults_data["verticalRadius"] = faults["verticalRadius"]
        faults_data["horizontalRadius"] = faults["horizontalRadius"]
        faults_data["colour"] = faults["colour"]
        faults_data["centreEasting"] = faults["centreX"]
        faults_data["centreNorthing"] = faults["centreY"]
        faults_data["centreAltitude"] = faults["centreZ"]
        faults_data["avgSlipDirEasting"] = faults["avgSlipDirX"]
        faults_data["avgSlipDirNorthing"] = faults["avgSlipDirY"]
        faults_data["avgSlipDirAltitude"] = faults["avgSlipDirZ"]
        faults_data["avgNormalEasting"] = faults["avgNormalX"]
        faults_data["avgNormalNorthing"] = faults["avgNormalY"]
        faults_data["avgNormalAltitude"] = faults["avgNormalZ"]
        LPF.Set(self.loop_filename, "faultLog", data=faults_data)

        # Save structural information
        observations = numpy.zeros(len(self.structure_samples), LPF.stratigraphicObservationType)
        observations["layer"] = "s0"
        observations["layerId"] = self.structure_samples["layerID"]
        observations["easting"] = self.structure_samples["X"]
        observations["northing"] = self.structure_samples["Y"]
        observations["altitude"] = self.structure_samples["Z"]
        observations["dipDir"] = self.structure_samples["DIPDIR"]
        observations["dip"] = self.structure_samples["DIP"]
        observations["dipPolarity"] = self.structure_samples["OVERTURNED"]
        LPF.Set(self.loop_filename, "stratigraphicObservations", data=observations)

        if self.topology.fault_fault_relationships is not None:
            ff_relationships = self.deformation_history.get_fault_relationships_with_ids(
                self.topology.fault_fault_relationships
            )
            relationships = numpy.zeros(len(ff_relationships), LPF.eventRelationshipType)

            relationships["eventId1"] = ff_relationships["eventId1"]
            relationships["eventId2"] = ff_relationships["eventId2"]
            relationships["bidirectional"] = True
            relationships["angle"] = ff_relationships["Angle"]
            relationships["type"] = LPF.EventRelationshipType.FAULT_FAULT_ABUT
            logger.info("Adding fault relationships to projectfile")
            logger.info(f"Fault relationships: {relationships}")
            LPF.Set(self.loop_filename, "eventRelationships", data=relationships)

    @beartype.beartype
    def draw_geology_map(self, points: pandas.DataFrame = None, overlay: str = ""):
        """
        Plots the geology map with optional points or specific data

        Args:
            points (pandas.DataFrame, optional):
                A dataframe to overlay on the geology map (must contains "X" and "Y" columns). Defaults to None.
            overlay (str, optional):
                Layer of points to overlay (options are "contacts", "basal_contacts", "orientations", "faults"). Defaults to "".
        """
        colour_lookup = (
            self.stratigraphic_column.stratigraphicUnits[["name", "colour"]]
            .set_index("name")
            .to_dict()["colour"]
        )
        geol = self.map_data.get_map_data(Datatype.GEOLOGY).copy()
        geol["colour"] = geol.apply(lambda row: colour_lookup[row.UNITNAME], axis=1)
        geol["colour_rgba"] = geol.apply(lambda row: hex_to_rgb(row["colour"]), axis=1)
        if points is None and overlay == "":
            geol.plot(color=geol["colour_rgba"])
            return
        else:
            base = geol.plot(color=geol["colour_rgba"])
        if overlay != "":
            if overlay == "basal_contacts":
                self.contact_extractor.basal_contacts[self.contact_extractor.basal_contacts["type"] == "BASAL"].plot(
                    ax=base
                )

                return
            elif overlay == "contacts":
                points = self.map_data.sampled_contacts
            elif overlay == "orientations":
                points = self.structure_samples
            elif overlay == "faults":
                points = self.fault_samples
            else:
                print(f"Invalid overlay option {overlay}")
                return
        gdf = geopandas.GeoDataFrame(
            points, geometry=geopandas.points_from_xy(points["X"], points["Y"], crs=geol.crs)
        )
        gdf.plot(ax=base, marker="o", color="red", markersize=5)

    @beartype.beartype
    def save_mapdata_to_files(self, save_path: Union[pathlib.Path,str], extension: str = ".shp.zip"):
        """
        Saves the map data frames to csv files

        Args:
            save_path (str, optional):
                The path to save the file to. Defaults to ".".
            extension (str, optional):
                An alternate extension to save the GeoDataFrame in. Defaults to ".csv".
        """
        
        save_path=pathlib.Path(save_path)
        if not save_path.exists():
            os.makedirs(save_path)
        self.map_data.save_all_map_data(save_path, extension)

    @beartype.beartype
    def save_geotiff_raster(
        self, filename: str = "test.tif", projection: str = "", pixel_size: int = 25
    ):
        """
        Saves the geology map to a geotiff

        Args:
            filename (str, optional):
                The filename of the geotiff file to save to. Defaults to "test.tif".
            projection (str, optional):
                A string of the format "EPSG:3857" that is the projection to output. Defaults to the project working projection
            pixel_size (int, optional):
                The size of a pixel in metres for the geotiff. Defaults to 25
        """
        colour_lookup = (
            self.stratigraphic_column.stratigraphicUnits[["name", "colour"]]
            .set_index("name")
            .to_dict()["colour"]
        )
        geol = self.map_data.get_map_data(Datatype.GEOLOGY).copy()
        geol["colour"] = geol.apply(lambda row: colour_lookup[row.UNITNAME], axis=1)
        geol["colour_red"] = geol.apply(lambda row: int(row["colour"][1:3], 16), axis=1)
        geol["colour_green"] = geol.apply(lambda row: int(row["colour"][3:5], 16), axis=1)
        geol["colour_blue"] = geol.apply(lambda row: int(row["colour"][5:7], 16), axis=1)
        source_ds = gdal.OpenEx(geol.to_json())
        source_layer = source_ds.GetLayer()
        x_min, x_max, y_min, y_max = source_layer.GetExtent()

        # Create the destination data source
        x_res = int((x_max - x_min) / pixel_size)
        y_res = int((y_max - y_min) / pixel_size)
        target_ds = gdal.GetDriverByName("GTiff").Create(filename, x_res, y_res, 4, gdal.GDT_Byte)
        target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(0)

        # Rasterize
        gdal.RasterizeLayer(target_ds, [1], source_layer, options=["ATTRIBUTE=colour_red"])
        gdal.RasterizeLayer(target_ds, [2], source_layer, options=["ATTRIBUTE=colour_green"])
        gdal.RasterizeLayer(target_ds, [3], source_layer, options=["ATTRIBUTE=colour_blue"])
        gdal.RasterizeLayer(target_ds, [4], source_layer, burn_values=[255])
        if re.search("^epsg:[0-9]+$", projection.lower()):
            print("Projection is :", projection)
            target_ds.SetProjection(projection)
        else:
            print("CRS is:", geol.crs.to_string())
            target_ds.SetProjection(geol.crs.to_string())
        target_ds = None
        source_layer = None
        source_ds = None