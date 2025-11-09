"""
Geological Scenario Builder

High-level declarative API for building geological models from topology specifications.
Provides a fluent interface for defining geological relationships before computational modeling.
"""

from typing import List, Optional, Dict
import pandas as pd
import numpy as np
from ...utils import getLogger
from .model_graph import GeologicalTopologyGraph, GeologicalObjectType, StratigraphicColumnView
from .geological_model import GeologicalModel
from .geological_observations import ObservationCollection

logger = getLogger(__name__)


class GeologicalScenario:
    """
    Declarative API for building geological models from topology specifications.
    
    This class provides a high-level interface for defining geological scenarios
    through topology relationships, which can then be compiled into a computational
    GeologicalModel. It acts as a bridge between pure topology (graph) and 
    the computational model.
    
    Examples
    --------
    >>> scenario = GeologicalScenario(origin, maximum)
    >>> scenario.add_conformable_sequence(['Unit A', 'Unit B', 'Unit C'])
    >>> scenario.add_unconformity('Unit D', 'Unit C', type='erode')
    >>> scenario.add_fault('Fault 1', displacement=100, cuts=['Unit A', 'Unit B'])
    >>> scenario.validate()
    >>> model = scenario.build(data=my_data)
    """
    
    def __init__(self, origin: np.ndarray, maximum: np.ndarray):
        """
        Initialize a geological scenario.
        
        Parameters
        ----------
        origin : np.ndarray
            Origin point of the model domain
        maximum : np.ndarray
            Maximum extent of the model domain
        """
        self.topology = GeologicalTopologyGraph()
        self.model: Optional[GeologicalModel] = None
        self._observations: Dict[str, ObservationCollection] = {}  # Feature-specific observations
        self._fault_parameters: Dict[str, Dict] = {}
        self._feature_parameters: Dict[str, Dict] = {}
        self.origin = np.array(origin)
        self.maximum = np.array(maximum)
        
    def add_unit(self, name: str, thickness: Optional[float] = None, 
                 observations: Optional[ObservationCollection] = None, **kwargs) -> 'GeologicalScenario':
        """
        Add a stratigraphic unit to the scenario.
        
        Parameters
        ----------
        name : str
            Name of the unit
        thickness : float, optional
            Thickness of the unit
        observations : ObservationCollection, optional
            Geological observations for this unit
        **kwargs
            Additional attributes for the unit
            
        Returns
        -------
        self : GeologicalScenario
            For method chaining
        """
        attributes = {'thickness': thickness, **kwargs}
        # Use name as object_id for easy lookup
        self.topology.add_geological_object(name, 'unit', object_id=name, attributes=attributes)
        
        if observations is not None:
            self._observations[name] = observations
            
        logger.info(f"Added unit '{name}' to scenario")
        return self
    
    def add_units(self, names: List[str], **kwargs) -> 'GeologicalScenario':
        """
        Add multiple units at once.
        
        Parameters
        ----------
        names : List[str]
            List of unit names
        **kwargs
            Common attributes for all units
            
        Returns
        -------
        self : GeologicalScenario
            For method chaining
        """
        for name in names:
            self.add_unit(name, **kwargs)
        return self
    
    def add_conformable_sequence(self, unit_names: List[str]) -> 'GeologicalScenario':
        """
        Define a conformable stratigraphic sequence.
        
        Parameters
        ----------
        unit_names : List[str]
            Unit names in order from oldest to youngest
            
        Returns
        -------
        self : GeologicalScenario
            For method chaining
        """
        # Ensure all units exist
        for name in unit_names:
            if name not in self.topology:
                self.add_unit(name)
        
        # Add conformable relationships
        for i in range(len(unit_names) - 1):
            older = unit_names[i]
            younger = unit_names[i + 1]
            self.topology.add_relationship(younger, older, 'conformable_overlies')
            logger.info(f"Added conformable relationship: {younger} overlies {older}")
            
        return self
    
    def add_unconformity(self, younger_unit: str, older_unit: str, 
                        type: str = 'erode') -> 'GeologicalScenario':
        """
        Define an unconformity between units.
        
        Parameters
        ----------
        younger_unit : str
            Name of the overlying (younger) unit
        older_unit : str
            Name of the underlying (older) unit
        type : str, optional
            Type of unconformity: 'erode' or 'onlap'
            
        Returns
        -------
        self : GeologicalScenario
            For method chaining
        """
        if type not in ['erode', 'onlap']:
            raise ValueError(f"Unconformity type must be 'erode' or 'onlap', got '{type}'")
        
        # Ensure units exist
        if younger_unit not in self.topology:
            self.add_unit(younger_unit)
        if older_unit not in self.topology:
            self.add_unit(older_unit)
        
        rel_type = f'{type}_unconformably_overlies'
        self.topology.add_relationship(younger_unit, older_unit, rel_type)
        logger.info(f"Added {type} unconformity: {younger_unit} unconformably overlies {older_unit}")
        
        return self
    
    def add_fault(self, name: str, displacement: float, 
                  cuts: Optional[List[str]] = None, 
                  observations: Optional[ObservationCollection] = None,
                  **kwargs) -> 'GeologicalScenario':
        """
        Add a fault to the scenario.
        
        Parameters
        ----------
        name : str
            Name of the fault
        displacement : float
            Fault displacement magnitude
        cuts : List[str], optional
            List of units that this fault cuts
        observations : ObservationCollection, optional
            Geological observations for this fault
        **kwargs
            Additional fault parameters (e.g., major_axis, minor_axis)
            
        Returns
        -------
        self : GeologicalScenario
            For method chaining
        """
        attributes = {'displacement': displacement, **kwargs}
        # Use name as object_id for easy lookup
        self.topology.add_geological_object(name, 'fault', object_id=name, attributes=attributes)
        self._fault_parameters[name] = kwargs
        
        if observations is not None:
            self._observations[name] = observations
        
        # Add cutting relationships
        if cuts:
            for unit in cuts:
                if unit not in self.topology:
                    logger.warning(f"Unit '{unit}' not in topology, adding it")
                    self.add_unit(unit)
                self.topology.add_relationship(name, unit, 'cuts')
                logger.info(f"Added relationship: {name} cuts {unit}")
        
        logger.info(f"Added fault '{name}' with displacement {displacement}")
        return self
    
    def add_fault_network(self, faults_oldest_to_youngest: List[str]) -> 'GeologicalScenario':
        """
        Define fault network with abutting relationships.
        
        Parameters
        ----------
        faults_oldest_to_youngest : List[str]
            Fault names in order from oldest to youngest
            
        Returns
        -------
        self : GeologicalScenario
            For method chaining
        """
        for i in range(len(faults_oldest_to_youngest)):
            older_fault = faults_oldest_to_youngest[i]
            for j in range(i + 1, len(faults_oldest_to_youngest)):
                younger_fault = faults_oldest_to_youngest[j]
                self.topology.add_relationship(younger_fault, older_fault, 'abuts')
                logger.info(f"Added relationship: {younger_fault} abuts {older_fault}")
        
        return self
    
    def add_fold(self, name: str, folds: Optional[List[str]] = None,
                 **kwargs) -> 'GeologicalScenario':
        """
        Add a fold to the scenario.
        
        Parameters
        ----------
        name : str
            Name of the fold
        folds : List[str], optional
            List of features that this fold affects
        **kwargs
            Additional fold parameters
            
        Returns
        -------
        self : GeologicalScenario
            For method chaining
        """
        attributes = {**kwargs}
        # Use name as object_id for easy lookup
        self.topology.add_geological_object(name, 'fold', object_id=name, attributes=attributes)
        
        if folds:
            for feature in folds:
                if feature not in self.topology:
                    logger.warning(f"Feature '{feature}' not in topology")
                    continue
                self.topology.add_relationship(name, feature, 'folds')
                logger.info(f"Added relationship: {name} folds {feature}")
        
        return self
    
    def set_feature_parameters(self, feature_name: str, **kwargs) -> 'GeologicalScenario':
        """
        Set interpolation and building parameters for a feature.
        
        Parameters
        ----------
        feature_name : str
            Name of the feature
        **kwargs
            Parameters for interpolation (nelements, interpolatortype, etc.)
            
        Returns
        -------
        self : GeologicalScenario
            For method chaining
        """
        self._feature_parameters[feature_name] = kwargs
        return self
    
    def observations(self, feature_name: str) -> ObservationCollection:
        """
        Get or create an observation collection for a feature.
        
        This provides a fluent interface for adding observations:
        
        Examples
        --------
        >>> scenario.observations('unit1')\
        ...     .add_contact([100, 200, 0])\
        ...     .add_orientation([100, 200, 0], strike=45, dip=30)\
        ...     .add_above_point([150, 250, 10])
        
        Parameters
        ----------
        feature_name : str
            Name of the geological feature
            
        Returns
        -------
        ObservationCollection
            Collection for adding observations
        """
        if feature_name not in self._observations:
            self._observations[feature_name] = ObservationCollection(feature_name)
        return self._observations[feature_name]
    
    def add_observations_from_dataframe(self, feature_name: str, 
                                       data: pd.DataFrame) -> 'GeologicalScenario':
        """
        Add observations from a traditional LoopStructural DataFrame.
        
        Parameters
        ----------
        feature_name : str
            Name of the feature
        data : pd.DataFrame
            DataFrame with standard LoopStructural columns
            
        Returns
        -------
        self : GeologicalScenario
            For method chaining
        """
        # Store the dataframe directly - will be used during build
        if feature_name not in self._observations:
            self._observations[feature_name] = ObservationCollection(feature_name)
        
        # We'll merge this with any observations already added
        # by storing the raw dataframe and merging during build
        self._observations[feature_name]._raw_dataframe = data
        
        return self
    
    def get_all_observations_dataframe(self) -> pd.DataFrame:
        """
        Get all observations as a combined LoopStructural DataFrame.
        
        Returns
        -------
        pd.DataFrame
            Combined observations for all features
        """
        if not self._observations:
            return pd.DataFrame()
        
        dfs = []
        for feature_name, obs_collection in self._observations.items():
            df = obs_collection.to_dataframe()
            # Also include any raw dataframe
            if hasattr(obs_collection, '_raw_dataframe'):
                dfs.append(obs_collection._raw_dataframe)
            if not df.empty:
                dfs.append(df)
        
        if dfs:
            return pd.concat(dfs, ignore_index=True)
        return pd.DataFrame()
    
    def validate(self) -> List[str]:
        """
        Validate the scenario before building.
        
        Returns
        -------
        warnings : List[str]
            List of validation warnings
            
        Raises
        ------
        ValueError
            If the scenario contains cycles or other critical errors
        """
        logger.info("Validating geological scenario...")
        
        # Check topology validity
        warnings = self.topology.validate_topology()
        
        # Check for cycles in stratigraphic relationships
        view = StratigraphicColumnView(self.topology)
        try:
            groups = view.identify_conformable_groups()
            logger.info(f"Identified {len(groups)} conformable group(s)")
        except ValueError as e:
            raise ValueError(f"Invalid scenario - circular dependencies detected: {e}")
        
        if warnings:
            logger.warning(f"Validation warnings: {warnings}")
        else:
            logger.info("Scenario validation passed")
        
        return warnings
    
    def build(self, data: Optional[pd.DataFrame] = None, 
              validate: bool = True, **kwargs) -> GeologicalModel:
        """
        Build the computational model from the scenario.
        
        This method converts the topology graph into a GeologicalModel with
        interpolated features, following the topological relationships.
        
        Parameters
        ----------
        data : pd.DataFrame, optional
            Additional observational data for the model. Will be combined
            with observations defined in the scenario.
        validate : bool, optional
            Whether to validate before building (default True)
        **kwargs
            Additional parameters for model building
            
        Returns
        -------
        model : GeologicalModel
            The built geological model
        """
        if validate:
            self.validate()
        
        logger.info("Building geological model from scenario...")
        
        # Create the model
        self.model = GeologicalModel(self.origin, self.maximum)
        
        # Transfer topology graph to model
        self.model.topology = self.topology
        
        # Prepare combined data
        scenario_data = self.get_all_observations_dataframe()
        
        if data is not None:
            # Combine scenario observations with additional data
            if not scenario_data.empty:
                combined_data = pd.concat([scenario_data, data], ignore_index=True)
            else:
                combined_data = data
        else:
            combined_data = scenario_data
        
        # Set data on model
        if not combined_data.empty:
            self.model.data = combined_data
            logger.info(f"Loaded {len(combined_data)} observation points")
        
        # Get build order from topology
        order_result = self.topology.topological_sort_by_dependencies()
        build_order = order_result['order']
        
        if order_result.get('cycles'):
            logger.error(f"Cycles detected in topology: {order_result['cycles']}")
            raise ValueError("Cannot build model with cyclic dependencies")
        
        logger.info(f"Building {len(build_order)} objects in topological order")
        
        # Build features in topological order
        for obj_id in build_order:
            obj = self.topology.get_object(obj_id)
            
            if obj.object_type == GeologicalObjectType.FAULT:
                self._build_fault(obj, **kwargs)
            elif obj.object_type == GeologicalObjectType.UNIT:
                self._build_unit(obj, **kwargs)
            elif obj.object_type == GeologicalObjectType.FOLIATION:
                self._build_foliation(obj, **kwargs)
            # Add other types as needed
        
        logger.info(f"Model built successfully with {len(self.model.features)} features")
        return self.model
    
    def _build_fault(self, fault_obj, **kwargs):
        """Build a fault feature from topology object."""
        fault_name = fault_obj.name
        params = self._fault_parameters.get(fault_name, {})
        feature_params = self._feature_parameters.get(fault_name, {})
        
        # Get displacement from object attributes or parameters
        displacement = fault_obj.attributes.get('displacement') or params.get('displacement', 100.0)
        
        # Get observations if available
        obs_collection = self._observations.get(fault_name)
        data = obs_collection.to_dataframe() if obs_collection else None
        
        # Skip if no data available
        if data is None or data.empty:
            logger.warning(f"No observations for fault '{fault_name}', skipping...")
            return
        
        # Merge parameters
        all_params = {**params, **feature_params, **kwargs}
        
        logger.info(f"Building fault '{fault_name}' with displacement {displacement}")
        
        try:
            self.model.create_and_add_fault(
                fault_name,
                displacement=displacement,
                data=data,
                **all_params
            )
        except Exception as e:
            logger.error(f"Failed to build fault '{fault_name}': {e}")
            raise
    
    def _build_unit(self, unit_obj, **kwargs):
        """Build a foliation feature from a unit topology object."""
        unit_name = unit_obj.name
        feature_params = self._feature_parameters.get(unit_name, {})
        
        # Get observations if available
        obs_collection = self._observations.get(unit_name)
        data = obs_collection.to_dataframe() if obs_collection else None
        
        # Skip if no data available
        if data is None or data.empty:
            logger.warning(f"No observations for unit '{unit_name}', skipping...")
            return
        
        # Determine which faults cut this unit
        faults = self._get_cutting_faults(unit_name)
        
        # Merge parameters
        all_params = {**feature_params, **kwargs}
        
        logger.info(f"Building unit '{unit_name}' (cut by {len(faults)} faults)")
        
        try:
            self.model.create_and_add_foliation(
                unit_name,
                data=data,
                faults=faults if faults else None,
                **all_params
            )
        except Exception as e:
            logger.error(f"Failed to build unit '{unit_name}': {e}")
            raise
    
    def _build_foliation(self, foliation_obj, **kwargs):
        """Build a foliation feature from topology object."""
        # Similar to _build_unit but with foliation-specific logic
        self._build_unit(foliation_obj, **kwargs)
    
    def _get_cutting_faults(self, unit_name: str) -> List:
        """Get faults that cut a specific unit."""
        rels = self.topology.get_relationships(
            to_object_id=unit_name,
            relationship_type='cuts'
        )
        
        faults = []
        for rel in rels:
            if rel.from_object in self.model:
                faults.append(self.model[rel.from_object])
        
        return faults
    
    def get_stratigraphic_groups(self) -> List[List[str]]:
        """
        Get conformable groups from the scenario.
        
        Returns
        -------
        groups : List[List[str]]
            List of groups, each group is a list of unit names
        """
        view = StratigraphicColumnView(self.topology)
        groups = view.identify_conformable_groups()
        return [[unit.name for unit in group] for group in groups]
    
    def plot_scenario(self, **kwargs):
        """
        Visualize the scenario topology graph.
        
        Parameters
        ----------
        **kwargs
            Arguments passed to GeologicalTopologyGraph.plot()
        """
        self.topology.plot(**kwargs)
    
    def to_dict(self) -> Dict:
        """
        Export scenario to dictionary format.
        
        Returns
        -------
        dict
            Dictionary representation of the scenario
        """
        return {
            'origin': self.origin.tolist(),
            'maximum': self.maximum.tolist(),
            'topology': self.topology.to_dict(),
            'fault_parameters': self._fault_parameters,
            'feature_parameters': self._feature_parameters,
        }
    
    @classmethod
    def from_dict(cls, data: Dict) -> 'GeologicalScenario':
        """
        Create scenario from dictionary format.
        
        Parameters
        ----------
        data : dict
            Dictionary representation of a scenario
            
        Returns
        -------
        scenario : GeologicalScenario
            Reconstructed scenario
        """
        scenario = cls(
            origin=np.array(data['origin']),
            maximum=np.array(data['maximum'])
        )
        
        # Reconstruct topology
        topo_data = data['topology']
        for obj_id, obj_info in topo_data['objects'].items():
            scenario.topology.add_geological_object(
                name=obj_info['name'],
                object_type=obj_info['type'],
                object_id=obj_id
            )
        
        for rel in topo_data['relationships']:
            scenario.topology.add_relationship(
                from_object_id=rel['from'],
                to_object_id=rel['to'],
                relationship_type=rel['type']
            )
        
        scenario._fault_parameters = data.get('fault_parameters', {})
        scenario._feature_parameters = data.get('feature_parameters', {})
        
        return scenario
