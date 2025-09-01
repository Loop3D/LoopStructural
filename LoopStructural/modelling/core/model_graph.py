"""
Geological Model Topology Graph

This module implements a graph-based data structure for managing topological 
relationships between geological objects including faults, foliations, and 
unconformities.
"""

from typing import Dict, List, Set, Optional, Tuple, Union, TYPE_CHECKING
from enum import Enum
import numpy as np
from dataclasses import dataclass, field
from collections import defaultdict, deque
import logging

from ...utils import getLogger

if TYPE_CHECKING:
    from typing import Callable

logger = getLogger(__name__)


class GeologicalObjectType(Enum):
    """Enumeration of geological object types supported by the topology graph."""
    FAULT = "fault"
    FOLIATION = "foliation"
    FOLD = "fold"

class RelationshipType(Enum):
    """Types of topological relationships between geological objects."""
    # Fault relationships
    CUTS = "cuts"                    # Fault cuts another feature
    IS_CUT_BY = "is_cut_by"         # Feature is cut by fault
    ABUTS = "abuts"                  # Fault abuts another fault 
    IS_ABUTTED_BY = "is_abutted_by" # Fault is abutted by another fault
    
    # Temporal relationships
    OLDER_THAN = "older_than"       # Feature is older than another
    YOUNGER_THAN = "younger_than"   # Feature is younger than another
    
    # Unconformity relationships
    ONLAP_UNCONFORMABLY_OVERLIES = "onlap_unconformably_overlies"
    ONLAP_UNCONFORMABLY_UNDERLIES = "onlap_unconformably_underlies"
    
    ERODE_UNCONFORMABLY_OVERLIES = "erode_unconformably_overlies"
    ERODE_UNCONFORMABLY_UNDERLIES = "erode_unconformably_underlies"

    # Structural relationships
    FOLDS = "folds"                 # Feature folds another
    IS_FOLDED_BY = "is_folded_by"   # Feature is folded by another


@dataclass
class GeologicalObject:
    """
    Represents a geological object in the topology graph.
    
    Attributes
    ----------
    id : str
        Unique identifier for the object
    name : str
        Human-readable name
    object_type : GeologicalObjectType
        Type of geological object
    
    """
    id: str
    name: str
    object_type: GeologicalObjectType
    
    def __post_init__(self):
        """Validate object after initialization."""
        if not isinstance(self.object_type, GeologicalObjectType):
            raise ValueError(f"Invalid object type: {self.object_type}")
    
    def __hash__(self):
        return hash(self.id)
    
    def __eq__(self, other):
        if not isinstance(other, GeologicalObject):
            return False
        return self.id == other.id


@dataclass
class TopologicalRelationship:
    """
    Represents a relationship between two geological objects.
    
    Attributes
    ----------
    from_object : str
        ID of the source object
    to_object : str
        ID of the target object
    relationship_type : RelationshipType
        Type of relationship

    """
    from_object: str
    to_object: str
    relationship_type: RelationshipType
    
    


class GeologicalTopologyGraph:
    """
    Graph-based data structure for managing geological object topology.
    
    This class maintains a directed graph where nodes represent geological objects
    (faults, foliations, unconformities) and edges represent topological 
    relationships between them.
    """
    
    def __init__(self):
        """Initialize an empty geological topology graph."""
        self._objects: Dict[str, GeologicalObject] = {}
        self._relationships: Dict[str, List[TopologicalRelationship]] = defaultdict(list)
        self._reverse_relationships: Dict[str, List[TopologicalRelationship]] = defaultdict(list)
        self._object_counter = 0
    
    def add_geological_object(self, 
                            name: str, 
                            object_type: Union[GeologicalObjectType, str],
                            object_id: Optional[str] = None,
                            ) -> GeologicalObject:
        """
        Add a geological object to the graph.
        
        Parameters
        ----------
        name : str
            Human-readable name for the object
        object_type : GeologicalObjectType or str
            Type of geological object
        object_id : str, optional
            Unique identifier. If None, generates automatically
        
            
        Returns
        -------
        GeologicalObject
            The created geological object
            
        Raises
        ------
        ValueError
            If object_id already exists or object_type is invalid
        """
        # Convert string type to enum if necessary
        if isinstance(object_type, str):
            try:
                object_type = GeologicalObjectType(object_type.lower())
            except ValueError:
                raise ValueError(f"Invalid object type: {object_type}")
        
        # Generate ID if not provided
        if object_id is None:
            object_id = f"{object_type.value}_{self._object_counter}"
            self._object_counter += 1
        
        # Check for duplicate IDs
        if object_id in self._objects:
            raise ValueError(f"Object with ID '{object_id}' already exists")
        
        # Create and store the object
        geo_object = GeologicalObject(
            id=object_id,
            name=name,
            object_type=object_type,
            
        )
        
        self._objects[object_id] = geo_object
        logger.info(f"Added {object_type.value} '{name}' with ID '{object_id}'")
        
        return geo_object
    
    def add_relationship(self,
                        from_object_id: str,
                        to_object_id: str,
                        relationship_type: Union[RelationshipType, str], **properties) -> TopologicalRelationship:
        """
        Add a topological relationship between two objects.
        
        Parameters
        ----------
        from_object_id : str
            ID of the source object
        to_object_id : str
            ID of the target object
        relationship_type : RelationshipType or str
            Type of relationship
        **properties
            Additional relationship properties
            
        Returns
        -------
        TopologicalRelationship
            The created relationship
            
        Raises
        ------
        ValueError
            If objects don't exist or relationship type is invalid
        """
        # Validate objects exist
        if from_object_id not in self._objects:
            raise ValueError(f"Source object '{from_object_id}' not found")
        if to_object_id not in self._objects:
            raise ValueError(f"Target object '{to_object_id}' not found")
        
        # Convert string type to enum if necessary
        if isinstance(relationship_type, str):
            try:
                relationship_type = RelationshipType(relationship_type.lower())
            except ValueError:
                raise ValueError(f"Invalid relationship type: {relationship_type}")
        
        # Create relationship
        relationship = TopologicalRelationship(
            from_object=from_object_id,
            to_object=to_object_id,
            relationship_type=relationship_type,
            
        )
        
        # Store in both forward and reverse indices
        self._relationships[from_object_id].append(relationship)
        self._reverse_relationships[to_object_id].append(relationship)
        
        # Automatically add reciprocal relationships where appropriate
        self._add_reciprocal_relationship(relationship)
        
        logger.info(f"Added relationship: {from_object_id} {relationship_type.value} {to_object_id}")
        
        return relationship
    
    def _add_reciprocal_relationship(self, relationship: TopologicalRelationship):
        """Add reciprocal relationships automatically."""
        reciprocal_map = {
            RelationshipType.CUTS: RelationshipType.IS_CUT_BY,
            RelationshipType.IS_CUT_BY: RelationshipType.CUTS,
            RelationshipType.OLDER_THAN: RelationshipType.YOUNGER_THAN,
            RelationshipType.YOUNGER_THAN: RelationshipType.OLDER_THAN,
            RelationshipType.ONLAP_UNCONFORMABLY_OVERLIES: RelationshipType.ONLAP_UNCONFORMABLY_UNDERLIES,
            RelationshipType.FOLDS: RelationshipType.IS_FOLDED_BY,
            RelationshipType.IS_FOLDED_BY: RelationshipType.FOLDS,
            RelationshipType.ABUTS: RelationshipType.IS_ABUTTED_BY,
            RelationshipType.ERODE_UNCONFORMABLY_OVERLIES: RelationshipType.ERODE_UNCONFORMABLY_UNDERLIES,
        }
        
        reciprocal_type = reciprocal_map.get(relationship.relationship_type)
        if reciprocal_type:
            # Check if reciprocal already exists
            existing = self.get_relationships(
                relationship.to_object, 
                relationship.from_object,
                reciprocal_type
            )
            
            if not existing:
                reciprocal = TopologicalRelationship(
                    from_object=relationship.to_object,
                    to_object=relationship.from_object,
                    relationship_type=reciprocal_type,
                )
                self._relationships[reciprocal.from_object].append(reciprocal)
                self._reverse_relationships[reciprocal.to_object].append(reciprocal)
    
    def get_object(self, object_id: str) -> Optional[GeologicalObject]:
        """Get a geological object by ID."""
        return self._objects.get(object_id)
    
    def get_objects_by_type(self, object_type: Union[GeologicalObjectType, str]) -> List[GeologicalObject]:
        """Get all objects of a specific type."""
        if isinstance(object_type, str):
            object_type = GeologicalObjectType(object_type.lower())
        
        return [obj for obj in self._objects.values() if obj.object_type == object_type]
    
    def get_relationships(self,
                         from_object_id: Optional[str] = None,
                         to_object_id: Optional[str] = None,
                         relationship_type: Optional[Union[RelationshipType, str]] = None) -> List[TopologicalRelationship]:
        """
        Get relationships matching the specified criteria.
        
        Parameters
        ----------
        from_object_id : str, optional
            Filter by source object ID
        to_object_id : str, optional
            Filter by target object ID
        relationship_type : RelationshipType or str, optional
            Filter by relationship type
            
        Returns
        -------
        List[TopologicalRelationship]
            Matching relationships
        """
        if isinstance(relationship_type, str):
            relationship_type = RelationshipType(relationship_type.lower())
        
        relationships = []
        
        if from_object_id and to_object_id:
            # Specific relationship between two objects
            candidates = self._relationships.get(from_object_id, [])
            relationships = [r for r in candidates if r.to_object == to_object_id]
        elif from_object_id:
            # All relationships from a specific object
            relationships = self._relationships.get(from_object_id, [])
        elif to_object_id:
            # All relationships to a specific object
            relationships = self._reverse_relationships.get(to_object_id, [])
        else:
            # All relationships
            for rel_list in self._relationships.values():
                relationships.extend(rel_list)
        
        # Filter by relationship type if specified
        if relationship_type:
            relationships = [r for r in relationships if r.relationship_type == relationship_type]
        
        return relationships
    
    def get_fault_network(self) -> Dict[str, List[str]]:
        """
        Get the fault network showing which faults cut which features.
        
        Returns
        -------
        Dict[str, List[str]]
            Dictionary mapping fault IDs to lists of feature IDs they cut
        """
        network = {}
        faults = self.get_objects_by_type(GeologicalObjectType.FAULT)
        
        for fault in faults:
            cut_relationships = self.get_relationships(
                from_object_id=fault.id,
                relationship_type=RelationshipType.CUTS
            )
            network[fault.id] = [rel.to_object for rel in cut_relationships]
        
        return network
    
    def get_unconformity_sequence(self) -> List[Tuple[str, str]]:
        """
        Get unconformity relationships in the model.
        
        Returns
        -------
        List[Tuple[str, str]]
            List of (upper_unit, lower_unit) tuples representing unconformities
        """
        unconformities = []
        relationships = self.get_relationships(
            relationship_type=RelationshipType.UNCONFORMABLY_OVERLIES
        )
        
        for rel in relationships:
            unconformities.append((rel.from_object, rel.to_object))
        
        return unconformities
    
    
    
    def validate_topology(self) -> List[str]:
        """
        Validate the geological topology for inconsistencies.
        
        Returns
        -------
        List[str]
            List of validation warnings/errors
        """
        warnings = []
                
        # Check for faults cutting themselves
        for fault in self.get_objects_by_type(GeologicalObjectType.FAULT):
            self_cuts = self.get_relationships(fault.id, fault.id, RelationshipType.CUTS)
            if self_cuts:
                warnings.append(f"Fault {fault.name} cuts itself")
        
        # Check for reciprocal relationships
        all_relationships = self.get_relationships()
        for rel in all_relationships:
            reciprocal_type_map = {
                RelationshipType.CUTS: RelationshipType.IS_CUT_BY,
                RelationshipType.OLDER_THAN: RelationshipType.YOUNGER_THAN,
                RelationshipType.ERODE_UNCONFORMABLY_OVERLIES: RelationshipType.ERODE_UNCONFORMABLY_UNDERLIES,
                RelationshipType.ERODE_UNCONFORMABLY_UNDERLIES: RelationshipType.ERODE_UNCONFORMABLY_OVERLIES,
                RelationshipType.FOLDS: RelationshipType.IS_FOLDED_BY,
            }
            
            expected_reciprocal = reciprocal_type_map.get(rel.relationship_type)
            if expected_reciprocal:
                reciprocal = self.get_relationships(
                    rel.to_object, rel.from_object, expected_reciprocal
                )
                if not reciprocal:
                    warnings.append(
                        f"Missing reciprocal relationship: {rel.to_object} "
                        f"{expected_reciprocal.value} {rel.from_object}"
                    )
        
        return warnings
    
    def remove_object(self, object_id: str) -> bool:
        """
        Remove a geological object and all its relationships.
        
        Parameters
        ----------
        object_id : str
            ID of the object to remove
            
        Returns
        -------
        bool
            True if object was removed, False if not found
        """
        if object_id not in self._objects:
            return False
        
        # Remove all relationships involving this object
        self._relationships.pop(object_id, None)
        self._reverse_relationships.pop(object_id, None)
        
        # Remove references from other objects' relationship lists
        for rel_list in self._relationships.values():
            rel_list[:] = [r for r in rel_list if r.to_object != object_id]
        
        for rel_list in self._reverse_relationships.values():
            rel_list[:] = [r for r in rel_list if r.from_object != object_id]
        
        # Remove the object itself
        del self._objects[object_id]
        
        logger.info(f"Removed object '{object_id}' and all its relationships")
        return True
    
    def to_dict(self) -> Dict:
        """
        Export the graph to a dictionary representation.
        
        Returns
        -------
        Dict
            Dictionary representation of the graph
        """
        return {
            'objects': {
                obj_id: {
                    'name': obj.name,
                    'type': obj.object_type.value,
                }
                for obj_id, obj in self._objects.items()
            },
            'relationships': [
                {
                    'from': rel.from_object,
                    'to': rel.to_object,
                    'type': rel.relationship_type.value,
                }
                for rel_list in self._relationships.values()
                for rel in rel_list
            ]
        }
    
    def __len__(self) -> int:
        """Return the number of objects in the graph."""
        return len(self._objects)
    
    def __contains__(self, object_id: str) -> bool:
        """Check if an object exists in the graph."""
        return object_id in self._objects
    
    def __str__(self) -> str:
        """String representation of the graph."""
        return (f"GeologicalTopologyGraph: {len(self._objects)} objects, "
                f"{sum(len(rels) for rels in self._relationships.values())} relationships")


    def build_dependency_graph(self) -> Dict[str, Set[str]]:
            """
            Build a dependency graph (adjacency list) where an edge u -> v means
            object u must be evaluated before object v. This is derived from
            topological relationships (cuts, older/younger, unconformities, folds).

            Returns
            -------
            Dict[str, Set[str]]
                Mapping from object_id to set of dependent object_ids (successors)
            """
            deps: Dict[str, Set[str]] = {obj_id: set() for obj_id in self._objects.keys()}

            # Helper to add edge u -> v meaning u must precede v
            def add_edge(u: str, v: str):
                if u not in deps:
                    deps[u] = set()
                deps[u].add(v)
                # ensure v exists in the map
                if v not in deps:
                    deps[v] = set()

            # Iterate relationships and translate to evaluation-order edges
            for rel_list in self._relationships.values():
                for rel in rel_list:
                    src = rel.from_object
                    tgt = rel.to_object
                    rtype = rel.relationship_type

                    # Interpret relationships to evaluation dependencies (younger/affecting -> older/affected)
                    if rtype == RelationshipType.CUTS:
                        # A CUTS B => A is younger (affects B), so A should be evaluated before B
                        add_edge(src, tgt)
                    elif rtype == RelationshipType.IS_CUT_BY:
                        # A IS_CUT_BY B => B cuts A, so B before A
                        add_edge(tgt, src)
                    elif rtype == RelationshipType.OLDER_THAN:
                        # A OLDER_THAN B => B is younger, so B should be evaluated before A
                        add_edge(tgt, src)
                    elif rtype == RelationshipType.YOUNGER_THAN:
                        # A YOUNGER_THAN B => A is younger, so A before B
                        add_edge(src, tgt)
                    elif rtype == RelationshipType.ERODE_UNCONFORMABLY_OVERLIES:
                        # A erodes_unconformably_overlies B => A is younger (overlying), so A before B
                        add_edge(src, tgt)
                    elif rtype == RelationshipType.ERODE_UNCONFORMABLY_UNDERLIES:
                        # A erodes_unconformably_underlies B => A is older, so B (younger) before A
                        add_edge(tgt, src)
                    elif rtype == RelationshipType.ONLAP_UNCONFORMABLY_OVERLIES:
                        # A onlap_unconformably_overlies B => A is younger (overlying), so A before B
                        add_edge(tgt, src)
                    elif rtype == RelationshipType.ONLAP_UNCONFORMABLY_UNDERLIES:
                        # A onlap_unconformably_underlies B => A is older, so B (younger) before A
                        add_edge(src, tgt)
                    elif rtype == RelationshipType.FOLDS:
                        # A folds B => A is younger (folds B), so A before B
                        add_edge(src, tgt)
                    elif rtype == RelationshipType.IS_FOLDED_BY:
                        # A is_folded_by B => B folds A, so B before A
                        add_edge(tgt, src)
                    else:
                        # For unknown types, do not add edges
                        continue

            return deps

    def topological_sort_by_dependencies(self) -> Dict[str, Union[List[str], List[List[str]]]]:
        """
        Return a topological ordering of object IDs based solely on dependency
        relationships (no ages). If cycles are present, return a partial order
        and report the cycles.

        Returns
        -------
        Dict[str, Union[List[str], List[List[str]]]]
            {'order': [obj_ids in evaluation order], 'cycles': [list of cycles]}
        """
        deps = self.build_dependency_graph()

        # Compute in-degree for Kahn's algorithm
        indeg: Dict[str, int] = {n: 0 for n in deps.keys()}
        for u, succs in deps.items():
            for v in succs:
                indeg[v] = indeg.get(v, 0) + 1

        # Start with nodes that have zero in-degree (no prerequisites)
        q = deque([n for n, d in indeg.items() if d == 0])
        order: List[str] = []

        while q:
            n = q.popleft()
            order.append(n)
            for m in list(deps.get(n, [])):
                indeg[m] -= 1
                if indeg[m] == 0:
                    q.append(m)

        # If order does not contain all nodes, there is a cycle
        cycles: List[List[str]] = []
        if len(order) != len(deps):
            # Identify strongly connected components / cycles simplistically
            remaining = set(deps.keys()) - set(order)
            # Try to find simple cycles by DFS
            visited = set()

            def dfs_cycle(node, stack, seen):
                seen.add(node)
                stack.append(node)
                for succ in deps.get(node, []):
                    if succ in stack:
                        # capture cycle
                        idx = stack.index(succ)
                        cycles.append(stack[idx:].copy())
                    elif succ not in seen:
                        dfs_cycle(succ, stack, seen)
                stack.pop()

            for n in remaining:
                if n not in visited:
                    dfs_cycle(n, [], visited)

        return {'order': order, 'cycles': cycles}

    def get_dependencies(self, object_id: str) -> List[str]:
        """
        Return the list of objects that must be evaluated before the given object
        (direct and transitive prerequisites), based on relationships.
        """
        deps = self.build_dependency_graph()
        # Build reverse adjacency for easy traversal: edge u->v means u before v
        rev = {k: set() for k in deps.keys()}
        for u, succs in deps.items():
            for v in succs:
                rev[v].add(u)

        # BFS/DFS from object_id over reverse edges to collect prerequisites
        result = []
        stack = [object_id]
        seen = set()
        while stack:
            n = stack.pop()
            for p in rev.get(n, []):
                if p not in seen:
                    seen.add(p)
                    result.append(p)
                    stack.append(p)
        return result

    def get_dependents(self, object_id: str) -> List[str]:
        """
        Return the list of objects that depend on the given object (direct and
        transitive successors), based on relationships.
        """
        deps = self.build_dependency_graph()
        result = []
        stack = [object_id]
        seen = set()
        while stack:
            n = stack.pop()
            for s in deps.get(n, []):
                if s not in seen:
                    seen.add(s)
                    result.append(s)
                    stack.append(s)
        return result

    def build_region_masks(self, xyz: np.ndarray, *, claim_overlaps: bool = True, model=None) -> Dict[str, np.ndarray]:
        """
        Compute boolean region masks for each topology object given points.

        Masks are computed in the topological evaluation order produced by
        :meth:`topological_sort_by_dependencies`. Objects earlier in the order
        (younger / affecting features) will claim points first when
        `claim_overlaps=True`, so they override older features where regions
        overlap.

        Parameters
        ----------
        xyz : np.ndarray
            Points at which to evaluate region callables (N,3).
        claim_overlaps : bool, optional
            If True, once a point is claimed by an earlier object it will not
            be assigned to later objects.

        Returns
        -------
        Dict[str, np.ndarray]
            Mapping from topology object name to boolean mask (length N).
        """
        N = int(np.asarray(xyz).shape[0])
        masks: Dict[str, np.ndarray] = {}

        if xyz is None or N == 0:
            return masks

        # Determine evaluation order (younger/affecting first)
        topo_order = self.topological_sort_by_dependencies().get('order', [])

        # Keep track of already claimed points
        claimed = np.zeros(N, dtype=bool) if claim_overlaps else None

        for obj_id in topo_order:
            topo_obj = self._objects.get(obj_id)
            if topo_obj is None:
                continue

            mask = np.ones(N, dtype=bool)
            # Attempt to derive mask from unconformity relationships.
            rels: List[TopologicalRelationship] = []
            # relationships where this object is the target
            rels.extend(self.get_relationships(to_object_id=obj_id))
            rels.extend(self.get_relationships(from_object_id=obj_id))
            for rel in rels:
                rtype = rel.relationship_type
                if rtype in (
                    RelationshipType.ERODE_UNCONFORMABLY_OVERLIES,
                    RelationshipType.ERODE_UNCONFORMABLY_UNDERLIES,
                    RelationshipType.ONLAP_UNCONFORMABLY_OVERLIES,
                    RelationshipType.ONLAP_UNCONFORMABLY_UNDERLIES,
                ):
                    # determine the feature id to sample scalar field from
                    if rtype in (
                        RelationshipType.ERODE_UNCONFORMABLY_OVERLIES,
                        RelationshipType.ERODE_UNCONFORMABLY_UNDERLIES,
                    ):
                        # erode -> use younger feature's geometry (value > 0)
                        if rtype == RelationshipType.ERODE_UNCONFORMABLY_OVERLIES:
                            source_id = rel.from_object
                            symbol = '>'
                        else:
                            source_id = rel.to_object
                            symbol = '<'
                    else:
                        # onlap -> use older feature's geometry (value > 0)
                        if rtype == RelationshipType.ONLAP_UNCONFORMABLY_OVERLIES:
                            source_id = rel.from_object
                            symbol = '<'
                        else:
                            source_id = rel.to_object
                            symbol = '>'
                    # If a model is provided, evaluate the scalar field for the source feature.
                    if model is not None:
                        try:
                            
                            vals = model.evaluate_feature_value(source_id, xyz, scale=False,use_regions=False)
                            mask = np.asarray(vals > 0, dtype=bool) if symbol == '>' else np.asarray(vals < 0, dtype=bool)
                            print(f"[DEBUG] Topology '{topo_obj.name}' using scalar field from '{source_id}' with symbol '{symbol}': true_count={np.sum(mask)}, false_count={np.sum(~mask)}")
                        except Exception:
                            logger.exception(
                                f"Failed to evaluate scalar field for topology source feature '{source_id}'"
                            )
                            mask = np.ones(N, dtype=bool)
                    else:
                        # No model available to evaluate scalar field -> empty mask
                        mask = np.ones(N, dtype=bool)

                    # We add all unconformities to a feature                    

            # If still no mask, default to empty mask
            if mask is None:
                mask = np.ones(N, dtype=bool)

            if claim_overlaps and rtype in (
                RelationshipType.ERODE_UNCONFORMABLY_OVERLIES,
                RelationshipType.ERODE_UNCONFORMABLY_UNDERLIES,
                RelationshipType.ONLAP_UNCONFORMABLY_OVERLIES,
                RelationshipType.ONLAP_UNCONFORMABLY_UNDERLIES,
            ):
                print(f"[DEBUG] Topology '{topo_obj.name}' before claiming: true_count={np.sum(mask)}, false_count={np.sum(~mask)}")
                mask = np.logical_and(mask, ~claimed)
                claimed |= mask
                print(f"[DEBUG] Topology '{topo_obj.name}' after claiming: true_count={np.sum(mask)}, false_count={np.sum(~mask)}")

            masks[topo_obj.name] = mask

        return masks