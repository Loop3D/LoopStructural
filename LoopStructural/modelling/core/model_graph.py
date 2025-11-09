"""
Geological Model Topology Graph

This module implements a graph-based data structure for managing topological 
relationships between geological objects including faults, foliations, and 
unconformities.
"""

from typing import Dict, List, Set, Optional, Tuple, Union, TYPE_CHECKING, Any
from enum import Enum
import numpy as np
from dataclasses import dataclass, field
from collections import defaultdict, deque

from ...utils import getLogger

if TYPE_CHECKING:
    pass

logger = getLogger(__name__)


class GeologicalObjectType(Enum):
    """Enumeration of geological object types supported by the topology graph."""
    FAULT = "fault"
    FOLIATION = "foliation"
    FOLD = "fold"
    UNCONFORMITY = "unconformity"
    UNIT = "unit"

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

    CONFORMABLE_OVERLIES = "conformable_overlies"
    CONFORMABLE_UNDERLIES = "conformable_underlies"
    # Structural relationships
    FOLDS = "folds"                 # Feature folds another
    IS_FOLDED_BY = "is_folded_by"   # Feature is folded by another

RECIPROCAL_MAP = {
            RelationshipType.CUTS: RelationshipType.IS_CUT_BY,
            RelationshipType.IS_CUT_BY: RelationshipType.CUTS,
            RelationshipType.OLDER_THAN: RelationshipType.YOUNGER_THAN,
            RelationshipType.YOUNGER_THAN: RelationshipType.OLDER_THAN,
            RelationshipType.ONLAP_UNCONFORMABLY_OVERLIES: RelationshipType.ONLAP_UNCONFORMABLY_UNDERLIES,
            RelationshipType.FOLDS: RelationshipType.IS_FOLDED_BY,
            RelationshipType.IS_FOLDED_BY: RelationshipType.FOLDS,
            RelationshipType.ABUTS: RelationshipType.IS_ABUTTED_BY,
            RelationshipType.ERODE_UNCONFORMABLY_OVERLIES: RelationshipType.ERODE_UNCONFORMABLY_UNDERLIES,
            RelationshipType.CONFORMABLE_OVERLIES: RelationshipType.CONFORMABLE_UNDERLIES
        }
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
    attributes: Dict[str, Any] = field(default_factory=dict)
    
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
                            attributes: Optional[Dict[str, Any]] = None) -> GeologicalObject:
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
        attributes : dict, optional
            Additional attributes for the geological object
            
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
        
        attributes = attributes or {}
        
        # Choose subclass based on object_type
        if object_type == GeologicalObjectType.UNIT:
            geo_object = UnitNode(id=object_id, name=name, object_type=object_type, attributes=attributes)
        elif object_type == GeologicalObjectType.FAULT:
            geo_object = FaultNode(id=object_id, name=name, object_type=object_type, attributes=attributes)
        elif object_type == GeologicalObjectType.FOLIATION:
            geo_object = FoliationNode(id=object_id, name=name, object_type=object_type, attributes=attributes)
        else:
            geo_object = GeologicalObject(id=object_id, name=name, object_type=object_type, attributes=attributes)
        
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
        # Do NOT add reciprocal relationships automatically anymore
        logger.info(f"Added relationship: {from_object_id} {relationship_type.value} {to_object_id}")
        return relationship
    
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
                         relationship_type: Optional[Union[RelationshipType, str]] = None,
                         include_reciprocal: bool = False) -> List[TopologicalRelationship]:
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
        include_reciprocal : bool, optional
            Whether to include reciprocal relationships (one level only)
            
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
        # Optionally include reciprocals (calculated on the fly)
        if include_reciprocal and relationship_type in RECIPROCAL_MAP:
            reciprocal_type = RECIPROCAL_MAP[relationship_type]
            reciprocal_rels = []
            if from_object_id and to_object_id:
                reciprocal_rels = [r for r in self._relationships.get(to_object_id, []) if r.to_object == from_object_id and r.relationship_type == reciprocal_type]
            elif from_object_id:
                reciprocal_rels = [r for r in self._reverse_relationships.get(from_object_id, []) if r.relationship_type == reciprocal_type]
            elif to_object_id:
                reciprocal_rels = [r for r in self._relationships.get(to_object_id, []) if r.relationship_type == reciprocal_type]
            else:
                for rel_list in self._relationships.values():
                    reciprocal_rels.extend([r for r in rel_list if r.relationship_type == reciprocal_type])
            relationships.extend(reciprocal_rels)
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
    
    def to_networkx(self):
        """
        Convert the geological topology graph to a NetworkX DiGraph for visualization/analysis.
        Returns
        -------
        nx.DiGraph
            Directed graph with node/edge attributes.
        """
        try:
            import networkx as nx
        except ImportError:
            raise ImportError("networkx is required for to_networkx(). Install with pip install networkx.")
        G = nx.DiGraph()
        for obj_id, obj in self._objects.items():
            G.add_node(obj_id, label=obj.name, type=obj.object_type.value, **obj.attributes)
        # Only add one direction for reciprocal relationships
        seen = set()
        for rel_list in self._relationships.values():
            for rel in rel_list:
                key = (rel.from_object, rel.to_object, rel.relationship_type)
                reciprocal_type = RECIPROCAL_MAP.get(rel.relationship_type)
                reciprocal_key = (rel.to_object, rel.from_object, reciprocal_type)
                if reciprocal_type and reciprocal_key in seen:
                    continue  # skip if reciprocal already added
                G.add_edge(rel.from_object, rel.to_object, type=rel.relationship_type.value)
                seen.add(key)
        return G

    def plot(self, with_labels=True, node_color_map=None, figsize=(10, 7)):
        """
        Visualize the topology graph using networkx and matplotlib.
        Parameters
        ----------
        with_labels : bool
            Whether to show node labels.
        node_color_map : dict, optional
            Mapping from object type to color.
        figsize : tuple
            Figure size for matplotlib.
        """
        try:
            import networkx as nx
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("networkx and matplotlib are required for plot(). Install with pip install networkx matplotlib.")
        G = self.to_networkx()
        plt.figure(figsize=figsize)
        pos = nx.spring_layout(G, seed=42)
        # Color nodes by type
        node_types = nx.get_node_attributes(G, 'type')
        if node_color_map is None:
            default_colors = {
                'unit': '#8dd3c7',
                'fault': '#fb8072',
                'foliation': '#80b1d3',
                'fold': '#bebada',
                'unconformity': '#fdb462',
            }
            node_color_map = default_colors
        node_colors = [node_color_map.get(node_types.get(n, ''), '#cccccc') for n in G.nodes]
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=700, alpha=0.9)
        nx.draw_networkx_edges(G, pos, arrows=True, arrowstyle='-|>', arrowsize=20)
        if with_labels:
            labels = nx.get_node_attributes(G, 'label')
            nx.draw_networkx_labels(G, pos, labels=labels, font_size=10)
        # Draw edge labels (relationship types)
        edge_labels = nx.get_edge_attributes(G, 'type')
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='gray', font_size=8)
        plt.axis('off')
        plt.title('Geological Topology Graph')
        plt.tight_layout()
        plt.show()
    
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
        indeg: Dict[str, int] = dict.fromkeys(deps.keys(), 0)
        for _u, succs in deps.items():
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

@dataclass
class UnitNode(GeologicalObject):
    thickness: Optional[float] = None
    data: Any = None
    def __post_init__(self):
        super().__post_init__()
        if self.thickness is None:
            self.thickness = self.attributes.get('thickness')
        if self.data is None:
            self.data = self.attributes.get('data')

@dataclass
class FaultNode(GeologicalObject):
    displacement: Optional[float] = None
    length: Optional[float] = None
    major_axis: Optional[float] = None
    minor_axis: Optional[float] = None

    def __post_init__(self):
        super().__post_init__()
        if self.displacement is None:
            self.displacement = self.attributes.get('displacement')
        if self.length is None:
            self.length = self.attributes.get('length')
        if self.major_axis is None:
            self.major_axis = self.attributes.get('major_axis')
        if self.minor_axis is None:
            self.minor_axis = self.attributes.get('minor_axis')

@dataclass
class FoliationNode(GeologicalObject):
    vector_field: Any = None
    def __post_init__(self):
        super().__post_init__()
        if self.vector_field is None:
            self.vector_field = self.attributes.get('vector_field')

class StratigraphicColumnView:
    """
    A view of the stratigraphic column derived from the GeologicalTopologyGraph.
    Provides methods to access units, groups, and their order as defined by the graph.
    """
    def __init__(self, topology_graph: 'GeologicalTopologyGraph'):
        self.graph = topology_graph

    def get_units(self) -> List[GeologicalObject]:
        """Return all unit objects in the topology graph."""
        return self.graph.get_objects_by_type(GeologicalObjectType.UNIT)

    def get_groups(self) -> List[GeologicalObject]:
        """Return all group objects in the topology graph (if groups are present)."""
        # If groups are not explicitly a GeologicalObjectType, this can be adapted.
        return self.graph.get_objects_by_type('group') if hasattr(GeologicalObjectType, 'GROUP') else []

    def get_units_in_group(self, group_name: str) -> List[GeologicalObject]:
        """Return all units belonging to a given group (by name)."""
        group_obj = next((g for g in self.get_groups() if g.name == group_name), None)
        if not group_obj:
            return []
        units = []
        for unit in self.get_units():
            rels = self.graph.get_relationships(from_object_id=unit.id, to_object_id=group_obj.id, relationship_type='belongs_to_group')
            if rels:
                units.append(unit)
        return units

    def get_ordered_units(self) -> List[GeologicalObject]:
        """Return units in topological (stratigraphic) order (oldest to youngest)."""
        # Use topological sort, filter to units only, and reverse for oldest to youngest
        order = self.graph.topological_sort_by_dependencies().get('order', [])
        units = [self.graph.get_object(obj_id) for obj_id in order if self.graph.get_object(obj_id) and self.graph.get_object(obj_id).object_type == GeologicalObjectType.UNIT]
        return list(reversed(units))  # Oldest first

    def get_unconformities(self) -> List[Tuple[str, str]]:
        """Return unconformity relationships as (upper_unit, lower_unit) tuples."""
        unconformities = []
        for rel in self.graph.get_relationships():
            if rel.relationship_type in [
                RelationshipType.ERODE_UNCONFORMABLY_OVERLIES,
                RelationshipType.ONLAP_UNCONFORMABLY_OVERLIES,
            ]:
                unconformities.append((rel.from_object, rel.to_object))
        return unconformities

    def identify_conformable_groups(self) -> List[List[GeologicalObject]]:
        """
        Identify groups of units by traversing the graph. Units belong to the same
        group if they are connected by conformable relationships (no unconformity).
        Groups are separated by unconformities.
        
        Returns
        -------
        List[List[GeologicalObject]]
            List of groups, where each group is a list of units in stratigraphic
            order (oldest to youngest within that group). Groups themselves are
            ordered from oldest to youngest.
            
        Raises
        ------
        ValueError
            If circular edges (cycles) are detected in the unit relationships.
            
        Examples
        --------
        >>> graph = GeologicalTopologyGraph()
        >>> u1 = graph.add_geological_object('unit1', 'unit')
        >>> u2 = graph.add_geological_object('unit2', 'unit')
        >>> u3 = graph.add_geological_object('unit3', 'unit')
        >>> graph.add_relationship(u2.id, u1.id, 'conformable_overlies')
        >>> graph.add_relationship(u3.id, u2.id, 'erode_unconformably_overlies')
        >>> column = StratigraphicColumnView(graph)
        >>> groups = column.identify_conformable_groups()
        >>> # Returns [[unit1, unit2], [unit3]] - two groups separated by unconformity
        """
        # First check for cycles in the entire graph
        self._check_for_cycles()
        
        # Get all units
        all_units = self.get_units()
        if not all_units:
            return []
        
        # Build adjacency information for conformable relationships only
        conformable_types = {
            RelationshipType.CONFORMABLE_OVERLIES,
            RelationshipType.CONFORMABLE_UNDERLIES,
            RelationshipType.OLDER_THAN,
            RelationshipType.YOUNGER_THAN,
        }
        
        _unconformity_types = {
            RelationshipType.ERODE_UNCONFORMABLY_OVERLIES,
            RelationshipType.ERODE_UNCONFORMABLY_UNDERLIES,
            RelationshipType.ONLAP_UNCONFORMABLY_OVERLIES,
            RelationshipType.ONLAP_UNCONFORMABLY_UNDERLIES,
        }
        
        # Build a conformable-only adjacency map (undirected for grouping purposes)
        conformable_neighbors: Dict[str, Set[str]] = defaultdict(set)
        unit_ids = {unit.id for unit in all_units}
        
        for unit_id in unit_ids:
            # Get all relationships involving this unit
            outgoing = self.graph.get_relationships(from_object_id=unit_id)
            incoming = self.graph.get_relationships(to_object_id=unit_id)
            
            for rel in outgoing:
                if rel.to_object in unit_ids:  # Only consider unit-to-unit relationships
                    if rel.relationship_type in conformable_types:
                        conformable_neighbors[unit_id].add(rel.to_object)
                    # Unconformities break groups - don't add to neighbors
            
            for rel in incoming:
                if rel.from_object in unit_ids:  # Only consider unit-to-unit relationships
                    if rel.relationship_type in conformable_types:
                        conformable_neighbors[unit_id].add(rel.from_object)
                    # Unconformities break groups - don't add to neighbors
        
        # Find connected components using DFS
        visited = set()
        groups = []
        
        for unit in all_units:
            if unit.id not in visited:
                # Start a new group with DFS
                group_ids = set()
                stack = [unit.id]
                
                while stack:
                    current_id = stack.pop()
                    if current_id in visited:
                        continue
                    
                    visited.add(current_id)
                    group_ids.add(current_id)
                    
                    # Add conformable neighbors to explore
                    for neighbor_id in conformable_neighbors[current_id]:
                        if neighbor_id not in visited:
                            stack.append(neighbor_id)
                
                # Convert IDs back to objects
                group = [self.graph.get_object(uid) for uid in group_ids]
                groups.append(group)
        
        # Sort units within each group (oldest to youngest)
        for i, group in enumerate(groups):
            groups[i] = self._sort_units_within_group(group)
        
        # Sort groups themselves (oldest to youngest)
        groups = self._sort_groups(groups)
        
        return groups
    
    def _check_for_cycles(self):
        """
        Check for cycles in the unit relationships using DFS.
        
        Raises
        ------
        ValueError
            If a cycle is detected in the graph.
        """
        all_units = self.get_units()
        unit_ids = {unit.id for unit in all_units}
        
        # Build directed adjacency for all unit relationships
        adj: Dict[str, List[str]] = defaultdict(list)
        
        for unit_id in unit_ids:
            rels = self.graph.get_relationships(from_object_id=unit_id)
            for rel in rels:
                if rel.to_object in unit_ids:
                    # Add directed edge based on relationship type
                    if rel.relationship_type in [
                        RelationshipType.CONFORMABLE_OVERLIES,
                        RelationshipType.ERODE_UNCONFORMABLY_OVERLIES,
                        RelationshipType.ONLAP_UNCONFORMABLY_OVERLIES,
                        RelationshipType.YOUNGER_THAN,
                    ]:
                        # from is younger than to
                        adj[unit_id].append(rel.to_object)
                    elif rel.relationship_type in [
                        RelationshipType.CONFORMABLE_UNDERLIES,
                        RelationshipType.ERODE_UNCONFORMABLY_UNDERLIES,
                        RelationshipType.ONLAP_UNCONFORMABLY_UNDERLIES,
                        RelationshipType.OLDER_THAN,
                    ]:
                        # from is older than to
                        adj[rel.to_object].append(unit_id)
        
        # DFS cycle detection with coloring
        WHITE, GRAY, BLACK = 0, 1, 2
        color = dict.fromkeys(unit_ids, WHITE)
        
        def dfs_visit(node_id: str, path: List[str]):
            color[node_id] = GRAY
            path.append(node_id)
            
            for neighbor_id in adj[node_id]:
                if color[neighbor_id] == GRAY:
                    # Back edge found - cycle detected
                    cycle_start = path.index(neighbor_id)
                    cycle = path[cycle_start:] + [neighbor_id]
                    cycle_names = [self.graph.get_object(uid).name for uid in cycle]
                    raise ValueError(
                        f"Circular dependency detected in stratigraphic relationships: "
                        f"{' -> '.join(cycle_names)}"
                    )
                elif color[neighbor_id] == WHITE:
                    dfs_visit(neighbor_id, path)
            
            path.pop()
            color[node_id] = BLACK
        
        for unit_id in unit_ids:
            if color[unit_id] == WHITE:
                dfs_visit(unit_id, [])
    
    def _sort_units_within_group(self, units: List[GeologicalObject]) -> List[GeologicalObject]:
        """
        Sort units within a conformable group from oldest to youngest.
        
        Parameters
        ----------
        units : List[GeologicalObject]
            Units in the same conformable group
            
        Returns
        -------
        List[GeologicalObject]
            Sorted units (oldest to youngest)
        """
        if len(units) <= 1:
            return units
        
        unit_ids = {u.id for u in units}
        
        # Build adjacency for relationships within this group
        adj: Dict[str, List[str]] = defaultdict(list)
        in_degree = dict.fromkeys(unit_ids, 0)
        
        for unit_id in unit_ids:
            rels = self.graph.get_relationships(from_object_id=unit_id)
            for rel in rels:
                if rel.to_object not in unit_ids:
                    continue
                    
                # Determine direction based on relationship type
                if rel.relationship_type in [
                    RelationshipType.CONFORMABLE_OVERLIES,
                    RelationshipType.YOUNGER_THAN,
                ]:
                    # from is younger, to is older
                    adj[rel.to_object].append(unit_id)
                    in_degree[unit_id] += 1
                elif rel.relationship_type in [
                    RelationshipType.CONFORMABLE_UNDERLIES,
                    RelationshipType.OLDER_THAN,
                ]:
                    # from is older, to is younger
                    adj[unit_id].append(rel.to_object)
                    in_degree[rel.to_object] += 1
        
        # Topological sort (Kahn's algorithm) - oldest to youngest
        queue = deque([uid for uid in unit_ids if in_degree[uid] == 0])
        sorted_ids = []
        
        while queue:
            current = queue.popleft()
            sorted_ids.append(current)
            
            for neighbor in adj[current]:
                in_degree[neighbor] -= 1
                if in_degree[neighbor] == 0:
                    queue.append(neighbor)
        
        # If not all units were sorted, there might be isolated units
        for uid in unit_ids:
            if uid not in sorted_ids:
                sorted_ids.append(uid)
        
        return [self.graph.get_object(uid) for uid in sorted_ids]
    
    def _sort_groups(self, groups: List[List[GeologicalObject]]) -> List[List[GeologicalObject]]:
        """
        Sort groups from oldest to youngest based on unconformity relationships.
        
        Parameters
        ----------
        groups : List[List[GeologicalObject]]
            List of conformable groups
            
        Returns
        -------
        List[List[GeologicalObject]]
            Sorted groups (oldest to youngest)
        """
        if len(groups) <= 1:
            return groups
        
        # Create a mapping from unit ID to group index
        unit_to_group = {}
        for i, group in enumerate(groups):
            for unit in group:
                unit_to_group[unit.id] = i
        
        # Build inter-group adjacency based on unconformity relationships
        group_adj: Dict[int, Set[int]] = defaultdict(set)
        group_in_degree = dict.fromkeys(range(len(groups)), 0)
        
        unconformity_types = {
            RelationshipType.ERODE_UNCONFORMABLY_OVERLIES,
            RelationshipType.ONLAP_UNCONFORMABLY_OVERLIES,
            RelationshipType.ERODE_UNCONFORMABLY_UNDERLIES,
            RelationshipType.ONLAP_UNCONFORMABLY_UNDERLIES,
        }
        
        for rel in self.graph.get_relationships():
            if rel.relationship_type not in unconformity_types:
                continue
            
            from_group = unit_to_group.get(rel.from_object)
            to_group = unit_to_group.get(rel.to_object)
            
            if from_group is None or to_group is None or from_group == to_group:
                continue
            
            # Determine direction
            if rel.relationship_type in [
                RelationshipType.ERODE_UNCONFORMABLY_OVERLIES,
                RelationshipType.ONLAP_UNCONFORMABLY_OVERLIES,
            ]:
                # from is younger, to is older
                if to_group not in group_adj[from_group]:
                    group_adj[to_group].add(from_group)
                    group_in_degree[from_group] += 1
            else:  # UNDERLIES
                # from is older, to is younger
                if from_group not in group_adj[to_group]:
                    group_adj[from_group].add(to_group)
                    group_in_degree[to_group] += 1
        
        # Topological sort for groups (oldest to youngest)
        queue = deque([i for i in range(len(groups)) if group_in_degree[i] == 0])
        sorted_indices = []
        
        while queue:
            current = queue.popleft()
            sorted_indices.append(current)
            
            for neighbor in group_adj[current]:
                group_in_degree[neighbor] -= 1
                if group_in_degree[neighbor] == 0:
                    queue.append(neighbor)
        
        # Handle any remaining groups (isolated)
        for i in range(len(groups)):
            if i not in sorted_indices:
                sorted_indices.append(i)
        
        return [groups[i] for i in sorted_indices]

    def __str__(self):
        units = self.get_ordered_units()
        return f"StratigraphicColumnView: {len(units)} units (oldest to youngest): {[u.name for u in units]}"