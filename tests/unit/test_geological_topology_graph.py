"""
Unit tests for GeologicalTopologyGraph and StratigraphicColumnView classes.
"""

import pytest
import numpy as np
from LoopStructural.modelling.core.model_graph import (
    GeologicalTopologyGraph,
    StratigraphicColumnView,
    GeologicalObject,
    GeologicalObjectType,
    RelationshipType,
    TopologicalRelationship,
)


class TestGeologicalTopologyGraph:
    """Test suite for GeologicalTopologyGraph class."""

    def test_graph_initialization(self):
        """Test that a graph can be initialized."""
        graph = GeologicalTopologyGraph()
        assert len(graph) == 0
        assert str(graph) == "GeologicalTopologyGraph: 0 objects, 0 relationships"

    def test_add_geological_object(self):
        """Test adding geological objects to the graph."""
        graph = GeologicalTopologyGraph()
        
        # Add a unit
        unit = graph.add_geological_object('unit1', 'unit')
        assert unit.name == 'unit1'
        assert unit.object_type == GeologicalObjectType.UNIT
        assert len(graph) == 1
        assert 'unit_0' in graph
        
        # Add a fault
        fault = graph.add_geological_object('fault1', GeologicalObjectType.FAULT)
        assert fault.object_type == GeologicalObjectType.FAULT
        assert len(graph) == 2

    def test_add_object_with_custom_id(self):
        """Test adding objects with custom IDs."""
        graph = GeologicalTopologyGraph()
        unit = graph.add_geological_object('unit1', 'unit', object_id='custom_id')
        assert unit.id == 'custom_id'
        assert 'custom_id' in graph

    def test_duplicate_id_raises_error(self):
        """Test that duplicate IDs raise an error."""
        graph = GeologicalTopologyGraph()
        graph.add_geological_object('unit1', 'unit', object_id='test_id')
        
        with pytest.raises(ValueError, match="already exists"):
            graph.add_geological_object('unit2', 'unit', object_id='test_id')

    def test_add_relationship(self):
        """Test adding relationships between objects."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        
        rel = graph.add_relationship(u2.id, u1.id, 'conformable_overlies')
        assert rel.from_object == u2.id
        assert rel.to_object == u1.id
        assert rel.relationship_type == RelationshipType.CONFORMABLE_OVERLIES

    def test_add_relationship_with_invalid_objects(self):
        """Test that adding relationships with invalid objects raises error."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        
        with pytest.raises(ValueError, match="not found"):
            graph.add_relationship(u1.id, 'nonexistent', 'conformable_overlies')

    def test_get_object(self):
        """Test retrieving objects by ID."""
        graph = GeologicalTopologyGraph()
        unit = graph.add_geological_object('unit1', 'unit')
        
        retrieved = graph.get_object(unit.id)
        assert retrieved == unit
        assert graph.get_object('nonexistent') is None

    def test_get_objects_by_type(self):
        """Test filtering objects by type."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        f1 = graph.add_geological_object('fault1', 'fault')
        
        units = graph.get_objects_by_type('unit')
        assert len(units) == 2
        assert u1 in units
        assert u2 in units
        
        faults = graph.get_objects_by_type(GeologicalObjectType.FAULT)
        assert len(faults) == 1
        assert f1 in faults

    def test_get_relationships(self):
        """Test retrieving relationships with various filters."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        u3 = graph.add_geological_object('unit3', 'unit')
        
        graph.add_relationship(u2.id, u1.id, 'conformable_overlies')
        graph.add_relationship(u3.id, u2.id, 'erode_unconformably_overlies')
        
        # Get all relationships
        all_rels = graph.get_relationships()
        assert len(all_rels) == 2
        
        # Filter by from_object
        rels_from_u2 = graph.get_relationships(from_object_id=u2.id)
        assert len(rels_from_u2) == 1
        
        # Filter by relationship type
        unconformities = graph.get_relationships(
            relationship_type=RelationshipType.ERODE_UNCONFORMABLY_OVERLIES
        )
        assert len(unconformities) == 1

    def test_remove_object(self):
        """Test removing objects from the graph."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        graph.add_relationship(u2.id, u1.id, 'conformable_overlies')
        
        assert len(graph) == 2
        result = graph.remove_object(u1.id)
        assert result is True
        assert len(graph) == 1
        assert u1.id not in graph
        
        # Check that relationships are removed
        rels = graph.get_relationships()
        assert len(rels) == 0

    def test_validate_topology(self):
        """Test topology validation."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        
        # Add a relationship without reciprocal
        graph.add_relationship(u2.id, u1.id, 'cuts')
        
        warnings = graph.validate_topology()
        # Should detect missing reciprocal relationship
        assert len(warnings) > 0
        assert any('reciprocal' in w.lower() for w in warnings)

    def test_to_dict(self):
        """Test exporting graph to dictionary."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        graph.add_relationship(u2.id, u1.id, 'conformable_overlies')
        
        result = graph.to_dict()
        assert 'objects' in result
        assert 'relationships' in result
        assert len(result['objects']) == 2
        assert len(result['relationships']) == 1


class TestStratigraphicColumnView:
    """Test suite for StratigraphicColumnView class."""

    def test_view_initialization(self):
        """Test that a view can be initialized with a graph."""
        graph = GeologicalTopologyGraph()
        view = StratigraphicColumnView(graph)
        assert view.graph == graph

    def test_get_units(self):
        """Test retrieving all units from the view."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        f1 = graph.add_geological_object('fault1', 'fault')
        
        view = StratigraphicColumnView(graph)
        units = view.get_units()
        
        assert len(units) == 2
        assert u1 in units
        assert u2 in units
        assert f1 not in units

    def test_get_unconformities(self):
        """Test retrieving unconformity relationships."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        u3 = graph.add_geological_object('unit3', 'unit')
        
        graph.add_relationship(u2.id, u1.id, 'conformable_overlies')
        graph.add_relationship(u3.id, u2.id, 'erode_unconformably_overlies')
        
        view = StratigraphicColumnView(graph)
        unconformities = view.get_unconformities()
        
        assert len(unconformities) == 1
        assert unconformities[0] == (u3.id, u2.id)

    def test_identify_conformable_groups_single_group(self):
        """Test identifying conformable groups with a single group."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        u3 = graph.add_geological_object('unit3', 'unit')
        
        # All conformable
        graph.add_relationship(u2.id, u1.id, 'conformable_overlies')
        graph.add_relationship(u3.id, u2.id, 'conformable_overlies')
        
        view = StratigraphicColumnView(graph)
        groups = view.identify_conformable_groups()
        
        assert len(groups) == 1
        assert len(groups[0]) == 3
        # Check order (oldest to youngest)
        assert groups[0][0].name == 'unit1'
        assert groups[0][1].name == 'unit2'
        assert groups[0][2].name == 'unit3'

    def test_identify_conformable_groups_multiple_groups(self):
        """Test identifying conformable groups with unconformities."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        u3 = graph.add_geological_object('unit3', 'unit')
        u4 = graph.add_geological_object('unit4', 'unit')
        
        # Group 1: unit1 and unit2 (conformable)
        graph.add_relationship(u2.id, u1.id, 'conformable_overlies')
        
        # Unconformity between group 1 and group 2
        graph.add_relationship(u3.id, u2.id, 'erode_unconformably_overlies')
        
        # Group 2: unit3 and unit4 (conformable)
        graph.add_relationship(u4.id, u3.id, 'conformable_overlies')
        
        view = StratigraphicColumnView(graph)
        groups = view.identify_conformable_groups()
        
        assert len(groups) == 2
        # Group 1 (oldest)
        assert len(groups[0]) == 2
        assert groups[0][0].name == 'unit1'
        assert groups[0][1].name == 'unit2'
        # Group 2 (youngest)
        assert len(groups[1]) == 2
        assert groups[1][0].name == 'unit3'
        assert groups[1][1].name == 'unit4'

    def test_identify_conformable_groups_onlap_unconformity(self):
        """Test identifying groups with onlap unconformity."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        u3 = graph.add_geological_object('unit3', 'unit')
        
        graph.add_relationship(u2.id, u1.id, 'conformable_overlies')
        graph.add_relationship(u3.id, u2.id, 'onlap_unconformably_overlies')
        
        view = StratigraphicColumnView(graph)
        groups = view.identify_conformable_groups()
        
        assert len(groups) == 2
        assert len(groups[0]) == 2  # unit1, unit2
        assert len(groups[1]) == 1  # unit3

    def test_identify_conformable_groups_isolated_units(self):
        """Test identifying groups with isolated units."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        u3 = graph.add_geological_object('unit3', 'unit')
        
        # No relationships between units
        
        view = StratigraphicColumnView(graph)
        groups = view.identify_conformable_groups()
        
        assert len(groups) == 3
        # Each unit in its own group
        for group in groups:
            assert len(group) == 1

    def test_identify_conformable_groups_with_faults(self):
        """Test that faults don't interfere with conformable grouping."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        f1 = graph.add_geological_object('fault1', 'fault')
        
        graph.add_relationship(u2.id, u1.id, 'conformable_overlies')
        graph.add_relationship(f1.id, u1.id, 'cuts')
        graph.add_relationship(f1.id, u2.id, 'cuts')
        
        view = StratigraphicColumnView(graph)
        groups = view.identify_conformable_groups()
        
        # Should still be one group (faults don't break conformable sequences)
        assert len(groups) == 1
        assert len(groups[0]) == 2

    def test_cycle_detection_direct_cycle(self):
        """Test that direct cycles are detected and raise an error."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        u3 = graph.add_geological_object('unit3', 'unit')
        
        # Create a cycle: u1 -> u2 -> u3 -> u1
        graph.add_relationship(u2.id, u1.id, 'conformable_overlies')
        graph.add_relationship(u3.id, u2.id, 'conformable_overlies')
        graph.add_relationship(u1.id, u3.id, 'conformable_overlies')
        
        view = StratigraphicColumnView(graph)
        
        with pytest.raises(ValueError, match="Circular dependency detected"):
            view.identify_conformable_groups()

    def test_cycle_detection_mixed_relationships(self):
        """Test cycle detection with mixed relationship types."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        u3 = graph.add_geological_object('unit3', 'unit')
        
        # Create a cycle with different relationship types
        graph.add_relationship(u2.id, u1.id, 'younger_than')
        graph.add_relationship(u3.id, u2.id, 'erode_unconformably_overlies')
        graph.add_relationship(u1.id, u3.id, 'conformable_overlies')
        
        view = StratigraphicColumnView(graph)
        
        with pytest.raises(ValueError, match="Circular dependency detected"):
            view.identify_conformable_groups()

    def test_no_cycle_with_valid_graph(self):
        """Test that valid graphs don't trigger false cycle detection."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        u3 = graph.add_geological_object('unit3', 'unit')
        u4 = graph.add_geological_object('unit4', 'unit')
        
        # Valid DAG structure
        graph.add_relationship(u2.id, u1.id, 'conformable_overlies')
        graph.add_relationship(u3.id, u1.id, 'conformable_overlies')
        graph.add_relationship(u4.id, u2.id, 'conformable_overlies')
        graph.add_relationship(u4.id, u3.id, 'conformable_overlies')
        
        view = StratigraphicColumnView(graph)
        
        # Should not raise an error
        groups = view.identify_conformable_groups()
        assert len(groups) == 1  # All connected conformably

    def test_empty_graph(self):
        """Test handling of empty graph."""
        graph = GeologicalTopologyGraph()
        view = StratigraphicColumnView(graph)
        
        groups = view.identify_conformable_groups()
        assert len(groups) == 0

    def test_single_unit(self):
        """Test handling of graph with single unit."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        
        view = StratigraphicColumnView(graph)
        groups = view.identify_conformable_groups()
        
        assert len(groups) == 1
        assert len(groups[0]) == 1
        assert groups[0][0] == u1

    def test_complex_stratigraphy(self):
        """Test a complex stratigraphic scenario."""
        graph = GeologicalTopologyGraph()
        
        # Create 7 units with complex relationships
        u1 = graph.add_geological_object('unit1', 'unit')  # Oldest
        u2 = graph.add_geological_object('unit2', 'unit')
        u3 = graph.add_geological_object('unit3', 'unit')
        u4 = graph.add_geological_object('unit4', 'unit')
        u5 = graph.add_geological_object('unit5', 'unit')
        u6 = graph.add_geological_object('unit6', 'unit')
        u7 = graph.add_geological_object('unit7', 'unit')  # Youngest
        
        # Group 1: u1, u2, u3 (conformable)
        graph.add_relationship(u2.id, u1.id, 'conformable_overlies')
        graph.add_relationship(u3.id, u2.id, 'conformable_overlies')
        
        # Unconformity 1
        graph.add_relationship(u4.id, u3.id, 'erode_unconformably_overlies')
        
        # Group 2: u4, u5 (conformable)
        graph.add_relationship(u5.id, u4.id, 'conformable_overlies')
        
        # Unconformity 2
        graph.add_relationship(u6.id, u5.id, 'onlap_unconformably_overlies')
        
        # Group 3: u6, u7 (conformable)
        graph.add_relationship(u7.id, u6.id, 'conformable_overlies')
        
        view = StratigraphicColumnView(graph)
        groups = view.identify_conformable_groups()
        
        assert len(groups) == 3
        
        # Check Group 1 (oldest)
        assert len(groups[0]) == 3
        assert [u.name for u in groups[0]] == ['unit1', 'unit2', 'unit3']
        
        # Check Group 2 (middle)
        assert len(groups[1]) == 2
        assert [u.name for u in groups[1]] == ['unit4', 'unit5']
        
        # Check Group 3 (youngest)
        assert len(groups[2]) == 2
        assert [u.name for u in groups[2]] == ['unit6', 'unit7']

    def test_underlies_relationships(self):
        """Test that underlies relationships work correctly."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        u3 = graph.add_geological_object('unit3', 'unit')
        
        # Use underlies instead of overlies
        graph.add_relationship(u1.id, u2.id, 'conformable_underlies')
        graph.add_relationship(u2.id, u3.id, 'conformable_underlies')
        
        view = StratigraphicColumnView(graph)
        groups = view.identify_conformable_groups()
        
        assert len(groups) == 1
        assert len(groups[0]) == 3
        # Check correct order (oldest to youngest)
        assert groups[0][0].name == 'unit1'
        assert groups[0][1].name == 'unit2'
        assert groups[0][2].name == 'unit3'

    def test_older_younger_relationships(self):
        """Test that older/younger relationships work correctly."""
        graph = GeologicalTopologyGraph()
        u1 = graph.add_geological_object('unit1', 'unit')
        u2 = graph.add_geological_object('unit2', 'unit')
        u3 = graph.add_geological_object('unit3', 'unit')
        
        graph.add_relationship(u1.id, u2.id, 'older_than')
        graph.add_relationship(u2.id, u3.id, 'older_than')
        
        view = StratigraphicColumnView(graph)
        groups = view.identify_conformable_groups()
        
        assert len(groups) == 1
        assert len(groups[0]) == 3
        assert groups[0][0].name == 'unit1'
        assert groups[0][1].name == 'unit2'
        assert groups[0][2].name == 'unit3'


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
