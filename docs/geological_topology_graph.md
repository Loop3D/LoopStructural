# Geological Topology Graph User Guide

## Overview

The Geological Topology Graph is a powerful data structure for managing and analyzing topological relationships between geological features in LoopStructural. It provides explicit representation of how geological objects (units, faults, unconformities, folds) interact with each other through time and space.

## Key Concepts

### Geological Objects

The topology graph supports several types of geological objects:

- **Units**: Stratigraphic units or rock formations
- **Faults**: Structural discontinuities that displace rock units
- **Foliations**: Planar fabric elements in rocks
- **Folds**: Structural features formed by deformation
- **Unconformities**: Surfaces representing gaps in the geological record

### Relationship Types

The graph supports various topological relationships:

#### Fault Relationships
- `CUTS` / `IS_CUT_BY`: Fault cutting relationships
- `ABUTS` / `IS_ABUTTED_BY`: Fault termination relationships

#### Temporal Relationships
- `OLDER_THAN` / `YOUNGER_THAN`: Relative age relationships
- `CONFORMABLE_OVERLIES` / `CONFORMABLE_UNDERLIES`: Conformable stratigraphic contacts

#### Unconformity Relationships
- `ERODE_UNCONFORMABLY_OVERLIES` / `ERODE_UNCONFORMABLY_UNDERLIES`: Erosional unconformities
- `ONLAP_UNCONFORMABLY_OVERLIES` / `ONLAP_UNCONFORMABLY_UNDERLIES`: Onlap unconformities

#### Structural Relationships
- `FOLDS` / `IS_FOLDED_BY`: Folding relationships

## Getting Started

### Basic Usage

```python
from LoopStructural.modelling.core.model_graph import GeologicalTopologyGraph

# Create a new topology graph
graph = GeologicalTopologyGraph()

# Add geological objects
unit1 = graph.add_geological_object('Unit A', 'unit')
unit2 = graph.add_geological_object('Unit B', 'unit')
unit3 = graph.add_geological_object('Unit C', 'unit')
fault1 = graph.add_geological_object('Fault 1', 'fault')

# Define relationships
graph.add_relationship(unit2.id, unit1.id, 'conformable_overlies')
graph.add_relationship(unit3.id, unit2.id, 'erode_unconformably_overlies')
graph.add_relationship(fault1.id, unit1.id, 'cuts')
graph.add_relationship(fault1.id, unit2.id, 'cuts')
```

### Using Custom Object IDs

You can specify custom IDs for better control:

```python
unit = graph.add_geological_object(
    name='Sandstone Formation',
    object_type='unit',
    object_id='sandstone_01',
    attributes={'thickness': 100.0, 'color': 'yellow'}
)
```

### Adding Attributes

Objects can store additional metadata:

```python
# Add a fault with displacement information
fault = graph.add_geological_object(
    name='Main Fault',
    object_type='fault',
    attributes={
        'displacement': 50.0,
        'strike': 045,
        'dip': 60
    }
)
```

## Working with Relationships

### Querying Relationships

```python
# Get all relationships from an object
rels = graph.get_relationships(from_object_id=unit1.id)

# Get all relationships to an object
rels = graph.get_relationships(to_object_id=unit1.id)

# Filter by relationship type
cuts_rels = graph.get_relationships(
    relationship_type='cuts'
)

# Get specific relationship between two objects
rel = graph.get_relationships(
    from_object_id=fault1.id,
    to_object_id=unit1.id,
    relationship_type='cuts'
)
```

### Reciprocal Relationships

Some relationships have reciprocals that can be queried:

```python
# Get relationships including reciprocals
rels = graph.get_relationships(
    from_object_id=unit1.id,
    include_reciprocal=True
)
```

## Advanced Features

### Dependency Analysis

The graph can determine evaluation order based on geological relationships:

```python
# Get topological sort order
result = graph.topological_sort_by_dependencies()

evaluation_order = result['order']  # List of object IDs
cycles = result['cycles']  # Any detected cycles

# Get dependencies for a specific object
deps = graph.get_dependencies('unit_1')  # Objects that must be evaluated first
dependents = graph.get_dependents('unit_1')  # Objects that depend on this one
```

### Region Masking

Build boolean masks to determine which points belong to each geological unit:

```python
import numpy as np

# Define evaluation points
xyz = np.array([
    [0, 0, 0],
    [10, 10, 10],
    # ... more points
])

# Build region masks based on topology
masks = graph.build_region_masks(xyz, model=geological_model)

# Access mask for specific unit
unit_mask = masks['Unit A']  # Boolean array
```

### Validation

Check for topological inconsistencies:

```python
warnings = graph.validate_topology()

for warning in warnings:
    print(f"Warning: {warning}")
```

## Stratigraphic Column View

The `StratigraphicColumnView` provides specialized tools for working with stratigraphic sequences:

```python
from LoopStructural.modelling.core.model_graph import StratigraphicColumnView

# Create a view
column = StratigraphicColumnView(graph)

# Get all units
units = column.get_units()

# Get units in stratigraphic order (oldest to youngest)
ordered_units = column.get_ordered_units()

# Get unconformity relationships
unconformities = column.get_unconformities()
```

### Identifying Conformable Groups

Group units based on conformable relationships, with groups separated by unconformities:

```python
# Identify conformable groups
groups = column.identify_conformable_groups()

# Each group is a list of units (oldest to youngest)
for i, group in enumerate(groups, 1):
    print(f"Group {i}:")
    for unit in group:
        print(f"  - {unit.name}")
```

**Important**: This method will raise a `ValueError` if circular dependencies are detected:

```python
try:
    groups = column.identify_conformable_groups()
except ValueError as e:
    print(f"Circular dependency detected: {e}")
```

## Visualization

### Plotting the Graph

```python
# Visualize the topology graph
graph.plot(
    with_labels=True,
    node_color_map={
        'unit': '#8dd3c7',
        'fault': '#fb8072',
        'foliation': '#80b1d3',
        'unconformity': '#fdb462'
    },
    figsize=(12, 8)
)
```

### NetworkX Integration

Export to NetworkX for advanced graph algorithms:

```python
import networkx as nx

# Convert to NetworkX DiGraph
G = graph.to_networkx()

# Use NetworkX algorithms
shortest_path = nx.shortest_path(G, source='unit_0', target='unit_2')
is_acyclic = nx.is_directed_acyclic_graph(G)
```

## Integration with GeologicalModel

The topology graph is integrated into the `GeologicalModel` class:

```python
from LoopStructural import GeologicalModel

# Create a geological model
model = GeologicalModel(origin, maximum)

# The model has a topology graph
model.topology.add_geological_object('Unit A', 'unit')

# Add features to the model (automatically adds to topology)
model.create_and_add_foliation('strati_1')

# Define relationships
model.topology.add_relationship('fault_1', 'strati_1', 'cuts')

# Use topology to determine interpolation order
order = model.topology.topological_sort_by_dependencies()
```

## Complete Example: Building a Complex Stratigraphic Sequence

```python
from LoopStructural.modelling.core.model_graph import (
    GeologicalTopologyGraph,
    StratigraphicColumnView
)

# Initialize graph
graph = GeologicalTopologyGraph()

# Add basement units (Group 1)
basement = graph.add_geological_object('Basement', 'unit')
unit1 = graph.add_geological_object('Lower Sandstone', 'unit')
unit2 = graph.add_geological_object('Middle Shale', 'unit')

graph.add_relationship(unit1.id, basement.id, 'conformable_overlies')
graph.add_relationship(unit2.id, unit1.id, 'conformable_overlies')

# Add first unconformity
unit3 = graph.add_geological_object('Upper Sandstone', 'unit')
graph.add_relationship(unit3.id, unit2.id, 'erode_unconformably_overlies')

# Add post-unconformity sequence (Group 2)
unit4 = graph.add_geological_object('Limestone', 'unit')
unit5 = graph.add_geological_object('Coal', 'unit')

graph.add_relationship(unit4.id, unit3.id, 'conformable_overlies')
graph.add_relationship(unit5.id, unit4.id, 'conformable_overlies')

# Add faults
main_fault = graph.add_geological_object(
    'Main Fault',
    'fault',
    attributes={'displacement': 100.0}
)
minor_fault = graph.add_geological_object(
    'Minor Fault',
    'fault',
    attributes={'displacement': 20.0}
)

# Main fault cuts all units
for unit in [basement, unit1, unit2, unit3, unit4, unit5]:
    graph.add_relationship(main_fault.id, unit.id, 'cuts')

# Minor fault only affects younger units
for unit in [unit3, unit4, unit5]:
    graph.add_relationship(minor_fault.id, unit.id, 'cuts')

# Minor fault abuts main fault
graph.add_relationship(minor_fault.id, main_fault.id, 'abuts')

# Analyze the stratigraphy
column = StratigraphicColumnView(graph)

print("Conformable Groups:")
groups = column.identify_conformable_groups()
for i, group in enumerate(groups, 1):
    print(f"\nGroup {i}:")
    for unit in group:
        print(f"  - {unit.name}")

print("\nUnconformities:")
for upper, lower in column.get_unconformities():
    upper_obj = graph.get_object(upper)
    lower_obj = graph.get_object(lower)
    print(f"  {upper_obj.name} unconformably overlies {lower_obj.name}")

print("\nFault Network:")
for fault in graph.get_objects_by_type('fault'):
    cuts_rels = graph.get_relationships(from_object_id=fault.id, relationship_type='cuts')
    cut_units = [graph.get_object(r.to_object).name for r in cuts_rels]
    print(f"  {fault.name} cuts: {', '.join(cut_units)}")

# Visualize
graph.plot(figsize=(14, 10))
```

## Best Practices

### 1. Consistent Relationship Directions

Always use consistent relationship directions:
- Use `OVERLIES` for younger-over-older relationships
- Use `CUTS` for fault-to-feature relationships
- The graph will handle reciprocal relationships automatically

### 2. Validate Topology

Always validate your topology after building:

```python
warnings = graph.validate_topology()
if warnings:
    for w in warnings:
        print(f"Warning: {w}")
```

### 3. Check for Cycles

Before using dependency-based operations, check for cycles:

```python
try:
    groups = column.identify_conformable_groups()
except ValueError as e:
    print(f"Fix circular dependencies: {e}")
```

### 4. Use Attributes for Metadata

Store geological metadata in object attributes:

```python
unit = graph.add_geological_object(
    'Sandstone',
    'unit',
    attributes={
        'thickness': 50.0,
        'porosity': 0.25,
        'age_ma': 150.0,
        'lithology': 'quartz sandstone'
    }
)
```

### 5. Document Complex Relationships

For complex models, document your topology decisions:

```python
# Document why certain relationships exist
graph.add_relationship(
    fault1.id,
    fault2.id,
    'abuts',
    # Note: Could also add a 'note' attribute if needed
)
```

## Common Patterns

### Pattern 1: Simple Layer Cake Stratigraphy

```python
units = []
for i, name in enumerate(['Unit A', 'Unit B', 'Unit C']):
    unit = graph.add_geological_object(name, 'unit')
    units.append(unit)
    if i > 0:
        graph.add_relationship(unit.id, units[i-1].id, 'conformable_overlies')
```

### Pattern 2: Fault Network

```python
# Create faults
fault1 = graph.add_geological_object('Fault 1', 'fault')
fault2 = graph.add_geological_object('Fault 2', 'fault')
fault3 = graph.add_geological_object('Fault 3', 'fault')

# Fault 1 is oldest (cuts no other faults)
# Fault 2 is younger (abuts Fault 1)
graph.add_relationship(fault2.id, fault1.id, 'abuts')

# Fault 3 is youngest (abuts both)
graph.add_relationship(fault3.id, fault1.id, 'abuts')
graph.add_relationship(fault3.id, fault2.id, 'abuts')
```

### Pattern 3: Multiple Unconformities

```python
# Group 1
u1 = graph.add_geological_object('U1', 'unit')
u2 = graph.add_geological_object('U2', 'unit')
graph.add_relationship(u2.id, u1.id, 'conformable_overlies')

# Unconformity 1
u3 = graph.add_geological_object('U3', 'unit')
graph.add_relationship(u3.id, u2.id, 'erode_unconformably_overlies')

# Group 2
u4 = graph.add_geological_object('U4', 'unit')
graph.add_relationship(u4.id, u3.id, 'conformable_overlies')

# Unconformity 2
u5 = graph.add_geological_object('U5', 'unit')
graph.add_relationship(u5.id, u4.id, 'onlap_unconformably_overlies')
```

## Export and Serialization

### Export to Dictionary

```python
# Export to dictionary format
data = graph.to_dict()

# Save to JSON
import json
with open('topology.json', 'w') as f:
    json.dump(data, f, indent=2)
```

### Export to NetworkX

```python
# Export for external analysis
G = graph.to_networkx()

# Save as various formats
import networkx as nx
nx.write_gexf(G, 'topology.gexf')
nx.write_graphml(G, 'topology.graphml')
```

## Troubleshooting

### Issue: Circular Dependency Error

**Problem**: `ValueError: Circular dependency detected...`

**Solution**: Check your relationships for logical inconsistencies:
```python
# Find the cycle in the error message
# Example: "unit1 -> unit2 -> unit3 -> unit1"
# Fix by removing or reversing one relationship
```

### Issue: Missing Objects in Groups

**Problem**: Some units don't appear in conformable groups

**Solution**: Check that all units have appropriate relationships:
```python
# Isolated units will form their own groups
units = column.get_units()
for unit in units:
    rels = graph.get_relationships(from_object_id=unit.id)
    if not rels:
        print(f"{unit.name} has no relationships")
```

### Issue: Unexpected Group Separation

**Problem**: Units expected to be conformable are in different groups

**Solution**: Check for unconformity relationships:
```python
unconformities = column.get_unconformities()
print("Unconformities:", unconformities)
```

## API Reference

For detailed API documentation, see:
- `GeologicalTopologyGraph` class documentation
- `StratigraphicColumnView` class documentation
- `GeologicalObject` dataclass documentation
- `RelationshipType` enum documentation

## Examples

See the following example scripts:
- `examples/5_graph/model_graph.py` - Basic topology graph usage
- `examples/5_graph/topology_graph_example.py` - Advanced examples

## Further Reading

- [Geological Relationships and Time](./geological_relationships.md)
- [Model Building Guide](./model_building.md)
- [API Documentation](./API.rst)
