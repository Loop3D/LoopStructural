# Integration Strategy: Topology Graph in GeologicalModel

## Overview

This document outlines the strategy for integrating the `GeologicalTopologyGraph` as the backbone of the `GeologicalModel` class, while providing both low-level automatic tracking and high-level declarative APIs.

## Architecture

### Three-Tier Approach

```
┌─────────────────────────────────────────────┐
│  High-Level: GeologicalScenario            │  ← User-friendly declarative API
│  (Define geology → Build model)            │
└──────────────┬──────────────────────────────┘
               │
┌──────────────▼──────────────────────────────┐
│  Mid-Level: GeologicalModel                │  ← Existing computational model
│  (Features + automatic topology tracking)  │
└──────────────┬──────────────────────────────┘
               │
┌──────────────▼──────────────────────────────┐
│  Low-Level: GeologicalTopologyGraph        │  ← Pure topology & relationships
│  (Graph algorithms & validation)           │
└─────────────────────────────────────────────┘
```

## 1. Low-Level Integration (GeologicalTopologyGraph)

**Status: ✓ Partially Complete**

The topology graph serves as the internal backbone for relationship management.

### What's Already Done

```python
class GeologicalModel:
    def __init__(self, ...):
        self.topology = GeologicalTopologyGraph()  # ✓
    
    def _add_feature(self, feature, ...):
        # Automatically add to topology ✓
        self.topology.add_geological_object(feature.name, obj_type, ...)
```

### What Needs Enhancement

#### A. Automatic Relationship Tracking

When features are added, automatically establish relationships:

```python
def _add_feature(self, feature, index=None):
    # ...existing code...
    
    # Automatically add cutting relationships for faults
    if feature.type == FeatureType.FAULT:
        # Fault cuts all units below it (before first unconformity)
        for f in reversed(self.features):
            if f.type == FeatureType.UNCONFORMITY:
                break
            if f.type == FeatureType.INTERPOLATED:
                self.topology.add_relationship(
                    feature.name, f.name, 'cuts'
                )
```

#### B. Use Topology for Build Order

Replace implicit ordering with topology-driven dependency resolution:

```python
def update(self, verbose=False, progressbar=True):
    # Get evaluation order from topology
    order_result = self.topology.topological_sort_by_dependencies()
    build_order = order_result['order']
    
    if order_result.get('cycles'):
        raise ValueError(f"Cyclic dependencies: {order_result['cycles']}")
    
    # Build in topological order
    for feature_id in build_order:
        feature = self[feature_id]
        if feature:
            feature.builder.up_to_date()
```

#### C. Region Masking Integration

Use topology for automatic region definition:

```python
def evaluate_model(self, xyz, scale=True):
    if scale:
        xyz = self.scale(xyz, inplace=False)
    
    # Get topology-driven region masks
    topology_masks = self.topology.build_region_masks(
        xyz, 
        claim_overlaps=True,  # Younger features override older
        model=self
    )
    
    # Apply masks when assigning stratigraphic IDs
    for group in self.stratigraphic_column.get_groups():
        feature_mask = topology_masks.get(group.name, None)
        # ...use mask in evaluation...
```

#### D. Querying Relationships

Add convenience methods:

```python
def get_faults_cutting(self, feature_name):
    """Get all faults that cut a feature."""
    rels = self.topology.get_relationships(
        to_object_id=feature_name,
        relationship_type='cuts'
    )
    return [self[rel.from_object] for rel in rels]

def get_features_affected_by_fault(self, fault_name):
    """Get all features cut by a fault."""
    rels = self.topology.get_relationships(
        from_object_id=fault_name,
        relationship_type='cuts'
    )
    return [self[rel.to_object] for rel in rels]

def get_unconformities(self):
    """Get all unconformity surfaces in the model."""
    rels = self.topology.get_relationships(
        relationship_type='erode_unconformably_overlies'
    )
    return [(self[rel.from_object], self[rel.to_object]) for rel in rels]
```

## 2. Mid-Level: Enhanced GeologicalModel

**Status: ⚠️ Needs Implementation**

The model should expose topology operations naturally.

### A. Explicit Relationship Definition

Allow users to explicitly define relationships:

```python
# Current (implicit)
model.create_and_add_fault('F1', displacement=100)  # cuts everything

# New (explicit)
model.create_and_add_fault('F1', displacement=100, cuts=['Unit1', 'Unit2'])

# Or post-hoc
model.add_relationship('F1', 'Unit1', 'cuts')
```

### B. Topology-Aware Methods

```python
def add_relationship(self, from_feature, to_feature, relationship_type):
    """Add a topological relationship between features."""
    self.topology.add_relationship(from_feature, to_feature, relationship_type)

def validate_topology(self):
    """Validate model topology for inconsistencies."""
    return self.topology.validate_topology()

def get_build_order(self):
    """Get the order in which features should be interpolated."""
    result = self.topology.topological_sort_by_dependencies()
    return result['order']

def visualize_topology(self, **kwargs):
    """Visualize the model's topology graph."""
    self.topology.plot(**kwargs)
```

### C. Stratigraphic Column Integration

Link stratigraphic column with topology:

```python
@property
def stratigraphic_view(self):
    """Get a stratigraphic column view of the topology."""
    return StratigraphicColumnView(self.topology)

def get_conformable_groups(self):
    """Get conformable stratigraphic groups."""
    return self.stratigraphic_view.identify_conformable_groups()
```

## 3. High-Level: GeologicalScenario API

**Status: ✓ Implemented (see geological_scenario.py)**

A declarative API for defining geology before building the model.

### Usage Example

```python
from LoopStructural import GeologicalScenario

# Define the geological scenario
scenario = GeologicalScenario(origin, maximum)

# Define stratigraphy
scenario.add_conformable_sequence([
    'Basement',
    'Sandstone_Lower',
    'Shale_Middle',
    'Sandstone_Upper'
])

# Add unconformity
scenario.add_unconformity('Limestone', 'Sandstone_Upper', type='erode')
scenario.add_conformable_sequence(['Limestone', 'Coal'])

# Add faults
scenario.add_fault(
    'Main_Fault',
    displacement=100,
    cuts=['Basement', 'Sandstone_Lower', 'Shale_Middle']
)

scenario.add_fault(
    'Minor_Fault',
    displacement=30,
    cuts=['Sandstone_Upper', 'Limestone', 'Coal']
)

# Define fault network
scenario.add_fault_network(['Main_Fault', 'Minor_Fault'])

# Validate scenario
warnings = scenario.validate()

# Build computational model
model = scenario.build(data=observation_data)

# Now use model normally
model.update()
model.evaluate_model(xyz)
```

### Benefits

1. **Separation of Concerns**: Topology definition separate from interpolation
2. **Validation Before Building**: Catch errors early
3. **Reproducibility**: Easy to save/load scenarios
4. **Clarity**: Geological intent is explicit
5. **Flexibility**: Can modify topology before building

## 4. Migration Path

### Phase 1: Foundation (Current)
- ✓ Topology graph class implemented
- ✓ Basic integration in GeologicalModel
- ✓ Stratigraphic column view

### Phase 2: Deep Integration (Recommended Next Steps)

1. **Enhance `_add_feature`** to automatically create relationships
2. **Use topology in `update()`** for build order
3. **Implement region masking** from topology
4. **Add relationship query methods**

Example implementation:

```python
def _add_feature(self, feature, index=None):
    # ...existing code to add feature...
    
    # Add to topology
    obj_type = self._determine_object_type(feature)
    self.topology.add_geological_object(
        name=feature.name,
        object_type=obj_type,
        object_id=feature.name,
        attributes={'displacement': getattr(feature, 'displacement', None)}
    )
    
    # Automatically establish relationships
    self._auto_establish_relationships(feature)

def _auto_establish_relationships(self, feature):
    """Automatically establish topological relationships."""
    if feature.type == FeatureType.FAULT:
        # Fault cuts all units below it until unconformity
        for f in reversed(self.features[:-1]):  # Exclude the fault itself
            if f.type == FeatureType.UNCONFORMITY:
                break
            if f.type in [FeatureType.INTERPOLATED, FeatureType.FOLDED]:
                self.topology.add_relationship(
                    feature.name, f.name, 'cuts'
                )
    
    elif feature.type == FeatureType.UNCONFORMITY:
        # Establish unconformity relationships
        base_feature = feature.feature  # The geological feature it's based on
        if base_feature and base_feature.name in self.topology:
            # Find the last conformable unit
            for f in reversed(self.features[:-1]):
                if f.type == FeatureType.INTERPOLATED:
                    self.topology.add_relationship(
                        feature.name, f.name, 
                        'erode_unconformably_overlies'
                    )
                    break
```

### Phase 3: High-Level API

1. **Implement GeologicalScenario class** (✓ Done)
2. **Add scenario import/export**
3. **Create example workflows**
4. **Write comprehensive documentation**

### Phase 4: Advanced Features

1. **Topology-driven uncertainty propagation**
2. **Scenario comparison tools**
3. **Automatic geological rule validation**
4. **Interactive topology editing**

## 5. Backward Compatibility

### Strategy: Gradual Enhancement

```python
# Old way still works
model = GeologicalModel(origin, maximum)
model.create_and_add_foliation('strati')
model.create_and_add_fault('fault1', 100)
# Topology tracked automatically in background

# New way (explicit)
model.add_relationship('fault1', 'strati', 'cuts')

# Advanced way (scenario)
scenario = GeologicalScenario(origin, maximum)
scenario.add_unit('strati')
scenario.add_fault('fault1', 100, cuts=['strati'])
model = scenario.build()
```

### No Breaking Changes

- All existing code continues to work
- Topology is built automatically
- New features are opt-in
- Default behavior unchanged

## 6. Testing Strategy

### Unit Tests

```python
def test_automatic_relationship_creation():
    model = GeologicalModel(origin, maximum)
    model.create_and_add_foliation('unit1')
    model.create_and_add_fault('fault1', 100)
    
    # Check relationship was created
    rels = model.topology.get_relationships(
        from_object_id='fault1',
        to_object_id='unit1',
        relationship_type='cuts'
    )
    assert len(rels) == 1

def test_topology_driven_build_order():
    model = GeologicalModel(origin, maximum)
    # Add features in random order
    model.create_and_add_fault('fault1', 100)
    model.create_and_add_foliation('unit1')
    model.create_and_add_foliation('unit2')
    
    # Topology should determine correct order
    order = model.get_build_order()
    assert order.index('unit1') < order.index('fault1')
    assert order.index('unit2') < order.index('fault1')
```

### Integration Tests

```python
def test_scenario_to_model_workflow():
    scenario = GeologicalScenario(origin, maximum)
    scenario.add_conformable_sequence(['u1', 'u2', 'u3'])
    scenario.add_fault('f1', 100, cuts=['u1', 'u2'])
    
    model = scenario.build(data=test_data)
    model.update()
    
    # Model should work correctly
    result = model.evaluate_model(test_points)
    assert result.shape[0] == len(test_points)
```

## 7. Documentation Updates

### User Guide

1. **Topology Basics**: Introduction to geological relationships
2. **Using Topology in Models**: How topology affects model building
3. **Scenario-Based Modeling**: High-level workflow
4. **Migration Guide**: Converting existing code
5. **Best Practices**: When to use each approach

### API Documentation

- Document all new methods
- Add examples to docstrings
- Create topology module documentation
- Update existing docstrings to mention topology

### Tutorials

1. **Simple Model with Topology**: Basic workflow
2. **Complex Fault Networks**: Using topology for faults
3. **Unconformity Modeling**: Topology-driven regions
4. **Scenario-Based Workflow**: End-to-end example
5. **Advanced: Custom Relationships**: Extending the system

## 8. Example: Complete Integration

```python
# ============================================
# Low-Level: GeologicalModel (automatic)
# ============================================
from LoopStructural import GeologicalModel

model = GeologicalModel(origin, maximum)
model.data = observation_data

# Features auto-register in topology
model.create_and_add_foliation('unit1')
model.create_and_add_foliation('unit2')
model.create_and_add_fault('fault1', displacement=100)

# Query topology
print(model.get_faults_cutting('unit1'))  # ['fault1']
print(model.validate_topology())  # []

# Topology drives interpolation order
model.update()  # Uses topological sort internally

# ============================================
# Mid-Level: Explicit relationships
# ============================================
model = GeologicalModel(origin, maximum)

# Define features
model.create_and_add_foliation('unit1')
model.create_and_add_foliation('unit2')

# Explicitly define relationships
model.add_relationship('unit2', 'unit1', 'conformable_overlies')
model.create_and_add_fault('fault1', 100, cuts=['unit1'])

# Visualize topology
model.visualize_topology()

# ============================================
# High-Level: Scenario-based
# ============================================
from LoopStructural import GeologicalScenario

scenario = GeologicalScenario(origin, maximum)

# Define geology declaratively
scenario.add_conformable_sequence(['unit1', 'unit2', 'unit3'])
scenario.add_unconformity('unit4', 'unit3', type='erode')
scenario.add_fault('fault1', displacement=100, cuts=['unit1', 'unit2', 'unit3'])
scenario.add_fault('fault2', displacement=50, cuts=['unit2', 'unit3', 'unit4'])
scenario.add_fault_network(['fault1', 'fault2'])

# Validate before building
scenario.validate()

# Build model
model = scenario.build(data=observation_data)
model.update()
```

## Summary

### Recommended Implementation Priority

1. **High Priority** (Do First):
   - Enhance `_add_feature` for automatic relationship tracking
   - Implement topology-driven `update()` order
   - Add convenience query methods
   - Complete `GeologicalScenario` class

2. **Medium Priority** (Do Next):
   - Region masking from topology
   - Export/import for scenarios
   - Comprehensive testing
   - Documentation updates

3. **Low Priority** (Future):
   - Advanced validation rules
   - Interactive topology editing
   - Scenario comparison tools
   - Uncertainty propagation

### Key Advantages

✓ **Backward Compatible**: Existing code works unchanged  
✓ **Progressive Enhancement**: Can adopt features gradually  
✓ **Multiple APIs**: Choose complexity level appropriate for task  
✓ **Explicit Relationships**: No more implicit assumptions  
✓ **Validation**: Catch errors before expensive interpolation  
✓ **Reproducibility**: Scenarios are easy to save and share  
✓ **Flexibility**: Can work at any level of abstraction  

The three-tier approach gives users flexibility while making the topology graph the single source of truth for all geological relationships.
