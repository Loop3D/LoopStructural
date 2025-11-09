# Geological Observation System - Summary

## Overview

The geological observation system provides an intuitive, geologically-meaningful way to define observational data for geological features in LoopStructural. Instead of working directly with DataFrames, you can use type-safe, feature-specific observation methods.

## Key Features

1. **Geologically Intuitive API**: Methods named after geological concepts (contacts, orientations, traces, etc.)
2. **Type-Safe**: Each observation type has specific required and optional parameters
3. **Feature-Specific**: Observations are attached to specific geological objects (units, faults, folds)
4. **Fluent Interface**: Method chaining for concise code
5. **Automatic Conversion**: Observations are automatically converted to LoopStructural's internal format

## Observation Types

### For Stratigraphic Units

```python
scenario.observations('unit_name')\
    .add_contact([x, y, z])                          # Point on contact surface
    .add_orientation([x, y, z], strike=45, dip=30)   # Strike/dip measurement
    .add_orientation([x, y, z], gradient=[gx, gy, gz]) # Gradient vector
    .add_inside_point([x, y, z])                     # Known to be inside unit
    .add_outside_point([x, y, z])                    # Known to be outside unit
    .add_above_point([x, y, z])                      # Above the contact
    .add_below_point([x, y, z])                      # Below the contact
    .set_thickness(25.0)                             # Unit thickness
```

### For Faults

```python
scenario.observations('fault_name')\
    .add_fault_trace([x, y, z])                      # Point on fault trace
    .add_fault_orientation([x, y, z], strike=90, dip=60)  # Fault plane orientation
    .add_hangingwall_point([x, y, z])                # In hanging wall
    .add_footwall_point([x, y, z])                   # In footwall
    .add_slip_vector([x, y, z], slip_vector=[0, 0, -1])  # Slip direction
    .add_displacement([x, y, z], displacement=100)   # Displacement measurement
```

## Complete Example

```python
import numpy as np
from LoopStructural.modelling.core.geological_scenario import GeologicalScenario

# Create scenario
origin = np.array([0, 0, 0])
maximum = np.array([1000, 1000, 500])
scenario = GeologicalScenario(origin, maximum)

# Define stratigraphy
scenario.add_conformable_sequence(['basement', 'sandstone', 'shale'])

# Add observations for basement
scenario.observations('basement')\
    .add_contact([100, 100, 50], comment="Drillhole 1")\
    .add_contact([200, 200, 45], comment="Drillhole 2")\
    .add_contact([300, 300, 48], comment="Drillhole 3")\
    .add_orientation([150, 150, 47], strike=10, dip=5, weight=2.0)

# Add observations for sandstone
scenario.observations('sandstone')\
    .add_contact([100, 100, 75])\
    .add_contact([200, 200, 72])\
    .add_orientation([150, 150, 73], strike=12, dip=8)\
    .add_above_point([150, 150, 85], comment="Known above contact")\
    .add_below_point([150, 150, 65], comment="Known below contact")

# Add fault with observations
scenario.add_fault('fault_1', displacement=100, cuts=['basement', 'sandstone'])

scenario.observations('fault_1')\
    .add_fault_trace([150, 100, 0], comment="Surface trace 1")\
    .add_fault_trace([250, 200, 0], comment="Surface trace 2")\
    .add_fault_orientation([200, 150, 50], strike=90, dip=60)\
    .add_hangingwall_point([210, 150, 50])\
    .add_footwall_point([190, 150, 50])

# Validate and build
scenario.validate()
model = scenario.build()
```

## How It Works

### 1. Observation Storage

Observations are stored in `ObservationCollection` objects, one per feature:

```python
collection = scenario.observations('unit1')  # Gets or creates collection
print(len(collection))  # Number of observations
print(collection)  # ObservationCollection('unit1', 10 observations)
```

### 2. Automatic Conversion

When building the model, observations are automatically converted to LoopStructural's DataFrame format:

```python
# Internal conversion during build
df = scenario.get_all_observations_dataframe()  # Combines all observations
# df has columns: feature_name, X, Y, Z, val, nx, ny, nz, gx, gy, gz, tx, ty, tz, w, coord, polarity
```

### 3. Coordinate Types

Different observation types map to different coordinate constraints:

- **coord=0**: Interface constraint (contacts, fault traces)
- **coord=1**: Gradient constraint (orientations, normals)
- **coord=2**: Inequality constraint (inside/outside, above/below, hangingwall/footwall)

### 4. Integration with Scenario Building

During `scenario.build()`:

1. All observations are collected and converted to a DataFrame
2. The DataFrame is passed to the GeologicalModel
3. Features are built in topological order
4. Each feature uses its specific observations from the combined dataset

## Advantages Over Traditional Approach

### Before (Traditional DataFrame):

```python
import pandas as pd

data = pd.DataFrame({
    'feature_name': ['unit1', 'unit1', 'unit1'],
    'X': [100, 200, 150],
    'Y': [100, 200, 150],
    'Z': [50, 45, 47],
    'coord': [0, 0, 1],
    'val': [np.nan, np.nan, np.nan],
    'gx': [np.nan, np.nan, 0.1],
    'gy': [np.nan, np.nan, 0.2],
    'gz': [np.nan, np.nan, 0.9],
    'w': [1.0, 1.0, 2.0],
    # ... many more columns
})

model = GeologicalModel(origin, maximum)
model.data = data
model.create_and_add_foliation('unit1')
```

### Now (Observation System):

```python
scenario = GeologicalScenario(origin, maximum)
scenario.add_unit('unit1')

scenario.observations('unit1')\
    .add_contact([100, 100, 50])\
    .add_contact([200, 200, 45])\
    .add_orientation([150, 150, 47], gradient=[0.1, 0.2, 0.9], weight=2.0)

model = scenario.build()
```

**Benefits:**
- ✅ Geologically intuitive method names
- ✅ Type-safe with validation
- ✅ No need to remember DataFrame column names
- ✅ Clear separation of data types (contacts vs orientations vs constraints)
- ✅ Easier to read and maintain
- ✅ Less error-prone

## Backwards Compatibility

You can still use traditional DataFrames if needed:

```python
# Add traditional DataFrame data
scenario.add_observations_from_dataframe('unit1', traditional_df)

# Or combine with scenario observations during build
model = scenario.build(data=additional_traditional_df)
```

## Validation

Observations are validated when created:

```python
# This will raise an error - location must be 3D
scenario.observations('unit1').add_contact([100, 200])  # ❌ Error

# This will raise an error - orientation needs strike/dip or gradient
scenario.observations('unit1').add_orientation([100, 200, 50])  # ❌ Error

# This is correct
scenario.observations('unit1').add_orientation([100, 200, 50], strike=45, dip=30)  # ✅
```

## Best Practices

1. **Add observations as you define geology**:
   ```python
   scenario.add_unit('sandstone')
   scenario.observations('sandstone').add_contact([100, 100, 50])
   ```

2. **Use method chaining for multiple observations**:
   ```python
   scenario.observations('sandstone')\
       .add_contact([100, 100, 50])\
       .add_contact([200, 200, 52])\
       .add_orientation([150, 150, 51], strike=45, dip=30)
   ```

3. **Add comments for traceability**:
   ```python
   scenario.observations('fault1')\
       .add_fault_trace([150, 100, 0], comment="GPS point A")\
       .add_fault_orientation([200, 150, 50], strike=90, dip=60, 
                            comment="Outcrop measurement site B")
   ```

4. **Use weights to indicate data quality**:
   ```python
   scenario.observations('unit1')\
       .add_orientation([100, 100, 50], strike=45, dip=30, weight=2.0,
                       comment="High-quality outcrop measurement")\
       .add_orientation([200, 200, 55], strike=50, dip=28, weight=0.5,
                       comment="Uncertain measurement from poor exposure")
   ```

## Next Steps

See the complete examples:
- `examples/5_graph/observations_example.py` - Comprehensive observation system demo
- `examples/5_graph/scenario_with_observations.py` - Full scenario building workflow
- `examples/5_graph/three_tier_example.py` - Three-tier API comparison

For API documentation, see:
- `docs/geological_observations.md` - Detailed API reference
- `docs/geological_topology_graph.md` - Topology system documentation
- `docs/topology_integration_strategy.md` - Integration strategy overview
