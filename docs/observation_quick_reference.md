# Geological Observation System - Quick Reference

## Creating a Scenario with Observations

```python
from LoopStructural.modelling.core.geological_scenario import GeologicalScenario
import numpy as np

# 1. Create scenario
scenario = GeologicalScenario(
    origin=np.array([0, 0, 0]),
    maximum=np.array([1000, 1000, 500])
)

# 2. Define geology
scenario.add_conformable_sequence(['unit1', 'unit2', 'unit3'])
scenario.add_fault('fault1', displacement=100, cuts=['unit1', 'unit2'])

# 3. Add observations
scenario.observations('unit1')\
    .add_contact([100, 100, 50])\
    .add_orientation([100, 100, 50], strike=45, dip=30)

scenario.observations('fault1')\
    .add_fault_trace([150, 100, 0])\
    .add_fault_orientation([150, 100, 50], strike=90, dip=60)

# 4. Build model
model = scenario.build()
```

## Unit Observation Methods

| Method | Purpose | Required Parameters | Optional Parameters |
|--------|---------|---------------------|---------------------|
| `.add_contact(location, ...)` | Point on contact surface | `location` (xyz) | `weight`, `scalar_value`, `comment` |
| `.add_orientation(location, ...)` | Strike/dip or gradient | `location` + (`strike`+`dip` OR `gradient` OR `tangent`) | `dip_direction`, `polarity`, `weight`, `comment` |
| `.add_inside_point(location, ...)` | Point inside unit | `location` | `weight`, `comment` |
| `.add_outside_point(location, ...)` | Point outside unit | `location` | `weight`, `comment` |
| `.add_above_point(location, ...)` | Point above contact | `location` | `weight`, `comment` |
| `.add_below_point(location, ...)` | Point below contact | `location` | `weight`, `comment` |
| `.set_thickness(thickness, ...)` | Unit thickness | `thickness` | `location` |

## Fault Observation Methods

| Method | Purpose | Required Parameters | Optional Parameters |
|--------|---------|---------------------|---------------------|
| `.add_fault_trace(location, ...)` | Point on fault trace | `location` (xyz) | `trace_direction`, `weight`, `comment` |
| `.add_fault_orientation(location, ...)` | Fault plane orientation | `location` + (`strike`+`dip` OR `normal_vector`) | `dip_direction`, `weight`, `comment` |
| `.add_displacement(location, ...)` | Displacement measurement | `location`, `displacement` | `direction`, `weight`, `comment` |
| `.add_hangingwall_point(location, ...)` | Point in hanging wall | `location` | `weight`, `comment` |
| `.add_footwall_point(location, ...)` | Point in footwall | `location` | `weight`, `comment` |
| `.add_slip_vector(location, ...)` | Slip direction | `location`, `slip_vector` | `weight`, `comment` |

## Common Patterns

### Pattern 1: Multiple Contact Points
```python
contact_points = [[100, 100, 50], [200, 200, 52], [300, 300, 48]]
obs = scenario.observations('unit1')
for point in contact_points:
    obs.add_contact(point)
```

### Pattern 2: Strike/Dip from Field Measurements
```python
measurements = [
    {'location': [100, 100, 50], 'strike': 45, 'dip': 30, 'quality': 'good'},
    {'location': [200, 200, 55], 'strike': 50, 'dip': 28, 'quality': 'fair'},
]
obs = scenario.observations('unit1')
for m in measurements:
    weight = 2.0 if m['quality'] == 'good' else 1.0
    obs.add_orientation(m['location'], strike=m['strike'], dip=m['dip'], 
                       weight=weight, comment=m['quality'])
```

### Pattern 3: Gradient Vectors from Seismic
```python
gradients = [
    {'xyz': [100, 100, 50], 'gradient': [0.1, 0.2, 0.9]},
    {'xyz': [200, 200, 55], 'gradient': [0.08, 0.15, 0.95]},
]
obs = scenario.observations('unit1')
for g in gradients:
    obs.add_orientation(g['xyz'], gradient=g['gradient'], 
                       comment="Seismic interpretation")
```

### Pattern 4: Fault Trace Mapping
```python
trace_points = [[100, 50, 0], [150, 100, 0], [200, 150, 0]]
fault_obs = scenario.observations('fault1')
for i, point in enumerate(trace_points, 1):
    fault_obs.add_fault_trace(point, comment=f"GPS point {i}")
```

### Pattern 5: Drillhole Data
```python
# Drillhole intersects multiple units
drillhole = [
    {'unit': 'shale', 'depth': 100, 'location': [500, 500, 400]},
    {'unit': 'sandstone', 'depth': 150, 'location': [500, 500, 350]},
    {'unit': 'basement', 'depth': 200, 'location': [500, 500, 300]},
]

for entry in drillhole:
    scenario.observations(entry['unit'])\
        .add_contact(entry['location'], 
                    comment=f"Drillhole at {entry['depth']}m")
```

## Orientation Input Formats

### Format 1: Strike and Dip
```python
.add_orientation([x, y, z], strike=45, dip=30)
```

### Format 2: Dip Direction and Dip
```python
.add_orientation([x, y, z], dip_direction=135, dip=30)
```

### Format 3: Gradient Vector
```python
.add_orientation([x, y, z], gradient=[gx, gy, gz])
```

### Format 4: Tangent Vector
```python
.add_orientation([x, y, z], tangent=[tx, ty, tz])
```

## Checking Observations

```python
# Get observation collection
obs = scenario.observations('unit1')

# Count observations
print(f"Number of observations: {len(obs)}")

# Convert to DataFrame to inspect
df = obs.to_dataframe()
print(df.head())

# Get all observations as combined DataFrame
all_obs = scenario.get_all_observations_dataframe()
print(f"Total observations: {len(all_obs)}")
print(f"Features with observations: {all_obs['feature_name'].unique()}")
```

## Error Handling

```python
# Will raise ValueError - location must be 3D [x, y, z]
try:
    scenario.observations('unit1').add_contact([100, 200])
except ValueError as e:
    print(f"Error: {e}")

# Will raise ValueError - orientation needs strike/dip or gradient
try:
    scenario.observations('unit1').add_orientation([100, 200, 50])
except ValueError as e:
    print(f"Error: {e}")

# Will raise ValueError - slip_vector is required
try:
    scenario.observations('fault1').add_slip_vector([100, 200, 50])
except ValueError as e:
    print(f"Error: {e}")
```

## Combining with Traditional Data

```python
# Option 1: Add traditional DataFrame to specific feature
import pandas as pd
traditional_data = pd.DataFrame({...})
scenario.add_observations_from_dataframe('unit1', traditional_data)

# Option 2: Provide additional data during build
extra_data = pd.DataFrame({...})
model = scenario.build(data=extra_data)

# Both scenario observations and extra_data will be combined
```

## Workflow Summary

1. **Define Geology** → `scenario.add_conformable_sequence([...])`
2. **Add Observations** → `scenario.observations('name').add_...()`
3. **Validate** → `scenario.validate()`
4. **Build** → `model = scenario.build()`
5. **Interpolate** → `model.update()`
6. **Evaluate** → `model.evaluate_model(points)`

## Tips

- Use `comment` parameter to document data sources
- Use `weight` to indicate data quality/certainty
- Add observations incrementally as you define geology
- Check observation counts before building: `len(scenario.observations('unit1'))`
- Every feature that will be built needs at least some observations
- Use `.validate()` before `.build()` to catch topology errors early
