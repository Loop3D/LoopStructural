# Geological Observations System

## Overview

The `geological_observations` module provides an intuitive, geologically-meaningful way to define observations for geological models. Instead of manually constructing DataFrames with coordinate types and polarity flags, you can use geological terminology like "contact", "orientation", "hanging wall", etc.

## Key Benefits

### 🎯 **Geological Intuition**
Use terms geologists understand:
- `add_contact()` instead of `coord=0, val=nan`
- `add_orientation()` instead of `coord=1, gx=..., gy=..., gz=...`
- `add_hangingwall_point()` instead of `coord=2, val=1.0`

### 🔒 **Type Safety**
Each observation type is a specific class with validation:
- `ContactObservation`
- `OrientationObservation`
- `FaultTraceObservation`
- etc.

### 📝 **Metadata Support**
Every observation can have:
- `weight`: Confidence/quality weighting
- `comment`: Data source, notes

### 🔄 **Automatic Conversion**
Observations automatically convert to LoopStructural's DataFrame format with correct `coord`, `polarity`, and constraint values.

## Basic Usage

### Creating Observations for a Unit

```python
from LoopStructural.modelling.core.geological_scenario import GeologicalScenario

scenario = GeologicalScenario(origin, maximum)

# Add a unit
scenario.add_unit('sandstone')

# Add observations using fluent interface
scenario.observations('sandstone')\
    .add_contact([100, 200, 50], comment="Drillhole DH-001")\
    .add_contact([150, 250, 48], comment="Drillhole DH-002")\
    .add_orientation([125, 225, 49], strike=45, dip=12, comment="Outcrop")\
    .add_above_point([125, 225, 60], comment="Surface sample")\
    .set_thickness(25.0)
```

## Observation Types

### Unit/Foliation Observations

#### Contact Points

Points on the unit boundary/contact surface:

```python
obs.add_contact(
    location=[100, 200, 50],
    weight=1.0,
    scalar_value=None,  # Optional known scalar field value
    comment="Drillhole intersection"
)
```

**LoopStructural equivalent**: `coord=0`, interface constraint

#### Orientation (Strike/Dip)

Orientation measurements from outcrops or oriented core:

```python
obs.add_orientation(
    location=[150, 250, 55],
    strike=45,       # Strike in degrees
    dip=15,          # Dip in degrees
    dip_direction=135,  # Optional, alternative to strike
    weight=2.0,
    comment="Outcrop measurement"
)
```

**LoopStructural equivalent**: `coord=1`, gradient constraint

#### Orientation (Gradient Vector)

Direct gradient vector (normal to surface):

```python
obs.add_orientation(
    location=[150, 250, 55],
    gradient=[0.05, 0.1, 0.9],  # [gx, gy, gz]
    polarity=1.0,   # +1 or -1
    weight=1.0
)
```

**LoopStructural equivalent**: `coord=1`, `gx`, `gy`, `gz` set

#### Orientation (Tangent Vector)

Tangent vector (parallel to surface):

```python
obs.add_orientation(
    location=[150, 250, 55],
    tangent=[1.0, 0.0, -0.05],  # [tx, ty, tz]
    weight=1.0
)
```

**LoopStructural equivalent**: `coord=1`, `tx`, `ty`, `tz` set

#### Inside/Outside Points

Points known to be inside or outside a unit:

```python
# Point inside the unit
obs.add_inside_point([200, 300, 60], comment="Core sample confirmed sandstone")

# Point outside the unit  
obs.add_outside_point([200, 300, 40], comment="Different lithology")
```

**LoopStructural equivalent**: `coord=2`, inequality constraint, `val=1.0` (inside) or `val=-1.0` (outside)

#### Above/Below Points

Points known to be above or below a contact surface:

```python
# Point above the surface
obs.add_above_point([200, 300, 70], comment="Overlying unit")

# Point below the surface
obs.add_below_point([200, 300, 30], comment="Underlying unit")
```

**LoopStructural equivalent**: `coord=2`, inequality constraint, `val=1.0` (above) or `val=-1.0` (below)

#### Thickness

Unit thickness (not location-specific):

```python
obs.set_thickness(25.0, location=[200, 300, 50])  # location optional
```

### Fault Observations

#### Fault Trace

Points on fault trace (surface intersection):

```python
obs.add_fault_trace(
    location=[150, 200, 0],
    trace_direction=[1.0, 1.0, 0.0],  # Optional direction along trace
    weight=2.0,
    comment="Mapped fault trace"
)
```

**LoopStructural equivalent**: `coord=0`, `val=0.0` (on fault surface)

#### Fault Orientation

Fault plane orientation:

```python
obs.add_fault_orientation(
    location=[150, 200, 50],
    strike=90,
    dip=60,
    comment="Fault plane orientation"
)

# Or with direct normal vector
obs.add_fault_orientation(
    location=[150, 200, 50],
    normal_vector=[0.0, 0.5, 0.866],  # Fault plane normal
    comment="Calculated from kinematic analysis"
)
```

**LoopStructural equivalent**: `coord=1`, gradient constraint for fault surface

#### Displacement

Fault displacement measurement:

```python
obs.add_displacement(
    location=[150, 200, 50],
    displacement=100.0,  # Magnitude
    direction=[0.0, 0.0, -1.0],  # Optional slip direction
    comment="Offset marker bed"
)
```

#### Hanging Wall / Footwall

Points in hanging wall or footwall blocks:

```python
# Hanging wall point
obs.add_hangingwall_point([170, 200, 50], comment="Confirmed from drillhole")

# Footwall point
obs.add_footwall_point([130, 200, 50], comment="Confirmed from drillhole")
```

**LoopStructural equivalent**: `coord=2`, inequality constraint

#### Slip Vector

Slip direction observation (e.g., from slickensides):

```python
obs.add_slip_vector(
    location=[150, 200, 50],
    slip_vector=[0.0, 0.0, -1.0],  # Pure dip slip
    comment="Slickensides measurement"
)
```

**LoopStructural equivalent**: `coord=1`, direction constraint

## Usage Patterns

### Pattern 1: Fluent Interface (Method Chaining)

```python
scenario.observations('unit1')\
    .add_contact([100, 100, 50])\
    .add_contact([200, 150, 48])\
    .add_orientation([150, 125, 49], strike=45, dip=12)\
    .add_above_point([150, 125, 60])\
    .set_thickness(20.0)
```

### Pattern 2: Get Collection Then Add

```python
obs = scenario.observations('unit1')

for point in contact_points:
    obs.add_contact(point, weight=1.0)

for orient in orientations:
    obs.add_orientation(orient['loc'], orient['strike'], orient['dip'])
```

### Pattern 3: Programmatic Generation

```python
obs = scenario.observations('synthetic_layer')

# Generate observations along a mathematical surface
for i in range(100):
    x = np.random.uniform(0, 1000)
    y = np.random.uniform(0, 1000)
    z = 50 + 0.1*x + 0.05*y + np.random.normal(0, 2)
    
    obs.add_contact([x, y, z])
```

### Pattern 4: Mix with Traditional Data

```python
# Load traditional LoopStructural DataFrame
traditional_data = pd.read_csv('observations.csv')

# Add to scenario
scenario.add_observations_from_dataframe('unit1', traditional_data)

# Still add new-style observations
scenario.observations('unit1')\
    .add_orientation([300, 300, 55], strike=48, dip=15)
```

## Complete Example

```python
from LoopStructural.modelling.core.geological_scenario import GeologicalScenario
import numpy as np

# Create scenario
scenario = GeologicalScenario(
    origin=np.array([0, 0, 0]),
    maximum=np.array([1000, 1000, 500])
)

# Add stratigraphic units
scenario.add_conformable_sequence(['basement', 'sandstone', 'shale'])

# Add observations for sandstone
scenario.observations('sandstone')\
    .add_contact([100, 100, 50], comment="DH-001")\
    .add_contact([200, 150, 48], comment="DH-002")\
    .add_contact([300, 200, 49], comment="DH-003")\
    .add_orientation([150, 125, 49], strike=45, dip=12, weight=2.0)\
    .add_orientation([250, 175, 48], strike=48, dip=15, weight=2.0)\
    .add_above_point([150, 125, 70])\
    .add_below_point([150, 125, 30])\
    .set_thickness(25.0)

# Add observations for shale  
scenario.observations('shale')\
    .add_contact([100, 100, 75])\
    .add_contact([200, 150, 72])\
    .add_orientation([150, 125, 73], strike=50, dip=18)\
    .add_inside_point([150, 125, 73], comment="Core confirmed shale")

# Add fault with observations
scenario.add_fault('main_fault', displacement=100, cuts=['sandstone', 'shale'])

scenario.observations('main_fault')\
    .add_fault_trace([150, 100, 0], comment="Surface trace")\
    .add_fault_trace([200, 200, 0])\
    .add_fault_trace([250, 300, 0])\
    .add_fault_orientation([200, 200, 50], strike=90, dip=60)\
    .add_hangingwall_point([220, 200, 50])\
    .add_footwall_point([180, 200, 50])\
    .add_slip_vector([200, 200, 50], slip_vector=[0.0, 0.0, -1.0])

# Build model with all observations
model = scenario.build()
model.update()
```

## API Reference

### ObservationCollection

The main class for managing observations for a feature.

#### Methods

- `add_contact(location, weight, scalar_value, comment)`
- `add_orientation(location, strike, dip, dip_direction, gradient, tangent, polarity, weight, comment)`
- `add_inside_point(location, weight, comment)`
- `add_outside_point(location, weight, comment)`
- `add_above_point(location, weight, comment)`
- `add_below_point(location, weight, comment)`
- `set_thickness(thickness, location)`
- `add_fault_trace(location, trace_direction, weight, comment)`
- `add_fault_orientation(location, strike, dip, dip_direction, normal_vector, weight, comment)`
- `add_displacement(location, displacement, direction, weight, comment)`
- `add_hangingwall_point(location, weight, comment)`
- `add_footwall_point(location, weight, comment)`
- `add_slip_vector(location, slip_vector, weight, comment)`
- `to_dataframe()` - Convert to LoopStructural DataFrame format

### Observation Classes

- `ContactObservation`
- `OrientationObservation`
- `InsideOutsideObservation`
- `AboveBelowObservation`
- `ThicknessObservation`
- `FaultTraceObservation`
- `FaultOrientationObservation`
- `DisplacementObservation`
- `HangingwallFootwallObservation`
- `SlipVectorObservation`

Each class has:
- `location`: np.ndarray [x, y, z]
- `obs_type`: ObservationType enum
- `weight`: float (default 1.0)
- `comment`: Optional str

Plus type-specific attributes (e.g., `strike`, `dip`, `gradient`, etc.)

## Best Practices

### 1. Use Geological Terms

```python
# Good - geologically intuitive
obs.add_contact([100, 200, 50])
obs.add_orientation([100, 200, 50], strike=45, dip=30)

# Avoid - requires knowing LoopStructural internals
df = pd.DataFrame({'coord': [0, 1], 'val': [np.nan, np.nan], ...})
```

### 2. Add Comments for Traceability

```python
obs.add_contact([100, 200, 50], comment="Drillhole DH-001, 50m depth")
obs.add_orientation([150, 250, 55], strike=45, dip=30, comment="Outcrop measurement, GPS: ...")
```

### 3. Weight by Data Quality

```python
# High quality outcrop measurement
obs.add_orientation([...], strike=45, dip=30, weight=3.0, comment="Excellent exposure")

# Lower quality interpretation
obs.add_orientation([...], strike=50, dip=25, weight=0.5, comment="Estimated from photos")
```

### 4. Use Appropriate Observation Types

```python
# For unit boundaries
obs.add_contact([...])

# For known presence/absence
obs.add_inside_point([...])  # Core confirms this unit
obs.add_outside_point([...])  # Core confirms NOT this unit

# For relative position
obs.add_above_point([...])  # Above the contact
obs.add_below_point([...])  # Below the contact
```

## See Also

- [Geological Scenario Guide](./geological_scenario.md)
- [Topology Graph Guide](./geological_topology_graph.md)
- [Model Building Guide](./model_building.md)
- Example: `examples/5_graph/observations_example.py`
