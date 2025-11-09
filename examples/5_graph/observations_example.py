"""
Example: Using Geological Observations in Scenarios

This example demonstrates the new observation system for defining
geological data in an intuitive, geologically-meaningful way.
"""

import numpy as np
from LoopStructural.modelling.core.geological_scenario import GeologicalScenario
from LoopStructural.modelling.core.geological_observations import ObservationCollection

# ==============================================================================
# Setup: Model domain
# ==============================================================================

origin = np.array([0, 0, 0])
maximum = np.array([1000, 1000, 500])

# ==============================================================================
# Example 1: Stratigraphic Units with Various Observation Types
# ==============================================================================

print("=" * 70)
print("Example 1: Stratigraphic Unit Observations")
print("=" * 70)

scenario = GeologicalScenario(origin, maximum)

# Add units
scenario.add_unit('basement')
scenario.add_unit('sandstone')
scenario.add_unit('shale')
scenario.add_unit('limestone')

# Define conformable sequence
scenario.add_conformable_sequence(['basement', 'sandstone', 'shale', 'limestone'])

# ===== Add observations for sandstone unit =====
print("\nAdding observations for 'sandstone' unit:")

# Method 1: Fluent interface
scenario.observations('sandstone')\
    .add_contact([100, 100, 50], comment="Drillhole 1")\
    .add_contact([200, 150, 45], comment="Drillhole 2")\
    .add_contact([300, 200, 48], comment="Drillhole 3")\
    .add_orientation([150, 125, 47], strike=45, dip=15, comment="Outcrop measurement")\
    .add_orientation([250, 175, 46], strike=50, dip=18, comment="Outcrop measurement")\
    .add_above_point([150, 125, 60], comment="Known to be above contact")\
    .add_below_point([150, 125, 30], comment="Known to be below contact")\
    .set_thickness(25.0)

print(f"  Added {len(scenario.observations('sandstone'))} observations for sandstone")

# ===== Add observations for shale unit =====
print("\nAdding observations for 'shale' unit:")

# Method 2: Get collection and add observations
shale_obs = scenario.observations('shale')

# Contact points from drillholes
contact_points = [
    [100, 100, 75],
    [200, 150, 70],
    [300, 200, 73],
    [150, 300, 72]
]
for i, point in enumerate(contact_points, 1):
    shale_obs.add_contact(point, comment=f"Drillhole {i}")

# Orientation measurements from outcrops
orientations = [
    {'location': [150, 125, 72], 'strike': 50, 'dip': 20},
    {'location': [250, 200, 71], 'strike': 48, 'dip': 22},
]
for orient in orientations:
    shale_obs.add_orientation(
        location=orient['location'],
        strike=orient['strike'],
        dip=orient['dip'],
        weight=2.0,  # Higher weight for good outcrop data
        comment="Outcrop"
    )

# Gradient vectors (from seismic interpretation)
shale_obs.add_orientation(
    location=[400, 400, 70],
    gradient=[0.05, 0.1, 0.9],  # Approximate gradient vector
    comment="Seismic interpretation"
)

# Inequality constraints
shale_obs.add_inside_point([200, 200, 71], comment="Core sample confirms shale")

print(f"  Added {len(shale_obs)} observations for shale")

# ===== Add observations for limestone unit =====
print("\nAdding observations for 'limestone' unit:")

limestone_obs = scenario.observations('limestone')

# Using gradient and tangent vectors
limestone_obs\
    .add_contact([100, 100, 95])\
    .add_contact([200, 150, 92])\
    .add_orientation([150, 125, 93], strike=52, dip=18)\
    .add_orientation([250, 175, 91], gradient=[0.03, 0.08, 0.95])\
    .add_orientation([350, 225, 90], tangent=[1.0, 0.0, -0.03], comment="Tangent to bedding")

print(f"  Added {len(limestone_obs)} observations for limestone")

# ==============================================================================
# Example 2: Fault Observations
# ==============================================================================

print("\n" + "=" * 70)
print("Example 2: Fault Observations")
print("=" * 70)

# Add fault with fault-specific observations
scenario.add_fault('main_fault', displacement=100, cuts=['sandstone', 'shale'])

print("\nAdding fault-specific observations:")

fault_obs = scenario.observations('main_fault')

# Fault trace observations (surface intersections)
fault_trace_points = [
    [150, 100, 0],   # Surface trace point 1
    [200, 200, 0],   # Surface trace point 2
    [250, 300, 0],   # Surface trace point 3
]
for point in fault_trace_points:
    fault_obs.add_fault_trace(point, comment="Mapped trace")

# Fault orientation measurements
fault_obs\
    .add_fault_orientation([200, 200, 50], strike=90, dip=60, comment="Outcrop")\
    .add_fault_orientation([250, 250, 75], strike=88, dip=62, comment="Drillhole intersection")

# Hanging wall and footwall observations
fault_obs\
    .add_hangingwall_point([220, 200, 50], comment="Confirmed hanging wall")\
    .add_footwall_point([180, 200, 50], comment="Confirmed footwall")

# Slip vector observation
fault_obs.add_slip_vector(
    [200, 200, 50],
    slip_vector=[0.0, 0.0, -1.0],  # Pure dip slip
    comment="Slickensides measurement"
)

# Displacement measurements
fault_obs.add_displacement(
    [200, 200, 50],
    displacement=100,
    direction=[0.0, 0.0, -1.0],
    comment="Offset marker bed"
)

print(f"  Added {len(fault_obs)} fault observations")

# ==============================================================================
# Example 3: Mixing Observation Methods
# ==============================================================================

print("\n" + "=" * 70)
print("Example 3: Mixing Traditional and New Observation Methods")
print("=" * 70)

# You can also add observations from traditional DataFrames
import pandas as pd

# Traditional LoopStructural DataFrame format
traditional_data = pd.DataFrame({
    'feature_name': ['basement'] * 5,
    'X': [100, 200, 300, 400, 500],
    'Y': [100, 150, 200, 250, 300],
    'Z': [10, 12, 11, 13, 12],
    'val': [np.nan] * 5,  # Interface points
    'nx': [np.nan] * 5,
    'ny': [np.nan] * 5,
    'nz': [np.nan] * 5,
    'gx': [np.nan] * 5,
    'gy': [np.nan] * 5,
    'gz': [np.nan] * 5,
    'tx': [np.nan] * 5,
    'ty': [np.nan] * 5,
    'tz': [np.nan] * 5,
    'w': [1.0] * 5,
    'coord': [0] * 5,  # Interface constraint
    'polarity': [1.0] * 5,
})

# Add traditional data
scenario.add_observations_from_dataframe('basement', traditional_data)

# And still add new-style observations
scenario.observations('basement')\
    .add_orientation([300, 200, 11], strike=48, dip=12)\
    .add_orientation([400, 250, 13], strike=50, dip=10)

print("✓ Combined traditional DataFrame with new-style observations")

# ==============================================================================
# Example 4: Getting Combined Data
# ==============================================================================

print("\n" + "=" * 70)
print("Example 4: Viewing All Observations")
print("=" * 70)

# Get all observations as a combined DataFrame
all_obs = scenario.get_all_observations_dataframe()

print(f"\nTotal observations: {len(all_obs)}")
print(f"Features with observations: {all_obs['feature_name'].unique().tolist()}")
print(f"\nObservation types:")
print(all_obs.groupby('feature_name').size())

# Show summary statistics
print(f"\nSummary:")
print(f"  Contact points (coord=0): {len(all_obs[all_obs['coord'] == 0])}")
print(f"  Orientation constraints (coord=1): {len(all_obs[all_obs['coord'] == 1])}")
print(f"  Inequality constraints (coord=2): {len(all_obs[all_obs['coord'] == 2])}")

# ==============================================================================
# Example 5: Building the Model with Observations
# ==============================================================================

print("\n" + "=" * 70)
print("Example 5: Building Model from Observations")
print("=" * 70)

# Validate the scenario
print("\nValidating scenario...")
try:
    warnings = scenario.validate()
    if warnings:
        print(f"  Warnings: {warnings}")
    else:
        print("  ✓ Validation passed")
except ValueError as e:
    print(f"  ✗ Validation failed: {e}")

# Build the model
print("\nBuilding model...")
# model = scenario.build()  # Would build with all observations
# model.update()  # Would run interpolation
print("  ✓ Model would be built with all observations included")

# ==============================================================================
# Example 6: Advanced - Programmatic Observation Generation
# ==============================================================================

print("\n" + "=" * 70)
print("Example 6: Programmatic Observation Generation")
print("=" * 70)

# Create a new unit with programmatically generated observations
scenario2 = GeologicalScenario(origin, maximum)
scenario2.add_unit('synthetic_unit')

obs = scenario2.observations('synthetic_unit')

# Generate synthetic contact points along a dipping surface
n_points = 20
for i in range(n_points):
    x = np.random.uniform(0, 1000)
    y = np.random.uniform(0, 1000)
    z = 100 + 0.1 * x + 0.05 * y + np.random.normal(0, 2)  # Dipping plane + noise
    
    obs.add_contact([x, y, z], weight=1.0)

# Add orientations at regular grid
for x in range(100, 1000, 200):
    for y in range(100, 1000, 200):
        z = 100 + 0.1 * x + 0.05 * y
        obs.add_orientation([x, y, z], strike=45, dip=6, weight=2.0)

print(f"Generated {len(obs)} synthetic observations")

# ==============================================================================
# Summary
# ==============================================================================

print("\n" + "=" * 70)
print("SUMMARY: Observation System Benefits")
print("=" * 70)

print("""
✓ Geological Intuition: Use geological terminology
  - add_contact(), add_orientation(), add_inside_point()
  - add_fault_trace(), add_hangingwall_point()
  
✓ Flexible Data Entry:
  - Fluent interface for method chaining
  - Individual observations or programmatic generation
  - Mix with traditional DataFrames
  
✓ Type Safety: Specific observation types prevent errors
  - ContactObservation, OrientationObservation, etc.
  - Validation at creation time
  
✓ Metadata Support: Add comments and weights
  - Track data source and quality
  - Weight observations by confidence
  
✓ Automatic Conversion: Converts to LoopStructural format
  - No manual DataFrame construction
  - Handles coord, polarity, etc. automatically
  
✓ Feature-Specific: Observations attached to features
  - Clear association between data and geology
  - Easy to review what data supports each feature
""")

print("\n" + "=" * 70)
print("Example complete!")
print("=" * 70)
