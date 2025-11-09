"""
Example: Three-Tier Topology Integration

This example demonstrates the three different ways to use the topology graph
in LoopStructural, from low-level automatic tracking to high-level declarative API.
"""

import numpy as np
from LoopStructural import GeologicalModel
from LoopStructural.modelling.core.geological_scenario import GeologicalScenario

# ==============================================================================
# Setup: Common data and domain
# ==============================================================================

origin = np.array([0, 0, 0])
maximum = np.array([1000, 1000, 500])

# Simulated observation data (in practice, load from file)
# ... would have columns: feature_name, X, Y, Z, nx, ny, nz, val, etc.

# ==============================================================================
# TIER 1: Low-Level - Automatic Topology Tracking
# ==============================================================================
# Use the standard GeologicalModel API. Topology is tracked automatically
# in the background. Good for simple models and backward compatibility.
# ==============================================================================

print("=" * 70)
print("TIER 1: Low-Level Automatic Tracking")
print("=" * 70)

model_auto = GeologicalModel(origin, maximum)

# Add features as usual - topology tracked automatically
unit1 = model_auto.create_and_add_foliation('unit1', interpolatortype='FDI')
unit2 = model_auto.create_and_add_foliation('unit2', interpolatortype='FDI')
fault1 = model_auto.create_and_add_fault('fault1', displacement=100)

# Topology is automatically maintained
print(f"\nObjects in topology: {len(model_auto.topology)}")
print(f"Relationships: {len(model_auto.topology.get_relationships())}")

# Query the automatically created topology
faults_cutting_unit1 = model_auto.topology.get_relationships(
    to_object_id='unit1',
    relationship_type='cuts'
)
print(f"Faults cutting unit1: {[r.from_object for r in faults_cutting_unit1]}")

# Validate topology
warnings = model_auto.topology.validate_topology()
print(f"Validation warnings: {len(warnings)}")

# ==============================================================================
# TIER 2: Mid-Level - Explicit Relationship Management
# ==============================================================================
# Use GeologicalModel with explicit relationship definition. Gives more
# control over topology while still using the familiar model API.
# ==============================================================================

print("\n" + "=" * 70)
print("TIER 2: Mid-Level Explicit Relationships")
print("=" * 70)

model_explicit = GeologicalModel(origin, maximum)

# Create features
model_explicit.create_and_add_foliation('basement')
model_explicit.create_and_add_foliation('sandstone')
model_explicit.create_and_add_foliation('shale')
model_explicit.create_and_add_foliation('limestone')

# Explicitly define conformable relationships
model_explicit.topology.add_relationship('sandstone', 'basement', 'conformable_overlies')
model_explicit.topology.add_relationship('shale', 'sandstone', 'conformable_overlies')

# Add unconformity
model_explicit.topology.add_relationship('limestone', 'shale', 'erode_unconformably_overlies')

# Add fault with explicit cutting relationships
model_explicit.create_and_add_fault('main_fault', displacement=150)
model_explicit.topology.add_relationship('main_fault', 'basement', 'cuts')
model_explicit.topology.add_relationship('main_fault', 'sandstone', 'cuts')
model_explicit.topology.add_relationship('main_fault', 'shale', 'cuts')
# Note: fault doesn't cut limestone (post-fault)

# Query the topology
from LoopStructural.modelling.core.model_graph import StratigraphicColumnView

column_view = StratigraphicColumnView(model_explicit.topology)

print(f"\nUnits: {[u.name for u in column_view.get_units()]}")

# Get conformable groups (separated by unconformities)
try:
    groups = column_view.identify_conformable_groups()
    print("\nConformable Groups:")
    for i, group in enumerate(groups, 1):
        print(f"  Group {i}: {[u.name for u in group]}")
except ValueError as e:
    print(f"Error identifying groups: {e}")

# Get unconformities
unconformities = column_view.get_unconformities()
print("\nUnconformities:")
for upper, lower in unconformities:
    upper_obj = model_explicit.topology.get_object(upper)
    lower_obj = model_explicit.topology.get_object(lower)
    print(f"  {upper_obj.name} unconformably overlies {lower_obj.name}")

# Visualize topology
print("\nVisualizing topology graph...")
# model_explicit.topology.plot()  # Uncomment to show plot

# ==============================================================================
# TIER 3: High-Level - Scenario-Based Declarative API
# ==============================================================================
# Use GeologicalScenario for fully declarative model definition.
# Define the geology first, then build the computational model.
# Best for complex models, reproducibility, and validation.
# ==============================================================================

print("\n" + "=" * 70)
print("TIER 3: High-Level Scenario-Based API")
print("=" * 70)

# Create a scenario
scenario = GeologicalScenario(origin, maximum)

# Define conformable sequences
print("\nDefining stratigraphy...")
scenario.add_conformable_sequence([
    'basement',
    'lower_sandstone',
    'middle_shale',
    'upper_sandstone'
])

# Add first unconformity
scenario.add_unconformity('post_unc_limestone', 'upper_sandstone', type='erode')

# Add post-unconformity sequence
scenario.add_conformable_sequence([
    'post_unc_limestone',
    'coal_seam',
    'top_sandstone'
])

# Add second unconformity
scenario.add_unconformity('recent_sediment', 'top_sandstone', type='onlap')

# Add faults
print("Defining fault network...")
scenario.add_fault(
    'fault_1',
    displacement=200,
    cuts=['basement', 'lower_sandstone', 'middle_shale', 'upper_sandstone'],
    major_axis=500,
    minor_axis=200
)

scenario.add_fault(
    'fault_2',
    displacement=100,
    cuts=['middle_shale', 'upper_sandstone', 'post_unc_limestone'],
    major_axis=300,
    minor_axis=150
)

scenario.add_fault(
    'fault_3',
    displacement=50,
    cuts=['post_unc_limestone', 'coal_seam', 'top_sandstone'],
    major_axis=200,
    minor_axis=100
)

# Define fault network (temporal relationships)
scenario.add_fault_network(['fault_1', 'fault_2', 'fault_3'])

# Add some observations to make the scenario buildable
print("Adding observations to features...")

# Add observations for basement
scenario.observations('basement')\
    .add_contact([100, 100, 50])\
    .add_contact([500, 500, 45])\
    .add_orientation([300, 300, 47], strike=0, dip=10)

# Add observations for lower_sandstone
scenario.observations('lower_sandstone')\
    .add_contact([100, 100, 100])\
    .add_contact([500, 500, 95])\
    .add_orientation([300, 300, 97], strike=0, dip=10)

# Add fault observations
scenario.observations('fault_1')\
    .add_fault_trace([200, 200, 0])\
    .add_fault_trace([400, 400, 0])\
    .add_fault_orientation([300, 300, 100], strike=90, dip=60)

# Set interpolation parameters for specific features
scenario.set_feature_parameters('coal_seam', nelements=5000, interpolatortype='DFI')

# Validate the scenario before building
print("\nValidating scenario...")
try:
    warnings = scenario.validate()
    if warnings:
        print(f"Validation warnings: {warnings}")
    else:
        print("✓ Scenario is valid!")
except ValueError as e:
    print(f"✗ Validation failed: {e}")
    exit(1)

# Get stratigraphic information before building
groups = scenario.get_stratigraphic_groups()
print("\nConformable Groups in Scenario:")
for i, group in enumerate(groups, 1):
    print(f"  Group {i}: {group}")

# Build the computational model from the scenario
print("\nBuilding computational model from scenario...")
print("Note: Only features with observations will be built")
model_scenario = scenario.build()  # Build with attached observations

# The built model now has the topology embedded
print(f"Built model with {len(model_scenario.features)} features")
# model_scenario.update()  # Interpolate features
# result = model_scenario.evaluate_model(evaluation_points)

# Visualize the scenario
print("\nVisualizing scenario topology...")
# scenario.plot_scenario(figsize=(14, 10))  # Uncomment to show plot

# Export scenario for reproducibility
scenario_dict = scenario.to_dict()
print(f"\nScenario exported with {len(scenario_dict['topology']['objects'])} objects")

# ==============================================================================
# Comparison: When to Use Each Tier
# ==============================================================================

print("\n" + "=" * 70)
print("WHEN TO USE EACH TIER")
print("=" * 70)

print("""
TIER 1 (Low-Level Automatic):
  ✓ Simple models with obvious relationships
  ✓ Quick prototyping
  ✓ Backward compatibility with existing code
  ✓ When you trust default relationship inference
  ✗ Requires data in traditional DataFrame format
  
TIER 2 (Mid-Level Explicit):
  ✓ Models needing specific relationship control
  ✓ Complex fault networks
  ✓ When you want to query topology during building
  ✓ Debugging and validation during development
  ✗ Still requires data in traditional DataFrame format
  
TIER 3 (High-Level Scenario):
  ✓ Complex models with many relationships
  ✓ When you want validation before computation
  ✓ Reproducible model definitions
  ✓ Collaborative work (save/share scenarios)
  ✓ Multiple alternative scenarios
  ✓ When geology definition should be separate from computation
  ✓ Geologically-intuitive observation API
  ✓ Can attach observations directly to features

NOTE: All tiers require observational data for features to be built.
TIER 3 provides the most intuitive way to add observations using
the scenario.observations(feature_name) API.
""")

# ==============================================================================
# Advanced: Combining Tiers
# ==============================================================================

print("\n" + "=" * 70)
print("ADVANCED: Combining Tiers")
print("=" * 70)

# You can mix approaches as needed
scenario_mixed = GeologicalScenario(origin, maximum)

# Define basic stratigraphy with high-level API
scenario_mixed.add_conformable_sequence(['unit_a', 'unit_b', 'unit_c'])

# Build the model
# model_mixed = scenario_mixed.build()

# Then add additional features with mid-level API
# model_mixed.create_and_add_fault('late_fault', displacement=75)
# model_mixed.topology.add_relationship('late_fault', 'unit_c', 'cuts')

print("""
You can:
1. Start with Scenario to define main structure
2. Build the model
3. Add additional features with standard API
4. Query and modify topology as needed
5. Re-validate before final interpolation

This gives maximum flexibility!
""")

print("\n" + "=" * 70)
print("Example complete!")
print("=" * 70)
