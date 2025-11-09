"""
Example: Building a Geological Model with Scenario and Observations

This example shows how to use GeologicalScenario with the observation
system to build a complete geological model with proper data.
"""

import numpy as np
from LoopStructural.modelling.core.geological_scenario import GeologicalScenario

# ==============================================================================
# Setup: Model domain
# ==============================================================================

origin = np.array([0, 0, 0])
maximum = np.array([1000, 1000, 500])

# ==============================================================================
# Create a Geological Scenario with Complete Observations
# ==============================================================================

print("=" * 70)
print("Building Geological Model from Scenario with Observations")
print("=" * 70)

# Create the scenario
scenario = GeologicalScenario(origin, maximum)

# ==============================================================================
# Define Stratigraphy
# ==============================================================================

print("\n1. Defining stratigraphic relationships...")

# Add conformable sequence
scenario.add_conformable_sequence([
    'basement',
    'sandstone',
    'shale',
    'limestone'
])

print("   ✓ Added 4 units in conformable sequence")

# ==============================================================================
# Add Observations for Each Unit
# ==============================================================================

print("\n2. Adding observations for stratigraphic units...")

# Basement observations
print("   - basement: adding contact points and orientations")
scenario.observations('basement')\
    .add_contact([100, 100, 50], comment="Drillhole 1")\
    .add_contact([300, 300, 48], comment="Drillhole 2")\
    .add_contact([500, 500, 45], comment="Drillhole 3")\
    .add_contact([700, 700, 47], comment="Drillhole 4")\
    .add_orientation([400, 400, 47], strike=10, dip=5, weight=2.0, comment="Outcrop")\
    .add_orientation([600, 600, 46], strike=15, dip=8, weight=2.0, comment="Outcrop")

# Sandstone observations
print("   - sandstone: adding contact points, orientations, and constraints")
scenario.observations('sandstone')\
    .add_contact([100, 100, 100], comment="Drillhole 1")\
    .add_contact([300, 300, 98], comment="Drillhole 2")\
    .add_contact([500, 500, 95], comment="Drillhole 3")\
    .add_contact([700, 700, 97], comment="Drillhole 4")\
    .add_orientation([400, 400, 97], strike=10, dip=5, weight=2.0)\
    .add_orientation([600, 600, 96], strike=15, dip=8, weight=2.0)\
    .add_above_point([400, 400, 110], comment="Known above sandstone")\
    .add_below_point([400, 400, 80], comment="Known below sandstone")

# Shale observations
print("   - shale: adding contact points and orientations")
scenario.observations('shale')\
    .add_contact([100, 100, 150], comment="Drillhole 1")\
    .add_contact([300, 300, 148], comment="Drillhole 2")\
    .add_contact([500, 500, 145], comment="Drillhole 3")\
    .add_contact([700, 700, 147], comment="Drillhole 4")\
    .add_orientation([400, 400, 147], strike=12, dip=6, weight=2.0)\
    .add_orientation([600, 600, 146], strike=14, dip=7, weight=2.0)

# Limestone observations
print("   - limestone: adding contact points and orientations")
scenario.observations('limestone')\
    .add_contact([100, 100, 200], comment="Drillhole 1")\
    .add_contact([300, 300, 198], comment="Drillhole 2")\
    .add_contact([500, 500, 195], comment="Drillhole 3")\
    .add_contact([700, 700, 197], comment="Drillhole 4")\
    .add_orientation([400, 400, 197], strike=12, dip=6, weight=2.0)\
    .add_orientation([600, 600, 196], strike=14, dip=7, weight=2.0)

# ==============================================================================
# Add Faults with Observations
# ==============================================================================

print("\n3. Adding faults and fault observations...")

# Add first fault
scenario.add_fault(
    'fault_1',
    displacement=100,
    cuts=['basement', 'sandstone', 'shale']
)

print("   - fault_1: adding trace and orientation observations")
scenario.observations('fault_1')\
    .add_fault_trace([200, 200, 0], comment="Surface trace")\
    .add_fault_trace([200, 400, 0], comment="Surface trace")\
    .add_fault_trace([200, 600, 0], comment="Surface trace")\
    .add_fault_orientation([200, 400, 100], strike=0, dip=70, weight=2.0, comment="Outcrop")\
    .add_fault_orientation([200, 400, 150], strike=5, dip=68, comment="Drillhole")\
    .add_hangingwall_point([220, 400, 100], comment="Confirmed hanging wall")\
    .add_footwall_point([180, 400, 100], comment="Confirmed footwall")

# Add second fault (post-dates first fault, doesn't cut limestone)
scenario.add_fault(
    'fault_2',
    displacement=50,
    cuts=['sandstone', 'shale']
)

print("   - fault_2: adding trace and orientation observations")
scenario.observations('fault_2')\
    .add_fault_trace([600, 200, 0], comment="Surface trace")\
    .add_fault_trace([600, 400, 0], comment="Surface trace")\
    .add_fault_trace([600, 600, 0], comment="Surface trace")\
    .add_fault_orientation([600, 400, 100], strike=0, dip=60, weight=2.0)\
    .add_fault_orientation([600, 400, 150], strike=2, dip=62)

# Define fault network
scenario.add_fault_network(['fault_1', 'fault_2'])

# ==============================================================================
# Validate and Build
# ==============================================================================

print("\n4. Validating scenario...")
try:
    warnings = scenario.validate()
    if warnings:
        print(f"   ⚠ Validation warnings: {warnings}")
    else:
        print("   ✓ Scenario is valid!")
except ValueError as e:
    print(f"   ✗ Validation failed: {e}")
    exit(1)

# Check observation counts
print("\n5. Observation summary:")
total_obs = 0
for feature_name, obs_collection in scenario._observations.items():
    count = len(obs_collection)
    total_obs += count
    print(f"   - {feature_name}: {count} observations")
print(f"   Total: {total_obs} observations")

# Build the model
print("\n6. Building computational model...")
try:
    model = scenario.build(validate=False)  # Already validated
    print(f"   ✓ Model built with {len(model.features)} features")
    
    # List built features
    print("\n   Built features:")
    for feature in model.features:
        print(f"     - {feature.name} ({feature.type})")
    
except Exception as e:
    print(f"   ✗ Build failed: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

# ==============================================================================
# Display Model Information
# ==============================================================================

print("\n7. Model information:")
print(f"   - Number of features: {len(model.features)}")
print(f"   - Topology objects: {len(model.topology)}")

# Get topology information
print("\n8. Topology relationships:")
relationships = model.topology.get_relationships()
print(f"   Total relationships: {len(relationships)}")
for rel in relationships[:10]:  # Show first 10
    from_obj = model.topology.get_object(rel.from_object)
    to_obj = model.topology.get_object(rel.to_object)
    print(f"   - {from_obj.name} {rel.relationship_type} {to_obj.name}")
if len(relationships) > 10:
    print(f"   ... and {len(relationships) - 10} more")

# ==============================================================================
# Export Scenario
# ==============================================================================

print("\n9. Exporting scenario for reproducibility...")
scenario_dict = scenario.to_dict()
print(f"   ✓ Scenario exported with {len(scenario_dict['topology']['objects'])} objects")

print("\n" + "=" * 70)
print("Complete! The model is ready for interpolation and evaluation.")
print("=" * 70)

print("""
Next steps:
1. Call model.update() to interpolate all features
2. Use model.evaluate_model(points) to evaluate at specific points
3. Visualize with model.visualisation.plot_model()
4. Export results using model export methods

The scenario can be saved and reloaded:
- Save: import json; json.dump(scenario_dict, open('scenario.json', 'w'))
- Load: scenario = GeologicalScenario.from_dict(json.load(open('scenario.json')))
""")
