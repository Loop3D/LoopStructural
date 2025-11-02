from LoopStructural.modelling.core.model_graph import GeologicalTopologyGraph, StratigraphicColumnView

graph = GeologicalTopologyGraph()

# Create units in different conformable groups
u1 = graph.add_geological_object('unit1','unit')
u2 = graph.add_geological_object('unit2','unit')
u3 = graph.add_geological_object('unit3','unit')
u4 = graph.add_geological_object('unit4','unit')
u5 = graph.add_geological_object('unit5','unit')

# Group 1: unit1 and unit2 (conformable)
graph.add_relationship(u2.id, u1.id, 'conformable_overlies')

# Unconformity between group 1 and group 2
graph.add_relationship(u3.id, u2.id, 'erode_unconformably_overlies')

# Group 2: unit3 and unit4 (conformable)
graph.add_relationship(u4.id, u3.id, 'conformable_overlies')

# Another unconformity
graph.add_relationship(u5.id, u4.id, 'onlap_unconformably_overlies')

# Add a fault that cuts multiple units
f1 = graph.add_geological_object('fault_1','fault')
graph.add_relationship(f1.id, u1.id, 'cuts')
graph.add_relationship(f1.id, u2.id, 'cuts')

# Test the stratigraphic column view
column = StratigraphicColumnView(graph)

print("=" * 60)
print("Conformable Groups (oldest to youngest):")
print("=" * 60)
try:
    groups = column.identify_conformable_groups()
    for i, group in enumerate(groups, 1):
        print(f"\nGroup {i}:")
        for unit in group:
            print(f"  - {unit.name}")
except ValueError as e:
    print(f"Error: {e}")

print("\n" + "=" * 60)
print("Graph Visualization:")
print("=" * 60)
graph.plot()

# Test cycle detection
print("\n" + "=" * 60)
print("Testing cycle detection:")
print("=" * 60)
graph2 = GeologicalTopologyGraph()
u1 = graph2.add_geological_object('unit1','unit')
u2 = graph2.add_geological_object('unit2','unit')
u3 = graph2.add_geological_object('unit3','unit')

# Create a cycle
graph2.add_relationship(u1.id, u2.id, 'conformable_overlies')
graph2.add_relationship(u2.id, u3.id, 'conformable_overlies')
graph2.add_relationship(u3.id, u1.id, 'conformable_overlies')

column2 = StratigraphicColumnView(graph2)
try:
    groups2 = column2.identify_conformable_groups()
    print("No cycle detected (this shouldn't happen!)")
except ValueError as e:
    print(f"Cycle correctly detected: {e}")
