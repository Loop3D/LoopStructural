"""
Analysis
========

Various tools for analysing loopstructural models, including calculating fault intersections and fault toplogies
"""
from ._fault_displacement import displacement_missfit
from ._fault_intersection import calculate_fault_intersections
from ._topology import calculate_fault_topology_matrix