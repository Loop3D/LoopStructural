Debugging LoopStructural
========================
Surfaces are not displaying
---------------------------
- check to see what the scalar field looks linked
- look at the values of the nodes of the implicit function supports `feature.interpolator.c` if all of the nodes are 0 the interpolator has not solved.
- check that there are sufficient data points to constrain a surface. There needs to be at least two unique value constraints or a single value constraint and a gradient norm constraint.
c

Surface geometry is erroneous
-----------------------------
- Try using a different solver
- Adding more elements
- Try using only orientation constraints
- Try using only value constraints

