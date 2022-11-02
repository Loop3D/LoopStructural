Structural Frames
=================

Structural frames form the basis of many of the geological features that are built using LoopStructural.
A structural frame is a curvilinear coordinate system with three coordinates:
* Major structrural feature
* Structural direction
* Additional/Intermediate direction

These coordinates correspond to the different elements of geological structures that are observed and recorded by geologists in the field.
For example, when modelling faults the major structural feature is the fault surface and the structural direction is the fault slip direction.
For folds, the major structural is the axial foliation and the structural direction is the fold axis.

Each coordinate of the structural frame are represented by three implicit functions. 
The major structural feature is interpolated first as this is the field usually associated with more observations.
The structural direction is interpolated using any associated observations, or conceptual knowledge (for example the expected fault slip direction or fold axis) and an additonal constraint to enforce orthogonality between the structural direction and the already interpolated major structural feature. 
