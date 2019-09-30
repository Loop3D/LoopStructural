
LoopStructural Tutorials
=============

Geological Features
-------------------

LoopStructural uses an abstract representation of geological objects/surfaces
inside the model using the **GeologicalFeature**. A geological feature
can be something that is interpolated using an interpolator, defined by
a constant value or by a function.

A **GeologicalFeature** has the following attribute: \* name - a unique
identified for the feature within the model \* support - the object that
stores the property values associated with the feature \* data - a
container holding all of the data that is associated with this feature
\* builder - the object that built the feature.

