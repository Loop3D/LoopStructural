Background
==========
LoopStructural is an opensource Python 3.9+ library for 3D geological modelling where geological objects are represented by an implicit function. 
Where the implicit function represents the distance or pseudodistance to a reference horizion. 
There is no analytical function that can be used to represent the geometry of geological objects, so the implicit function has to be approximated. 
The implicit function is approximated from the observations of the surface
* location of the geological feature can be used to define the value of the implicit function at a location
* orientation of the geological feature (e.g. strike and dip of a surface) can be used to constrain the gradient of the implicit function

Complex features such as folds, faults and unconformities need additional information to be incorporated into the model. 
Using LoopStructural the overprinting relationships between folds can be incorporated by applying specific fold constraints. 
Faults require the fault displacement and fault slip direction vector to be known to incorporate these directly into the model.  
Unconformities are incorporated by specifying the value of an existing implicit function that defines the unconformity surface. 

Intepreting a map into model input
----------------------------------
To build a geological model the user needs to convert their observations into the appropriate format for LoopStructural.


Automatic intepretation into a model using map2loop
----------------------------------------------------
 



