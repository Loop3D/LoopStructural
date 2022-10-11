LoopStructural Philosophy
==========================
Loopstructural is an open source 3D geological modelling library with minimal dependencies - *numpy, scipy and scikit-learn for model building and evaluation*. 
For visualisation *LavaVu and scikit-image* are required. 
 
The main entry point for using LoopStructural is to create a GeologicalModel. 
The GeologicalModel manages the dataset, model area, creating of different elements inside the GeologicalModel.
Within the GeologicalModel, different elements in the model are represented by GeologicalFeatures. 

GeologicalModel 
~~~~~~~~~~~~~~~
* Manages the order of the geological features
* Creates geological features, addings 
* Rescales data/region
* Relates geological features - e.g. faults are added to stratigraphy

GeologicalFeature
~~~~~~~~~~~~~~~~~~
A base GeologicalFeature is something that can be evaluated anywhere in the model space and return the value of the scalar field
and the orientation of the gradient of the scalar field.
* Interprets scalar field for geology 

GeologicalFeatureBuilder
~~~~~~~~~~~~~~~~~~~~~~~~

StructuralFrame
~~~~~~~~~~~~~~~
* Curvilinear coordinate system made up of 3 geological features that relate to the same feature being modelled
* 

StructuralFrameBuilder
~~~~~~~~~~~~~~~~~~~~~~

GeologicalInterpolator
~~~~~~~~~~~~~~~~~~~~~~


DiscreteInterpolator
~~~~~~~~~~~~~~~~~~~~


PiecewiseLinearInterpolator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


FiniteDifferenceInterpolator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


