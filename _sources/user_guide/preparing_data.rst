From geological map to model
----------------------------
LoopStructural is an implicit geological modelling library with different interpolation algorithms (discrete and data supported).
A geological model consists of multiple implicit functions that each define a distance or pseudo-distance within the geological object.
The interaction between the implicit function is defined by the relative timing of the feature and the type of geological feature being modelled.
Each implicit function is approximated to fit the observations of the geological feature. 
The geological observations need to be converted into mathematical constraints for the function defining the value of the implicit function and the gradient of the implicit function.

This guide details the `modelling.input` module which provides the tools for converting from a geological dataset into a 3D modelling.

Input dataset
-------------

To build a geological model of the location and distribution of stratigraphic units the input dataset must include:  
1. Location of the interface between stratigraphic units
2. Orientation of stratigraphy
3. Order of stratigraphic units
4. Stratigraphic thicknesses
5. Location of fault traces
6. Orientation of fault
7. Fault properties
8. Fault intersections
9. Fault stratigraphy intersections
10. Visualisation colours 


Requirements for modelling stratigraphy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

contacts
********

To build a model showing the distribution and location of stratigraphy the location of the basal interface between stratigraphic layers needs to recorded. 
This can be included as a DataFrame like object with the following columns:
* `X` eastings of the contact location
* `Y` northings of the contact Location
* `Z` altitude of the contact location
* `formation` a unique name defining the interface

orientations
************

Structural observations can be included using a DataFrame like object with the following columns:
* `X` eastings of the contact location
* `Y` northings of the contact Location
* `Z` altitude of the contact location
* `strike` using the right hand thumb rule
* `dip` angle between horizontal and contact

stratigraphic_order
*******************

The stratigraphic order defines the relative age of the stratigraphic interfaces. 
This can be included using a nested list containing the formation names. 
The lists are to be formated so that the first level of the list contains the youngest conformable layers.
The next level contains the order of the formations within the conformable group, with them order youngest to oldest.
In the example below, there are two groups the youngest group contains a,b,c and the older group contains d,e.

.. code-block::
    [[a,b,c],['d','e]]

Note: If the stratigraphy can be interpreted as conformable, it is advisable to include all stratigraphy inside a single group. 

thicknesses
***********


I
should be incolumn 

ProcessInputData
--------------------

The minimum requirements for building a geological model are:
* contact location observations
* contact orientation observations
* stratigraphic order
