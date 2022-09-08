LoopStructural design
=====================

**Authors:** Lachlan Grose, Laurent Ailleres, Gautier Laurent, Mark Jessell, Roy Thomson
 
**Last updated:** 1/11/2021

Overview
~~~~~~~~

LoopStructural is a modular 3D geological modelling library with a flexible API for building 3D geological models. It provides access to implicit interpolators as well as geological feature/events.


Context
~~~~~~~

The objective of geological modelling is to predict the geometry and distribution of geological objects in the subsurface. Different geological objects are represented by approximating an implicit function where isovalues/level sets of the implicit function define the geometry of geological features. The interaction between different geological objects is defined by the geological feature type (fault, fold, unconformity). 

Goals and Non-Goals
~~~~~~~~~~~~~~~~~~~

The aim of LoopStructural are:

* to provide an easy to use python API for building 3D geological models from data.
* a platform for developing 3D modelling ideas, remove the requirement for boiler plate code.
* integrate more geological observation types and knowledge into the modelling workflow


Existing Solutions
~~~~~~~~~~~~~~~~~~

There are a number of commercial and open source 3D modelling packages available. 

* Gocad/Skua provides a 3D modelling environment and an implicit interpolation algorithm
* GemPy is a python based open source 3D modelling library that implements co-kriging and event management.
* Geomodeller is a commerical package that implements the co-kriging algorithm and event management.
* Petrel is a commercial package that implements implicit modelling.
* Leapfrog is a commercial package that implements implicit modelling.
* StructuralLab is a plugin for Gocad that implements discrete implicit modelling.

Technical Architecture 
~~~~~~~~~~~~~~~~~~~~~~

**Module design** LoopStructural is divided into submodules which contain specific independent elements of the code. 
The submodules separate clear differences in the functionality of the code.
The main modules for the library are:

1. interpolators 
2. modelling 
3. visualisation
4. utils

**Coding paradigm** LoopStructural is written in an object oriented paradigm.
The four pillars of object oriented coding should be adhered to:

* **Encapsulation**: Attributes in python are designed to be accessible both within the class and from outside of the class. In LoopStructural, the where an attribute of a class should be private the ``@property`` decorator should be used.
* **Inheritance**: Classes should inheret attributes and methods from parent classes. 
* **Polymorphism**: Child classes, where a method is different from a parent class should override methods from parent classes.

**Dependency management** The library is written in python and is dependent on numpy, scipy. Care should be taken to avoid incorporating of uncessary dependencies. For dependencies that are not part of the standard scientific python toolkit (numfocus projects) the dependencies should be optional and not be required. 

**Documentation** The documentation is written for sphinx and should be formatted using the numpy style. 

**Tests** The tests are written using pytest and should be run using the command ``pytest``. 
There should be unit tests for all classes and methods.
Integration tests should be included to test the interaction between classes.

**CI/CD** The code is automatically versioned using release-please. 
When creating a commit use the conventional commit style e.g. ``fix: bug fix``, ``feat: new feature``, ``chore: code change``.
A new patch pull request will be created when a bug fix is applied and passes the tests. 
A new feature pull request will be created when a new feature is added and passes the tests.
When the new release is merged the library will be deployed to pypi, loop3d conda channel and dockerhub. 
The documentation is automatically built when a new release is created and will be uploaded.
The entire CI/CD workflow is run using github actions.


1. interpolators
----------------
The interpolation module contains different interpolators and associated supports used by LoopStructural.
The purpose of the interpolator is to build a 2D/3D scalar field from constraints (e.g. orientation, contact, fold).
The base of the interpolation module is the **GeologicalInterpolator** class.
The **GeologicalInterpolator** class contains the entry point for all interpolators and the functions for accessing constraints.
All implementations of interpolators are subclassed from the **GeologicalInterpolator**. 
This module should include all interpolators and associated classes (e.g. mesh, grid, operators).
Future update: there should be a separate API so that the interpolators can be used without a GeologicalModel. 


2. modelling
------------

The modelling is the main entry point to LoopStructural and includes the **GeologicalModel** class.
A **GeologicalModel** consists of **GeologicalFeatures** which define a particular geological phenomenon.
For example, a **GeologicalFeature** can be a stratigraphic sequence, a fault, a foliation, an unconformity, an intrusion.
The **GeologicalFeature** is something that can be evaluated within the model for a value, or gradient. 

In LoopStructural the **GeologicalFeature** is built by a **GeologicalFeatureBuilder**.
The **GeologicalFeatureBuilder** is a class that manages the assembly of the interpolation constraints.
The builder is in charge of keeping track of whether the interpolator has been run and whether the data/constraints have changed.
The builder could also be used to build a **GeologicalFeature** without an interpolator.


3. visualisation
----------------

The visualisation module is in charge of creating 2D/3D visualisations for the GeolgicalModel and GeologicalFeatures.

3D visualisation
****************

The 3D visualisation module for LoopStructural has been written using an abstract base class.
**BaseModelPlotter** contains all of the methods for plotting the various different attributes for a GeologicalModel.
A visualisation module can be implemented by simply subclassing the **BaseModelPlotter** and implementing the following methods:

* `_add_surface(self,vertices,face,name)` to add a triangular surface to the plot.
* `_add_points(self,points,name,value=None)` to add xyz points with an optional value attribute to colour by. 
* `_add_vector_marker(self, location, vector, name symbol_type=arrow)` to add a vector marker to the plot.


4. utils
---------

The utils module is a container for utility functions that are used throughout the library.





