LoopStructural design
=====================

**Authors:** Lachlan Grose, Laurent Ailleres, Gautier Laurent, Mark Jessell, Roy Thomson
 
**Last updated:** 1/11/2021

Overview
~~~~~~~~

LoopStructural is a modular 3D geological modelling library with a flexible API for building 3D geological models. 


Context
~~~~~~~

The objective of geological modelling is to predict the geometry and distribution of geological objects in the subsurface. 
LoopStructural provides a Explain implicit modelling

Goals and Non-Goals
~~~~~~~~~~~~~~~~~~~

The aim of LoopStructural are:
* to provide an easy to use python API for building 3D geological models from data.
* a platform for developing 3D modelling ideas
*  

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

* **Encapsulation**: Attributes in python are designed to be accessible both within the class and from outside of the class. In LoopStructural, where an attribute of a class should be private the ``@property`` decorator should be used.
* **Inheritance**: Classes should inheret attributes and methods from parent classes. 
* **Polymorphism**: Child classes, where a method is different from a parent class should override methods from parent classes.

**Dependency management** The library is written in python and is dependent on numpy, scipy.

**Documentation** The documentation is written for sphinx and should be formatted using the numpy style. 

**Tests** The tests are written using pytest and should be run using the command ``pytest``. 
There should be unit tests for all classes and methods.
Integration tests should be included to test the interaction between classes.

**CI/CD** The code is automatically versioned using release-please. 
When creating a commit use the conventional commit style e.g. ``fix: bug fix``, ``feat: new feature``, ``chore: code change``.
A new patch pull request will be created when a bug fix is applied and passes the tests. 
A new feature pull request will be created when a new feature is added and passes the tests.
When the new release is merged the library will be deployed to pypi, loop3d conda channel and dockerhub. 
The documentation is automatically built when a new release is created and will be uploaded to `loop3d.github.io/LoopStructural`_.


1. interpolators
----------------
The interpolation module contains the tools necessary for approximating a function given the constraints. 
LoopStructural uses a generic design for the interpolation algorithm implemented in the base class **GeologicalInterpolator**.
**GeologicalInterpolator** implements all of the functions required to interact with the interpolator (adding data, solving, evaluating).

Different interpolation algorithms can be implemented by inhertiting from the **GeologicalInterpolator**.

2. modelling
------------

The modelling module contains the main 

3. visualisation
----------------

4. utils
---------


Versioning
-----------
LoopStructural uses schemantic versioning to 



Continuous integration
----------------------

LoopStructural uses github actions to automatically run tests and checks for every 
Unit tests
**********

Integration tests
******************

Linting
*******

