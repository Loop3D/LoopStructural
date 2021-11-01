LoopStructural design
=====================

**Authors: ** Lachlan Grose, Laurent Ailleres, Gautier Laurent, Mark Jessell, Roy Thomson. 
**Last updated:** 1/11/2021

Overview
~~~~~~~~

LoopStructural is a modular 3D geological modelling library with a flexible API for building 3D geological models. 


Context
~~~~~~~

Explain implicit modelling

Goals and Non-Goals
~~~~~~~~~~~~~~~~~~~

The aim of LoopStructural are:
* to provide an easy to use python API for building 3D geological models from data.
* a platform for developing 3D modelling ideas
*  

Existing Solution
~~~~~~~~~~~~~~~~~

Technical Architecture 
~~~~~~~~~~~~~~~~~~~~~~

LoopStructural is divided into submodules which contain specific elements of the code. 
The submodules separate clear differences in the functionality of the code.
Where possible, object oriented coding is used to define 
The main modules for the library are:

1. interpolators 
2. modelling 
3. visualisation
4. utils

1. interpolators
----------------
LoopStructural uses a generic design for the interpolation algorithm implemented in the base class **GeologicalInterpolator**.
**GeologicalInterpolator** implements all of the functions required to interact with the interpolator (adding data, solving, evaluating).

Different interpolation algorithms can be implemented by inhertiting from the **GeologicalInterpolator**.

2. modelling
------------

The modelling module contains the main 

3. visualisation
----------------

4. utils


Versioning
-----------
LoopStructural uses 



Continuous integration
----------------------

LoopStructural uses github actions to automatically run tests and checks for every 
Unit tests
**********

Integration tests
******************

Linting
*******

