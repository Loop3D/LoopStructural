.. LoopStructural documentation master file, created by
   sphinx-quickstart on Thu May 28 15:51:42 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Loop Structural
===============
.. image:: images/image823.png

Loop structural is an open source 3D modelling library providing access to numerous interpolation algorithms with a
high level and easy to use API for creating geological models. The library has been written for the Loop platform by
Lachlan Grose at Monash University.

Loop is an open source 3D probabilistic geological and geophysical modelling platform,
initiated by Geoscience Australia and the OneGeology consortium. The project is funded by Australian territory,
State and Federal Geological Surveys, the Australian Research Council and the MinEx Collaborative Research Centre.

LoopStructural is the forward geological modelling engine for the loop and contains the classes and functions to manage the various
elements in a 3D geological model. Including stratigraphy, unconformities, fault and folds. LoopStructural contains
three native interpolation algorithms:

1. Piecewise-linear interpolation Frank 2007
2. Discrete Fold Interpolator Laurent 2016
3. Finite Difference Interpolation

and a wrapper for the generalised radial basis functions provided by Surfe.

.. toctree::
   :caption: Getting Started
   :hidden:

   installation
   index
   API
   auto_examples/index


.. automodule:: LoopStructural

