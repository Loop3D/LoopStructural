.. LoopStructural documentation master file, created by
   sphinx-quickstart on Thu May 28 15:51:42 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Loop Structural
===============
.. image:: ./images/image823.png

.. important:: **Upgrading from 1.5.x to 1.6.x**

   Version 1.6.x of LoopStructural has been released and there are a number of changes to how the library is used. The main changes are:
   
   - The `LavaVuModelViewer` class has been removed and the model visualisation has been ported to `pyvista` using a LoopStructural 
     wrapped `loopstructuralvisualisation`. Many of the functions have been re-named and can be combined with the `pyvista.Plotter` methods for advanced
     visualisation. pyvista does not currently allow for interactive plots while using colab notebooks, so `jupyter_backend='static'` is needed.
   - Datastructures for representing model input/output have been introduced including the `LoopStructural.datatypes.StructuredGrid`,
     :class:`LoopStructural.datatypes.ValuePoints`, :class:`LoopStructural.datatypes.VectorPoints` and :class:`LoopStructural.datatypes.Surface`. These Datastructures
     are intended to improve the accessibility of the modelling and allow for easier export into other software packages.
     All data objects have a `.save` method that allows the user to save the points into various formats including `json`, `csv`, `vtk` and `geoh5`.
     In the future exporters for gocad and other formats will be added. All objects also have a `.vtk()` method to return a `pyvista.PolyData` object.
   - The solver options for the interpolators have been removed and only :meth:`scipy.sparse.linalg.cg` and :meth:`scipy.sparse.linalg.lsmr` are supported. 
     To use an external solver the `external_solver` argument can be provided to the `interpolator.solve` method, usually passed through
     the `model.create_and_add*` methods. To specify specific arguments to the solver these can be passed through `solver_kwargs` dict.
   - Improvements have been made to how unconformities are handled and their interaction with faults
   - A framework for including inequalities into the interpolation has been added but there is no solver that can be used for solving this problem.
   - A convenience class :class:`LoopStructural.LoopInterpolator` was added for building an implicit function without requiring the :class:`LoopStructural.GeologicalModel`
   **Please report any bugs or issues on the github repository.**
 
Overview
========

LoopStructural is an open source 3D modelling library providing access to multiple interpolation schemes  with a
high level and easy to use API for creating geological models. The library has been written for the Loop platform by
Lachlan Grose at Monash University.

Loop is an open source 3D probabilistic geological and geophysical modelling platform,
initiated by Geoscience Australia and the OneGeology consortium. The project is funded by Australian territory,
State and Federal Geological Surveys, the Australian Research Council and the MinEx Collaborative Research Centre.



   
.. toctree::
   :hidden:

   getting_started/index
   _auto_examples/index
   user_guide/index

   
   .. toctree::
   :caption: LoopStructural API
   :hidden:

   API


