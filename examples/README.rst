Examples
========
These examples show how to build, visualise and export implicit
geological models with LoopStructural, roughly in the order you would
learn them:

1. **Basics** - loading data, building a model, adding stratigraphy,
   unconformities and faults, and visualising and exporting the result.
2. **Modelling folds** - constraining folded surfaces with fold frames,
   including refolded (multiply-deformed) folds.
3. **Modelling faults** - fault networks, custom displacement profiles,
   and updating fault geometry after a model has been built.
4. **Advanced use** - building models directly from geological map data,
   logging, controlling data weighting, and comparing interpolators.

Each example is a standalone, runnable Python script. Most examples load
one of the sample datasets bundled in :code:`LoopStructural.datasets`, so
no external data is required to follow along.

Visualisation in these examples uses :code:`Loop3DView` from the
`loopstructuralvisualisation <https://github.com/Loop3D/LoopStructuralVisualisation>`_
package, a PyVista-based 3D viewer - install it (and matplotlib, used for
2D plots) with :code:`pip install loopstructural[visualisation]`.
