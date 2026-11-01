Using map2loop output
=====================

`map2loop <https://github.com/Loop3D/map2loop>`_ extracts the observations needed for a 3D model
(stratigraphic contacts, orientations, fault traces and their relationships) from a geological map.
LoopStructural turns that output into a model in two steps: a *processor* gathers the observations
and the order of the geological features, and the model is built from the processor.

.. code-block:: python

    from LoopStructural import GeologicalModel
    from LoopStructural.modelling.input import Map2LoopProcessor

    processor = Map2LoopProcessor("path/to/map2loop/project")
    model = GeologicalModel.from_processor(processor)
    model.update()      # run the interpolation

:class:`~LoopStructural.modelling.input.Map2LoopProcessor` reads the directory layout written by
the original map2loop (the ``output`` and ``tmp`` folders). For other sources, or to control the
inputs yourself, build a :class:`~LoopStructural.modelling.input.ProcessInputData` from data frames
and use :meth:`GeologicalModel.from_processor` in the same way; see
:doc:`/_auto_examples/4_advanced/plot_1_model_from_geological_map` for a worked example.

``from_processor`` creates the model for the extent given by the processor (``origin`` to ``maximum``)
and attaches its data. It then adds the faults of the processor's fault network, including splay and
abutting relationships between faults, and adds each stratigraphic group as a foliation affected by the
faults listed for it in the fault-stratigraphy table. The base of each interpolated group is added as an
unconformity, and the processor's stratigraphic column is attached to the model. The model is returned
without having been solved; call ``model.update()``.
