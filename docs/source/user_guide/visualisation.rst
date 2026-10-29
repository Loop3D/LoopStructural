Visualisation
=============

LoopStructural visualises models with :class:`LoopStructural.visualisation.Loop3DView`, a thin
subclass of the `PyVista <https://docs.pyvista.org>`_ ``Plotter`` that knows how to draw
geological features, faults, data and whole models. It comes from the separate
``loopstructuralvisualisation`` package; install it with

.. code-block:: bash

    pip install loopstructural[visualisation]

Because ``Loop3DView`` *is* a PyVista plotter, everything PyVista can do (``add_mesh``,
``add_points``, ``add_arrows``, ``camera_position``, ``screenshot``, ...) is available as well.

.. code-block:: python

    from LoopStructural.visualisation import Loop3DView

    view = Loop3DView(model)       # the model is used for bounds and defaults
    view.plot_surface(feature, value=0)
    view.display()                 # opens the viewer (inline in a notebook)

Some methods need a bounding box to know where to evaluate the feature. By default this is
the bounding box of the model passed to ``Loop3DView``; ``plot_surface``,
``plot_scalar_field``, ``plot_vector_field``, ``plot_fault`` and ``plot_fault_displacements``
also accept a ``bounding_box=`` argument to override it.

Visualising geological features
-------------------------------

There are four main ways to look at a geological feature.

1. A scalar field: the value of the implicit function on a volume.

   .. code-block:: python

       view.plot_scalar_field(feature, cmap="viridis")

2. Surfaces: isosurfaces of the implicit function. ``value`` is one isovalue, a list of
   isovalues, or a positive integer number of regularly spaced surfaces between the minimum and
   maximum of the feature (``0`` is read as the isovalue 0.0, not a count). If ``value`` is
   omitted the surface at the middle of the feature's range is drawn.

   .. code-block:: python

       view.plot_surface(feature, value=[0, 100, 300], cmap="rainbow")

3. A vector field: the gradient of the implicit function, drawn as arrows.

   .. code-block:: python

       view.plot_vector_field(feature, scale=50)

4. Sliders: pass ``slicer=True`` to ``plot_scalar_field`` or ``plot_surface`` to add an
   interactive slicing widget.

The methods share a few arguments:

.. list-table:: Common arguments
    :widths: 25 75
    :header-rows: 1

    * - Parameter
      - Meaning
    * - ``vmin``, ``vmax``
      - limits of the colour map
    * - ``cmap``
      - a matplotlib colour map, or its name e.g. ``"rainbow"``
    * - ``opacity``
      - transparency of the object, between 0 and 1
    * - ``name``
      - name of the object in the viewer. If the name is already used a numeric suffix is added
        (``name_1``, ``name_2``, ...) so nothing is replaced
    * - ``pyvista_kwargs``
      - a dictionary passed on to :meth:`pyvista.Plotter.add_mesh`
    * - ``bounding_box``
      - the region to evaluate the feature in, by default the model's bounding box

Resolution
~~~~~~~~~~

The scalar field and isosurfaces are evaluated on a regular grid over the bounding box.
To change its resolution either set the number of points in each direction or the approximate
number of cells:

.. code-block:: python

    model.bounding_box.nsteps = np.array([nx, ny, nz])
    print(model.bounding_box.nsteps)

    model.bounding_box.nelements = 1e6    # updates nsteps for roughly cubic cells
    print(model.bounding_box.nsteps)

Visualising data
----------------

The data attached to a feature can be drawn with the feature:

.. code-block:: python

    view.plot_data(feature, value=True, vector=True, scale=50)

``value`` draws the value constraints as points and ``vector`` draws the orientation
constraints as arrows. Arbitrary points and arrows can be added with the PyVista methods:

.. code-block:: python

    view.add_points(xyz, name="my points")
    view.add_points(xyz, scalars=values, name="my values")
    view.add_arrows(xyz, directions, name="my vectors")

.. warning::

    Pass a distinct ``name`` for each object. ``Loop3DView`` does not replace an existing
    object when a name is reused: the new object is given a numeric suffix (``name_1``,
    ``name_2``, ...) and both stay in the scene.

Visualising the geological model
--------------------------------

When a model has a stratigraphic column (and, optionally, faults), the whole model can be
drawn at once. A block model colours each cell by the ``id`` of its stratigraphic unit:

.. code-block:: python

    view.plot_block_model()

or draw a surface for each stratigraphic unit and, optionally, each fault:

.. code-block:: python

    view.plot_model_surfaces(strati=True, faults=True)

Faults can also be drawn individually, with the slip vector and displacement volume:

.. code-block:: python

    view.plot_fault(model["fault_name"], surface=True, slip_vector=True)

Custom visualisation
--------------------

A :class:`LoopStructural.modelling.features.LambdaGeologicalFeature` wraps any function of
location, so arbitrary quantities can be drawn with the same methods.

.. code-block:: python

    from LoopStructural.modelling.features import LambdaGeologicalFeature

    def x_function(xyz):
        return xyz[:, 0]

    custom_feature = LambdaGeologicalFeature(x_function, name="x feature")

    # a lambda feature has no model, so say where to evaluate it
    bb = model.bounding_box
    view.plot_surface(custom_feature, value=0.5 * (bb.origin[0] + bb.maximum[0]), bounding_box=bb)
    view.plot_scalar_field(custom_feature, bounding_box=bb)

The function is evaluated at the locations the visualisation needs, normally between the
origin and maximum of the bounding box. ``min()`` and ``max()`` on a feature are methods, not
settable attributes, and a feature without a model (such as a lambda feature) returns 0 for
both. To control the colour range pass ``vmin`` and ``vmax`` to the plotting method instead.
If you use the lambda feature as ``paint_with``, also pass ``vmin`` and ``vmax`` since the
default range comes from ``min()`` and ``max()``.

Using the grid directly
-----------------------

``feature.scalar_field()`` evaluates a feature on the bounding box's regular grid and returns a
:class:`LoopStructural.geometry.StructuredGrid`. Its ``properties`` hold the values at the nodes and
``cell_properties`` hold values at the cell centres. ``vtk()`` returns a PyVista
grid and ``save()`` writes it to VTK, GeoH5, OMF or GOCAD files.

.. code-block:: python

    grid = feature.scalar_field()
    grid.plot()                         # quick look with PyVista
    grid.save("scalar_field.vtk")
