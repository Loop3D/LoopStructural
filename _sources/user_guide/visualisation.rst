Visualisation
=================
LoopStructural uses the powerful lavaVu visualisaton library for visualising the 3D models as surfaces and scalar fields.
LavaVu provides an interative viewer inside jupyter notebooks and can use matplotlib colourmaps as well as export the model to webgl.
The LavaVuModelViewer provides an easy to use interfaces for adding GeologicalFeatures to the viewer, the base lavaVu object can still be used.
Some of the functions in the viewer require a bounding box to be able to draw isosurfaces or scalar fields.
By default the viewer is initialised with a GeologicalModel and the origin and maximum of the model are used to set the viewers extents.
It is also possible to modify this manually by changing the bounding_box attribute.

Visualising GeologicalFeatures
------------------------------

There are four main ways a GeologicalFeature can be visualised.
1. Scalar field showing the value of the implicit function on a block representing the model domain

.. code-block::

    view.plot_scalar_field(geological_feature)

2. Surfaces extracted from isovalues of the implicit function

.. code-block::

    view.plot_surface(geological_feature,value=0)

3. Vector field representing the gradient of the implicit function

.. code-block::

    view.plot_surface(geological_feature)

4. Cross section showing the value of the implicit function at a particular location inside the model

.. code-block::

    view.add_section(geological_feature,axis='x',value=0)

When adding to the viewer there are a number of common arguments that can be specified to change the behaviour:

  .. list-table:: Viewer kwargs
      :widths: 25 75
      :header-rows: 1

      * - Parameter
        - Options
      * - vmin
        - minimum value for colourmap
      * - vmax
        - maximum value for colourmap
      * - cmap
        - matplotlib colourmap or string with name of matplotlib colourmap e.g. 'rainbow'
      * - pyvista_kwargs
        - any arguments to be passed to the :meth:`pyvista.Plotter.add_mesh` method
      
The default behaviour of the viewer can also be modified.
To change the resolution of the scalar field or isosurfacing, either set the number of steps in x,y,z or the total number of elements of tht bounding box

.. code-block::
    
    model.bounding_box.nsteps = np.array([nx,ny,nz])
    print(model.bounding_box.nsteps)
    model.bounding_box.nelements = 1e6
    print(model.bounding_box.nsteps)


    
Visualising data
----------------

The input data for the geological model can also be visualised within the viewer. 
To add all of the data attached to a GeologicalFeature

.. code-block::

    view.plot_data(geological_feature,value=True,vector=True,scale=10,pyvista_kwargs={})

Alternatively, points can be added into the model using the:

.. code-block::

    view.add_points(xyz)

or points with an associated value:

.. code-block::
    view.add_points(xyz,scalars=np.array(values))

or points with an associated vector:

.. code-block::

    view.add_arrows(xyz,direction)

**Warning, you cannot have two objects in the viewer with the same name**

Visualising geological model
----------------------------

Where the GeologicalModel has an assocaited stratigraphic column including a list of faults in the model, this can be used for automatically visualising the model.
To add a block model showing the stratigraphic unit distribution, where the value is linked to the *id* in the stratigraphic column.

.. code-block::

    view.plot_block_model()

alternatively, the surfaces can be added to the visualisation.

.. code-block::

    view.plot_model_surfaces(faults=True,strati=True)

Where the faults flag allows you to disable visualising the fault surface. 
This function will add a surface for every stratigraphic unit in the stratigraphic column.

Custom visualisation
--------------------

It is also possible to pass custom functions to the viewer using the LambdaGeologicalFeature.
The LambdaGeologicalFeature is a class allowing for a function describing the value at locations or vector field.

.. code-block::

    from LoopStructural.modelling.features import LambdaGeologicalFeature
    def x_function(xyz):
        return xyz[:,0]
    custom_feature = LambdaGeologicalFeature(x_function,name='x feature')
    view.plot_surface(custom_feature,value)
    view.plot_scalar_field(custom_feature)

The function passed to the LambdaGeologicalFeature can be as simple or complicated as required. 
It will be evaluated for the locations within the model that the visualisation requires, usually between the origin and maximum of the geoloical model.
By default the min and max values are 0, however these can be overwritten by setting the attribute min and max to the required values.




