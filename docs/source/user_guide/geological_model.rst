Geological Model
================
The GeologicalModel is the main entry point for creating a LoopStructural model. 
GeologicalModel's define the bounding box of the model, determine how the different scalar fields interact and apply the time aware modelling.

Creating a model manually
-------------------------
A GeologicalModel can be created using the default constructor.
.. code-block::
    model = GeologicalModel(origin,maximum)

Additional arguments that can be provided are:
  .. list-table:: Finite Difference Interpolator (FDI)
      :widths: 25 75
      :header-rows: 1
      
      * - Argument 
        - Description
      * - rescale
        - Boolean (True), whether to rescale the model so between 0,1
      * - nsteps
        - 
      * - reuse_supports
        - Boolean (False), whether to reuse the inteprolation supports for similar objects
      * - logfile
        - String (None), file to store the loopstructural log to
      * - loglevel
        - logging level to use ('info'), 'warn','error','debug'

A dataset then needs to be added to the model, the data set is a pandas dataframe.
Within this method the data will be rescaled to the model coordinates.
The dataframe will be checked to make sure the correct columns exists.
The model creates a copy of the dataframe meaning that the original data frame is not modified.

.. code-block::

    model.set_model_data(data)

Adding a geological feature
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Once the model has been initialised and the data has been linked to the model, 
we can add a geological feature to the model.
There are a number of different methods that can be used to add a geological feature
to the model.

Adding a foliation
~~~~~~~~~~~~~~~~~~
We can add a foliation to the model using the *create_and_add_foliation* method for a GeologicalModel object.
The only required argument is the feature_name, which is the identifier to the datasets.



.. code-block::

    model.create_and_add_foliation(feature_name,
                                    interpolatortype = 'PLI',
                                    nelements = 1e5,
                                    buffer = 0.2,
                                    solver = 'cg',
                                    damp = True)
                                    
The other parameters are used to tune the behaviour of LoopStructural. 
Different interpolation schemes can be used, interpolation_options_ for more information 

Adding a fault
~~~~~~~~~~~~~~
Faults are modelled by building a structural frame to define the fault geometry.
The structural frame has three coordinates:

* Coordinate 0: Fault surface 
* Coordinate 1: Fault slip direction
* Coordinate 2: Fault extent

Observations can be linked to particular features using the *coord* column in the data frame.

When the fault is modelled the fault surface is defined by the 0 isosurface of the scalar field.
The positive values of the scalar field represent the hanging wall of the fault and the negative values the footwall.

The displacement of the fault can be provided as either a constant displacement for the hanging wall, or a function of 
the fault frame coordinates scaled by the maximum displacement.


The fault can be modelled with either a constant displacement on the hanging wall. 

* displacement magnitude or displacement function
* fault surface is isovalue of 0 for scalar field
* fault is incorporated by restoring the model area and data points
* fault needs to be given slip direction
* fault extent is estimated to be the map expression of the fault if no data is given  

.. code-block::

    model.create_and_add_fault(feature_name, displacement)

Adding an uncomformity
~~~~~~~~~~~~~~~~~~~~~~
There are two ways an unconformity can be added into the model. 
The first approach creates an unconformity as a new geological feature assigning the unconformity surface to the 0 isovalue.
Using this approach does not ensure that the unconformity is conformable to a stratigraphic series. 

.. code-block::

    model.create_and_add_unconformity(feature_name)

The second approach uses an existing GeologicalFeature and definies a specific value of this scalar field to represent the unconformity.

.. code-block::

    model.add_unconformity(feature, value)

The unconformity acts as a mask clipping any older stratigraphic units by the geometry of the unconformity.

Onlap unconformities can be added in the same way using the 

.. code-block::

    model.add_onlap_unconformity(feature, value)

This will clip any younger stratigraphic units to be clipped by the value of the unconformity.

Setting a stratigraphic column
------------------------------
The GeologicalModel can have a stratigraphic column assigned to it. 
The stratigraphic column performs a map between the GeologicalFeatures representing stratigraphy and a unique identification for each stratigraphic unit.
The stratigraphic column can also be used to specify colouring for faults.

The stratigraphic column is a dictionary and has the form:

.. code-block::

    stratrigraphic column = {'group': #must be the GeologicalFeature name
                {'series1': # name that will appear on legend
                            {'min':0., 'max':10.,'id':0,'colour':}
                }
        }
    model.set_stratigraphic_column(stratigraphic_column)

When used to evaluate the block model or to evaluate the model surfaces the name of the first entry in the dictionary is checked for whether it contains fault.
If it is a fault, this scalar field is not used for evaluating the stratigraphy.


Using the GeologicalModel
-------------------------
The GeologicalModel has a number of helper functions allowing you to easily access different 
aspects of the model.

A regular grid inside the model bounding box can be retrieved:

.. code-block::

    regular_grid = model.regular_grid(nsteps=(50,50,25),shuffle=True,rescale=False)

The parameters nsteps define how many points in x, y and z. 
Shuffle defines whether the points should be ordered by axis or random (note this is useful when visualising vector data)
as regular sampling becomes obvious when slicing numpy arrays.
The rescale parameter defines whether the returned points should be in model coordinates or real world coordinates, it is 
in model coordinates by default.

A GeologicalFeature can be extracted from the model either by name

.. code-block::

    myfeature = model.create_and_add_feature('myfeature')
    myfeature_also = model['myfeature']
    myfeature_by_index = model.features[0]
    myfeature_by_name = model.get_feature_by_name('myfeature')

Evaluating a scalar field
~~~~~~~~~~~~~~~~~~~~~~~~~

The scalar field for a GeologicalFeature can be evaluated using the evaluate_feature_value method:

.. code-block::

    regular_grid = model.regular_grid(shuffle=False,rescale=True)
    sf = model.evaluate_feature_value('myfeature',regular_grid,scale=True)

The scale parameter determines whether the points are in model coordinates or real world. 
In this example, the points were returned in real world and therefore need scaling before evaluating in 
the model.
A numpy array of N dimensions will be returned where N is the number of points.
If the value cannot be evaulated it will be np.nan.

Evaluating scalar field gradient
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The gradient scalar field for a GeologicalFeature can be evaluated using the evaluate_feature_gradient method:

.. code-block::

    regular_grid = model.regular_grid(shuffle=False,rescale=False)
    vf = model.evaluate_feature_gradient('myfeature',regular_grid,scale=False)

Evaluating fault displacement magnitude
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fault displacement magnitude can be evaulated on a pointset. 
This will evaluate the fault function for every fault in the model and can be used
for analysing models.

.. code-block::

    regular_grid = model.regular_grid(shuffle=False,rescale=False)
    fault_displacement = model.evaluate_fault_displacements('myfeature',regular_grid,scale=False)

Evaluating lithology id
~~~~~~~~~~~~~~~~~~~~~~~~~

The model can be evaluated for the lithology id defined in the stratigraphic column.

.. code-block::

    regular_grid = model.regular_grid(shuffle=False,rescale=False)
    litho = model.evaluate_model('myfeature',regular_grid,scale=False)
