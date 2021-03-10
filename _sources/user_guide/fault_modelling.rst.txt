Fault modelling
===============

Fault can be added into LoopStructural using the :meth:`LoopStructural.GeologicalModel.create_and_add_fault` function.

.. code-block::
    
    model.create_and_add_fault(feature_name, 
                               displacement,
                               fault_slip_vector=None,
                               fault_center = None, 
                               fault_extent = None, 
                               fault_influence = None, 
                               fault_vectical_radius = None, 
                               faultfunction = None,
                               **kwargs
                               )

In version 1.1 the additional arguments for the fault geometry were added. 
Additional parameters can be specified to the interpolator including:

        
  .. list-table:: Fault parameters
      :widths: 25 75
      :header-rows: 1

      * - Parameter
        - Options
      * - feature_name
        - identifying string in the model dataframe
      * - displacement
        - magnitude of the displacement, or maximum displacement, 
      * - fault_slip_vector
        - numpy array describing the fault slip
      * - fault_center
        - the center of the fault 
      * - fault_extent 
        - the radius of the fault ellipsoid parallel to fault surface and orthogonal to slip direction
      * - fault_influence
        - radius of ellipsoid orthogonal to fault
      * - fault_vertical_radius
        - fault ellipsoid radius in slip direction
      * - faultfunction
        - function describing the fault slip in the fault frame coordinates


        
Fault function
--------------

The fault function describes the displacement magnitude of the fault relative to the fault frame coordinate.
A basic fault function is provided in LoopStructural. 