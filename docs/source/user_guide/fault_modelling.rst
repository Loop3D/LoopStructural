Fault modelling
===============
LoopStructural incorporates the fault displacement into the implicit representation of the faulted feature. 
The fault displacement is defined by the product of the fault slip vector and the displacement magnitude. 
A StructuralFrame is constructed for a fault where the three coordinates represent distance to the fault surface, distance in the fault slip direction and distance along the fault extent. 
The fault is integrated into the implicit model by restoring the observations of the faulted surface prior to interpolating the faulted surface. 

To model faults in LoopStructural the fault slip direction and fault displacement magnitude are necessary input into the implicit model.
Faults can either be incorporated where the assumption is the fault displacement is constant on the hanging wall of the fault, or the displacement decreases away from the fault centre. 


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


If only the required parameters are specified (`feature_name`, `displacement`) the fault is added using a constant displacement.
The displacement is incorporated by multiplying the normalised gradient of the fault slip direction field by the displacement magnitude on the hanging wall of the fault.
The hanging wall of the fault is identified by the positive component of the fault surface field. 

If `faultfunction` is used the displacement used is provided by faultfunction(coord_0,coord_1,coord_2) for every location in the model.
A generic faultfunction can be used by passing `faultfunction=`BaseFault`.
