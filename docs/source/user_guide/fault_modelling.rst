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
                               fault_major_axis = None, 
                               fault_minor_axis = None, 
                               fault_intermediate_axis = None, 
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

How fault modelling works
=========================
Faults are modelled in LoopStructural by applying a kinematic operator to the geological feature that is being modelled. 
The kinematic operator applies the kinematics of the fault in reverse to effectively restore observations or locations to their position prior to faulting. 
For example to interpolated a folded surface we first construct a fault frame using the observations of the fault surface and slip direction. 
The displacement of the fault can be constrained using a function of the fault coordinates.

.. image-sg:: /images/fault_frame_figure.png
   :alt: figure showing fault frame 
   :srcset: /images/fault_frame_figure.png
   :class: sphx-glr-single-img

The fault can then be added to the older features that are faulted.
Before interpolating a geological feature the associated faults, stored in the :attr:`LoopStructural.modelling.features.BaseFeature.faults` attribute of t are applied to the data points constraining the interpolation.
The geological feature can be built using the restored data. 
When the feature is evaluated, the locations being evaluated are first past through the list of faults using the :meth:`LoopStructural.modelling.features.BaseFeature._apply_faults`.
The `apply_faults` method should be called whenever :meth:`LoopStructural.modelling.features.BaseFeature.evaluate_value` or :meth:`LoopStructural.modelling.features.BaseFeature.evaluate_gradient` are overriden.
