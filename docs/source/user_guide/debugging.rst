Debugging LoopStructural
========================

Interpolation
-------------
Surfaces are not displaying
~~~~~~~~~~~~~~~~~~~~~~~~~~~
- check to see what the scalar field looks linked
- look at the values of the nodes of the implicit function supports `feature.interpolator.c` if all of the nodes are 0 the interpolator has not solved.
- Does the scalar field show different colours? `view.add_scalar_field(feature)`
- If the scalar field shows the right geometry, are the surface values correct?


Surface geometry is erroneous
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Try using a different solver - sometimes :code:`pyamg` does not work as it tries to pick a pattern in the interpolation matrix indicating the mesh geometry.  
- If the solution doesn't fit your data points try adding more elements
- Try using only orientation constraints
- Try using only value constraints
- Change gradient norm constraints to gradient direction constraints (rename :code:`nx,ny,nz` to :code:`gx,gy,gz`)
- Try removing some of the data points
- Try changing to using interface constraints instead of value, to do this rename column :code:`val` to code:`interface` but only use :code:`FDI`. 
Alternatively use surfe with `method=interface`. Using either of these methods the value of the scalar field will change.
- The tolerance of the conjugate gradient solver is scaled to the size of the model assuming that the scalar field is interpolated to have a magnitude of the gradient
norm to be close to 1. The tolerance is saved in `feature.builder.build_arguments['tol']` but can be overwritten using the argument in code:`create_and_add` `tol=val`. 
Smaller values mean the solution will be a closer fit to the data + regularisation. 
- Try changing the relative weighting of the data points and regularisation.
 * :code:`regularisation` - changes the weighting of the constant gradient or the finite difference second derivative. 1.0 is equivalent to 0.1 constant gradient.
 * :code:`npw=1.0` changes the weighting of the gradient norm constraints
 * :code:`gpw=1.0` changes the weighting of gradient direction constraints
 * :code:`cpw=1.0` changes the weighting of the control points.

- A good strategy for debugging a model is to lower the weights of the data points until a over smoothed solution is found. Then increase the weights until the model starts to overfit.
- 

Interpolator hasn't solved
~~~~~~~~~~~~~~~~~~~~~~~~~~
- check that there are sufficient data points to constrain a surface. There needs to be at least two unique value constraints or a single value constraint and a gradient norm constraint.
- Try using the default solver, :code:`solver="cg"` 

Interpolator is taking a long time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Check how many elements you are using, :code:`1e5` should take under 5 minutes on a modern laptop
- How many faults are associated with the feature? :code:`len(feature.faults)`, you can remove the faults by setting this to an empty list

.. code-block:
    faults = feature.faults 
    feature.faults = []
    feature.update()

 

Folds
-----
When using the :code:`create_and_add_folded_foliation` or :code:`create_and_add_folded_fold_frame` LoopStructural the quality of the results can be dependent on 
the correct parameterisation.

S-Plot doesn't fit the observations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Is the fold axis correct? If you expect a constant fold axis check that the calculated fold axis is consistent with your interpretation.
- Is the wavelength obvious in the s-variogram? If it is, check the guessed wavelength is reasonable Check the s 
- Try specifying the wavelength using :code:`limb_wl` or :code:`axis_wl` 
- If the data points in the S-Plot do not look correct, check the data points in the model by plotting the axial foliation and the data. The fold geometry
should be visible.
- If the fold is non-cylindrical, try building a cylindrical fold by removing some data
- Check the polarity of the datapoints, the fold modelling requires the polarity of the vectors is constant. Check uising vector plots :code:`view.add_data(folded_foliation,vectors=True)`

Faults
------

No displacement
~~~~~~~~~~~~~~~
- Check that the displacement is large enough that it should be visible. E.g. compare displacement to model bounding box
- Has the fault frame interpolated correctly? Coordinates 0 and 1 should both have values.

.. code-block:
    view.add_scalar_field(fault[0])
    view.add_scalar_field(fault[1])
    view.add_scalar_field(fault[2])

- Increase the resolution of the visualisation mesh :code:`view.nelements=1e6`
- Is the fault parallel to the feature being faulted?
- Has the fault been added to the feature being faulted?

.. code-block:
    print([f.name for f in faulted_feature.faults])

- Is the fault displacement vector correct? Add the vector field to the visualisation

.. code-block:
    view = LavaVuModelViewer(model)
    view.add_vector_field(model['fault'][1],locations=model.regular_grid()[::200]) #random 200 locations
    view.add_isosurface(model['fault'][0],0)
    view.interactive()
