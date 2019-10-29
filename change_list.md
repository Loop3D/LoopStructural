* Dockerfile updated to use requirements.txt and no longer mounts notebooks directory as editable
* GeologicalModel class implemented where multiple different geological features can be built from a high level api. Currently supports unconformities, faults and stratigraphy NOT folds
* Modified base discrete interpolator so that solvers are separate functions and building the interpolation matrix is performed from a separate function. Inequaltiies are added for all interpolators and external solvers can be used by passing a function(A,B).
* Changed gradient constraints to be true gradients not norm and magnitude. 
* NPoints are vectors that constrain the normal to the scalar field
* Vector data is currently automatically added as normal constraints if the system is not constrained by value constraints. Gradient otherwise. For folds there is a hack where the data can be added as gradeint not normal by running add_data_to_interpolator(constrained=True) because the normal constraint from the fold constraints is enough to constrain the system - and can conflict with normal constraints.
## StructuralFrameBuilder
* StructuralFrameBuilder now containts three GeologicalFeatureInterpolators - instead of accessing the interpolators by `structural_frame_builder.interpolator[i]` you should now use `structural_frame_builder.builders[i].interpolator`
* All calculations using the fold frame now take a GeologicalFeatureInterpolator as arguments rather than the numpy arrays. This allows for both gradient and normal points to be included. The fold limb and fold axis rotation angle functions now also return the fold frame values.
* Using itype='gx' etc is depreciated and will be removed in next release. Use coord=0 instead.
## Faults 
* Faults are now added to a GeologicalFeature using `.add_fault(FaultSegment)` rather than creating a `FaultedGeologicalFeature` this means you won't have to chain together multiple FaultedGeologicalFeatures.
## GeologicalFeatures
* Moving towards having lambda functions as regions rather than precomputing the region on the support. This means that regions can be calculated indpendent from the support where the interpolation occurs.
* multiple regions can be added and the boolean addition of all regions will be used (np.logical_and)
* Evaluate_value and evaluate_gradient now check whether the locations are within the regions 

 
## Visualisation
* Slice now colours isosurfaces by different colours rather than applying colourmap to the whole surface!
* Can add a section to the model now using add_section
