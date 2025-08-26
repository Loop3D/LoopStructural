from LoopStructural.modelling.features.fold import FoldEvent, FoldFrame
from LoopStructural.modelling.features.builders import FoldedFeatureBuilder, StructuralFrameBuilder
def add_fold_to_feature(feature, fold_frame,**kwargs):
    if not isinstance(fold_frame, FoldFrame):
        raise ValueError("fold_frame must be a FoldFrame instance")
        
    fold = FoldEvent(fold_frame, name=f"Fold_{feature.name}", invert_norm=kwargs.get('invert_fold_norm', False))

    builder = FoldedFeatureBuilder.from_feature_builder(
        feature.builder,
        fold, 
        **kwargs
    )
    feature = builder.feature
    feature.fold = fold
    return feature

def convert_feature_to_structural_frame(feature, **kwargs):
    """
    Convert a geological feature to a structural frame by adding the feature to the frame

    Parameters
    ----------
    feature : GeologicalFeature
        the geological feature to convert
    
    Returns
    -------
    StructuralFrame
        the updated structural frame with the feature added
    """
    builder = feature.builder

    new_builder = StructuralFrameBuilder.from_feature_builder(
        builder,
        **kwargs
    )
    return new_builder.frame
    