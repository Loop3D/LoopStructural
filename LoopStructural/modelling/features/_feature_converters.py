from LoopStructural.modelling.features.fold import FoldEvent
from LoopStructural.modelling.features.builders import FoldedFeatureBuilder
def add_fold_to_feature(feature, fold_frame,**kwargs):
    fold = FoldEvent(fold_frame, name=f"Fold_{feature.name}", invert_norm=kwargs.get('invert_fold_norm', False))

    builder = FoldedFeatureBuilder.from_feature_builder(
        feature.builder,
        fold, 
        **kwargs
    )
    feature = builder.feature
    feature.fold = fold
    return feature
