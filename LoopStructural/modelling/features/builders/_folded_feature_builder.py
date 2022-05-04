from LoopStructural.modelling.features.builders import GeologicalFeatureBuilder
from LoopStructural.modelling.features.fold import FoldRotationAngle
import numpy as np

from LoopStructural.utils import getLogger, InterpolatorError

logger = getLogger(__name__)


class FoldedFeatureBuilder(GeologicalFeatureBuilder):
    def __init__(
        self, interpolator, fold, fold_weights={}, name="Feature", region=None, **kwargs
    ):
        """Builder for creating a geological feature using fold constraints

        Parameters
        ----------
        interpolator : GeologicalInterpolator
            the interpolator to add the fold constraints to
        fold : FoldEvent
            a fold event object that contains the geometry of the fold
        fold_weights : dict, optional
            interpolation weights for the fold, by default {}
        name : str, optional
            name of the geological feature, by default "Feature"
        region : _type_, optional
            _description_, by default None
        """
        GeologicalFeatureBuilder.__init__(
            self, interpolator, name=name, region=region, **kwargs
        )
        self.fold = fold
        self.fold_weights = fold_weights
        self.kwargs = kwargs
        self.svario = True

    def set_fold_axis(self):
        """calculates the fold axis/ fold axis rotation and adds this to the fold"""
        kwargs = self.kwargs
        fold_axis = kwargs.get("fold_axis", None)
        if fold_axis is not None:
            fold_axis = np.array(fold_axis)
            if len(fold_axis.shape) == 1:
                self.fold.fold_axis = fold_axis

        if "av_fold_axis" in kwargs:
            l2 = self.fold.foldframe.calculate_intersection_lineation(self)
            self.fold.fold_axis = np.mean(l2, axis=0)
        if self.fold.fold_axis is None:
            if self.fold.foldframe[1].is_valid() == False:
                raise InterpolatorError("Fold frame direction coordinate is not valid")
            far, fad = self.fold.foldframe.calculate_fold_axis_rotation(self)
            fold_axis_rotation = FoldRotationAngle(far, fad, svario=self.svario)
            a_wl = kwargs.get("axis_wl", None)
            if "axis_function" in kwargs:
                # allow predefined function to be used
                fold_axis_rotation.set_function(kwargs["axis_function"])
            else:
                fold_axis_rotation.fit_fourier_series(wl=a_wl)
            self.fold.fold_axis_rotation = fold_axis_rotation

    def set_fold_limb_rotation(self):
        """Calculates the limb rotation of the fold and adds it to the fold object"""
        kwargs = self.kwargs
        # give option of passing own fold limb rotation function
        flr, fld = self.fold.foldframe.calculate_fold_limb_rotation(
            self, self.fold.get_fold_axis_orientation
        )
        fold_limb_rotation = FoldRotationAngle(flr, fld, svario=self.svario)
        l_wl = kwargs.get("limb_wl", None)
        if "limb_function" in kwargs:
            # allow for predefined functions to be used
            fold_limb_rotation.set_function(kwargs["limb_function"])
        else:

            fold_limb_rotation.fit_fourier_series(wl=l_wl, **kwargs)
        self.fold.fold_limb_rotation = fold_limb_rotation

    def build(self, data_region=None, constrained=None, **kwargs):
        """the main function to run the interpolation and set up the parameters

        Parameters
        ----------
        data_region : [type], optional
            [description], by default None
        """
        # add the data to the interpolator and force constraints to be
        # gradient not norm, to prevent issues with fold norm constraint
        # TODO folding norm constraint should be minimising the difference in norm
        # not setting the norm

        # Use norm constraints if the fold normalisation weight is 0.
        if constrained is None:
            if "fold_normalisation" in kwargs:
                if kwargs["fold_normalisation"] == 0.0:
                    constrained = False
                else:
                    constrained = True
        self.add_data_to_interpolator(constrained=constrained)
        if self.fold.foldframe[0].is_valid() == False:
            raise InterpolatorError("Fold frame main coordinate is not valid")
        self.set_fold_axis()
        self.set_fold_limb_rotation()
        logger.info("Adding fold to {}".format(self.name))
        self.interpolator.fold = self.fold
        # if we have fold weights use those, otherwise just use default
        self.interpolator.add_fold_constraints(**self.fold_weights)
        if "cgw" not in kwargs:
            # try adding very small cg
            kwargs["cgw"] = 0.0
        # now the fold is set up run the standard interpolation
        GeologicalFeatureBuilder.build(self, data_region=data_region, **kwargs)
