from ....modelling.features.builders import GeologicalFeatureBuilder
from ....modelling.features.fold.fold_function import FoldRotationType, get_fold_rotation_profile
from ....modelling.features import FeatureType
import numpy as np

from ....utils import getLogger, InterpolatorError
from ....datatypes import BoundingBox

logger = getLogger(__name__)


class FoldedFeatureBuilder(GeologicalFeatureBuilder):
    def __init__(
        self,
        interpolatortype: str,
        bounding_box: BoundingBox,
        fold,
        nelements: int = 1000,
        fold_weights={},
        name="Feature",
        region=None,
        svario=True,
        axis_profile_type=FoldRotationType.FOURIER_SERIES,
        limb_profile_type=FoldRotationType.FOURIER_SERIES,
        **kwargs,
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
        # create the feature builder, this intialises the interpolator
        GeologicalFeatureBuilder.__init__(
            self,
            interpolatortype=interpolatortype,
            bounding_box=bounding_box,
            nelements=nelements,
            name=name,
            region=region,
            **kwargs,
        )
        self._feature.type = FeatureType.FOLDED
        # link the interpolator to the fold object
        self.interpolator.fold = fold
        self.fold = fold
        self.fold_weights = fold_weights
        self.update_build_arguments(kwargs)
        # self.kwargs = kwargs
        self.svario = svario
        self.axis_profile_type = axis_profile_type
        self.limb_profile_type = limb_profile_type
    @classmethod
    def from_feature_builder(cls, feature_builder, fold, **kwargs):
        """Create a FoldedFeatureBuilder from an existing feature builder"""
        if not isinstance(feature_builder, GeologicalFeatureBuilder):
            logger.error(f'Feature builder is {type(feature_builder)} not GeologicalFeatureBuilder')
            raise TypeError("feature_builder must be an instance of GeologicalFeatureBuilder")
        builder = cls(
            interpolatortype='DFI',
            bounding_box=feature_builder.model.bounding_box,
            fold=fold,
            nelements=feature_builder.interpolator.n_elements,
            name=feature_builder.name,
            **kwargs
        )
        builder.data = feature_builder.data
        return builder
    @property
    def fold_axis_rotation(self):
        if self.fold.fold_axis_rotation is None:
            self.set_fold_axis()
        return self.fold.fold_axis_rotation

    @property
    def fold_limb_rotation(self):
        _axis = self.fold.fold_axis  # get axis to make sure its initialised
        if self.fold.fold_limb_rotation is None:
            self.set_fold_limb_rotation()
        return self.fold.fold_limb_rotation

    def set_fold_axis(self):
        """calculates the fold axis/ fold axis rotation and adds this to the fold"""
        kwargs = self.build_arguments
        fold_axis = kwargs.get("fold_axis", None)
        if fold_axis is not None:
            fold_axis = np.array(fold_axis)
            if len(fold_axis.shape) == 1:
                self.fold.fold_axis = fold_axis

        if "av_fold_axis" in kwargs:
            l2 = self.fold.foldframe.calculate_intersection_lineation(self)
            self.fold.fold_axis = np.mean(l2, axis=0)
        if self.fold.fold_axis is None:
            if not self.fold.foldframe[1].is_valid():
                raise InterpolatorError("Fold frame direction coordinate is not valid")
            far, fad = self.fold.foldframe.calculate_fold_axis_rotation(self)
            fold_axis_rotation = get_fold_rotation_profile(self.axis_profile_type, far, fad)
            if "axis_function" in kwargs:
                # allow predefined function to be used
                logger.error("axis_function is deprecated, use a specific fold rotation angle profile type")
            else:
                fold_axis_rotation.fit(params={'wavelength': kwargs.get("axis_wl", None)})
            self.fold.fold_axis_rotation = fold_axis_rotation
            fold_axis_rotation.add_observer(self)

    def set_fold_limb_rotation(self):
        """Calculates the limb rotation of the fold and adds it to the fold object"""
        kwargs = self.build_arguments
        # need to calculate the fold axis before the fold limb rotation angle
        if self.fold.fold_axis is None:
            self.set_fold_axis()
        # give option of passing own fold limb rotation function
        flr, fld = self.calculate_fold_limb_rotation_angle()

        fold_limb_rotation = get_fold_rotation_profile(self.limb_profile_type, flr, fld)
        if "limb_function" in kwargs:
            # allow for predefined functions to be used
            logger.error("limb_function is deprecated, use a specific fold rotation angle profile type")
        else:
            fold_limb_rotation.fit(params={'wavelength': kwargs.get("limb_wl", None)})

        self.fold.fold_limb_rotation = fold_limb_rotation
        fold_limb_rotation.add_observer(self)

    def calculate_fold_limb_rotation_angle(self):
        flr, fld = self.fold.foldframe.calculate_fold_limb_rotation(
            self, self.fold.get_fold_axis_orientation
        )
        return flr, fld

    # def
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
        if not self.fold.foldframe[0].is_valid():
            raise InterpolatorError("Fold frame main coordinate is not valid")
        if self.fold.fold_axis is None:
            self.set_fold_axis()
        if self.fold.fold_limb_rotation is None:
            self.set_fold_limb_rotation()
        logger.info("Adding fold to {}".format(self.name))
        self.interpolator.fold = self.fold
        # if we have fold weights use those, otherwise just use default
        # self.interpolator.add_fold_constraints(**self.fold_weights)
        # kwargs["fold_weights"] = self.fold_weights
        if "cgw" not in kwargs:
            # try adding very small cg
            kwargs["cgw"] = 0.0
        # now the fold is set up run the standard interpolation
        return super().build(data_region=data_region, **kwargs)
