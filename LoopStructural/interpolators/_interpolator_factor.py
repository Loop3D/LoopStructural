from LoopStructural.interpolators import interpolator_map, InterpolatorType
from LoopStructural.interpolators.supports._support_factory import SupportFactory


class InterpolatorFactory:
    @staticmethod
    def create_interpolator(interpolator_type, support, data, c, up_to_date, **kwargs):
        if interpolator_type is None:
            raise ValueError("No interpolator type specified")
        if type(interpolator_type) == str:
            interpolator_type = InterpolatorType._member_map_[
                interpolator_type
            ].numerator
        # TODO add a warning for all kwargs that are not used
        return interpolator_map[interpolator_type](support, data, c, up_to_date)

    @staticmethod
    def from_dict(d):
        d = d.copy()
        interpolator_type = d.pop("type", None)
        if interpolator_type is None:
            raise ValueError("No interpolator type specified")
        support = d.pop("support", None)
        if support is None:
            raise ValueError("No support specified")
        support = SupportFactory.from_dict(support)
        return InterpolatorFactory.create_interpolator(interpolator_type, support, **d)
