from LoopStructural.interpolators.supports import support_map, SupportType


class SupportFactory:
    @staticmethod
    def create_support(support_type, **kwargs):
        if support_type is None:
            raise ValueError("No support type specified")
        if isinstance(support_type, str):
            support_type = SupportType._member_map_[support_type].numerator
        return support_map[support_type](**kwargs)

    @staticmethod
    def from_dict(d):
        d = d.copy()
        support_type = d.pop("type", None)
        if support_type is None:
            raise ValueError("No support type specified")
        return SupportFactory.create_support(support_type, **d)

    @staticmethod
    def create_support_from_bbox(
        support_type, bounding_box, nelements, element_volume=None, buffer: float = 0.2
    ):
        if isinstance(support_type, str):
            support_type = SupportType._member_map_[support_type].numerator
        bbox = bounding_box.with_buffer(buffer=buffer)
        bbox.nelements = nelements

        return support_map[support_type](
            origin=bbox.origin, step_vector=bbox.step_vector, nsteps=bbox.nsteps
        )
