from ....modelling.features import BaseFeature
from ....utils import getLogger

logger = getLogger(__name__)


class FoldRotationAngleFeature(BaseFeature):
    """ """

    def __init__(
        self,
        fold_frame,
        rotation,
        name="fold_rotation_angle",
        model=None,
        faults=[],
        regions=[],
        builder=None,
    ):
        """

        Parameters
        ----------
        fold_frame
        rotation
        """
        BaseFeature.__init__(self, f"{name}_displacement", model, faults, regions, builder)
        self.fold_frame = fold_frame
        self.rotation = rotation

    def evaluate_value(self, location):
        """

        Parameters
        ----------
        location

        Returns
        -------

        """
        s1 = self.fold_frame.features[0].evaluate_value(location)
        r = self.rotation(s1)
        return r
    def copy(self, name = None):
        raise NotImplementedError("FoldRotationAngleFeature cannot be copied directly, copy the fold frame and rotation function separately")
    def evaluate_gradient(self, pos, ignore_regions=False):
        raise NotImplementedError("FoldRotationAngleFeature does not have a gradient")
    def get_data(self, value_map = None):
        raise NotImplementedError("FoldRotationAngleFeature does not have data associated with it directly, get data from the fold frame and rotation function separately")