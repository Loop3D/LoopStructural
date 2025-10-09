from ....modelling.features import BaseFeature, StructuralFrame
from typing import Optional
from ....utils import getLogger

logger = getLogger(__name__)


class FaultDisplacementFeature(BaseFeature):
    """Geological feature representing fault displacement.

    This class models the displacement associated with a fault surface using
    a fault frame and displacement function.

    Parameters
    ----------
    fault_frame : StructuralFrame
        The geometric frame describing the fault surface
    displacement : callable
        Function defining the fault displacement
    name : str, optional
        Name of the fault displacement feature, by default "fault_displacement"
    model : GeologicalModel, optional
        The geological model containing this feature, by default None
    faults : list, optional
        List of associated faults, by default []
    regions : list, optional
        List of regions where this feature applies, by default []
    builder : object, optional
        Builder object used to create this feature, by default None
    """

    def __init__(
        self,
        fault_frame,
        displacement,
        name="fault_displacement",
        model=None,
        faults=[],
        regions=[],
        builder=None,
    ):
        """Initialize the fault displacement feature.

        Parameters
        ----------
        fault_frame : StructuralFrame
            The geometric frame describing the fault surface
        displacement : callable
            Function defining the fault displacement
        name : str, optional
            Name of the fault displacement feature, by default "fault_displacement"
        model : GeologicalModel, optional
            The geological model containing this feature, by default None
        faults : list, optional
            List of associated faults, by default []
        regions : list, optional
            List of regions where this feature applies, by default []
        builder : object, optional
            Builder object used to create this feature, by default None
        """
        BaseFeature.__init__(self, f"{name}_displacement", model, faults, regions, builder)
        self.fault_frame = StructuralFrame(
            f"{fault_frame.name}_displacementframe",
            [fault_frame[0].copy(), fault_frame[1].copy(), fault_frame[2].copy()],
        )
        self.displacement = displacement

    def evaluate_value(self, location):
        """Return the value of the fault displacement at given locations.

        Parameters
        ----------
        location : np.ndarray
            Array of xyz coordinates where displacement should be evaluated

        Returns
        -------
        np.ndarray
            Fault displacement values at the given locations
        """
        fault_suface = self.fault_frame.features[0].evaluate_value(location)
        fault_displacement = self.fault_frame.features[1].evaluate_value(location)
        fault_strike = self.fault_frame.features[2].evaluate_value(location)
        d = self.displacement(fault_suface, fault_displacement, fault_strike)
        return d

    def evaluate_gradient(self, location):
        """Get the scaled displacement gradient at given locations.

        Parameters
        ----------
        location : np.ndarray
            Array of xyz coordinates where displacement gradient should be evaluated

        Returns
        -------
        np.ndarray
            Fault displacement gradient values at the given locations
        """
        fault_suface = self.fault_frame.features[0].evaluate_value(location)
        fault_displacement = self.fault_frame.features[1].evaluate_value(location)
        fault_strike = self.fault_frame.features[2].evaluate_value(location)
        d = self.displacement(fault_suface, fault_displacement, fault_strike)
        return d

    def evaluate_on_surface(self, location):
        """Evaluate displacement specifically on the fault surface.

        Parameters
        ----------
        location : np.ndarray
            Array of xyz coordinates on the fault surface

        Returns
        -------
        np.ndarray
            Fault displacement values on the surface

        Notes
        -----
        This method evaluates displacement only considering the fault displacement
        and strike components, not the fault surface component.
        """
        fault_displacement = self.fault_frame.features[1].evaluate_value(location)
        fault_strike = self.fault_frame.features[2].evaluate_value(location)
        d = self.displacement.evaluate(fault_displacement, fault_strike)
        return d

    def get_data(self, value_map: Optional[dict] = None):
        """Get data associated with this fault displacement feature.

        Parameters
        ----------
        value_map : dict, optional
            Optional mapping of values, by default None

        Notes
        -----
        This method is not yet implemented for fault displacement features.
        """
        pass

    def copy(self, name: Optional[str] = None):
        """Create a copy of this fault displacement feature.

        Parameters
        ----------
        name : str, optional
            Name for the copied feature, by default None

        Raises
        ------
        NotImplementedError
            This method is not yet implemented
        """
        raise NotImplementedError("Not implemented yet")
