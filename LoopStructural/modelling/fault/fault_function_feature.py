import logging

from LoopStructural.utils import getLogger
logger = getLogger(__name__)


class FaultDisplacementFeature:
    """

    """
    def __init__(self, fault_frame, displacement, name = 'fault_displacement'):
        """
        Geological feature representing the fault displacement

        Parameters
        ----------
        fault_frame - geometry of the fault
        displacement - function defining fault displacement
        """
        self.fault_frame = fault_frame
        self.displacement = displacement
        self.name = name +'_displacement'

    def evaluate_value(self, location):
        """
        Return the value of the fault displacement

        Parameters
        ----------
        location

        Returns
        -------

        """
        fault_suface = self.fault_frame.features[0].evaluate_value(location)
        fault_displacement = self.fault_frame.features[1].evaluate_value(
            location)
        fault_strike = self.fault_frame.features[2].evaluate_value(location)
        d = self.displacement(fault_suface, fault_displacement, fault_strike)
        return d

    def evaluate_gradient(self, location):
        """
        get the scaled displacement

        Parameters
        ----------
        location

        Returns
        -------

        """
        fault_suface = self.fault_frame.features[0].evaluate_value(location)
        fault_displacement = self.fault_frame.features[1].evaluate_value(
            location)
        fault_strike = self.fault_frame.features[2].evaluate_value(location)
        d = self.displacement(fault_suface, fault_displacement, fault_strike)
        return d

    def min(self):
        return -1

    def max(self):
        return 1

    def evaluate_on_surface(self, location):
        fault_displacement = self.fault_frame.features[1].evaluate_value(
            location)
        fault_strike = self.fault_frame.features[2].evaluate_value(location)
        d = self.displacement.evaluate(fault_displacement, fault_strike)
        return d

    def mean(self):
        return 0