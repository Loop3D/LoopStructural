class FaultDisplacementFeature:
    def __init__(self, fault_frame, displacement):
        self.fault_frame = fault_frame
        self.displacement = displacement

    def evaluate_value(self, location):
        fault_suface = self.fault_frame.features[0].evaluate_value(location)
        fault_displacement = self.fault_frame.features[1].evaluate_value(location)
        fault_strike = self.fault_frame.features[2].evaluate_value(location)
        d = self.displacement(fault_suface, fault_displacement, fault_strike)
        return d