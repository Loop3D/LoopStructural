class FoldRotationAngleFeature:
    def __init__(self, fold_frame, rotation):
        self.fold_frame = fold_frame
        self.rotation = rotation

    def evaluate_value(self, location):
        s1 = self.fold_frame.features[0].evaluate_value(location)
        r = self.rotation(s1)
        return r