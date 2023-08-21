class Region:
    def __init__(self, feature, value, sign):
        self.feature = feature
        self.value = value
        self.sign = sign

    def __call__(self, xyz):
        if self.sign:
            return self.feature.evaluate_value(xyz) > 0
        else:
            return self.feature.evaluate_value(xyz) < 0

    def to_json(self):
        return {
            "feature": self.feature.name,
            "value": self.value,
            "sign": self.sign,
        }
