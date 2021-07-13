class BaseFeatureProperties:
    def __init__(self):
        pass
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return self.__class__.__name__

    