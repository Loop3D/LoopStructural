import json
class LoopJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        """All jsonable loop objects should have a tojson method

        Parameters
        ----------
        obj : LoopStructuralObject
            An object from loopstructural

        Returns
        -------
        str
            string representing the json encoding
        """
        return obj.__tojson__()
        