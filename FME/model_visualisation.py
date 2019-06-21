import lavavu
class ModelPlotter():
    def __init__(self,interpolator,**kwargs):
        """
        ModelPlotter is a wrapper to link a geologicalinterpolator to a lavavu layout
        :param interpolator: the geological interpolator
        :param kwargs: kwargs for lava vu
        """
        self.lv = lavavu.Viewer(**kwargs)
        self.objects = {}
    def add_isosurface(self):

    def add_vector_field(self):

    def add_volume(self):

    def plot_gradient_constraints(self):

    def plot_value_constraints(self):