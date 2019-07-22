import lavavu
# class ModelViewer():
#     def __init__(self,modelsupport,**kwargs):
#         """
#         Visualisation object for visualising FME models.
#         :param modelsupport: the support the interpolator is stored on e.g. structured grid or unstructured mesh
#         :param kwargs: possible kwargs 'backend' lavavu, vista, matplotlib
#         """
#
#     def add_isosurface(self,**kwargs):
#
#     def add_vector_field(self,**kwargs):
#
#     def add_volume(self,**kwargs):
#
#     def plot_gradient_constraints(self,**kwargs):
#
#     def plot_value_constraints(self,**kwargs):
#
#     def save_state(self,**kwargs):
#
#     def load_state(self,**kwargs):


class LavaVuModelViewer:
    def __init__(self, **kwargs):
        """
        ModelPlotter is a wrapper to link a geologicalinterpolator to a lavavu layout
        :param interpolator: the geological interpolator
        :param kwargs: kwargs for lava vu
        """
        self.lv = lavavu.Viewer(**kwargs)
        self.objects = {}

    def plot_isosurface(self, geological_feature, **kwargs):
        mean_property_val = geological_feature.support.mean_property_value()
        min_property_val = geological_feature.support.min_property_value()
        max_property_val = geological_feature.support.max_property_value()
        slices = [mean_property_val]
        colour = 'red'

        if 'isovalue' in kwargs:
            slices = [kwargs['isovalue']]
        if 'slices' in kwargs:
            slices = kwargs['slices']
        if 'nslices' in kwargs:
            slices = np.linspace(min_property_val, max_property_val, kwargs['nslices'])
        if 'colour' in kwargs:
            colour = kwargs['colour']
        for isovalue in slices:
            if isovalue < min_property_val or isovalue > max_property_val:
                print("No surface to create for isovalue")
                isovalue = kwargs['isovalue']

            tris, nodes = geological_feature.support.slice(isovalue)
            # reg = np.zeros(self.properties[propertyname].shape).astype(bool)
            # reg[:] = True
            # if 'region' in kwargs:
            #     reg = self.regions[kwargs['region']]
            name = geological_feature.name + '_iso_%f' % isovalue
            if 'name' in kwargs:
                name = kwargs['name']
            surf = self.lv.triangles(name)
            surf.vertices(nodes)
            surf.indices(tris)
            surf.colours(colour)

    def lv_plot_vector_field(self, propertyname, lv, **kwargs):
        try:
            import lavavu
        except ImportError:
            print("Cannot import Lavavu")
        if 'colour' not in kwargs:
            kwargs['colour'] = 'black'
        vectorslicing = 100
        if 'vectorslicing' in kwargs:
            vectorslicing = kwargs['vectorslicing']
        vector = self.property_gradients[propertyname]
        # normalise
        vector /= np.linalg.norm(vector, axis=1)[:, None]
        vectorfield = lv.vectors(propertyname + "_grad", **kwargs)
        vectorfield.vertices(self.barycentre[::vector_slicing, :])
        vectorfield.vectors(vectors)
        return

