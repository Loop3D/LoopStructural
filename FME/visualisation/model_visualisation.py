import lavavu
import numpy as np
import meshpy
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
def surface_cutter(feature,isovalue,nodes,tris):
    values = feature.evaluate_value(nodes)
    property_bool = values > isovalue
    tri_type_index = np.sum(property_bool[tris]*np.array([1,2,4]),axis=1)



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
        mean_property_val = geological_feature.mean_property_value()
        min_property_val = geological_feature.min_property_value()
        max_property_val = geological_feature.max_property_value()
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
                continue #isovalue = kwargs['isovalue']

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
            if "normals" in kwargs:
                a = nodes[tris[:, 0], :] - nodes[tris[:, 1], :]
                b = nodes[tris[:, 0], :] - nodes[tris[:, 2], :]

                crosses = np.cross(a, b)
                crosses = crosses / (np.sum(crosses ** 2, axis=1) ** (0.5))[:, np.newaxis]
                tribc = np.mean(nodes[tris, :], axis=1)
                vec = self.lv.vectors(name+"grad", colour="black")
                vec.vertices(tribc)
                vec.vectors(crosses)
    def plot_structural_frame_isosurface(self, structural_frame, i, **kwargs):
        mean_property_val = structural_frame.supports[i].mean_property_value()
        min_property_val = structural_frame.supports[i].min_property_value()
        max_property_val = structural_frame.supports[i].max_property_value()
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

            tris, nodes = structural_frame.supports[i].slice(isovalue)
            # reg = np.zeros(self.properties[propertyname].shape).astype(bool)
            # reg[:] = True
            # if 'region' in kwargs:
            #     reg = self.regions[kwargs['region']]
            name = structural_frame.name + '_%i_iso_%f' % (i,isovalue)
            if 'name' in kwargs:
                name = kwargs['name']
            surf = self.lv.triangles(name)
            surf.vertices(nodes)
            surf.indices(tris)
            surf.colours(colour)

    def lv_plot_vector_field(self, propertyname, lv, **kwargs):
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

    def plot_points(self, points, name, col='red'):
        p = self.lv.points(name, pointsize=4, pointtype="sphere", colour=col)
        p.vertices(points)

    def plot_vector_data(self, position, vector, name, **kwargs):
        if 'colour' not in kwargs:
            kwargs['colour'] = 'black'
        # normalise
        vector /= np.linalg.norm(vector, axis=1)[:, None]
        vectorfield = self.lv.vectors(name, **kwargs)
        vectorfield.vertices(position)
        vectorfield.vectors(vector)
        return


