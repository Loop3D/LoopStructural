import lavavu
from lavavu.vutils import is_notebook
import numpy as np
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
        mean_property_val = geological_feature.mean()
        min_property_val = geological_feature.min()
        max_property_val = geological_feature.max()
        slices = [mean_property_val]
        colour = 'red'
        painter = None
        voxet = None
        tris = None
        nodes = None
        if 'isovalue' in kwargs:
            slices = [kwargs['isovalue']]
        if 'slices' in kwargs:
            slices = kwargs['slices']
        if 'nslices' in kwargs:
            slices = np.linspace(min_property_val, max_property_val, kwargs['nslices'])
        if 'colour' in kwargs:
            colour = kwargs['colour']
        if 'paint_with' in kwargs:
            painter = kwargs['paint_with']
        if 'voxet' in kwargs:
            voxet = kwargs['voxet']
        for isovalue in slices:
            print("Creating isosurface for %f"%isovalue)
            if isovalue < min_property_val or isovalue > max_property_val:
                print("No surface to create for isovalue")
                continue #isovalue = kwargs['isovalue']
            if voxet is None:
                tris, nodes = geological_feature.support.slice(isovalue)
            # if voxet:
            #
            if nodes.shape[0] == 0:
                continue

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
            if painter is None:
                surf.colours(colour)
            if painter is not None:
                # add a property to the surface nodes for visualisation
                surf.values(painter.evaluate_value(nodes),painter.name)
                surf["colourby"] = painter.name
                surf.colourmap(lavavu.cubehelix(100))#nodes.shape[0]))
            if "normals" in kwargs:
                a = nodes[tris[:, 0], :] - nodes[tris[:, 1], :]
                b = nodes[tris[:, 0], :] - nodes[tris[:, 2], :]

                crosses = np.cross(a, b)
                crosses = crosses / (np.sum(crosses ** 2, axis=1) ** (0.5))[:, np.newaxis]
                tribc = np.mean(nodes[tris, :], axis=1)
                vec = self.lv.vectors(name+"grad", colour="black")
                vec.vertices(tribc)
                vec.vectors(crosses)


    def plot_vector_field(self, geological_feature, locations, **kwargs):
        if 'colour' not in kwargs:
            kwargs['colour'] = 'black'
        vectorslicing = 100
        if 'vectorslicing' in kwargs:
            vectorslicing = kwargs['vectorslicing']
        vector = geological_feature.evaluate_gradient(locations)
        # normalise
        vector /= np.linalg.norm(vector, axis=1)[:, None]
        vectorfield = self.lv.vectors(geological_feature.name + "_grad", **kwargs)
        vectorfield.vertices(locations)
        vectorfield.vectors(vector)
        return

    def plot_points(self, points, name, col='red'):
        p = self.lv.points(name, pointsize=4, pointtype="sphere", colour=col)
        p.vertices(points)

    def plot_vector_data(self, position, vector, name, **kwargs):
        if 'colour' not in kwargs:
            kwargs['colour'] = 'black'
        # normalise
        if position.shape[0] > 0:
            vector /= np.linalg.norm(vector, axis=1)[:, None]
            vectorfield = self.lv.vectors(name, **kwargs)
            vectorfield.vertices(position)
            vectorfield.vectors(vector)
            return

    def plot_value_data(self, position, value, name, **kwargs):
        if "pointtype" not in kwargs:
            kwargs["pointtype"] = "sphere"
        if "pointsize" not in kwargs:
            kwargs["pointsize"] = 4
        # if "colour" not in kwargs:
        #     kwargs["colour"] = "red"
        cmap = "diverge"
        if "colourmap" in kwargs:
            cmap = kwargs["colourmap"]
        p = self.lv.points(name, **kwargs)
        p.vertices(position)
        p.values(value,"v")
        p["colourby"] = "v"
        p.colourmap(cmap)
    # def plot_feature_data(self,feature,**kwargs):
    #     grad = feature.support.interpolator.
    def interactive(self):
        if is_notebook():
            self.lv.control.Panel()
            self.lv.control.ObjectList()
            self.lv.control.show()
        if not is_notebook():
            self.lv.interactive()


