import lavavu
class ModelViewer():
    def __init__(self,modelsupport,**kwargs):
        """
        Visualisation object for visualising FME models.
        :param modelsupport: the support the interpolator is stored on e.g. structured grid or unstructured mesh
        :param kwargs: possible kwargs 'backend' lavavu, vista, matplotlib
        """

    def add_isosurface(self,**kwargs):

    def add_vector_field(self,**kwargs):

    def add_volume(self,**kwargs):

    def plot_gradient_constraints(self,**kwargs):

    def plot_value_constraints(self,**kwargs):

    def save_state(self,**kwargs):

    def load_state(self,**kwargs):

class LavaVuModelViewer()
    def __init__(self,interpolator,**kwargs):
        """
        ModelPlotter is a wrapper to link a geologicalinterpolator to a lavavu layout
        :param interpolator: the geological interpolator
        :param kwargs: kwargs for lava vu
        """
        self.lv = lavavu.Viewer(**kwargs)
        self.objects = {}



        def lv_plot_isosurface(self, propertyname, lv, **kwargs):
            # import lavavu in case its not imported
            try:
                import lavavu  # visualisation library
            except ImportError:
                print("Cannot import Lavavu")
                return
                # if no isovalue is specified just use the average
            property = self.properties[propertyname]
            slices = [np.mean(property)]
            colour = 'red'
            if 'isovalue' in kwargs:
                slices = [kwargs['isovalue']]
            if 'slices' in kwargs:
                slices = kwargs['slices']
            if 'nslices' in kwargs:
                slices = np.linspace(np.min(property), np.max(property), kwargs['nslices'])
            if 'colour' in kwargs:
                colour = kwargs['colour']
            for isovalue in slices:
                if isovalue < np.min(property) or isovalue > np.max(property):
                    print("No surface to create for isovalue")
                    isovalue = kwargs['isovalue']
                reg = np.zeros(self.properties[propertyname].shape).astype(bool)
                reg[:] = True
                if 'region' in kwargs:
                    reg = self.regions[kwargs['region']]
                name = propertyname + '_%f' % isovalue
                if 'name' in kwargs:
                    name = kwargs['name']

                tri, ntri = marching_tetra(isovalue, self.elements, self.nodes, reg, self.properties[propertyname])

                ##convert from memoryview to np array
                tri = np.array(tri)
                ntri = np.array(ntri)[0]
                ##create a triangle indices array and initialise to -1
                tris = np.zeros((ntri, 3)).astype(int)
                tris[:, :] = -1
                ##create a dict for storing nodes index where the key is the node as as a tuple.
                # A dict is preferable because it is very quick to check if a key exists
                # assemble arrays for unique vertex and triangles defined by vertex indices
                nodes = {}
                n = 0  # counter
                for i in range(ntri):
                    for j in range(3):
                        if tuple(tri[i, j, :]) in nodes:
                            tris[i, j] = nodes[tuple(tri[i, j, :])]
                        else:
                            nodes[tuple(tri[i, j, :])] = n
                            tris[i, j] = n
                            n += 1
                # find the normal vector to the faces using the vertex order
                a = tri[:ntri, 0, :] - tri[:ntri, 1, :]
                b = tri[:ntri, 0, :] - tri[:ntri, 2, :]

                crosses = np.cross(a, b)
                crosses = crosses / (np.sum(crosses ** 2, axis=1) ** (0.5))[:, np.newaxis]

                # get barycentre of faces and findproperty gradient from scalar field
                tribc = np.mean(tri[:ntri, :, :], axis=1)
                tribc = tribc[np.any(np.isnan(tribc), axis=1)]
                propertygrad = self.eval_gradient(tribc, prop='strati')
                propertygrad /= np.linalg.norm(propertygrad, axis=1)[:, None]

                # dot product between gradient and normal indicates if faces are incorrectly ordered
                dotproducts = (propertygrad * crosses[np.any(np.isnan(tribc), axis=1)]).sum(axis=1)
                #
                if dot product > 0 then adjust triangle indexing
                indices = (dotproducts > 0).nonzero()[0]
                tris[indices] = tris[indices, ::-1]
            nodes_np = np.zeros((n, 3))
            for v in nodes.keys():
                nodes_np[nodes[v], :] = np.array(v)
            surf = lv.triangles(name)
            surf.vertices(nodes_np)
            surf.indices(tris)
            surf.colours(colour)

        return lv

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
class VistaModelViewer()
    def __init__(self,**kwargs):
