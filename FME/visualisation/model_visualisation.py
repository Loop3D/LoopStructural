import lavavu
from lavavu.vutils import is_notebook
import numpy as np


class LavaVuModelViewer:
    def __init__(self, **kwargs):
        """
        LavaVuModelViewer


        A wrapper to plot FME object with lavavu

        Parameters
        ----------
        **kwargs : lavavu viewere kwargs

        Attributes
        ----------
        lv Lavavu.Viewer object
        objects : dictionary of objects that have been plotted
        """

        self.lv = lavavu.Viewer(**kwargs)
        self.objects = {}

    def plot_isosurface(self, geological_feature, **kwargs):
        """
        Plot the surface of a geological feature if no kwargs are given plots the
        average surface and colours it red.
        Parameters
        ----------
        geological_feature
        isovalue

        Returns
        -------

        """
        # update the feature to make sure its current
        if 'update' in kwargs:
            geological_feature.update()
        # get the stats to check what we are plotting
        mean_property_val = geological_feature.mean()
        min_property_val = geological_feature.min()
        max_property_val = geological_feature.max()

        # set default parameters
        slices = [mean_property_val]
        colour = 'red'
        painter = None
        voxet = None
        tris = None
        nodes = None

        # parse kwargs for parameters
        if 'isovalue' in kwargs:
            slices = [kwargs['isovalue']]
        if 'slices' in kwargs:
            slices = kwargs['slices']
        if 'nslices' in kwargs:
            var = max_property_val - min_property_val
            # buffer slices by 5%
            slices = np.linspace(min_property_val+var*0.05, max_property_val-var*0.05, kwargs['nslices'])
        if 'colour' in kwargs:
            colour = kwargs['colour']
        if 'paint_with' in kwargs:
            painter = kwargs['paint_with']

        # do isosurfacing of tetrahedral mesh using marching tetras
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
                cmap = lavavu.cubehelix(100)
                if 'cmap' in kwargs:
                    cmap = kwargs['cmap']
                surf.colourmap(cmap)#nodes.shape[0]))
            if "normals" in kwargs:
                a = nodes[tris[:, 0], :] - nodes[tris[:, 1], :]
                b = nodes[tris[:, 0], :] - nodes[tris[:, 2], :]

                crosses = np.cross(a, b)
                crosses = crosses / (np.sum(crosses ** 2, axis=1) ** (0.5))[:, np.newaxis]
                tribc = np.mean(nodes[tris, :], axis=1)
                vec = self.lv.vectors(name+"grad", colour="black")
                vec.vertices(tribc[::kwargs['normals'],:])
                vec.vectors(crosses[::kwargs['normals'],::])

    def plot_vector_field(self, geological_feature, locations, **kwargs):
        """
        Plot the gradient of a geological feature at given locations
        Parameters
        ----------
        geological_feature : Geological Feature to evaluate gradient
        locations : ((N,3)) array of evaluation locations
        kwargs : kwargs for lavavu vector

        Returns
        -------

        """
        vector = geological_feature.evaluate_gradient(locations)
        # normalise
        vector /= np.linalg.norm(vector, axis=1)[:, None]
        vectorfield = self.lv.vectors(geological_feature.name + "_grad", **kwargs)
        vectorfield.vertices(locations)
        vectorfield.vectors(vector)
        return

    def plot_points(self, points, name, **kwargs):
        """
        Plot points location in the lavavu viewer
        Parameters
        ----------
        points : numpy array of the points locations
        name :  string name of the object for lavavu
        **kwargs : lavavu points kwargs

        Returns
        -------

        """
        p = self.lv.points(name, **kwargs)
        p.vertices(points)

    def plot_vector_data(self, position, vector, name, **kwargs):
        """
        Plot point data with a vector component into the lavavu viewer
        Parameters
        ----------
        position : numpy array N,3 for xyz locations
        vector : numpy array of vector N,3
        name :  string name for the object in lavavu
        kwargs to pass to lavavu

        Returns
        -------

        """
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
        """
        Plot points data with a value component
        Parameters
        ----------
        position : numpy array N,3 for xyz locations
        value : N array of values
        name :  string name of the object for lavavu
        kwargs : kwargs to pass to lavavu

        Returns
        -------

        """
        if "pointtype" not in kwargs:
            kwargs["pointtype"] = "sphere"
        if "pointsize" not in kwargs:
            kwargs["pointsize"] = 4
        # set the colour map to diverge unless user decides otherwise
        cmap = "diverge"
        if "colourmap" in kwargs:
            cmap = kwargs["colourmap"]
        p = self.lv.points(name, **kwargs)
        p.vertices(position)
        p.values(value,"v")
        p["colourby"] = "v"
        p.colourmap(cmap)

    def interactive(self):
        """
        Runs the lavavu viewer as either a jupyter notebook
        inline interactive viewer or as a separate window
        Returns
        -------

        """
        if is_notebook():
            self.lv.control.Panel()
            self.lv.control.ObjectList()
            self.lv.control.show()
        if not is_notebook():
            self.lv.interactive()


