import lavavu
from lavavu.vutils import is_notebook
from LoopStructural.utils.helper import create_surface
import numpy as np
import scipy as sp
import logging
logger = logging.getLogger(__name__)

class LavaVuModelViewer:
    def __init__(self, **kwargs):
        """
        LavaVuModelViewer


        A wrapper to plot LoopStructural object with lavavu

        Parameters
        ----------
        **kwargs : lavavu viewer kwargs

        Attributes
        ----------
        lv Lavavu.Viewer object
        objects : dictionary of objects that have been plotted
        """

        self.lv = lavavu.Viewer(**kwargs)
        self.lv['orthographic'] = True
        self.objects = {}
        self.model = kwargs.get("model",None)


    def add_section(self, geological_feature, axis='x', value = None, boundary_points = None, nsteps = None, **kwargs):
        if boundary_points is None:
            return False #boundary_points = geological_feature.support.support.

        # nsteps = np.array(nsteps)
        # xx = None
        # yy = None
        # zz = None
        if axis == 'x':
            tri, yy, zz = create_surface(boundary_points[:, [1, 2]], nsteps[[1, 2]])
            xx = np.zeros(zz.shape)
            xx[:] = boundary_points[1, 2]
        if axis == 'y':
            tri, xx, zz = create_surface(boundary_points[:, [0, 2]], nsteps[[0, 2]])
            yy = np.zeros(xx.shape)
            yy[:] = boundary_points[1, 2]
        if axis == 'z':
            tri, xx, yy = create_surface(boundary_points[0:2, :], nsteps[0:2])
            zz = np.zeros(xx.shape)
            zz[:] = boundary_points[1, 2]
        name = kwargs.get('name',axis+'_slice')
        colour = kwargs.get('colour','red')
        points = np.zeros((len(xx), 3))  #
        points[:, 0] = xx
        points[:, 1] = yy
        points[:, 2] = zz

        surf = self.lv.triangles(name)
        surf.vertices(points)
        surf.indices(tri)
        if geological_feature is None:
            surf.colours(colour)
        if geological_feature is not None:
            if 'norm' in kwargs:
                surf.values(np.linalg.norm(geological_feature.evaluate_gradient(points), axis=1),
                            geological_feature.name)
            else:
                surf.values(geological_feature.evaluate_value(points), geological_feature.name)
            surf["colourby"] = geological_feature.name
        cmap = lavavu.cubehelix(100)
        if 'cmap' in kwargs:
            cmap = kwargs['cmap']
        surf.colourmap(cmap)

    def add_isosurface(self, geological_feature, **kwargs):
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
        if 'voxet' in kwargs:
            voxet = kwargs['voxet']
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
        # do isosurfacing of support using marching tetras/cubes
        for isovalue in slices:
            logger.debug("Creating isosurface for %f"%isovalue)
            if isovalue < min_property_val or isovalue > max_property_val:
                logger.debug("No surface to create for isovalue")
                continue #isovalue = kwargs['isovalue']
            if voxet is None:
                tris, nodes = geological_feature.slice(isovalue)
            if voxet:
                tris, nodes = geological_feature.slice(isovalue,
                                                       bounding_box=voxet['bounding_box'],
                                                       nsteps=voxet['nsteps'])
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
                # calculate the mode value, just to get the most common value
                val = np.zeros(nodes.shape[0])
                val[:] = isovalue
                surf.values(val, painter.name)
                surf["colourby"] = painter.name
                cmap = lavavu.cubehelix(100)
                if 'cmap' in kwargs:
                    cmap = kwargs['cmap']
                surf.colourmap(cmap,range= (geological_feature.min(),geological_feature.max()))#nodes.shape[0]))

            if "normals" in kwargs:
                a = nodes[tris[:, 0], :] - nodes[tris[:, 1], :]
                b = nodes[tris[:, 0], :] - nodes[tris[:, 2], :]

                crosses = np.cross(a, b)
                crosses = crosses / (np.sum(crosses ** 2, axis=1) ** (0.5))[:, np.newaxis]
                tribc = np.mean(nodes[tris, :], axis=1)
                vec = self.lv.vectors(name+"grad", colour="black")
                vec.vertices(tribc[::kwargs['normals'],:])
                vec.vectors(crosses[::kwargs['normals'],::])

    def add_scalar_field(self, boundary_points, nsteps, name, **kwargs):
        """

        Parameters
        ----------
        boundary_points
        dims
        name
        kwargs

        Returns
        -------

        """
        colour = 'red'
        painter = None
        cmap = lavavu.cubehelix(100)
        if 'cmap' in kwargs:
            cmap = kwargs['cmap']
            if isinstance(cmap,str):
                cmap = lavavu.matplotlib_colourmap(cmap)
        if 'colour' in kwargs:
            colour = kwargs['colour']
        if 'paint_with' in kwargs:
            painter = kwargs['paint_with']



        nsteps = np.array(nsteps)
        tri, xx, yy = create_surface(boundary_points[0:2, :], nsteps[0:2])

        zz = np.zeros(xx.shape)
        zz[:] = boundary_points[1, 2]

        tri = np.vstack([tri, tri + np.max(tri)+1])
        xx = np.hstack([xx, xx])
        yy = np.hstack([yy, yy])

        z = np.zeros(zz.shape)
        z[:] = boundary_points[0, 2]
        zz = np.hstack([zz, z])
        # y faces
        t, x, z = create_surface(boundary_points[:, [0, 2]], nsteps[[0, 2]])
        tri = np.vstack([tri, t + np.max(tri)+1])
        y = np.zeros(x.shape)
        y[:] = boundary_points[0, 1]
        xx = np.hstack([xx, x])
        zz = np.hstack([zz, z])
        yy = np.hstack([yy, y])

        tri = np.vstack([tri, t + np.max(tri)+1])
        y[:] = boundary_points[1, 1]
        xx = np.hstack([xx, x])
        zz = np.hstack([zz, z])
        yy = np.hstack([yy, y])

        # x faces
        t, y, z = create_surface(boundary_points[:, [1, 2]], nsteps[[1, 2]])
        tri = np.vstack([tri, t + np.max(tri)+1])
        x = np.zeros(y.shape)
        x[:] = boundary_points[0, 0]
        xx = np.hstack([xx, x])
        zz = np.hstack([zz, z])
        yy = np.hstack([yy, y])

        tri = np.vstack([tri, t + np.max(tri)+1])
        x[:] = boundary_points[1, 0]
        xx = np.hstack([xx, x])
        zz = np.hstack([zz, z])
        yy = np.hstack([yy, y])

        points = np.zeros((len(xx), 3))  #
        points[:, 0] = xx
        points[:, 1] = yy
        points[:, 2] = zz

        surf = self.lv.triangles(name)
        surf.vertices(points)
        surf.indices(tri)
        if painter is None:
            surf.colours(colour)
        if painter is not None:
            if 'norm' in kwargs:
                surf.values(np.linalg.norm(painter.evaluate_gradient(points),axis=1),painter.name)
            else:
                surf.values(painter.evaluate_value(points),painter.name)
            surf["colourby"] = painter.name
        cmap = lavavu.cubehelix(100)
        if 'cmap' in kwargs:
            cmap = kwargs['cmap']
        surf.colourmap(cmap)  # nodes.shape[0]))
        #return tri, xx, yy, zz

    def add_vector_field(self, geological_feature, locations, **kwargs):
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

    def add_data(self, feature, **kwargs):
        """
        Plot the data linked to the feature
        Parameters
        ----------
        feature
        kwargs

        Returns
        -------

        """
        name = feature.name
        add_grad = True
        add_value = True
        if 'name' in kwargs:
            name = kwargs['name']
            del kwargs['name']
        if 'grad' in kwargs:
            add_grad = kwargs['grad']
        if 'value' in kwargs:
            add_value = kwargs['value']
        grad = feature.support.interpolator.get_gradient_constraints()
        value = feature.support.interpolator.get_value_constraints()
        if grad.shape[0] > 0 and add_grad:
            self.add_vector_data(grad[:, :3], grad[:, 3:], name + "_grad_cp", **kwargs)
        if value.shape[0] > 0 and add_value:
            kwargs['range'] = [feature.min(), feature.max()]
            self.add_value_data(value[:, :3], value[:, 3], name + "_value_cp", **kwargs)
    def add_points(self, points, name, **kwargs):
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

    def add_vector_data(self, position, vector, name, **kwargs):
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

    def add_value_data(self, position, value, name, **kwargs):
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
        cmap = kwargs.get('cmap',"diverge")
        p = self.lv.points(name, **kwargs)
        p.vertices(position)
        p.values(value,"v")
        p["colourby"] = "v"
        if 'range' in kwargs:
            p.colourmap(cmap,range=kwargs['range'])
        else:
            p.colourmap(cmap)

    def add_fold(self, fold, locations):
        """
        Draw the vector components of the fold at the locations
        Parameters
        ----------
        fold - fold object
        locations - numpy array of xyz

        Returns
        -------

        """
        r2r, fold_axis, dgz = fold.get_deformed_orientation(locations)
        self.add_vector_data(locations, r2r, 'foldr2r', colour='red')
        self.add_vector_data(locations, fold_axis, 'fold_axis', colour='black')
        self.add_vector_data(locations, dgz, 'fold_norm', colour='green')

    def interactive(self, popout=False):
        """
        Runs the lavavu viewer as either a jupyter notebook
        inline interactive viewer or as a separate window
        Returns
        -------

        """
        if is_notebook() and popout is False:
            self.lv.control.Panel()
            self.lv.control.ObjectList()
            self.lv.control.show()
        if not is_notebook() or popout:
            self.lv.control.Panel()
            self.lv.control.ObjectList()
            self.lv.interactive()

    def set_viewer_rotation(self, rotation):
        """
        Set the viewer rotation given a list of rotations x,y,z
        Parameters
        ----------
        rotation numpy array of 3 rotation

        Returns
        -------

        """
        self.lv.rotate(rotation)

    def save(self,fname,**kwargs):
        """
        Calls lavavu.Viewer.image to save the viewer current state as an image
        Parameters
        ----------
        fname - file name string including relative path
        kwargs - optional kwargs to give to lavavu e.g. transparent, resolution

        Returns
        -------

        """
        self.lv.image(fname, **kwargs)


