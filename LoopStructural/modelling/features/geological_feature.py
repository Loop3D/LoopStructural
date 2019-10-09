from LoopStructural.supports.scalar_field import ScalarField
from LoopStructural.modelling.core.geological_points import GPoint, IPoint, TPoint
from skimage.measure import marching_cubes_lewiner as marching_cubes

import numpy as np

import logging
logger = logging.getLogger(__name__)


class GeologicalFeatureInterpolator:
    """
    A builder for a GeologicalFeature will link data to the interpolator
    and run the interpolation
    """
    def __init__(self, interpolator, **kwargs):
        """
        The interpolator to use to build the geological feature
        Parameters
        ----------
        interpolator - a GeologicalInterpolator
        kwargs - name of the feature, region to interpolate the feature
        """
        self.interpolator = interpolator
        self.name = "UnnamedFeature"
        if 'name' in kwargs:
            self.name = kwargs['name']
            self.interpolator.set_property_name(self.name)
        self.region = 'everywhere'
        if 'region' in kwargs:
            self.region = kwargs['region']
        self.data = []

    def update(self):
        pass

    def add_data(self, pos, strike = None, dip_dir = None, dip = None, dir = None,
                 val = None, plunge = None, plunge_dir = None,polarity = None):
        """
        Generic function to add data to a geological feature.
        Parameters
        ----------
        pos - required numpy array for position
        strike - optional strike
        dip_dir - optional dip_dir
        dip - optional
        dir - numpy array for vector
        val - value of constraint
        polarity - polarity of vector
        plunge - plunge value
        plunge_dir - plunge direction


        Returns
        -------

        """
        pass

    def add_strike_dip_and_value(self, pos, strike, dip, val, polarity = 1):
        """

        Parameters
        ----------
        pos
        strike
        dip
        val

        Returns
        -------

        """
        self.data.append(GPoint.from_strike_and_dip(pos, strike, dip, polarity))
        self.interpolator.add_data(self.data[-1])
        self.data.append(IPoint(pos, val))
        self.interpolator.add_data(self.data[-1])

    def add_point(self, pos, val):
        """

        Parameters
        ----------
        pos
        val

        Returns
        -------

        """
        self.data.append(IPoint(pos, val))
        self.interpolator.add_data(self.data[-1])

    def add_planar_constraint(self, pos, val):
        """

        Parameters
        ----------
        pos -
        val

        Returns
        -------

        """
        self.data.append(GPoint(pos, val))
        self.interpolator.add_data(self.data[-1])

    def add_strike_and_dip(self, pos, s, d, polarity=1):
        """

        Parameters
        ----------
        pos
        s
        d

        Returns
        -------

        """
        self.data.append(GPoint.from_strike_and_dip(pos, s, d, polarity))
        self.interpolator.add_data(self.data[-1])

    def add_plunge_and_plunge_dir(self,pos,plunge,plunge_dir):
        """

        Parameters
        ----------
        pos
        plunge
        plunge_dir

        Returns
        -------

        """
        self.data.append(GPoint.from_plunge_plunge_dir(pos,plunge,plunge_dir))
        self.interpolator.add_data(self.data[-1])

    def add_tangent_constraint(self, pos, val):
        """

        Parameters
        ----------
        pos
        val

        Returns
        -------

        """
        self.data.append(TPoint(pos, val))
        self.interpolator.add_data(self.data[-1])

    def add_tangent_constraint_angle(self, pos, s, d):
        """

        Parameters
        ----------
        pos
        s
        d

        Returns
        -------

        """
        self.data.append(TPoint(pos, s, d))
        self.interpolator.add_data(self.data[-1])

    def build(self, solver='cg', **kwargs):
        """
        Runs the interpolation
        Parameters
        ----------
        solver
        kwargs

        Returns
        -------

        """
        # for d in self.data:
        #     self.interpolator.add_data(d)
        # we can add a fold to the interpolator if the interpolator is a fold interpolator
        # pass the dict with weights as kwargs to the fold interpolator
        # self.interpolator.reset()
        self.interpolator.set_region(regionname=self.region)
        if "fold" in kwargs and "fold_weights" in kwargs:
            self.interpolator.update_fold(kwargs['fold'])
            self.interpolator.add_fold_constraints(**kwargs['fold_weights'])
            if kwargs['cgw'] == 0:
                kwargs['cg'] = False
            # kwargs['cg'] = False
        self.interpolator.setup_interpolator(**kwargs)
        self.interpolator.solve_system(solver=solver,**kwargs)
        return GeologicalFeature(self.name,
                                 ScalarField.from_interpolator(self.interpolator),
                                 builder=self, data=self.data, region=self.region)


class GeologicalFeature:
    """
    Geological feature is class that is used to represent a geometrical element in a geological
    model. For example foliations, fault planes, fold rotation angles etc. The feature has a support
    which 
    """
    def __init__(self, name, support, builder = None, data = None, region = None):
        """

        Parameters
        ----------
        name: string

        support
        builder
        data
        region
        """
        self.name = name
        self.support = support
        self.ndim = 1
        self.data = data
        self.builder = builder
        self.region = region
        if region is None:
            self.region = 'everywhere'

    def set_builder(self, builder):
        """

        Parameters
        ----------
        builder

        Returns
        -------

        """
        self.builder = builder

    def evaluate_value(self, evaluation_points):
        """

        Parameters
        ----------
        evaluation_points

        Returns
        -------

        """
        v = np.zeros(evaluation_points.shape[0])
        v[:] = np.nan
        return self.support.evaluate_value(evaluation_points)

    def evaluate_gradient(self, locations):
        """

        Parameters
        ----------
        locations

        Returns
        -------

        """
        return self.support.evaluate_gradient(locations)

    def mean(self):
        """
        Calculate average of the support values
        Returns
        -------

        """
        return np.nanmean(self.support.get_node_values())

    def min(self):
        """

        Returns
        -------

        """
        return np.nanmin(self.support.get_node_values())

    def max(self):
        """
        Calculate average of the support values
        Returns
        -------

        """
        return np.nanmax(self.support.get_node_values())

    def update(self):
        """
        Calculate average of the support values
        Returns
        -------

        """
        # re-run the interpolator and update the support.
        # this is a bit clumsy and not abstract, i think
        # if evaluating the property doesn't require the dictionary on
        # the nodes and actually just uses the interpolator values this would be
        # much better.
        self.support.interpolator.up_to_date = False
        self.support.interpolator.update()
        self.support.update_property(self.support.interpolator.c)

    def get_interpolator(self):
        """
        Get the interpolator used to build this feature
        Returns
        -------
        GeologicalInterpolator
        """
        return self.support.interpolator

    def get_node_values(self):
        """
        Get the node values of the support used to build this interpolator if the
        interpolator is a discrete interpolator
        Returns
        -------
        numpy array of values
        """
        return self.support.get_node_values()

    def slice(self, isovalue, bounding_box = None, nsteps = None):
        """
        Calculate an isosurface of a geological feature.
        Option to specify a new support to calculate the isosurface on
        Parameters
        ----------
        isovalue
        bounding_box
        nsteps

        Returns
        -------

        """
        if bounding_box is not None and nsteps is not None:
            x = np.linspace(bounding_box[0,0],bounding_box[1,0],nsteps[0])
            y = np.linspace(bounding_box[0,1],bounding_box[1,1],nsteps[1])
            z = np.linspace(bounding_box[1,2],bounding_box[0,2],nsteps[2])

            xx,yy,zz = np.meshgrid(x,y,z, indexing='ij')
            val = self.evaluate_value(np.array([xx.flatten(),yy.flatten(),zz.flatten()]).T)
            step_vector = np.array([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
            try:
                verts, faces, normals, values = marching_cubes(
                val.reshape(nsteps, order='C'),
                isovalue,
                spacing=step_vector)

                return faces, verts + np.array([bounding_box[0,0],bounding_box[0,1],bounding_box[1,2]])
            except ValueError:
                logger.warning("No surface to mesh, skipping")
                return np.zeros((3,1)).astype(int),np.zeros((3,1))
        else:
            try:
                return self.support.slice(isovalue,self.region)
            except RuntimeError:
                logger.warning("No surface to mesh, skipping")
                return np.zeros((3,1)).astype(int),np.zeros((3,1))


