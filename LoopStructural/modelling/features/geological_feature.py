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
        # everywhere region is just a lambda that returns true for all locations
        self.region = lambda pos : np.ones(pos.shape[0], dtype=bool)

        if 'region' in kwargs:
            print('region kwarg')
            self.region = kwargs['region']
        self.data = []
        self.data_original = []
        self.faults = []
        self.data_added = False

    def update(self):
        pass

    def add_fault(self, fault):
        """
        Add a fault to the geological feature builder
        Parameters
        ----------
        fault FaultSegment
            A faultsegment to add to the geological feature

        Returns
        -------

        """
        self.faults.append(fault)

    def add_data_from_data_frame(self, data_frame):
        """
        Extract data from a pandas dataframe with columns for
        Parameters
        ----------
        data_frame - pandas data frame

        Returns
        -------

        """
        if 'X' not in data_frame.columns or 'X' not in data_frame.columns or 'X' not in data_frame.columns:
            logger.error("No location in data frame")
            return
        for i, r in data_frame.iterrows():

            if np.isnan(r['X']) or np.isnan(r['X']) or np.isnan(r['X']):
                continue
            pos = r[['X', 'Y', 'Z']]
            if 'val' in data_frame.columns and ~np.isnan(r['val']):
                self.add_point(pos,r['val'])
            if 'strike' in data_frame.columns and  'dip' in data_frame.columns and  \
                    ~np.isnan(r['strike']) and ~np.isnan(r['dip']):
                polarity = 1
                if 'polarity' in data_frame.columns and ~np.isnan(r['polarity']):
                    polarity = r['polarity']
                self.add_strike_and_dip(pos, r['strike'], r['dip'], polarity=polarity)
            if 'azimuth' in data_frame.columns and 'dip' in data_frame.columns and \
                    ~np.isnan(r['azimuth']) and ~np.isnan(r['dip']):
                polarity = 1
                if 'polarity' in data_frame.columns and ~np.isnan(r['polarity']):
                    polarity = r['polarity']
                self.add_plunge_and_plunge_dir(pos,r['dip'],r['azimuth'],polarity=polarity)

            if 'nx' in data_frame.columns and 'ny' in data_frame.columns and 'nz' in data_frame.columns and \
                    ~np.isnan(r['nx']) and ~np.isnan(r['ny'])and ~np.isnan(r['nz']):
                 self.add_planar_constraint(r[['X','Y','Z']],r[['nx','ny','nz']])
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
        # self.interpolator.add_data(self.data[-1])
        self.data.append(IPoint(pos, val))
        # self.interpolator.add_data(self.data[-1])

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
        # self.interpolator.add_data(self.data[-1])

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
        # self.interpolator.add_data(self.data[-1])

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
        # self.interpolator.add_data(self.data[-1])

    def add_plunge_and_plunge_dir(self, pos, plunge, plunge_dir, polarity=1):
        """

        Parameters
        ----------
        pos
        plunge
        plunge_dir

        Returns
        -------

        """
        self.data.append(GPoint.from_plunge_plunge_dir(pos,plunge,plunge_dir,polarity))
        # self.interpolator.add_data(self.data[-1])

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
        # self.interpolator.add_data(self.data[-1])

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
        # self.interpolator.add_data(self.data[-1])

    def add_data_to_interpolator(self, constrained = False):
        """
        Iterates through the list of data and applies any faults active on the
        data in the order they are added
        Returns
        -------

        """
        # first move the data for the fault
        for f in self.faults:
            f.apply_to_data(self.data)
        # Now check whether there are enough constraints for the interpolator to be able to solve
        # we need at least 2 different value points or a single norm constraint. If there are not enough
        # try converting grad to norms, if still not enough send user an error
        vals = []
        for d in self.data:
            if d.type == "GPoint" and d.norm == True:
                constrained = True
                break
            if d.type == 'IPoint':
                vals.append(d.val)
        if len(np.unique(vals)) > 1:
            constrained = True
        if not constrained:
            for d in self.data:
                if d.type == "GPoint":
                    d.norm = True
                    logger.warning("Setting gradient points to norm constraints")
                    constrained = True
        if not constrained:
            logger.error("Not enough constraints for scalar field add more")
        for d in self.data:
            self.interpolator.add_data(d)

    def build(self, solver='pyamg', **kwargs):
        """
        Runs the interpolation
        Parameters
        ----------
        solver
        kwargs

        Returns
        -------

        """
        if not self.data_added:
            self.add_data_to_interpolator()
        self.interpolator.set_region(region=self.region)
        if "fold" in kwargs and "fold_weights" in kwargs:
            self.interpolator.update_fold(kwargs['fold'])
            self.interpolator.add_fold_constraints(**kwargs['fold_weights'])
            if 'cgw' not in kwargs:
                kwargs['cgw'] = 0.
            # if self.interpolator.interpolation_weights['cgw'] == 0:
            #     kwargs['cg'] = False
            # kwargs['cg'] = False
        self.interpolator.setup_interpolator(**kwargs)
        self.interpolator.solve_system(solver=solver,**kwargs)
        return GeologicalFeature(self.name,
                                 ScalarField.from_interpolator(self.interpolator),
                                 builder=self, data=self.data, region=self.region,
                                 faults = self.faults)


class GeologicalFeature:
    """
    Geological feature is class that is used to represent a geometrical element in a geological
    model. For example foliations, fault planes, fold rotation angles etc. The feature has a support
    which 
    """
    def __init__(self, name, support, builder = None, data = None, region = None, type = None, faults = []):
        """

        Parameters
        ----------
        name: string

        support
        builder
        data
        region

        Attributes
        ----------
        name - string
            should be a unique name for the geological feature
        support - a ScalarField
            holds the property values for the feature and links to the support geometry
        data
        regions - list of boolean functions defining whether the feature is active
        faults - list of faults that affect this feature
        """
        self.name = name
        self.support = support
        self.ndim = 1
        self.data = data
        self.builder = builder
        self.region = region
        self.regions = []
        self.type = type
        self.faults = faults
        if region is None:
            self.region = 'everywhere'

    def __str__(self):
        return self.name

    def add_region(self,region):
        """

        Parameters
        ----------
        region - boolean function(x,y,z)
                - returns true if inside region, false if outside
                can be passed as a lambda function e.g.
                lambda pos : feature.evaluate_value(pos) > 0
        Returns
        -------

        """
        self.regions.append(region)

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


        # check if the points are within the display region
        v = np.zeros(evaluation_points.shape[0])
        v[:] = np.nan
        mask = np.zeros(evaluation_points.shape[0]).astype(bool)
        mask[:] = True
        # check regions
        for r in self.regions:
            mask = np.logical_and(mask,r(evaluation_points))
        # apply faulting after working out which regions are visible
        for f in self.faults:
            evaluation_points = f.apply_to_points(evaluation_points)
        v[mask] = self.support.evaluate_value(evaluation_points[mask, :])
        return v#self.support.evaluate_value(evaluation_points)

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
            if isovalue > np.nanmax(val) or isovalue < np.nanmin(val):
                logger.warning("Isovalue doesn't exist inside bounding box")
                return np.zeros((3,1)).astype(int),np.zeros((3,1))
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


