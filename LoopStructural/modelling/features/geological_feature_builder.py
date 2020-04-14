import copy
import logging

import numpy as np

logger = logging.getLogger(__name__)

from LoopStructural.modelling.core.geological_points import GPoint, IPoint, \
    TPoint
from LoopStructural.modelling.features import GeologicalFeature
from LoopStructural.utils.helper import get_data_axis_aligned_bounding_box


class GeologicalFeatureInterpolator:
    def __init__(self, interpolator, name='Feature', region=None, **kwargs):
        """
        A builder for a GeologicalFeature will link data to the interpolator
        and run the interpolation

        Parameters
        ----------
        interpolator - a GeologicalInterpolator
        region : lambda function
            defining whether the location (xyz) should be included in the
        kwargs - name of the feature, region to interpolate the feature
        """
        self.interpolator = interpolator
        self.name = name
        self.interpolator.set_property_name(self.name)
        # everywhere region is just a lambda that returns true for all locations
        if region is None:
            self.region = lambda pos: np.ones(pos.shape[0], dtype=bool)
        else:
            self.region = region
        self.data = []
        self.data_original = []
        self.faults = []
        self.data_added = False
        self.interpolator.set_region(region=self.region)

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
        if 'X' not in data_frame.columns or 'Y' not in data_frame.columns or \
                'Z' not in data_frame.columns:
            logger.error("No location in data frame")
            return
        for i, r in data_frame.iterrows():

            if np.isnan(r['X']) or np.isnan(r['X']) or np.isnan(r['X']):
                continue
            pos = r[['X', 'Y', 'Z']]
            if 'val' in data_frame.columns and ~np.isnan(r['val']):
                self.add_point(pos, r['val'])
            if 'strike' in data_frame.columns and 'dip' in \
                    data_frame.columns and \
                    ~np.isnan(r['strike']) and ~np.isnan(r['dip']):
                polarity = 1
                if 'polarity' in data_frame.columns and ~np.isnan(
                        r['polarity']):
                    polarity = r['polarity']
                self.add_strike_and_dip(pos, r['strike'], r['dip'],
                                        polarity=polarity)
            if 'azimuth' in data_frame.columns and 'dip' in \
                    data_frame.columns and \
                    ~np.isnan(r['azimuth']) and ~np.isnan(r['dip']):
                polarity = 1
                if 'polarity' in data_frame.columns and ~np.isnan(
                        r['polarity']):
                    polarity = r['polarity']
                self.add_plunge_and_plunge_dir(pos, r['dip'], r['azimuth'],
                                               polarity=polarity)

            if 'nx' in data_frame.columns and 'ny' in data_frame.columns and \
                    'nz' in data_frame.columns and \
                    ~np.isnan(r['nx']) and ~np.isnan(r['ny']) and ~np.isnan(
                r['nz']):
                self.add_planar_constraint(r[['X', 'Y', 'Z']],
                                           r[['nx', 'ny', 'nz']])

    def add_data(self, pos, strike=None, dip_dir=None, dip=None, dir=None,
                 val=None, plunge=None, plunge_dir=None, polarity=None):
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

    def add_strike_dip_and_value(self, pos, strike, dip, val, polarity=1):
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
        self.data.append(
            GPoint.from_strike_and_dip(pos, strike, dip, polarity))
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

    def add_strike_and_dip(self, pos, s, d, polarity=1, weight=1.):
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
        self.data[0].weight = weight
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
        self.data.append(
            GPoint.from_plunge_plunge_dir(pos, plunge, plunge_dir, polarity))
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

    def add_orthogonal_feature(self, feature, w=1., region=None):
        self.interpolator.add_gradient_orthogonal_constraint(
            self.interpolator.support.barycentre(),
            feature.evaluate_gradient(self.interpolator.support.barycentre()),
            w=w
        )

    def add_data_to_interpolator(self, constrained=False, force_constrained=False, **kwargs):
        """
        Iterates through the list of data and applies any faults active on the
        data in the order they are added

        Returns
        -------

        """
        # first move the data for the fault
        logger.info("Adding %i faults to %s" % (len(self.faults), self.name))
        data = copy.deepcopy(self.data)
        # convert data locations to numpy array and then update
        locations = self.get_data_locations()
        for f in self.faults:
            locations = f.apply_to_points(locations)
        i = 0
        for d in data:
            d.pos = locations[i, :]
            i += 1

        # Now check whether there are enough constraints for the
        # interpolator to be able to solve
        # we need at least 2 different value points or a single norm
        # constraint. If there are not enough
        # try converting grad to norms, if still not enough send user an error
        vals = []
        for d in data:
            if d.type == "GPoint" and d.norm == True:
                constrained = True
                break
            if d.type == 'IPoint':
                vals.append(d.val)
        if len(np.unique(vals)) > 1:
            constrained = True
        if not constrained or force_constrained:
            for d in data:
                if d.type == "GPoint":
                    d.norm = True
                    logger.debug(
                        "Setting gradient points to norm constraints")
                    constrained = True
        if not constrained:
            logger.error("Not enough constraints for scalar field add more")
        for d in data:
            self.interpolator.add_data(d)

        self.data_added = True

    def get_value_constraints(self):
        """
        Get the value constraints for this geological feature

        Returns
        -------
        numpy array
        """
        points = np.zeros((len(self.data), 4))  # array
        c = 0
        for d in self.data:
            if d.type == 'IPoint':
                points[c, :3] = d.pos
                points[c, 4] = d.val
                c += 1
        return points[:c, :]

    def get_gradient_constraints(self):
        """
        Get the gradient direction constraints

        Returns
        -------
        numpy array
        """
        points = np.zeros((len(self.data), 6))  # array
        c = 0
        for d in self.data:
            if d.type == 'GPoint':
                points[c, :3] = d.pos
                points[c, 3:] = d.vec
                c += 1
        return points[:c, :]

    def get_tangent_constraints(self):
        """

        Returns
        -------
        numpy array
        """
        points = np.zeros((len(self.data), 6))  # array
        c = 0
        for d in self.data:
            if d.type == 'TPoint':
                points[c, :3] = d.pos
                points[c, 3:] = d.vec
                c += 1
        return points[:c, :]

    def get_norm_constraints(self):
        """
        Get the gradient norm constraints

        Returns
        -------
        numpy array
        """
        points = np.zeros((len(self.data), 6))  # array
        c = 0
        for d in self.data:
            if d.type == 'NPoint':
                points[c, :3] = d.pos
                points[c, 3:] = d.vec
                c += 1
        return points[:c, :]

    def get_data_locations(self):
        """
        Get only the location for all data points

        Returns
        -------

        """
        points = np.zeros((len(self.data), 3))  # array
        c = 0
        for d in self.data:
            points[c, :] = d.pos
            c += 1
        return points[:c, :]

    def update_data_locations(self, locations):
        i = 0
        for d in self.data:
            d.pos = locations[i, :]
            i += 1

    def build(self, fold=None, fold_weights=None, data_region=None, **kwargs):
        """
        Runs the interpolation and builds the geological feature

        Parameters
        ----------
        fold : FoldEvent
        fold_weights : dict
        data_region : double <1
            If not none adds a region around the data points to the interpolation
            with data_region as a buffer
        kwargs

        Returns
        -------

        """

        if data_region is not None:
            xyz = self.get_data_locations()
            bb, region = get_data_axis_aligned_bounding_box(xyz, data_region)
            self.interpolator.set_region(region=region)
        if not self.data_added:
            self.add_data_to_interpolator(**kwargs)

        # moving this to init because it needs to be done before constraints
        # are added?
        if fold is not None:
            logger.info("Adding fold to %s" % self.name)
            self.interpolator.fold = fold
            # if we have fold weights use those, otherwise just use default
            if fold_weights is None:
                self.interpolator.add_fold_constraints()
            else:
                self.interpolator.add_fold_constraints(fold_weights)
            if 'cgw' not in kwargs:
                kwargs['cgw'] = 0.

        self.interpolator.setup_interpolator(**kwargs)
        self.interpolator.solve_system(**kwargs)
        return GeologicalFeature(self.name,
                                 self.interpolator,
                                 builder=self, data=self.data,
                                 region=self.region,
                                 faults=self.faults
                                 )
