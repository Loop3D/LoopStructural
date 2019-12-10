import logging

import numpy as np

logger = logging.getLogger(__name__)

from LoopStructural.modelling.core.geological_points import GPoint, IPoint, \
    TPoint
from LoopStructural.modelling.features import GeologicalFeature
from LoopStructural.supports.scalar_field import ScalarField


class GeologicalFeatureInterpolator:
    def __init__(self, interpolator, name = 'Feature', region = None, **kwargs):
        """
        A builder for a GeologicalFeature will link data to the interpolator
        and run the interpolation

        Parameters
        ----------
        interpolator - a GeologicalInterpolator
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
        if 'X' not in data_frame.columns or 'X' not in data_frame.columns or \
                'X' not in data_frame.columns:
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

    def add_strike_and_dip(self, pos, s, d, polarity=1, weight =1.):
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
        self.data[0].weight=weight
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
            np.arange(0, self.interpolator.support.n_elements),
            feature.evaluate_gradient(self.interpolator.support.barycentre),
            w=w
        )

    def add_data_to_interpolator(self, constrained=False):
        """
        Iterates through the list of data and applies any faults active on the
        data in the order they are added

        Returns
        -------

        """
        # first move the data for the fault
        for f in self.faults:
            f.apply_to_data(self.data)
        # Now check whether there are enough constraints for the
        # interpolator to be able to solve
        # we need at least 2 different value points or a single norm
        # constraint. If there are not enough
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
                    logger.warning(
                        "Setting gradient points to norm constraints")
                    constrained = True
        if not constrained:
            logger.error("Not enough constraints for scalar field add more")
        for d in self.data:
            self.interpolator.add_data(d)

    def build(self, fold = None, fold_weights = None, **kwargs):
        """
        Runs the interpolation and builds the geological feature

        Parameters
        ----------
        solver
        kwargs

        Returns
        -------

        """
        if not self.data_added:
            self.add_data_to_interpolator()
        # moving this to init because it needs to be done before constraints
        # are added?
        if fold is not None:

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

                                 ScalarField.from_interpolator(
                                     self.interpolator),
                                 builder=self, data=self.data,
                                 region=self.region,
                                 faults=self.faults
                                 )
