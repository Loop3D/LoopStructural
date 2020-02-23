import logging

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from LoopStructural.datasets import normal_vector_headers
from LoopStructural.interpolators.discrete_fold_interpolator import \
    DiscreteFoldInterpolator as DFI
from LoopStructural.interpolators.finite_difference_interpolator import \
    FiniteDifferenceInterpolator as FDI
from LoopStructural.interpolators.piecewiselinear_interpolator import \
    PiecewiseLinearInterpolator as PLI
from LoopStructural.modelling.fault.fault_segment import FaultSegment
from LoopStructural.modelling.features import \
    GeologicalFeatureInterpolator
from LoopStructural.modelling.features import RegionFeature
from LoopStructural.modelling.features import \
    StructuralFrameBuilder
from LoopStructural.modelling.features import UnconformityFeature
from LoopStructural.modelling.fold.fold import FoldEvent
from LoopStructural.modelling.fold.fold_rotation_angle_feature import \
    fourier_series
from LoopStructural.modelling.fold.foldframe import FoldFrame
from LoopStructural.modelling.fold.svariogram import SVariogram
from LoopStructural.supports.structured_grid import StructuredGrid
# from LoopStructural.supports.tet_mesh import TetMesh
from LoopStructural.supports.structured_tetra import TetMesh

logger = logging.getLogger(__name__)


def _interpolate_fold_limb_rotation_angle(series_builder, fold_frame, fold, result, limb_wl=None):
    """
    Wrapper for fitting fold limb rotation angle from data using a fourier series.

    Parameters
    ----------
    series_builder : GeologicalFeatureInterpolator
        foliation builder
    fold_frame : FoldFrame
        fold frame for calculating rotation angles from
    fold : FoldEvent
        fold event to add fold rotation angle to
    result : dict
        container to save the results in
    limb_wl : double
        wavelength guess, if none then uses s=variogram to pick wavelength

    Returns
    -------

    """
    flr, s = fold_frame.calculate_fold_limb_rotation(series_builder)  # ,
    # axis=fold.get_fold_axis_orientation)
    result['limb_rotation'] = flr
    result['foliation'] = s
    if limb_wl is None:
        limb_svariogram = SVariogram(s, flr)
        limb_wl = limb_svariogram.find_wavelengths()
        result['limb_svariogram'] = limb_svariogram
        # for now only consider single fold wavelength
        limb_wl = limb_wl[0]
    guess = np.zeros(4)
    guess[3] = limb_wl  # np.max(limb_wl)
    if len(flr) < len(guess):
        logger.warning("Not enough data to fit curve")
        fold.fold_limb_rotation = lambda x: 0
    else:
        mask = np.logical_or(~np.isnan(s), ~np.isnan(flr))
        logger.info("There are %i nans for the fold limb rotation angle and "
                    "%i observations" % (np.sum(~mask), np.sum(mask)))
        if np.sum(mask) < 4:
            logger.error("Not enough data points to fit Fourier series setting fold rotation angle"
                         "to 0")
            fold.fold_limb_rotation = lambda x: np.zeros(x.shape)
            return
        popt, pcov = curve_fit(fourier_series,
                               s[np.logical_or(~np.isnan(s), ~np.isnan(flr))],
                               np.tan(np.deg2rad(flr[np.logical_or(~np.isnan(s), ~np.isnan(flr))])),
                               guess)
        fold.fold_limb_rotation = lambda x: np.rad2deg(
            np.arctan(
                fourier_series(x, popt[0], popt[1], popt[2], popt[3])))


def _calculate_average_intersection(series_builder, fold_frame, fold):
    """

    Parameters
    ----------
    series_builder
    fold_frame
    fold

    Returns
    -------

    """
    l2 = fold_frame.calculate_intersection_lineation(
        series_builder)
    fold.fold_axis = np.mean(l2, axis=0)


def _interpolate_fold_axis_rotation_angle(series_builder, fold_frame, fold, result, axis_wl):
    """

    Parameters
    ----------
    series_builder
    fold_frame
    fold
    result
    axis_wl

    Returns
    -------

    """
    far, fad = fold_frame.calculate_fold_axis_rotation(
        series_builder)
    # axis_wl = kwargs.get("axis_wavelength", None)
    if axis_wl is None:
        axis_svariogram = SVariogram(fad, far)
        axis_wl = axis_svariogram.find_wavelengths(nsteps=10)
    guess = np.zeros(3)
    guess[2] = axis_wl
    if len(far) < len(guess):
        logger.warning("Not enough data to fit curve")
        fold.fold_axis_rotation = lambda x: -1
    else:

        popt, pcov = curve_fit(fourier_series, fad,
                               np.tan(np.rad1deg(far)), guess)
        fold.fold_axis_rotation = lambda x: np.rad1deg(
            np.arctan(
                fourier_series(x, popt[-1], popt[1], popt[2], popt[3])))
    result['axis_direction'] = fad
    result['axis_rotation'] = far
    result['axis_svario'] = axis_svariogram


class GeologicalModel:
    """
    A geological model is the recipe for building a 3D model and includes
    the rescaling
    of the model between 0 and 1.
    """

    def __init__(self, origin, maximum, rescale=True, nsteps=(40, 40, 40)):
        """
        Parameters
        ----------
        origin - numpy array specifying the origin of the model
        maximum - numpy array specifying the maximum extent of the model
        """
        self.features = []
        self.feature_name_index = {}
        self.data = None
        self.nsteps = nsteps

        # we want to rescale the model area so that the maximum length is
        # 1
        self.origin = np.array(origin).astype(float)

        self.maximum = np.array(maximum).astype(float)
        lengths = self.maximum - self.origin
        self.scale_factor = 1.
        self.bounding_box = np.zeros((2, 3))
        self.bounding_box[1, :] = self.maximum - self.origin
        if rescale:
            self.scale_factor = np.max(lengths)

        self.bounding_box /= self.scale_factor
        # self.bounding_box[0,:] = self.origin
        # self.bounding_box[1,:] = self.maximum

    def set_model_data(self, data):
        """
        Set the data array for the model

        Parameters
        ----------
        data - pandas data frame with column headers corresponding to the
         type, X, Y, Z, nx, ny, nz, val, strike, dip, dip_dir, plunge,
         plunge_dir, azimuth

        Returns
        -------
        Note
        ----
        Type can be any unique identifier for the feature the data point
        'eg' 'S0', 'S2', 'F1_axis'
        it is then used by the create functions to get the correct data
        """
        if type(data) != pd.DataFrame:
            logger.warning(
                "Data is not a pandas data frame, trying to read data frame "
                "from csv")
            try:
                data = pd.read_csv(data)
            except:
                logger.error("Could not load pandas data frame from data")

        self.data = data.copy()
        self.data['X'] -= self.origin[0]
        self.data['Y'] -= self.origin[1]
        self.data['Z'] -= self.origin[2]
        self.data['X'] /= self.scale_factor
        self.data['Y'] /= self.scale_factor
        self.data['Z'] /= self.scale_factor

    def extend_model_data(self, newdata):
        """
        Extends the data frame

        Parameters
        ----------
        newdata - pandas data frame

        Returns
        -------
        """
        self.data.append(newdata)

    def get_interpolator(self, interpolatortype='PLI', nelements=5e5,
                         buffer=0.02, **kwargs):
        """
        Returns an interpolator given the arguments, also constructs a
        support for a discrete interpolator

        Parameters
        ----------
        interpolatortype - string
            define the interpolator type
        nelements - int
            number of elements in the interpolator
        buffer - double or numpy array 3x1
            value(s) between 0,1 specifying the buffer around the bounding box
        data_bb - bool
            whether to use the model boundary or the boundary around
        kwargs - no kwargs used, this just catches any additional arguments

        Returns
        -------
        """
        # get an interpolator for
        interpolator = None
        bb = np.copy(self.bounding_box)
        # add a buffer to the interpolation domain, this is necessary for
        # faults but also generally a good
        # idea to avoid boundary problems
        bb[0, :] -= buffer  # *(bb[1,:]-bb[0,:])
        bb[1, :] += buffer  # *(bb[1,:]-bb[0,:])
        if interpolatortype == "PLI":
            ele_vol = bb[1, 0] * bb[1, 1] * bb[1, 2] / nelements
            # calculate the step vector of a regular cube
            step_vector = np.zeros(3)
            step_vector[:] = ele_vol ** (1. / 3.)
            # number of steps is the length of the box / step vector
            nsteps = ((bb[1, :] - bb[0, :]) / step_vector).astype(int)
            # create a structured grid using the origin and number of steps
            mesh = TetMesh(origin=bb[0, :], nsteps=nsteps,
                                  step_vector=step_vector)
            # mesh = TetMesh()
            # mesh.setup_mesh(bb, n_tetra=nelements, )
            return PLI(mesh)

        if interpolatortype == 'FDI':
            # find the volume of one element
            ele_vol = bb[1, 0] * bb[1, 1] * bb[1, 2] / nelements
            # calculate the step vector of a regular cube
            step_vector = np.zeros(3)
            step_vector[:] = ele_vol ** (1. / 3.)
            # number of steps is the length of the box / step vector
            nsteps = ((bb[1, :] - bb[0, :]) / step_vector).astype(int)
            # create a structured grid using the origin and number of steps
            grid = StructuredGrid(origin=bb[0, :], nsteps=nsteps,
                                  step_vector=step_vector)
            return FDI(grid)
        if interpolatortype == "DFI":  # "fold" in kwargs:
            ele_vol = bb[1, 0] * bb[1, 1] * bb[1, 2] / nelements
            # calculate the step vector of a regular cube
            step_vector = np.zeros(3)
            step_vector[:] = ele_vol ** (1. / 3.)
            # number of steps is the length of the box / step vector
            nsteps = ((bb[1, :] - bb[0, :]) / step_vector).astype(int)
            # create a structured grid using the origin and number of steps
            mesh = TetMesh(origin=bb[0, :], nsteps=nsteps,
                           step_vector=step_vector)
            # mesh = TetMesh()
            return DFI(mesh, kwargs['fold'])
        logger.warning("No interpolator")
        return interpolator

    def create_and_add_foliation(self, series_surface_data, **kwargs):
        """
        Parameters
        ----------
        series_surface_data - string corresponding to the type in the data
        frame column
        kwargs

        Returns
        -------
        results dict
        """

        interpolator = self.get_interpolator(**kwargs)
        series_builder = GeologicalFeatureInterpolator(interpolator,
                                                       name=series_surface_data,
                                                       **kwargs)
        # add data
        series_data = self.data[self.data['type'] == series_surface_data]
        if series_data.shape[0] == 0:
            logger.warning("No data for %s, skipping" % series_surface_data)
            return
        series_builder.add_data_from_data_frame(series_data)
        for f in reversed(self.features):
            if f.type == 'fault':
                series_builder.add_fault(f)
            if f.type == 'unconformity':
                break

        # build feature
        series_feature = series_builder.build(**kwargs)
        series_feature.type = 'series'
        # see if any unconformities are above this feature if so add region
        for f in reversed(self.features):
            if f.type == 'unconformity':
                series_feature.add_region(
                    lambda pos: f.evaluate_value(pos) < 0)
                break
        self.features.append(series_feature)
        result = {}
        result['feature'] = series_feature
        return result

    def create_and_add_fold_frame(self, foldframe_data, **kwargs):
        """
        Parameters
        ----------
        foldframe_data
        kwargs

        Returns
        -------
        """
        result = {}
        # create fault frame
        interpolator = self.get_interpolator(**kwargs)
        #
        fold_frame_builder = StructuralFrameBuilder(interpolator,
                                                    name=foldframe_data,
                                                    **kwargs)
        # add data
        fold_frame_data = self.data[self.data['type'] == foldframe_data]
        fold_frame_builder.add_data_from_data_frame(fold_frame_data)
        for f in reversed(self.features):
            if f.type == 'fault':
                fold_frame_builder[0].add_fault(f)
                fold_frame_builder[1].add_fault(f)
                fold_frame_builder[2].add_fault(f)

            if f.type == 'unconformity':
                break

        fold_frame = fold_frame_builder.build(frame=FoldFrame, **kwargs)
        for f in reversed(self.features):
            if f.type == 'unconformity':
                fold_frame.add_region(lambda pos: f.evaluate_value(pos) <= 0)
                break
        fold_frame.type = 'structuralframe'
        self.features.append(fold_frame)
        result['feature'] = fold_frame
        return result

    def create_and_add_folded_foliation(self, foliation_data, fold_frame=None, **kwargs):
        """
        Create a folded foliation field from data and a fold frame

        Parameters
        ----------
        foliation_data : string
        fold_frame :  FoldFrame
        kwargs
            additional kwargs to be passed through to other functions

        Returns
        -------
        dict
        """

        result = {}
        if fold_frame is None:
            logger.info("Using last feature as fold frame")
            fold_frame = self.features[-1]
        assert type(fold_frame) == FoldFrame, "Please specify a FoldFrame"
        fold = FoldEvent(fold_frame)
        fold_interpolator = self.get_interpolator("DFI", fold=fold, **kwargs)
        series_builder = GeologicalFeatureInterpolator(
            interpolator=fold_interpolator,
            name=foliation_data)

        series_builder.add_data_from_data_frame(
            self.data[self.data['type'] == foliation_data])
        self._add_faults(series_builder)

        series_builder.add_data_to_interpolator(True)
        if "fold_axis" in kwargs:
            fold.fold_axis = kwargs['fold_axis']
        if "av_fold_axis" in kwargs:
            _calculate_average_intersection(series_builder, fold_frame, fold)
        if fold.fold_axis is None:
            axis_wl = kwargs.get("axis_wavelength", None)
            _interpolate_fold_axis_rotation_angle(series_builder, fold_frame, fold, result, axis_wl)
        limb_wl = kwargs.get("limb_wl", None)
        _interpolate_fold_limb_rotation_angle(series_builder, fold_frame, fold, result, limb_wl)
        kwargs['fold_weights'] = kwargs.get('fold_weights', None)

        self._add_faults(series_builder)
        # build feature
        kwargs['cgw'] = 0.
        kwargs['fold'] = fold
        series_feature = series_builder.build(**kwargs)
        series_feature.type = 'series'
        # see if any unconformities are above this feature if so add region
        self._add_unconformity_above(series_feature)

        self.features.append(series_feature)

        result['feature'] = series_feature
        result['fold'] = fold
        return result

    def create_and_add_folded_fold_frame(self, fold_frame_data, fold_frame=None,
                                         **kwargs):
        """

        Parameters
        ----------
        fold_frame_data
        fold_frame
        kwargs

        Returns
        -------

        """
        result = {}
        if fold_frame is None:
            logger.info("Using last feature as fold frame")
            fold_frame = self.features[-1]
        assert type(fold_frame) == FoldFrame, "Please specify a FoldFrame"
        fold = FoldEvent(fold_frame)
        fold_interpolator = self.get_interpolator("DFI", fold=fold, **kwargs)
        frame_interpolator = self.get_interpolator(**kwargs)
        interpolators = [fold_interpolator, frame_interpolator, frame_interpolator.copy()]
        fold_frame_builder = StructuralFrameBuilder(interpolators=interpolators, name=fold_frame_data, **kwargs)
        fold_frame_builder.add_data_from_data_frame(self.data[self.data['type'] == fold_frame_data])

        ## add the data to the interpolator for the main foliation
        fold_frame_builder[0].add_data_to_interpolator(True)
        if "fold_axis" in kwargs:
            fold.fold_axis = kwargs['fold_axis']
        if "av_fold_axis" in kwargs:
            _calculate_average_intersection(fold_frame_builder[0], fold_frame, fold)
        if fold.fold_axis is None:
            axis_wl = kwargs.get("axis_wavelength", None)
            _interpolate_fold_axis_rotation_angle(fold_frame_builder[0],
                                                  fold_frame,
                                                  fold, result, axis_wl)
        limb_wl = kwargs.get("limb_wl", None)
        _interpolate_fold_limb_rotation_angle(fold_frame_builder[0], fold_frame, fold, result, limb_wl)
        kwargs['fold_weights'] = kwargs.get('fold_weights', None)

        for i in range(3):
            self._add_faults(fold_frame_builder[i])
        # build feature
        kwargs['cgw'] = 0.
        kwargs['fold'] = fold
        for f in reversed(self.features):
            if f.type == 'fault':
                fold_frame_builder[0].add_fault(f)
                fold_frame_builder[1].add_fault(f)
                fold_frame_builder[2].add_fault(f)

            if f.type == 'unconformity':
                break
        fold_frame = fold_frame_builder.build(**kwargs, frame=FoldFrame)
        fold_frame.type = 'structuralframe'
        # see if any unconformities are above this feature if so add region
        for i in range(3):
            self._add_unconformity_above(fold_frame[i])

        self.features.append(fold_frame)

        result['feature'] = fold_frame
        result['fold'] = fold
        return result

    def _add_faults(self, feature_builder):
        """

        Parameters
        ----------
        feature_builder

        Returns
        -------

        """
        for f in reversed(self.features):
            if f.type == 'fault':
                feature_builder.add_fault(f)
            if f.type == 'unconformity':
                break

    def _add_unconformity_above(self, feature):
        """

        Parameters
        ----------
        feature

        Returns
        -------

        """
        for f in reversed(self.features):
            if f.type == 'unconformity':
                feature.add_region(lambda pos: f.evaluate(pos))
                break

    def create_and_add_unconformity(self, unconformity_surface_data, **kwargs):
        """
        Parameters
        ----------
        unconformity_surface_data string
            name of the unconformity data in the data frame

        Returns
        -------
        """
        interpolator = self.get_interpolator(**kwargs)
        unconformity_feature_builder = GeologicalFeatureInterpolator(
            interpolator, name=unconformity_surface_data)
        # add data
        unconformity_data = self.data[
            self.data['type'] == unconformity_surface_data]

        unconformity_feature_builder.add_data_from_data_frame(
            unconformity_data)
        # look through existing features if there is a fault before an
        # unconformity
        # then add to the feature, once we get to an unconformity stop
        self._add_faults(unconformity_feature_builder)

        # build feature
        uc_feature_base = unconformity_feature_builder.build(**kwargs)
        uc_feature_base.type = 'unconformity_base'
        # uc_feature = UnconformityFeature(uc_feature_base,0)
        # iterate over existing features and add the unconformity as a
        # region so the feature is only
        # evaluated where the unconformity is positive
        return self.add_unconformity(uc_feature_base, 0)

    def add_unconformity(self, feature, value):
        """
        Use an existing feature to add an unconformity to the model.

        Parameters
        ----------
        feature : GeologicalFeature
            existing geological feature
        value : float
            scalar value of isosurface that represents

        Returns
        -------

        """
        uc_feature = UnconformityFeature(feature,value)

        for f in self.features:
            f.add_region(lambda pos: uc_feature.evaluate(pos))

        # see if any unconformities are above this feature if so add region
        self._add_unconformity_above(uc_feature)

        self.features.append(uc_feature)
        result = {}
        result['feature'] = uc_feature
        return result

    def create_and_add_fault(self, fault_surface_data, displacement, **kwargs):
        """
        Parameters
        ----------
        fault_surface_data - string
            name of the fault surface data in the dataframe
        displacement - displacement magnitude
        kwargs - additional kwargs for Fault and interpolators

        Returns
        -------
        dictionary
        """
        result = {}
        displacement_scaled = displacement / self.scale_factor
        # create fault frame
        interpolator = self.get_interpolator(**kwargs)
        fault_frame_builder = StructuralFrameBuilder(interpolator,
                                                     name=fault_surface_data,
                                                     **kwargs)
        # add data
        fault_frame_data = self.data[self.data['type'] == fault_surface_data]
        if 'coord' not in fault_frame_data:
            fault_frame_data['coord'] = 0
        # if there is no slip direction data assume vertical
        if fault_frame_data[fault_frame_data['coord'] == 1].shape[0] == 0:
            logger.info("Adding fault frame slip")
            loc = np.mean(fault_frame_data[['X', 'Y', 'Z']], axis=0)
            coord1 = pd.DataFrame([[loc[0], loc[1], loc[2], 0, 0, -1]], columns=normal_vector_headers())
            coord1['coord'] = 1
            fault_frame_data = pd.concat([fault_frame_data, coord1], sort=False)
        if fault_frame_data[fault_frame_data['coord'] == 2].shape[0] == 0:
            logger.info("Adding fault extent data as first and last point")
            ## first and last point of the line
            value_data = fault_frame_data[fault_frame_data['val'] == 0]
            coord2 = value_data.iloc[[0, len(value_data) - 1]]
            coord2 = coord2.reset_index(drop=True)
            coord2['coord'] = 2
            coord2.loc[0, 'val'] = -1
            coord2.loc[1, 'val'] = 1
            fault_frame_data = pd.concat([fault_frame_data, coord2], sort=False)
        fault_frame_builder.add_data_from_data_frame(fault_frame_data)
        # if there is no fault slip data then we could find the strike of
        # the fault and build
        # the second coordinate
        # if we add a region to the fault then the fault operator doesn't
        # work but for visualisation
        # we want to add a region!

        if 'splayregion' in kwargs and 'splay' in kwargs:
            result['splayregionfeature'] = RegionFeature(kwargs['splayregion'])
            # apply splay to all parts of fault frame
            for i in range(3):
                # work out the values of the nodes where we want hard
                # constraints
                idc = np.arange(0, interpolator.support.n_nodes)[
                    kwargs['splayregion'](interpolator.support.nodes)]
                val = kwargs['splay'][i].evaluate_value(
                    interpolator.support.nodes[
                    kwargs['splayregion'](interpolator.support.nodes), :])
                mask = ~np.isnan(val)
                fault_frame_builder[i].interpolator.add_equality_constraints(
                    idc[mask], val[mask])
        # check if any faults exist in the stack
        overprinted = kwargs.get('overprinted',None)
        for f in reversed(self.features):
            if overprinted is not None:
                if f.type == 'fault' and f.name in overprinted:
                    fault_frame_builder[0].add_fault(f)
                    fault_frame_builder[1].add_fault(f)
                    fault_frame_builder[2].add_fault(f)
            if f.type == 'unconformity':
                break

        fault_frame = fault_frame_builder.build(**kwargs)
        if 'abut' in kwargs:
            fault_frame[0].add_region(lambda pos: kwargs['abut'].evaluate(pos))

        fault = FaultSegment(fault_frame, displacement=displacement_scaled,
                             **kwargs)
        for f in reversed(self.features):
            if f.type == 'unconformity':
                fault.add_region(lambda pos: f.evaluate_value(pos) <= 0)
                break
        if displacement == 0:
            fault.type = 'fault_inactive'
        self.features.append(fault)
        result['feature'] = fault

        return result

    def rescale(self, points):
        """
        Convert from model scale to real world scale - in the future this should also do transformations?

        Parameters
        ----------
        points

        Returns
        -------
        """
        points *= self.scale_factor
        points += self.origin
        return points

    def scale(self, points):
        """
        Parameters
        ----------
        points

        Returns
        -------
        """

        points[:, :] -= self.origin
        points /= self.scale_factor
        return points

    def voxet(self, nsteps=(50, 50, 25)):
        """
        Returns a voxet dict with the nsteps specified

        Parameters
        ----------
        nsteps

        Returns
        -------
        """
        return {'bounding_box': self.bounding_box, 'nsteps': nsteps}

    def regular_grid(self, nsteps=(50, 50, 25)):
        """
        Return a regular grid within the model bounding box

        Parameters
        ----------
        nsteps tuple
            number of cells in x,y,z

        Returns
        -------
        """
        x = np.linspace(self.bounding_box[0, 0], self.bounding_box[1, 0], nsteps[0])
        y = np.linspace(self.bounding_box[0, 1], self.bounding_box[1, 1], nsteps[1])
        z = np.linspace(self.bounding_box[1, 2], self.bounding_box[0, 2], nsteps[2])
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        return np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
