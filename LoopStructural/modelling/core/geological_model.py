from LoopStructural.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from LoopStructural.interpolators.discrete_fold_interpolator import DiscreteFoldInterpolator as DFI
from LoopStructural.supports.tet_mesh import TetMesh
from LoopStructural.modelling.features.geological_feature import GeologicalFeatureInterpolator
from LoopStructural.interpolators.finite_difference_interpolator import FiniteDifferenceInterpolator as FDI
from LoopStructural.supports.structured_grid import StructuredGrid
from LoopStructural.modelling.features.structural_frame import StructuralFrameBuilder
from LoopStructural.modelling.fault.fault_segment import FaultSegment
from LoopStructural.modelling.fold.foldframe import FoldFrame
from LoopStructural.modelling.fold.fold import FoldEvent
from LoopStructural.modelling.fold.svariogram import SVariogram
from LoopStructural.modelling.fold.fold_rotation_angle_feature import fourier_series
from LoopStructural.modelling.features import RegionFeature
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)


class GeologicalModel:
    """
    A geological model is the recipe for building a 3D model and includes the rescaling
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
         type, X, Y, Z, nx, ny, nz, val, strike, dip, dip_dir, plunge, plunge_dir, azimuth

        Returns
        -------

        Note
        ----
        Type can be any unique identifier for the feature the data point 'eg' 'S0', 'S2', 'F1_axis'
        it is then used by the create functions to get the correct data
        """
        if type(data) != pd.DataFrame:
            logger.warning("Data is not a pandas data frame, trying to read data frame from csv")
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

    def get_interpolator(self, interpolatortype='PLI', nelements=5e5, buffer=0.02, **kwargs):
        """
        Returns an interpolator given the arguments, also constructs a support for a discrete
        interpolator
        Parameters
        ----------
        interpolatortype - string
            define the interpolator type
        nelements - int
            number of elements in the interpolator
        buffer - double or numpy array 3x1
            value(s) between 0,1 specifying the buffer around the bounding box

        kwargs - no kwargs used, this just catches any additional arguments

        Returns
        -------

        """
        # get an interpolator for 
        interpolator = None
        bb = np.copy(self.bounding_box)
        # add a buffer to the interpolation domain, this is necessary for faults but also generally a good
        # idea to avoid boundary problems
        bb[0, :] -= buffer  # *(bb[1,:]-bb[0,:])
        bb[1, :] += buffer  # *(bb[1,:]-bb[0,:])
        if interpolatortype == "PLI":
            mesh = TetMesh()
            mesh.setup_mesh(bb, n_tetra=nelements, )
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
            grid = StructuredGrid(origin=bb[0, :], nsteps=nsteps, step_vector=step_vector)
            return FDI(grid)
        if interpolatortype == "DFI":  # "fold" in kwargs:
            mesh = TetMesh()
            mesh.setup_mesh(bb, n_tetra=nelements, )
            return DFI(mesh, kwargs['fold'])
        logger.warning("No interpolator")
        return interpolator

    def create_and_add_foliation(self, series_surface_data, **kwargs):
        """
        Parameters
        ----------
        series_surface_data - string corresponding to the type in the data frame column
        kwargs

        Returns
        -------
        results dict

        """
        interpolator = self.get_interpolator(**kwargs)
        series_builder = GeologicalFeatureInterpolator(interpolator, name=series_surface_data, **kwargs)
        # add data
        series_data = self.data[self.data['type'] == series_surface_data]
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
            if f.type is 'unconformity':
                series_feature.add_region(lambda pos: f.evaluate_value(pos) < 0)
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
        # create fault frame
        interpolator = self.get_interpolator(**kwargs)
        #
        fold_frame_builder = StructuralFrameBuilder(interpolator, name=foldframe_data, **kwargs)
        # add data
        fold_frame_data = self.data[self.data['type'] == foldframe_data]
        fold_frame_builder.add_data_from_data_frame(fold_frame_data)
        # if there is no fault slip data then we could find the strike of the fault and build
        # the second coordinate
        fold_frame = fold_frame_builder.build(frame=FoldFrame, **kwargs)
        for f in reversed(self.features):
            if f.type == 'unconformity':
                fold_frame.add_region(lambda pos: f.evaluate_value(pos) <= 0)
                break
        fold_frame.type = 'structuralframe'
        self.features.append(fold_frame)
        return fold_frame

    def create_and_add_folded_foliation(self, foliation_data, fold_frame, **kwargs):
        """
        Create a folded foliation field from data and a fold frame
        Parameters
        ----------
        foliation_data
        fold_frame
        kwargs

        Returns
        -------
        dict
        """
        result = {}

        fold = FoldEvent(fold_frame)
        fold_interpolator = self.get_interpolator("DFI", fold=fold, **kwargs)
        series_builder = GeologicalFeatureInterpolator(
            interpolator=fold_interpolator,
            name=foliation_data)

        series_builder.add_data_from_data_frame(self.data[self.data['type'] == foliation_data])
        series_builder.add_data_to_interpolator(True)
        if "fold_axis" in kwargs:
            fold.fold_axis = kwargs['fold_axis']
        if "av_fold_axis" in kwargs:
            pass
        if "fold_axis" not in kwargs:
            far, fad = fold_frame.calculate_fold_axis_rotation(
                series_builder)
            axis_wl = kwargs.get("axis_wavelength", None)
            if axis_wl is None:
                axis_svariogram = SVariogram(fad, far)
                axis_wl = axis_svariogram.find_wavelengths()
            guess = np.zeros(4)
            guess[3] = axis_wl
            if len(far) < len(guess):
                logger.warning("Not enough data to fit curve")
                fold.fold_axis_rotation = lambda x: 0
            else:
                popt, pcov = curve_fit(fourier_series, fad, np.tan(np.rad2deg(far)), guess)
                fold.fold_axis_rotation = lambda x: np.rad2deg(
                    np.arctan(fourier_series(x, popt[0], popt[1], popt[2], popt[3])))
            result['axis_direction'] = fad
            result['axis_rotation'] = far
            result['axis_svario'] = axis_svariogram
        flr, s = fold_frame.calculate_fold_limb_rotation(series_builder, axis=fold.get_fold_axis_orientation)
        result['limb_rotation'] = flr
        result['foliation'] = s
        limb_wl = kwargs.get("limb_wl", None)
        if limb_wl is None:
            limb_svariogram = SVariogram(s, flr)
            limb_wl = limb_svariogram.find_wavelengths()
        guess = np.zeros(4)
        guess[3] = limb_wl
        if len(flr) < len(guess):
            logger.warning("Not enough data to fit curve")
            fold.fold_limb_rotation = lambda x: 0
        else:
            popt, pcov = curve_fit(fourier_series, s, np.tan(np.deg2rad(flr)), guess)
            fold.fold_limb_rotation = lambda x: np.rad2deg(
                np.arctan(fourier_series(x, popt[0], popt[1], popt[2], popt[3])))
        fold_weights = kwargs.get('fold_weights', None)

        if fold_weights is None:
            fold_weights = {}
            fold_weights['fold_orientation'] = 10.  # reference values?
            fold_weights['fold_axis'] = 10.  # reference values?
            fold_weights['fold_normalisation'] = 1.  # reference values?
            fold_weights['fold_regularisation'] = .100  # reference values?
            kwargs['fold_weights'] = fold_weights

            for f in reversed(self.features):
                if f.type == 'fault':
                    series_builder.add_fault(f)
                if f.type == 'unconformity':
                    break
        # build feature
        kwargs['cgw'] = 0.
        kwargs['fold'] = fold
        series_feature = series_builder.build(**kwargs)
        series_feature.type = 'series'
        # see if any unconformities are above this feature if so add region
        for f in reversed(self.features):
            if f.type is 'unconformity':
                series_feature.add_region(lambda pos: f.evaluate_value(pos) < 0)
                break
        self.features.append(series_feature)
        result['feature'] = series_feature
        result['fold'] = fold
        return result

    def create_and_add_folded_fold_frame(self, foliation_data, fold_frame, **kwargs):

        pass

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
        unconformity_feature_builder = GeologicalFeatureInterpolator(interpolator, name=unconformity_surface_data)
        # add data
        unconformity_data = self.data[self.data['type'] == unconformity_surface_data]

        unconformity_feature_builder.add_data_from_data_frame(unconformity_data)
        # look through existing features if there is a fault before an unconformity
        # then add to the feature, once we get to an unconformity stop
        for f in reversed(self.features):
            if f.type == 'fault':
                unconformity_feature_builder.add_fault(f)
            if f.type == 'unconformity':
                break

        # build feature
        uc_feature = unconformity_feature_builder.build(**kwargs)
        uc_feature.type = 'unconformity'

        # iterate over existing features and add the unconformity as a region so the feature is only
        # evaluated where the unconformity is positive
        for f in self.features:
            f.add_region(lambda pos: uc_feature.evaluate_value(pos) >= 0)

        # see if any unconformities are above this feature if so add region
        for f in reversed(self.features):
            if f.type == 'unconformity':
                uc_feature.add_region(lambda pos: f.evaluate_value(pos) <= 0)
                break

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
        fault_frame_builder = StructuralFrameBuilder(interpolator, name=fault_surface_data, **kwargs)
        # add data
        fault_frame_data = self.data[self.data['type'] == fault_surface_data]
        fault_frame_builder.add_data_from_data_frame(fault_frame_data)
        # if there is no fault slip data then we could find the strike of the fault and build
        # the second coordinate
        # if we add a region to the fault then the fault operator doesn't work but for visualisation
        # we want to add a region!

        if 'splayregion' in kwargs and 'splay' in kwargs:
            result['splayregionfeature'] = RegionFeature(kwargs['splayregion'])
            # apply splay to all parts of fault frame
            for i in range(3):
                # work out the values of the nodes where we want hard constraints
                idc = np.arange(0, interpolator.support.n_nodes)[
                    kwargs['splayregion'](interpolator.support.nodes)]
                val = kwargs['splay'][i].evaluate_value(
                    interpolator.support.nodes[kwargs['splayregion'](interpolator.support.nodes), :])
                mask = ~np.isnan(val)
                fault_frame_builder[i].interpolator.add_equality_constraints(idc[mask], val[mask])
        # check if any faults exist in the stack

        for f in reversed(self.features):
            if f.type == 'fault':
                fault_frame_builder[0].add_fault(f)
                fault_frame_builder[1].add_fault(f)
                fault_frame_builder[2].add_fault(f)
            if f.type == 'unconformity':
                break

        fault_frame = fault_frame_builder.build(**kwargs)
        if 'abut' in kwargs:
            fault_frame[0].add_region(lambda pos: kwargs['abut'].evaluate(pos))

        fault = FaultSegment(fault_frame, displacement=displacement_scaled, **kwargs)
        for f in reversed(self.features):
            if f.type is 'unconformity':
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
