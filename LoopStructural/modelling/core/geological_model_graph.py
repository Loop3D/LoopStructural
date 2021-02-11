"""
Main entry point for creating a geological model
"""
import logging

import numpy as np
import pandas as pd
import networkx as nx

from LoopStructural.datasets import normal_vector_headers
from LoopStructural.interpolators.discrete_fold_interpolator import \
    DiscreteFoldInterpolator as DFI
from LoopStructural.interpolators.finite_difference_interpolator import \
    FiniteDifferenceInterpolator as FDI
from LoopStructural.interpolators.piecewiselinear_interpolator import \
    PiecewiseLinearInterpolator as PLI

try:
    from LoopStructural.interpolators.surfe_wrapper import \
        SurfeRBFInterpolator as Surfe

    surfe = True

except ImportError:
    surfe = False

from LoopStructural.interpolators.structured_grid import StructuredGrid
from LoopStructural.interpolators.structured_tetra import TetMesh
from LoopStructural.modelling.fault.fault_segment import FaultSegment
from LoopStructural.modelling.features import (GeologicalFeatureInterpolator,
                                               RegionFeature,
                                               StructuralFrameBuilder,
                                               UnconformityFeature)
from LoopStructural.modelling.fold import FoldRotationAngle
from LoopStructural.modelling.fold.fold import FoldEvent
from LoopStructural.modelling.fold.foldframe import FoldFrame
from LoopStructural.utils.exceptions import LoopBaseException
from LoopStructural.utils.helper import (all_heading, gradient_vec_names,
                                         strike_dip_vector)

from LoopStructural.utils import getLogger, log_to_file
logger = getLogger(__name__)
if not surfe:
    logger.warning("Cannot import Surfe")


def _calculate_average_intersection(series_builder, fold_frame, fold,
                                    **kwargs):
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


class GeologicalModel:
    """
    A geological model is the recipe for building a 3D model and  can include
    the rescaling of the model between 0 and 1.

    Attributes
    ----------
    features : list
        Contains all features youngest to oldest
    feature_name_index : dict
        maps feature name to the list index of the features
    data : pandas dataframe
        the dataframe used for building the geological model
    nsteps : tuple/np.array(3,dtype=int)
        the number of steps x,y,z to evaluate the model
    origin : tuple/np.array(3,dtype=doubles)
        the origin of the model box
    parameters : dict
        a dictionary tracking the parameters used to build the model
    scale_factor : double
        the scale factor used to rescale the model


    """
    def __init__(self, origin, maximum, rescale=True, nsteps=(40, 40, 40),
                 reuse_supports=False, logfile=None, loglevel='info'):
        """
        Parameters
        ----------
        origin : numpy array
            specifying the origin of the model
        maximum : numpy array
            specifying the maximum extent of the model
        rescale : bool
            whether to rescale the model to between 0/1

        Examples
        --------
        Demo data
        
        >>> from LoopStructural.datasets import load_claudius
        >>> from LoopStructural import GeologicalModel
        
        >>> data, bb = load_claudius()
        
        >>> model = GeologicalModel(bb[:,0],bb[:,1] 
        >>> model.set_model_data(data)
        >>> model.create_and_add_foliation('strati')
        
        >>> y = np.linspace(model.bounding_box[0, 1], model.bounding_box[1, 1],
                        nsteps[1])
        >>> z = np.linspace(model.bounding_box[1, 2], model.bounding_box[0, 2],
                        nsteps[2])
        >>> xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        >>> xyz = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        >>> model.evaluate_feature_value('strati',xyz,scale=False)


        """
        if logfile:
            self.logfile = logfile
            log_to_file(logfile,loglevel)
        
        logger.info('Initialising geological model')
        self.features = []
        self.feature_name_index = {}
        self.data = None
        self.nsteps = nsteps
        self._str = 'Instance of LoopStructural.GeologicalModel \n'
        self._str += '------------------------------------------ \n'
        # we want to rescale the model area so that the maximum length is
        # 1
        self.origin = np.array(origin).astype(float)
        originstr = 'Model origin: {} {} {}'.format(self.origin[0],self.origin[1],self.origin[2])
        logger.info(originstr)
        self._str+=originstr+'\n'
        self.maximum = np.array(maximum).astype(float)
        maximumstr = 'Model maximum: {} {} {}'.format(self.maximum[0],self.maximum[1],self.maximum[2])
        logger.info(maximumstr)
        self._str+=maximumstr+'\n'

        lengths = self.maximum - self.origin
        self.scale_factor = 1.
        self.bounding_box = np.zeros((2, 3))
        self.bounding_box[1, :] = self.maximum - self.origin
        if rescale:
            self.scale_factor = np.max(lengths)
            logger.info('Rescaling model using scale factor {}'.format(self.scale_factor))
        self._str+='Model rescale factor: {} \n'.format(self.scale_factor)
        self._str+='The model contains {} GeologicalFeatures \n'.format(len(self.features))
        self._str+=''
        self._str += '------------------------------------------ \n'
        self._str += ''
        self.bounding_box /= self.scale_factor
        self.support = {}
        self.reuse_supports = reuse_supports
        if self.reuse_supports:
            logger.warning("Supports are shared between geological features \n"
                                "this may cause unexpected behaviour and should only\n"
                                "be use by advanced users")
        logger.info('Reusing interpolation supports: {}'.format(self.reuse_supports))
        self.stratigraphic_column = None
        self.feature_graph = nx.DiGraph()
        self.parameters = {'features': [], 'model': {'bounding_box': self.origin.tolist() + self.maximum.tolist(),
                                                     'rescale': rescale,
                                                     'nsteps': nsteps,
                                                     'reuse_supports': reuse_supports}}
        
    def __str__(self):
        return self._str

    def _ipython_key_completions_(self):
        return self.feature_name_index.keys()
        
    @classmethod
    def from_map2loop_directory(cls, m2l_directory,**kwargs):
        """Alternate constructor for a geological model using m2l output

        Uses the information saved in the map2loop files to build a geological model.
        You can specify kwargs for building foliation using foliation_params and for 
        faults using fault_params.  faults is a flag that allows for the faults to be skipped.

        Parameters
        ----------
        m2l_directory : string
            path to map2loop directory

        Returns
        -------
        (GeologicalModel, dict)
            the created geological model and a dictionary of the map2loop data
        """
        from LoopStructural.utils import build_model, process_map2loop
        logger.info('LoopStructural model initialised from m2l directory: {}'.format(m2l_directory))
        m2lflags = kwargs.pop('m2lflags',{})
        m2l_data = process_map2loop(m2l_directory,m2lflags)
        return build_model(m2l_data,**kwargs), m2l_data

    @classmethod
    def from_file(cls, file):
        """Load a geological model from file

        Parameters
        ----------
        file : string
            path to the file

        Returns
        -------
        GeologicalModel
            the geological model object
        """        
        try:
            import dill as pickle
        except ImportError:
            logger.error("Cannot import from file, dill not installed")
            return None
        model = pickle.load(open(file,'rb'))
        if type(model) == GeologicalModel:
            logger.info('GeologicalModel initialised from file')
            return model
        else:
            logger.error('{} does not contain a geological model'.format(file))
            return None
        
    def __getitem__(self, feature_name):
        """Accessor for feature in features using feature_name_index

        Parameters
        ----------
        feature_name : string
            name of the feature to return
        """
        return self.get_feature_by_name(feature_name)
    
    def feature_names(self):
        return self.feature_name_index.keys()

    def fault_names(self):
        pass

    def check_inialisation(self):
        if self.data is None:
            logger.error("Data not associated with GeologicalModel. Run set_data")
            return False
        
    def to_file(self, file):
        """Save a model to a pickle file requires dill

        Parameters
        ----------
        file : string
            path to file location
        """        
        try:
            import dill as pickle
        except ImportError:
            logger.error("Cannot write to file, dill not installed \n"
                        "pip install dill")
            return
        try:
            logger.info('Writing GeologicalModel to: {}'.format(file))
            pickle.dump(self,open(file,'wb'))
        except pickle.PicklingError:
            logger.error('Error saving file')

    def _add_feature(self, feature):
        """
        Add a feature to the model stack

        Parameters
        ----------
        feature : GeologicalFeature
            the geological feature to add

        """

        if feature.name in self.feature_name_index:
            logger.info("Feature %s already exists at %i, overwriting" %
                        (feature.name, self.feature_name_index[feature.name]))
            self.features[self.feature_name_index[feature.name]] = feature
        else:
            self._str += 'GeologicalFeature: {} of type - {} \n'.format(feature.name,feature.type)
            self.features.append(feature)
            self.feature_name_index[feature.name] = len(self.features) - 1
            logger.info("Adding %s to model at location %i" % (
                feature.name, len(self.features)))
        # self._add_domain_fault_above(feature)
        # self._add_unconformity_above(feature)
        feature.set_model(self)
        
    def _add_faults(self, feature_builder, features=None):
        """Adds all existing faults to a geological feature builder 
        
        Parameters
        ----------
        feature_builder : GeologicalFeatureInterpolator/StructuralFrameBuilder
            The feature buider to add the faults to
        features : list, optional
            A specific list of features rather than all features in the model
        Returns
        -------

        """
        if features is None:
            features = self.features
        for f in reversed(features):
            if f.type == 'fault':
                feature_builder.add_fault(f)
            # if f.type == 'unconformity':
            #     break
                
    def data_for_feature(self,feature_name):
        return self.data.loc[self.data['feature_name'] == feature_name,:]
        
    def set_model_data(self, data):
        """
        Set the data array for the model

        Parameters
        ----------
        data : pandas data frame
            with column headers corresponding to the
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
        logger.info('Adding data to GeologicalModel with {} data points'.format(len(data)))
        self.data = data.copy()
        self.data['X'] -= self.origin[0]
        self.data['Y'] -= self.origin[1]
        self.data['Z'] -= self.origin[2]
        self.data['X'] /= self.scale_factor
        self.data['Y'] /= self.scale_factor
        self.data['Z'] /= self.scale_factor
        if 'type' in self.data:
            logger.warning("'type' is depreciated replace with 'feature_name' \n")
            self.data.rename(columns={'type':'feature_name'},inplace=True)
        for h in all_heading():
            if h not in self.data:
                self.data[h] = np.nan
                if h == 'w':
                    self.data[h] = 1.
                if h == 'coord':
                    self.data[h] = 0
        
        if 'strike' in self.data and 'dip' in self.data:
            logger.info('Converting strike and dip to vectors')
            mask = np.all(~np.isnan(self.data.loc[:, ['strike', 'dip']]),
                          axis=1)
            self.data.loc[mask, gradient_vec_names()] = strike_dip_vector(
                self.data.loc[mask, 'strike'], self.data.loc[mask, 'dip'])
            self.data.drop(['strike', 'dip'], axis=1, inplace=True)
        #     self.data.loc
        # if 'nx' in self.data and 'ny' in self.data and 'nz' in self.data:
        #     mask = np.all(~np.isnan(self.data.loc[:, ['nx', 'ny','nz']]),
        #                   axis=1)
        #     self.data.loc[mask,['nx', 'ny','nz']] /= self.scale_factor
    
    def set_stratigraphic_column(self, stratigraphic_column,cmap='tab20'):
        """
        Adds a stratigraphic column to the model

        Parameters
        ----------
        stratigraphic_column : dictionary
        cmap : matplotlib.cmap
        Returns
        -------

        Notes
        -----
        stratigraphic_column is a nested dictionary with the format
        {'group':
                {'series1':
                            {'min':0., 'max':10.,'id':0,'colour':}
                }
        }

        """
        # if the colour for a unit hasn't been specified we can just sample from 
        # a colour map e.g. tab20
        logger.info('Adding stratigraphic column to model')
        random_colour = True
        n_units=0
        for g in stratigraphic_column.keys():
            for u in stratigraphic_column[g].keys():
                if 'colour' in stratigraphic_column[g][u]:
                    random_colour = False
                    break
                n_units+=1
        if random_colour:
            import matplotlib.cm as cm
            cmap = cm.get_cmap(cmap,n_units)
            cmap_colours = cmap.colors
            ci = 0
            for g in stratigraphic_column.keys():
                for u in stratigraphic_column[g].keys():
                    stratigraphic_column[g][u]['colour'] = cmap_colours[ci,:]
        
        self.stratigraphic_column = stratigraphic_column

    def get_interpolator(self, interpolatortype='PLI', nelements=1e5,
                         buffer=0.2, element_volume = None, **kwargs):
        """
        Returns an interpolator given the arguments, also constructs a
        support for a discrete interpolator

        Parameters
        ----------
        interpolatortype : string
            define the interpolator type
        nelements : int
            number of elements in the interpolator
        buffer : double or numpy array 3x1
            value(s) between 0,1 specifying the buffer around the bounding box
        data_bb : bool
            whether to use the model boundary or the boundary around
        kwargs : no kwargs used, this just catches any additional arguments

        Returns
        -------
        """
        # get an interpolator for
        interpolator = None
        bb = np.copy(self.bounding_box)
        # add a buffer to the interpolation domain, this is necessary for
        # faults but also generally a good
        # idea to avoid boundary problems
        # buffer = bb[1, :]
        buffer = (bb[1,:]-bb[0,:])*buffer
        bb[0, :] -= buffer  # *(bb[1,:]-bb[0,:])
        bb[1, :] += buffer  # *(bb[1,:]-bb[0,:])
        box_vol = (bb[1, 0]-bb[0, 0]) * (bb[1, 1]-bb[0, 1]) * (bb[1, 2]-bb[0, 2])
        if interpolatortype == "PLI":
            if element_volume is None:
                nelements /= 5
                element_volume = box_vol / nelements
            # calculate the step vector of a regular cube
            step_vector = np.zeros(3)
            step_vector[:] = element_volume ** (1. / 3.)
            # step_vector /= np.array([1,1,2])
            # number of steps is the length of the box / step vector
            nsteps = np.ceil((bb[1, :] - bb[0, :]) / step_vector).astype(int)
            # create a structured grid using the origin and number of steps
            if self.reuse_supports:
                mesh_id = 'mesh_{}'.format(nelements)
                mesh = self.support.get(mesh_id,
                                        TetMesh(origin=bb[0, :], nsteps=nsteps,
                                                step_vector=step_vector))
                if mesh_id not in self.support:
                    self.support[mesh_id] = mesh
            else:
                mesh = TetMesh(origin=bb[0, :], nsteps=nsteps, step_vector=step_vector)
            logger.info("Creating regular tetrahedron mesh with %i elements \n"
                        "for modelling using PLI" % (mesh.ntetra))

            return PLI(mesh)

        if interpolatortype == 'FDI':

            # find the volume of one element
            if element_volume is None:
                element_volume = box_vol / nelements
            # calculate the step vector of a regular cube
            step_vector = np.zeros(3)
            step_vector[:] = element_volume ** (1. / 3.)
            # number of steps is the length of the box / step vector
            nsteps = np.ceil((bb[1, :] - bb[0, :]) / step_vector).astype(int)
            if np.any(np.less(nsteps, 3)):
                logger.error("Cannot create interpolator: number of steps is too small")
                return None
            # create a structured grid using the origin and number of steps
            if self.reuse_supports:
                grid_id = 'grid_{}'.format(nelements)
                grid = self.support.get(grid_id, StructuredGrid(origin=bb[0, :],
                                                            nsteps=nsteps,
                                                            step_vector=step_vector))
                if grid_id not in self.support:
                    self.support[grid_id] = grid
            else:
                grid = StructuredGrid(origin=bb[0, :], nsteps=nsteps,step_vector=step_vector)
            logger.info("Creating regular grid with %i elements \n"
                        "for modelling using FDI" % grid.n_elements)
            return FDI(grid)

        if interpolatortype == "DFI":  # "fold" in kwargs:
            if element_volume is None:
                nelements /= 5
                element_volume = box_vol / nelements
            # calculate the step vector of a regular cube
            step_vector = np.zeros(3)
            step_vector[:] = element_volume ** (1. / 3.)
            # number of steps is the length of the box / step vector
            nsteps = np.ceil((bb[1, :] - bb[0, :]) / step_vector).astype(int)
            # create a structured grid using the origin and number of steps
            mesh = kwargs.get('mesh', TetMesh(origin=bb[0, :], nsteps=nsteps,
                                              step_vector=step_vector))
            logger.info("Creating regular tetrahedron mesh with %i elements \n"
                        "for modelling using DFI" % mesh.ntetra)
            return DFI(mesh, kwargs['fold'])
        if interpolatortype == 'Surfe' or interpolatortype == 'surfe' and \
                surfe:
            method = kwargs.get('method', 'single_surface')
            logger.info("Using surfe interpolator")
            return Surfe(method)
        logger.warning("No interpolator")
        return interpolator

    def create_and_add_foliation(self, series_surface_data, **kwargs):
        """
        Parameters
        ----------
        series_surface_data : string
            corresponding to the feature_name in the data
        kwargs

        Returns
        -------
        feature : GeologicalFeature
            the created geological feature
        """
        if self.check_inialisation() == False:
            return False
        interpolator = self.get_interpolator(**kwargs)
        series_builder = GeologicalFeatureInterpolator(interpolator,
                                                       name=series_surface_data,
                                                       **kwargs)
        # add data
        series_data = self.data[self.data['feature_name'] == series_surface_data]
        if series_data.shape[0] == 0:
            logger.warning("No data for %s, skipping" % series_surface_data)
            return
        series_builder.add_data_from_data_frame(series_data)
        self._add_faults(series_builder)

        series_builder.build_arguments = kwargs
        self._add_feature(series_builder.feature)
        return series_builder.feature


    def _add_unconformity_above(self, feature):
        """

        Adds a region to the feature to prevent the value from being
        interpolated where the unconformities exists above e.g.
        if there is another feature above and the unconformity is at 0
        then the features added below (after) will only be visible where the
        uncomformity is <0

        Parameters
        ----------
        feature - GeologicalFeature

        Returns
        -------

        """
        for f in reversed(self.features):
            if f.type == 'unconformity':
                feature.add_region(lambda pos: f.evaluate(pos))
                break

    def _add_unconformity_below(self, feature):
        """
        Adds a region to the features that represents the
        unconformity so it is not evaluated below the unconformity

        Parameters
        ----------
        feature

        Returns
        -------

        """
        for f in self.features:
            if f.type == 'series' and feature.feature.name != f.name:
                f.add_region(lambda pos: ~feature.evaluate(pos))

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
        unconformity : GeologicalFeature
            unconformity feature  

        """
        self.parameters['features'].append({'feature_type': 'unconformity', 'feature': feature, 'value': value})
        uc_feature = UnconformityFeature(feature, value)

        # for f in self.features:
        #     f.add_region(lambda pos: uc_feature.evaluate(pos))

        # see if any unconformities are above this feature if so add region
        # self._add_unconformity_above(uc_feature)
        # self._add_unconformity_below(feature)#, uc_feature)
        self._add_feature(uc_feature)

        
        return uc_feature

    def add_onlap_unconformity(self, feature, value):
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
        unconformity_feature : GeologicalFeature
            the created unconformity

        """
        self.parameters['features'].append({'feature_type': 'onlap', 'feature': feature, 'value': value})

        uc_feature = UnconformityFeature(feature, value)

        # for f in self.features:
        #     f.add_region(lambda pos: uc_feature.evaluate(pos))

        # see if any unconformities are above this feature if so add region
        # self._add_unconformity_above(uc_feature)
        self._add_unconformity_below(uc_feature)  # , uc_feature)
        self._add_feature(uc_feature)


        return uc_feature

    def create_and_add_fault(self, fault_surface_data, displacement, renormalise=True, **kwargs):
        """
        Parameters
        ----------
        fault_surface_data : string
            name of the fault surface data in the dataframe
        displacement : displacement magnitude
        kwargs : additional kwargs for Fault and interpolators

        Returns
        -------
        fault : FaultSegment
            created fault
        """
        self.parameters['features'].append(
            {'feature_type': 'fault', 'feature_name': fault_surface_data, 'displacement': displacement, **kwargs})

        displacement_scaled = displacement / self.scale_factor
        # create fault frame
        interpolator = self.get_interpolator(**kwargs)
        fault_frame_builder = StructuralFrameBuilder(interpolator,
                                                     name=fault_surface_data,
                                                     **kwargs)
        # add data
        fault_frame_data = self.data[
            self.data['feature_name'] == fault_surface_data].copy()
        if 'coord' not in fault_frame_data:
            fault_frame_data['coord'] = 0
        vals = fault_frame_data['val']
        fault_frame_builder.add_data_from_data_frame(fault_frame_data)
        # if there is no fault slip data then we could find the strike of
        # the fault and build
        # the second coordinate
        # if we add a region to the fault then the fault operator doesn't
        # work but for visualisation
        # we want to add a region!
        # check if this fault overprint any existing faults exist in the stack

        fault_frame = fault_frame_builder.build(**kwargs)
        

        fault = FaultSegment(fault_frame, displacement=displacement_scaled,
                             **kwargs)
        # for f in reversed(self.features):
        #     if f.type == 'unconformity':
        #         fault.add_region(lambda pos: f.evaluate_value(pos) <= 0)
        #         break
        # if displacement == 0:
        #     fault.type = 'fault_inactive'
        # self._add_feature(fault)
        self._add_feature(fault)

        return fault

    def rescale(self, points, inplace=True):
        """
        Convert from model scale to real world scale - in the future this
        should also do transformations?

        Parameters
        ----------
        points : np.array((N,3),dtype=double)
        inplace : boolean
            whether to return a modified copy or modify the original array

        Returns
        -------
        points : np.array((N,3),dtype=double)

        """
        if inplace == False:
            points = points.copy()
        points *= self.scale_factor
        points += self.origin
        return points

    def scale(self, points, inplace=True):
        """ Take points in UTM coordinates and reproject
        into scaled model space

        Parameters
        ----------
        points : np.array((N,3),dtype=float)
            points to 
        inplace : bool, optional default = True
            whether to copy the points array or update the passed array
        Returns
        -------
        points : np.array((N,3),dtype=double)

        """
        if inplace==False:
            points = points.copy()

        points[:, :] -= self.origin
        points /= self.scale_factor
        return points


    def regular_grid(self, nsteps=(50, 50, 25), shuffle = True, rescale=False):
        """
        Return a regular grid within the model bounding box

        Parameters
        ----------
        nsteps : tuple
            number of cells in x,y,z

        Returns
        -------
        xyz : np.array((N,3),dtype=float)
            locations of points in regular grid
        """
        x = np.linspace(self.bounding_box[0, 0], self.bounding_box[1, 0],
                        nsteps[0])
        y = np.linspace(self.bounding_box[0, 1], self.bounding_box[1, 1],
                        nsteps[1])
        z = np.linspace(self.bounding_box[1, 2], self.bounding_box[0, 2],
                        nsteps[2])
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        locs = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        if shuffle:
            logger.info("Shuffling points")
            np.random.shuffle(locs)
        if rescale:
            locs = self.rescale(locs)
        return locs

    def evaluate_model(self, xyz, scale=True):
        """Evaluate the stratigraphic id at each location
        

        Parameters
        ----------
        xyz : np.array((N,3),dtype=float)
            locations
        scale : bool
            whether to rescale the xyz before evaluating model

        Returns
        -------
        stratigraphic_id : np.array(N,dtype=int)
            the stratigraphic index for locations
        
        Examples
        --------
        Evaluate on a voxet

        >>> x = np.linspace(model.bounding_box[0, 0], model.bounding_box[1, 0],
                        nsteps[0])
        >>> y = np.linspace(model.bounding_box[0, 1], model.bounding_box[1, 1],
                        nsteps[1])
        >>> z = np.linspace(model.bounding_box[1, 2], model.bounding_box[0, 2],
                        nsteps[2])
        >>> xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        >>> xyz = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        >>> model.evaluate_model(xyz)

        Evaluate on points defined by regular grid function
    
        >>> model.evaluate_model(model.regular_grid())


        Evaluate on a map
        
        >>> x = np.linspace(self.bounding_box[0, 0], self.bounding_box[1, 0],
                        nsteps[0])
        >>> y = np.linspace(self.bounding_box[0, 1], self.bounding_box[1, 1],
                        nsteps[1])
        >>> xx, yy = np.meshgrid(x, y, indexing='ij')
        >>> zz = np.zeros_like(yy)
        >>> xyz = np.array([xx.flatten(), yy.flatten(), zz.flatten()]).T
        >>> model.evaluate_model(xyz)
        
        """
        if scale:
            xyz = self.scale(xyz,inplace=False)
        strat_id = np.zeros(xyz.shape[0],dtype=int)
        for group in self.stratigraphic_column.keys():
            if group == 'faults':
                continue
            feature_id = self.feature_name_index.get(group, -1)
            if feature_id >= 0:
                feature = self.features[feature_id]
                vals = feature.evaluate_value(xyz)
                for series in self.stratigraphic_column[group].values():
                    strat_id[np.logical_and(vals < series.get('max',feature.max()), vals > series.get('min',feature.min()))] = series['id']
            if feature_id == -1:
                logger.error('Model does not contain {}'.format(group))
        return strat_id

    def evaluate_fault_displacements(self,points,scale=True):
        """Evaluate the fault displacement magnitude at each location
        

        Parameters
        ----------
        xyz : np.array((N,3),dtype=float)
            locations
        scale : bool
            whether to rescale the xyz before evaluating model

        Returns
        -------
        fault_displacement : np.array(N,dtype=float)
            the fault displacement magnitude
        """
        if scale:
            points = self.scale(points,inplace=False)
        vals = np.zeros(points.shape[0])
        for f in self.features:
            if f.type == 'fault':
                disp = f.displacementfeature.evaluate_value(points)
                vals[~np.isnan(disp)] += disp[~np.isnan(disp)]
        return vals*-self.scale_factor # convert from restoration magnutude to displacement


