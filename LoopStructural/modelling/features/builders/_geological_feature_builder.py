"""
Feature builder
"""

import numpy as np
import pandas as pd

from ....utils import getLogger


from ....interpolators import GeologicalInterpolator
from ....utils.helper import (
    xyz_names,
    val_name,
    normal_vec_names,
    weight_name,
    gradient_vec_names,
    tangent_vec_names,
    interface_name,
    inequality_name,
    pairs_name,
)
from ....modelling.features import GeologicalFeature
from ....modelling.features.builders import BaseBuilder
from ....utils.helper import (
    get_data_bounding_box_map as get_data_bounding_box,
)
from ....interpolators import DiscreteInterpolator
from ....interpolators import InterpolatorFactory

logger = getLogger(__name__)


class GeologicalFeatureBuilder(BaseBuilder):
    def __init__(
        self,
        interpolatortype: str,
        bounding_box,
        nelements: int = 1000,
        name="Feature",
        model=None,
        **kwargs,
    ):
        """
        Constructor for a GeologicalFeatureBuilder

        Parameters
        ----------
        interpolator : GeologicalInterpolator
            An empty GeologicalInterpolator
        region : lambda function
            defining whether the location (xyz) should be included in the
        kwargs - name of the feature, region to interpolate the feature
        """
        BaseBuilder.__init__(self, model, name)
        interpolator = InterpolatorFactory.create_interpolator(
            interpolatortype=interpolatortype,
            boundingbox=bounding_box,
            nelements=nelements,
            buffer=kwargs.get("buffer", 0.2),
        )
        
        if not issubclass(type(interpolator), GeologicalInterpolator):
            raise TypeError(
                "interpolator is {} and must be a GeologicalInterpolator".format(type(interpolator))
            )
        self._interpolator = interpolator
        self._up_to_date = self._interpolator.up_to_date
        header = (
            xyz_names()
            + val_name()
            + gradient_vec_names()
            + normal_vec_names()
            + tangent_vec_names()
            + weight_name()
        )
        self.data = pd.DataFrame(columns=header)
        self.data_added = False

        self._feature = GeologicalFeature(
            self._name,
            builder=self,
            regions=[],
            faults=self.faults,
        )
        self._orthogonal_features = {}
        self._equality_constraints = {}
        # add default parameters
        self.update_build_arguments({'cpw':1.0,'npw':1.0,'regularisation':1.,'nelements':self.interpolator.n_elements})
    def set_not_up_to_date(self, caller):
        logger.info(
            f"Setting {self.name} to not up to date from an instance of {caller.__class__.__name__}"
        )
        self._up_to_date = False
        self._interpolator.up_to_date = False

    @property
    def interpolator(self):
        return self._interpolator

    @interpolator.setter
    def interpolator(self, interpolator):
        if not issubclass(type(interpolator), GeologicalInterpolator):
            raise TypeError(
                "interpolator is {} and must be a GeologicalInterpolator".format(type(interpolator))
            )

    def add_data_from_data_frame(self, data_frame, overwrite=False):
        """
        Extract data from a pandas dataframe with columns for

        Parameters
        ----------
        data_frame : pd.DataFrame
            a dataframe containing the data to be added

        Returns
        -------

        """
        self.data = data_frame.copy()

    def add_orthogonal_feature(self, feature, w=1.0, region=None, step=1, B=0):
        """
        Add a constraint to the interpolator so that the gradient of an exisitng
        feature is orthogonal to the feature being built. E.g. dot product
        between gradients should be = 0

        Parameters
        ----------
        feature : GeologicalFeature
            feature which we want to be orthogonal to
        w :  double
            how much to weight in least squares sense
        region : unused
        step : int
            numpy slicing step size to see how many tetras to add

        Returns
        -------

        Notes
        -----
        The constraint can be applied to a random subset of the tetrahedral
        elements in the mesh
        """
        try:
            step = int(step)  # cast as int in case it was a float
        except ValueError:
            logger.error("Cannot cast {} as integer, setting step to 1".format(step))
            step = 1
        self._orthogonal_features[feature.name] = [feature, w, region, step, B]

        logger.info(f"Adding orthogonal constraint {feature.name} to {self.name}")
        self._up_to_date = False

    def add_data_to_interpolator(self, constrained=False, force_constrained=False, **kwargs):
        """
        Iterates through the list of data and applies any faults active on the
        data in the order they are added

        Parameters
        -----------
        constrained : boolean
        force_constrained : boolean

        Returns
        -------

        """
        logger.info('Adding data to interpolator for {}'.format(self.name))
        logger.info(f"Data shape: {self.data.shape}")
        logger.info(f'Constrained: {constrained}, force_constrained: {force_constrained}')
        if self.data_added:
            logger.info("Data already added to interpolator")
            return
        # first move the data for the fault
        logger.info(f"Adding {len(self.faults)} faults to {self.name}")
        data = self.data.copy()

        # convert data locations to numpy array and then update
        # for gradient data we need to calculate the rotation of the vector field
        # create a tetrahedron at each point, and assign the corner values as the values that would
        # fit the vector for the gradient. Then apply the fault to these points and then recalculate
        # the gradient using the tetrahedron.

        mask_gradient = np.all(~np.isnan(data.loc[:, gradient_vec_names()].to_numpy(float)), axis=1)
        mask_normal = np.all(~np.isnan(data.loc[:, normal_vec_names()].to_numpy(float)), axis=1)
        mask_tangent = np.all(~np.isnan(data.loc[:, tangent_vec_names()].to_numpy(float)), axis=1)
        mask_vector = mask_gradient | mask_normal | mask_tangent
        if np.sum(mask_vector) > 0:
            for f in self.faults:
                if np.sum(mask_normal) > 0:
                    data.loc[mask_normal, normal_vec_names()] = f.apply_to_vectors(
                        data.loc[mask_normal, xyz_names() + normal_vec_names()].to_numpy(float)
                    )
                if np.sum(mask_gradient) > 0:
                    data.loc[mask_gradient, gradient_vec_names()] = f.apply_to_vectors(
                        data.loc[mask_gradient, xyz_names() + gradient_vec_names()].to_numpy(float)
                    )
                if np.sum(mask_tangent) > 0:
                    data.loc[mask_tangent, tangent_vec_names()] = f.apply_to_vectors(
                        data.loc[mask_tangent, xyz_names() + tangent_vec_names()].to_numpy(float)
                    )

        for f in self.faults:
            data.loc[:, xyz_names()] = f.apply_to_points(data.loc[:, xyz_names()])

        # self.check_interpolation_geometry(data.loc[:,xyz_names()].to_numpy())
        # Now check whether there are enough constraints for the
        # interpolator to be able to solve
        # we need at least 2 different value points or a single norm
        # constraint. If there are not enough
        # try converting grad to norms, if still not enough send user an error
        if constrained:
            # Change normals to gradients
            mask = np.all(~np.isnan(data.loc[:, normal_vec_names()]), axis=1)
            if mask.sum() > 0:
                data.loc[mask, gradient_vec_names()] = data.loc[mask, normal_vec_names()].to_numpy(
                    float
                )
                data.loc[mask, normal_vec_names()] = np.nan
        if self.get_norm_constraints().shape[0] > 0:
            constrained = True

        if np.unique(self.get_value_constraints()[:, 3]).shape[0] > 1:
            constrained = True

        if not constrained or force_constrained:
            # change gradient constraints to normal vector constraints
            mask = np.all(~np.isnan(data.loc[:, gradient_vec_names()]), axis=1)
            if mask.sum() > 0:
                data.loc[mask, normal_vec_names()] = data.loc[mask, gradient_vec_names()].to_numpy(
                    float
                )
                data.loc[mask, gradient_vec_names()] = np.nan
                logger.info("Setting gradient points to norm constraints")
                constrained = True
                mask = np.all(~np.isnan(data.loc[:, normal_vec_names()].to_numpy(float)), axis=1)

        if not constrained:
            logger.error("Not enough constraints for scalar field add more")
        # self.interpolator.reset()
        mask = ~np.isnan(data.loc[:, val_name()].to_numpy(float))
        # add value constraints
        if mask.sum() > 0:
            value_data = data.loc[mask[:, 0], xyz_names() + val_name() + weight_name()].to_numpy(
                float
            )
            self.interpolator.set_value_constraints(value_data)

        # add gradient constraints
        mask = np.all(~np.isnan(data.loc[:, gradient_vec_names()].to_numpy(float)), axis=1)
        if mask.sum() > 0:
            gradient_data = data.loc[
                mask, xyz_names() + gradient_vec_names() + weight_name()
            ].to_numpy(float)
            self.interpolator.set_gradient_constraints(gradient_data)

        # add normal vector data
        mask = np.all(~np.isnan(data.loc[:, normal_vec_names()].to_numpy(float)), axis=1)
        if mask.sum() > 0:
            normal_data = data.loc[mask, xyz_names() + normal_vec_names() + weight_name()].to_numpy(
                float
            )
            self.interpolator.set_normal_constraints(normal_data)

        # add tangent data
        mask = np.all(~np.isnan(data.loc[:, tangent_vec_names()].to_numpy(float)), axis=1)
        if mask.sum() > 0:
            tangent_data = data.loc[
                mask, xyz_names() + tangent_vec_names() + weight_name()
            ].to_numpy(float)
            self.interpolator.set_tangent_constraints(tangent_data)

        # add interface constraints
        mask = np.all(~np.isnan(data.loc[:, interface_name()].to_numpy(float)), axis=1)
        if mask.sum() > 0:
            interface_data = data.loc[
                mask, xyz_names() + interface_name() + weight_name()
            ].to_numpy(float)
            self.interpolator.set_interface_constraints(interface_data)
        # add inequality constraints
        mask = np.all(~np.isnan(data.loc[:, inequality_name()].to_numpy(float)), axis=1)
        if mask.sum() > 0:
            inequality_data = data.loc[mask, xyz_names() + inequality_name()].to_numpy(float)
            self.interpolator.set_value_inequality_constraints(inequality_data)
        mask = np.all(~np.isnan(data.loc[:, pairs_name()].to_numpy(float)), axis=1)
        if mask.sum() > 0:
            pairs_data = data.loc[mask, xyz_names() + pairs_name()].to_numpy(float)
            self.interpolator.set_inequality_pairs_constraints(pairs_data)
        self.data_added = True
        self._up_to_date = False

    def install_gradient_constraint(self):
        if issubclass(type(self.interpolator), DiscreteInterpolator):
            for g in self._orthogonal_features.values():
                feature, w, region, step, B = g
                if w == 0:
                    continue
                logger.info(f"Adding gradient orthogonal constraint {feature.name} to {self.name}")
                logger.info(f'Evaluating gradient {feature.name} for {self.name}')
                # don't worry about masking for unconformities here
                vector = feature.evaluate_gradient(
                    self.interpolator.support.barycentre, ignore_regions=True
                )

                norm = np.linalg.norm(vector, axis=1)

                vector[norm > 0] /= norm[norm > 0, None]
                element_idx = np.arange(self.interpolator.support.n_elements)
                logger.info(f"Adding to least squares matrix: {self.name}")

                self.interpolator.add_gradient_orthogonal_constraints(
                    self.interpolator.support.barycentre[element_idx[::step], :],
                    vector[element_idx[::step], :],
                    w=w,
                    b=B,
                    name=f'{feature.name}_orthogonal',
                )

    def add_equality_constraints(self, feature, region, scalefactor=1.0):
        self._equality_constraints[feature.name] = [feature, region, scalefactor]
        logger.info(f'Adding equality constraints to {self.name}')
        self._up_to_date = False

    def install_equality_constraints(self):
        for e in self._equality_constraints.values():
            try:
                # assume all parts of structural frame have the same support
                support = self.interpolator.support

                # work out the values of the nodes where we want hard
                # constraints
                idc = np.arange(0, support.n_nodes)[e[1](support.nodes)]
                val = e[0].evaluate_value(support.nodes[e[1](support.nodes), :])
                mask = ~np.isnan(val)
                self.interpolator.add_equality_constraints(idc[mask], val[mask] * e[2])
            except BaseException as e:
                logger.error(f"Could not add equality for {self.name}")
                logger.error(f"Exception: {e}")

    def get_value_constraints(self):
        """
        Get the value constraints for this geological feature

        Returns
        -------
        np.array((N,4),dtype=double)
        """
        header = xyz_names() + val_name() + weight_name()
        mask = ~np.isnan(self.data.loc[:, val_name()].to_numpy(float))
        return self.data.loc[mask[:, 0], header].to_numpy(float)

    def get_gradient_constraints(self):
        """
        Get the gradient direction constraints

        Returns
        -------
        numpy array
        """
        mask = np.all(~np.isnan(self.data.loc[:, gradient_vec_names()].to_numpy(float)), axis=1)
        if mask.shape[0] > 0:
            return self.data.loc[mask, xyz_names() + gradient_vec_names() + weight_name()].to_numpy(
                float
            )
        else:
            return np.zeros((0, 7))

    def get_tangent_constraints(self):
        """

        Returns
        -------
        numpy array
        """
        mask = np.all(~np.isnan(self.data.loc[:, tangent_vec_names()].to_numpy(float)), axis=1)
        if mask.shape[0] > 0:
            return self.data.loc[mask, xyz_names() + tangent_vec_names() + weight_name()].to_numpy(
                float
            )
        else:
            return np.zeros((0, 7))

    def get_norm_constraints(self):
        """
        Get the gradient norm constraints

        Returns
        -------
        numpy array
        """
        mask = np.all(~np.isnan(self.data.loc[:, normal_vec_names()].to_numpy(float)), axis=1)
        if mask.shape[0] > 0:
            return self.data.loc[mask, xyz_names() + normal_vec_names() + weight_name()].to_numpy(
                float
            )
        else:
            return np.zeros((0, 7))

    def get_orientation_constraints(self):
        """
        Get the orientation constraints

        Returns
        -------
        numpy array
        """
        gradient_constraints = self.get_gradient_constraints()
        normal_constraints = self.get_norm_constraints()
        return np.vstack([gradient_constraints, normal_constraints])

    def get_interface_constraints(self):
        mask = np.all(~np.isnan(self.data.loc[:, interface_name()].to_numpy(float)), axis=1)
        if mask.shape[0] > 0:
            return self.data.loc[mask, xyz_names() + interface_name() + weight_name()].to_numpy(
                float
            )
        else:
            return np.zeros((0, 5))
    # def get_inequality_pair_constraints(self):
    #     mask = np.all(~np.isnan(self.data.loc[:,pairs_name()].to_numpy(float())),axis=1)
    #     if mask.shape[0] > 0:
    #         return self.data.loc[mask, xyz_names() + pairs_name()+weight_name()].to_numpy(float)
    # def get_inequality_constraints(self):
    #     mask = np.all(~np.isnan(self.data.loc[:,inequality_name()].to_numpy(float())),axis=1)
    #     if mask.shape[0] > 0:
    #         return self.data.loc[mask, xyz_names() + inequality_name()+weight_name()].to_numpy(float)
    def get_data_locations(self):
        """
        Get only the location for all data points

        Returns
        -------

        """
        return self.data.loc[:, xyz_names()].to_numpy(float)

    def set_interpolation_geometry(self, origin, maximum, rotation=None):
        """Update the interpolation support geometry to new bounding box

        Parameters
        ----------
        origin : np.array(3)
            origin vector
        maximum : np.array(3)
            maximum vector
        """
        logger.info(f"Setting mesh origin: {origin[0]} {origin[1]} {origin[2]} ")
        logger.info(f"Setting mesh maximum: {maximum[0]} {maximum[1]} {maximum[2]}")
        if np.any(np.isnan(origin)):
            logger.warning("Origin is NaN, not updating")
            return

        if np.any(np.isnan(maximum)):
            logger.warning("Maximum is NaN, not updating")
            return

        self.interpolator.support.origin = origin
        self.interpolator.support.maximum = maximum
        self.interpolator.support.rotation_xy = rotation
        self._up_to_date = False

        while self.interpolator.dof < 100:
            self.interpolator.support.step_vector = self.interpolator.support.step_vector * 0.9

    def check_interpolation_geometry(self, data, buffer=0.3):
        """Check the interpolation support geometry
        to data to make sure everything fits
        Apply the fault to the model grid to ensure that the support
        is big enough to capture the faulted feature.
        """
        if self.interpolator.support is not None:
            logger.info(f"Checking interpolation geometry for {self.name}")
            origin = self.interpolator.support.origin
            maximum = self.interpolator.support.maximum
            pts = self.model.bounding_box.with_buffer(buffer).regular_grid(local=True)
            for f in self.faults:
                pts = f.apply_to_points(pts)

            origin[origin > np.min(pts, axis=0)] = np.min(pts, axis=0)[origin > np.min(pts, axis=0)]
            maximum[maximum < np.max(pts, axis=0)] = np.max(pts, axis=0)[
                maximum < np.max(pts, axis=0)
            ]
            self.interpolator.support.origin = origin
            self.interpolator.support.maximum = maximum
            self.update_build_arguments({'nelements':self.interpolator.n_elements})
    def build(self, data_region=None, **kwargs):
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
        logger.info(f'Building {self.name}')
        # self.get_interpolator(**kwargs)
        for f in self.faults:
            f.builder.update()
        domain = kwargs.get("domain", None)
        if 'nelements' in kwargs:
            # if the number of elements has changed then update the interpolator
            logger.info(f'Interpolator has {self.interpolator.n_elements} elements')
            logger.info(f'Kwargs has {kwargs["nelements"]} elements')
            if self.interpolator.n_elements != kwargs['nelements']:
                logger.info('Setting nelements to {} for {}'.format(kwargs['nelements'], self.name))
                self.build_arguments['nelements'] = self.interpolator.set_nelements(kwargs['nelements'])
                logger.info(f'Interpolator nelements {self.interpolator.n_elements}')
                self._up_to_date = False
        if domain:
            self.check_interpolation_geometry(None)
        self.add_data_to_interpolator(**kwargs)
        if data_region is not None:
            xyz = self.interpolator.get_data_locations()
            bb, region = get_data_bounding_box(xyz, data_region)
            # if self.model.reuse_supports == False:
            if np.any(np.min(bb[0, :]) < self.interpolator.support.origin):
                neworigin = np.min([self.interpolator.support.origin, bb[0, :]], axis=0)
                logger.info(
                    f"Changing origin of support for {self.name} from \
                        {self.interpolator.support.origin[0]} \
                        {self.interpolator.support.origin[1]} \
                        {self.interpolator.support.origin[2]} \
                        to {neworigin[0]} {neworigin[1]} {neworigin[2]}"
                )

                self.interpolator.support.origin = neworigin
            if np.any(np.max(bb[0, :]) < self.interpolator.support.maximum):
                newmax = np.max([self.interpolator.support.maximum, bb[0, :]], axis=0)
                logger.info(
                    f"Changing origin of support for {self.name} from \
                        from {self.interpolator.support.maximum[0]} \
                        {self.interpolator.support.maximum[1]} \
                            {self.interpolator.support.maximum[2]} \
                                to {newmax[0]} {newmax[1]} {newmax[2]}"
                )

                self.interpolator.support.maximum = newmax
            self.interpolator.set_region(region=region)
        logger.info(f'Setting up interpolator for {self.name}')
        self.interpolator.setup_interpolator(**kwargs)
        logger.info(f'installing gradient constraints {self.name}')
        self.install_gradient_constraint()
        logger.info(f'installing equality constraints {self.name}')
        self.install_equality_constraints()
        logger.info(f'running interpolation for {self.name}')

        self.interpolator.solve_system(
            solver=kwargs.get('solver', 'cg'),
            tol=kwargs.get('tol', None),
            solver_kwargs=kwargs.get('solver_kwargs', {}),
        )
        logger.info(f'Finished building  {self.name}')
        self._up_to_date = True
        return self._feature
