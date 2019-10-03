import numpy as np
from LoopStructural.modelling.core.geological_points import IPoint, GPoint
from LoopStructural.modelling.features.geological_feature import GeologicalFeature
from LoopStructural.modelling.features.cross_product_geological_feature import CrossProductGeologicalFeature
from LoopStructural.supports.scalar_field import ScalarField


class StructuralFrame:
    """
    Structural frame is a curvilinear coordinate system defined by structural
    observations associated with a fault or fold.
    """
    def __init__(self, name, features):
        """

        Parameters
        ----------
        name - name of the structural frame
        features - list of features to build the frame with
        """
        self.name = name
        self.features = features
        self.data = None

    def __getitem__(self, item):
        """

        Parameters
        ----------
        item index of feature to access

        Returns
        -------
        the structural frame geological feature
        """
        return self.features[item]

    def get_feature(self, i):
        """
        Return the ith feature
        Parameters
        ----------
        i

        Returns
        -------

        """
        return self.features[i]

    def set_data(self, data):
        """
        Associate data with structural frame
        Parameters
        ----------
        data

        Returns
        -------

        """
        self.data = data

    def evaluate_value(self, evaluation_points, i=None):
        """
        Evaluate the value of the structural frame for the points.
        Can optionally only evaluate one coordinate
        Parameters
        ----------
        evaluation_points
        i

        Returns
        -------

        """
        if i is not None:
            self.features[i].support.evaluate_value(evaluation_points)
        return (self.features[0].support.evaluate_value(evaluation_points),
                self.features[1].support.evaluate_value(evaluation_points),
                self.features[2].support.evaluate_value(evaluation_points))

    def evaluate_gradient(self, evaluation_points, i=None):
        """
        Evaluate the gradient of the structural frame.
        Can optionally only evaluate the ith coordinate
        Parameters
        ----------
        evaluation_points
        i

        Returns
        -------

        """
        if i is not None:
            return self.features[i].support.evaluate_gradient(evaluation_points)
        return (self.features[0].support.evaluate_gradient(evaluation_points),
                self.features[1].support.evaluate_gradient(evaluation_points),
                self.features[2].support.evaluate_gradient(evaluation_points))



class StructuralFrameBuilder:
    """
    Class for representing a slip event of a fault
    """
    def __init__(self, interpolator, **kwargs):
        """

        Parameters
        ----------
        interpolator - a template interpolator for the frame
        kwargs
        """
        """
        mesh:  support for interpolation
        structural_feature: the geological feature that this frame describes
        data: the data that are builds this frame
        name: a name for this frame
        region: region where the interpolation occurs for this frame
        you can pass kwargs for the DSI interpolator
        """
        self.support = None
        self.fault_event = None
        self.name = 'Undefined'
        self.region = 'everywhere'
        if 'support' in kwargs:
            self.support = kwargs['support']
            del kwargs['support']
        if 'mesh' in kwargs:
            self.support = kwargs['mesh']
            del kwargs['mesh']
        if 'name' in kwargs:
            self.name = kwargs['name']

        self.data = [ [] , [] ,[]]
        # list of interpolators
        self.interpolators = []
        # Create the interpolation objects by copying the template
        self.interpolators.append(interpolator)
        self.interpolators.append(interpolator.copy())
        self.interpolators.append(interpolator.copy())
        self.interpolators[0].set_property_name(self.name+'_gx')
        self.interpolators[1].set_property_name(self.name+'_gy')
        self.interpolators[2].set_property_name(self.name+'_gz')
        for i in range(3):
            self.interpolators[i].set_region(regionname=self.region)

    def add_strike_dip_and_value(self, pos, strike, dip, val, polarity = 1, coord = None, itype= None):
        if itype == 'gx':
            coord = 0
        if itype == 'gy':
            coord = 1
        if itype == 'gz':
            coord = 2
        self.data[coord].append(GPoint.from_strike_and_dip(pos, strike, dip, polarity))
        self.data[coord].append(IPoint(pos, val))

    def add_point(self, pos, val,  coord = None, itype= None):
        if itype == 'gx':
            coord = 0
        if itype == 'gy':
            coord = 1
        if itype == 'gz':
            coord = 2
        self.data[coord].append(IPoint(pos, val))


    def add_planar_constraint(self, pos, val,  coord = None, itype= None):
        if itype == 'gx':
            coord = 0
        if itype == 'gy':
            coord = 1
        if itype == 'gz':
            coord = 2
        self.data[coord].append(GPoint(pos, val))

    def add_plunge_and_plunge_dir(self, pos, plunge, plunge_dir, polarity =1 , coord = None, itype= None):
        if itype == 'gx':
            coord = 0
        if itype == 'gy':
            coord = 1
        if itype == 'gz':
            coord = 2
        self.data[coord].append(GPoint.from_plunge_plunge_dir(pos,plunge,plunge_dir))

    def add_strike_and_dip(self, pos, strike, dip, polarity = 1, coord = None, itype= None):
        if itype == 'gx':
            coord = 0
        if itype == 'gy':
            coord = 1
        if itype == 'gz':
            coord = 2
        self.data[coord].append(GPoint.from_strike_and_dip(pos, strike, dip, polarity))
    def add_tangent_constraint(self, pos, val,coord = None, itype= None):
        # if itype == 'gx':
        #     self.data[0].append(TGPoint.from_plunge_plunge_dir(pos, plunge, plunge_dir))
        # if itype == 'gy':
        #     self.data[1].append(GPoint.from_plunge_plunge_dir(pos, plunge, plunge_dir))
        # if itype == 'gz':
        #     self.data[2].append(GPoint.from_plunge_plunge_dir(pos, plunge, plunge_dir))
        pass
    def add_tangent_constraint_angle(self,pos,s,d,itype):
        # if itype == 'gx':
        #     self.data[0].append(GPoint.from_plunge_plunge_dir(pos, plunge, plunge_dir))
        # if itype == 'gy':
        #     self.data[1].append(GPoint.from_plunge_plunge_dir(pos, plunge, plunge_dir))
        # if itype == 'gz':
        #     self.data[2].append(GPoint.from_plunge_plunge_dir(pos, plunge, plunge_dir))
        pass
    def build(self, solver='lsqr', frame=StructuralFrame, **kwargs):
        """
        Build the structural frame
        Parameters
        ----------
        solver solver to use
        frame - type of frame to build StructuralFrame or FoldFrame
        kwargs

        Returns
        -------

        """
        # reset the interpolators
        # for i in range(3):
        #     self.interpolators[i].reset()
        # all of the inteprolator weights should be set by accessing the interpolator weights
        # dictionary directly
        gxxgy = 1.
        gxxgz = 1.
        gyxgz = 1.
        if 'gxxgy' in kwargs:
            gxxgy = kwargs['gxxgy']
        if 'gxxgz' in kwargs:
            gxxgz = kwargs['gxxgz']
        if 'gyxgz' in kwargs:
            gyxgz = kwargs['gyxgz']

        for i in range(3):
            for d in self.data[i]:
                self.interpolators[i].add_data(d)
        # initialise features as none then where data exists build
        gx_feature = None
        gy_feature = None
        gz_feature = None
        if len(self.data[0]) > 0:
            print("Building structural frame coordinate 0")
            if "fold" in kwargs and "fold_weights" in kwargs:
                self.interpolators[0].update_fold(kwargs['fold'])
                self.interpolators[0].add_fold_constraints(**kwargs['fold_weights'])
                # TODO check this still works
                if 'cgw' in kwargs:
                    if kwargs['cgw'] == 0:
                        kwargs['cg'] = False
            self.interpolators[0].setup_interpolator()
            self.interpolators[0].solve_system(solver=solver,**kwargs)
            gx_feature = GeologicalFeature(self.name + '_gx',
                                      ScalarField.from_interpolator(self.interpolators[0]),
                                            data=self.data[0])
            if gx_feature is None:
                print("Not enough constraints for fold frame coordinate 0, \n"
                      "Add some more and try again.")
        if len(self.data[1]) > 0:
            print("Building structural frame coordinate 1")
            if gx_feature is not None:
                self.interpolators[1].add_gradient_orthogonal_constraint(
                    np.arange(0,self.support.n_elements),
                    gx_feature.evaluate_gradient(self.support.barycentre),
                    w=gxxgy)
            self.interpolators[1].setup_interpolator()

            self.interpolators[1].solve_system(solver=solver)
            gy_feature = GeologicalFeature(self.name + '_gy',
                                      ScalarField.from_interpolator(self.interpolators[1]),
                                            data=self.data[1])
            if gy_feature is None:
                print("Not enough constraints for fold frame coordinate 1, \n"
                      "Add some more and try again.")
        if len(self.data[2]) > 0:
            print("Building structural frame coordinate 2")
            if gy_feature is not None:
                self.interpolators[2].add_gradient_orthogonal_constraint(
                    np.arange(0, self.support.n_elements),
                    gy_feature.evaluate_gradient(self.support.barycentre),
                    w=gyxgz)
            if gx_feature is not None:
                self.interpolators[2].add_gradient_orthogonal_constraint(
                    np.arange(0, self.support.n_elements),
                    gx_feature.evaluate_gradient(self.support.barycentre),
                    w=gxxgz)

            self.interpolators[2].setup_interpolator()  # cgw=0.1)
            self.interpolators[2].solve_system(solver=solver)
            gz_feature = GeologicalFeature(self.name + '_gz',
                                      ScalarField.from_interpolator(self.interpolators[2]),
                                            data=self.data[2])
        if len(self.data[2]) == 0:
            if gy_feature is not None:
                print("Creating analytical structural frame coordinate 2")
                gz_feature = CrossProductGeologicalFeature(self.name + '_gz',  gy_feature, gx_feature)
            if gy_feature is None or gx_feature is None:
                print("Not enough constraints for fold frame coordinate 1, \n"
                      "Add some more and try again.")
        # use the frame argument to build a structural frame
        return frame(self.name, [gx_feature, gy_feature, gz_feature])
        
