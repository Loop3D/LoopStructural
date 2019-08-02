import numpy as np
from FME.modelling.geological_points import IPoint, GPoint, TPoint
from FME.modelling.features.geological_feature import GeologicalFeature
from FME.modelling.features.cross_product_geological_feature import CrossProductGeologicalFeature
from FME.modelling.scalar_field import TetrahedralMeshScalarField

class StructuralFrame:
    def __init__(self, name, features):
        self.name = name
        self.features = features
        self.data = None
    def get_feature(self,i):
        return self.features[i]
    def set_data(self, data):
        self.data = data
    def evaluate_value(self, evaluation_points, i=None):
        if i is not None:
            self.features[i].support.evaluate_value(evaluation_points)
        return (self.features[0].support.evaluate_value(evaluation_points),
                self.features[1].support.evaluate_value(evaluation_points),
                self.features[2].support.evaluate_value(evaluation_points))

    def evaluate_gradient(self, evaluation_points, i=None):
        if i is not None:
            return self.features[i].support.evaluate_gradient(evaluation_points)
        return (self.features[0].support.evaluate_gradient(evaluation_points),
                self.features[1].support.evaluate_gradient(evaluation_points),
                self.features[2].support.evaluate_gradient(evaluation_points))

    def get_data(self,itype=None,ptype=None):
        if itype == None and ptype== None:
            return self.data
        data = []
        for d in self.data:
            if d['type'] == itype:
                if ptype == None:
                    data.append(d)
                if type(d['data']) == ptype:
                    data.append(d)
        return data

    def get_values(self, i):
        return self.features[i].support.get_node_values()


class StructuralFrameBuilder:
    """
    Class for representing a slip event of a fault
    """
    def __init__(self, interpolator, **kwargs):
        """
        mesh:  support for interpolation
        structural_feature: the geological feature that this frame describes
        data: the data that are builds this frame
        name: a name for this frame
        region: region where the interpolation occurs for this frame
        you can pass kwargs for the DSI interpolator
        """
        self.mesh = None
        self.fault_event = None
        self.name = 'Undefined'
        self.region = 'everywhere'
        if 'mesh' in kwargs:
            self.mesh = kwargs['mesh']
            del kwargs['mesh']
        if 'fault_event' in kwargs:
            self.fault_event = kwargs['fault_event']
        if 'name' in kwargs:
            self.name = kwargs['name']
        if 'region' in kwargs:
            self.region = kwargs['region']
        self.data = []
        #dictionary to contain all of the interpolation objects for this frame
        self.interpolators = []
        #Create the interpolation objects you
        self.interpolators.append(interpolator)
        self.interpolators.append(interpolator.copy())
        self.interpolators.append(interpolator.copy())
        self.interpolators[0].set_property_name(self.name+'_gx')
        self.interpolators[1].set_property_name(self.name+'_gy')
        self.interpolators[2].set_property_name(self.name+'_gz')
        for i in range(3):
            self.interpolators[i].set_region(regionname=self.region)

    def add_strike_dip_and_value(self, pos, strike, dip, val, itype):
        self.data.append({'type':itype,'data':GPoint(pos,strike,dip)})
        self.data.append({'type':itype,'data':IPoint(pos,val)})

    def add_point(self, pos, val, itype):
        self.data.append({'type':itype,'data':IPoint(pos,val)})

    def add_planar_constraint(self, pos, val, itype):
        self.data.append({'type':itype,'data':GPoint(pos,val)})

    def add_strike_and_dip(self, pos, s, d, itype):
        self.data.append({'type':itype,'data':GPoint(pos,s,d)})

    def add_tangent_constraint(self, pos, val,itype):
        self.data.append({'type':itype,'data':TPoint(pos,val)})

    def add_tangent_constraint_angle(self,pos,s,d,itype):
        self.data.append({'type': itype,'data':TPoint(pos,s,d)})

    def build(self, solver='lsqr', frame = StructuralFrame, **kwargs):
        """
        Build the fault frame for this segment using the solver specified, default is scipy lsqr

        """
        gxxgy = 1.
        gxxgz = 1.
        gyxgz = 1.
        gxcg = 0.1
        gycg = 0.1
        gzcg = 0.1
        gxcp = 1
        gycp = 1
        gzcp = 1
        gxgcp = 1
        gygcp = 1
        gzgcp = 1
        if 'gxxgy' in kwargs:
            gxxgy = kwargs['gxxgy']
        if 'gxxgz' in kwargs:
            gxxgz = kwargs['gxxgz']
        if 'gyxgz' in kwargs:
            gyxgz = kwargs['gyxgz']
        if 'gxcg' in kwargs:
            gxcg = kwargs['gxcg']
        if 'gycg' in kwargs:
            gycg = kwargs['gycg']
        if 'gzcg' in kwargs:
            gzcg = kwargs['gzcg']
        if 'gxcp' in kwargs:
            gxcp = kwargs['gxcp']
        if 'gycp' in kwargs:
            gycp = kwargs['gycp']
        if 'gzcp' in kwargs:
            gzcp = kwargs['gzcp']
        if 'gxgcp' in kwargs:
            gxgcp = kwargs['gxgcp']
        if 'gygcp' in kwargs:
            gygcp = kwargs['gygcp']
        if 'gzgcp' in kwargs:
            gzgcp = kwargs['gzgcp']
        shape = 'rectangular'
        if 'shape' in kwargs:
            shape = kwargs['shape']
        gx = False
        gy = False
        gz = False
        for d in self.data:
            if d['type'] == 'gx':
                gx = True
                self.interpolators[0].add_data(d['data'])
            if d['type'] == 'gy':
                gy = True
                self.interpolators[1].add_data(d['data'])
            if d['type'] == 'gz':
                gz = True
                self.interpolators[2].add_data(d['data'])
        gx_feature = None
        gy_feature = None
        gz_feature = None
        if gx:
            print("Building gx")
            self.interpolators[0].setup_interpolator(cgw=gxcg, cpw=gxcp, gpw=gxgcp)
            self.interpolators[0].solve_system(solver=solver)
            gx_feature =  GeologicalFeature(self.name + '_gx',
                                      TetrahedralMeshScalarField.from_interpolator(self.interpolators[0]))
            # self.mesh.update_property(self.name + '_' + 'gx', self.interpolators[0].c)
        if gy:
            print("Building gy")
            self.interpolators[1].add_elements_gradient_orthogonal_constraint(
                np.arange(0,self.mesh.n_elements),
                gx_feature.evaluate_gradient(self.mesh.barycentre),
                w=gxxgy)
            self.interpolators[1].setup_interpolator(cgw=gycg, cpw=gycp, gpw=gygcp)

            self.interpolators[1].solve_system(solver=solver)
            gy_feature = GeologicalFeature(self.name + '_gy',
                                      TetrahedralMeshScalarField.from_interpolator(self.interpolators[1]))
        if gz:
            print("Building gz")
            self.interpolators[2].add_elements_gradient_orthogonal_constraint(
                np.arange(0, self.mesh.n_elements),
                gx_feature.evaluate_gradient(self.mesh.barycentre),
                w=gyxgz)
            self.interpolators[2].add_elements_gradient_orthogonal_constraint(
                np.arange(0, self.mesh.n_elements),
                gy_feature.evaluate_gradient(self.mesh.barycentre),
                w=gyxgz)
            self.interpolators[2].setup_interpolator(cgw=gzcg, cpw=gzcp, gpw=gzgcp)  # cgw=0.1)
            self.interpolators[2].solve_system(solver=solver)
            gz_feature = GeologicalFeature(self.name + '_gz',
                                      TetrahedralMeshScalarField.from_interpolator(self.interpolators[2]))
        if gz is False:
            print("Creating analytical gz")
            gz_feature = CrossProductGeologicalFeature(self.name + '_gz',gx_feature, gy_feature)

        return frame(self.name, [gx_feature, gy_feature, gz_feature])
        
