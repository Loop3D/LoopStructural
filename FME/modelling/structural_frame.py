import numpy as np
from FME.modelling.geological_points import IPoint, GPoint, TPoint


class StructuralFrame(GeologicalFeature):
    def __init__(self, age, name, supports):
        self.age = age
        self.name = name
        self.supports = supports
        self.ndim = 3

    def evaluate_value(self, evaluation_points):
        return (self.support[0].evaluate_value(evaluation_points),
                self.support[1].evaluate_value(evaluation_points),
                self.support[2].evaluate_value(evaluation_points))

    def evaluate_gradient(self, evaluation_points):
        return (self.support[0].evaluate_gradient(evaluation_points),
                self.support[1].evaluate_gradient(evaluation_points),
                self.support[2].evaluate_gradient(evaluation_points))

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
        self.interpolators = {}
        #If this frame overlaps with another frame, store that frame here and the condition
        #of overlap
        self.overlap = {}
        #Create the interpolation objects you
        self.interpolators['gx'] = interpolator(itype='gx',**kwargs)
        self.interpolators['gy'] = interpolator(itype='gy',**kwargs)
        self.interpolators['gz'] = interpolator(itype='gz',**kwargs)
    def add_strike_dip_and_value(self,pos,strike,dip,val,itype):
        self.data.append({'type':itype,'data':GPoint(pos,strike,dip)})
        self.data.append({'type':itype,'data':IPoint(pos,val)})
    def add_point(self,pos,val,itype):
        self.data.append({'type':itype,'data':IPoint(pos,val)})
    def add_planar_constraint(self,pos,val,itype):
        self.data.append({'type':itype,'data':GPoint(pos,val)})
    def add_strike_and_dip(self,pos,s,d,itype):
        self.data.append({'type':itype,'data':GPoint(pos,s,d)})
    def add_tangent_constraint(self,pos,val,itype):
        self.data.append({'type':itype,'data':TPoint(pos,val)})
    def add_tangent_constraint_angle(self,pos,s,d,itype):
        self.data.append({'type':itype,'data':TPoint(pos,s,d)})

    def build(self, solver='lsqr', **kwargs):
        """
        Build the fault frame for this segment using the solver specified, default is scipy lsqr

        """
        # determine region
        overlap = False
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
        gx = True
        gy = True
        gz = True
        if 'gx' in kwargs:
            gx = kwargs['gx']
        if 'gy' in kwargs:
            gy = kwargs['gy']
        if 'gz' in kwargs:
            gz = kwargs['gz']
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
        if 'segment' in self.overlap:
            overlapsegment = self.overlap['segment']
            overlapregion = self.overlap['region']
            overlap = True
        shape = 'rectangular'
        if 'shape' in kwargs:
            shape = kwargs['shape']

        for d in self.data:
            if d['type'] == 'gx':
                self.interpolators['gx'].add_data(d['data'])
            if d['type'] == 'gy':
                self.interpolators['gy'].add_data(d['data'])
            if d['type'] == 'gz':
                self.interpolators['gz'].add_data(d['data'])
        self.interpolators['gx'].setup_interpolator(cgw=gxcg, cpw=gxcp, gpw=gxgcp)
        self.interpolators['gx'].solve_system(solver=solver)
        self.interpolators['gy'].add_elements_gradient_orthogonal_constraint(np.arange(0, \
                                                                                       self.mesh.n_elements), \
                                                                             self.mesh.property_gradients[
                                                                                 self.name + '_' + 'gx'], w=gxxgy)
        self.interpolators['gy'].setup_interpolator(cgw=gycg, cpw=gycp, gpw=gygcp)  # cgw=0.1)#

        self.mesh.update_property(self.name + '_' + 'gx', self.interpolators['gx'].c)
        self.interpolators['gy'].solve_system(solver=solver)
        self.mesh.update_property(self.name + '_' + 'gy', self.interpolators['gy'].c)
        self.interpolators['gz'].add_elements_gradient_orthogonal_constraint(np.arange(0, self.mesh.n_elements) \
                                                                             , \
                                                                             self.mesh.property_gradients[
                                                                                 self.name + '_' + 'gx']
                                                                             , w=gxxgz)
        self.interpolators['gz'].add_elements_gradient_orthogonal_constraint(np.arange(0, self.mesh.n_elements) \
                                                                             , self.mesh.property_gradients[
                                                                                 self.name + '_' + 'gy'], w=gyxgz)

        self.interpolators['gz'].setup_interpolator(cgw=gzcg, cpw=gzcp, gpw=gzgcp)  # cgw=0.1)
        self.interpolators['gz'].solve_system(solver=solver)
        self.mesh.update_property(self.name + '_' + 'gz', self.interpolators['gz'].c)

        return StructuralFrame(0, 'name',
                               [self.interpolators['gx'].get_support(),
                                self.interpolators['gy'].get_support(),
                                self.interpolators['gz'].get_support()])
        
