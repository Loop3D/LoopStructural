import numpy as np
from FME.modelling.geological_points import IPoint,GPoint,TPoint


class StructuralFrame:
    """
    Class for representing a slip event of a fault
    """
    def __init__(self,interpolator,**kwargs):#mesh,fault_event,data,name,region):
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

    def buildFaultFrame(self, solver='lsqr', **kwargs):
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
        if overlap:
            self.mesh.indices = np.array(range(0, self.mesh.n_nodes))
            sharedindices = self.mesh.indices[overlapregion]
            node_values = self.mesh.properties[overlapsegment.interpolators['gx'].propertyname] \
                [self.mesh.indices[overlapregion]]
            As = np.zeros(node_values.shape)
            As[:] = 1
            rows = np.array(range(self.interpolators['gx'].c_, self.interpolators['gx'].c_ + len(sharedindices)))
            self.interpolators['gx'].A.extend(As.tolist())
            self.interpolators['gx'].col.extend(sharedindices.tolist())
            if shape == 'rectangular':
                self.interpolators['gx'].B.extend(node_values)
                self.interpolators['gx'].row.extend(rows.tolist())
                self.interpolators['gx'].c_ += +len(sharedindices)
            if shape == 'square':
                self.interpolators['gx'].B[overlapregion] += node_values
                self.interpolators['gx'].row.extend(sharedindices.tolist())

        if gx:
            self.interpolators['gx'].solve_system(solver=solver)
            self.mesh.update_property(self.name + '_' + 'gx', self.interpolators['gx'].c)
            self.interpolators['gy'].add_elements_gradient_orthogonal_constraint(np.arange(0, \
                                                                                           self.mesh.n_elements), \
                                                                                 self.mesh.property_gradients[
                                                                                     self.name + '_' + 'gx'], w=gxxgy)

        self.interpolators['gy'].setup_interpolator(cgw=gycg, cpw=gycp, gpw=gygcp)  # cgw=0.1)#
        if overlap:
            self.mesh.indices = np.array(range(0, self.mesh.n_nodes))
            sharedindices = self.mesh.indices[overlapregion]
            node_values = self.mesh.properties[overlapsegment.interpolators['gy'].propertyname] \
                [self.mesh.indices[overlapregion]]
            As = np.zeros(node_values.shape)
            As[:] = 1
            rows = np.array(range(self.interpolators['gy'].c_, self.interpolators['gy'].c_ + len(sharedindices)))
            self.interpolators['gy'].A.extend(As.tolist())
            self.interpolators['gy'].col.extend(sharedindices.tolist())

            if shape == 'rectangular':
                self.interpolators['gy'].B.extend(node_values)
                self.interpolators['gy'].row.extend(rows.tolist())
                self.interpolators['gy'].c_ += +len(sharedindices)
            if shape == 'square':
                self.interpolators['gy'].B[overlapregion] += node_values
                self.interpolators['gy'].row.extend(sharedindices.tolist())

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
    def get_gx(self,points,grad=False):
        if grad:
            return self.mesh.eval_gradient(points,self.name+'_gx',k=100)
        if not grad:
            return self.mesh.eval_interpolant(points,self.name+'_gx',k=100)
    def get_gy(self,points,grad=False):
        if grad:
            return self.mesh.eval_gradient(points,self.name+'_gy',k=100)
        if not grad:
            return self.mesh.eval_interpolant(points,self.name+'_gy',k=100)
    def get_gz(self,points,grad=False):
        #if self.name+'_gz' in self.mesh.properties:
        #    if grad:
        #        return self.mesh.eval_gradient(points,self.name+'_gz',k=100)
        #    if not grad:
        #        return self.mesh.eval_interpolant(points,self.name+'_gz',k=100)
        #else:
        if grad:
            #calculate dgz using cross product
            return np.cross(self.get_gy(points,grad),self.get_gx(points,grad),axisa=1,axisb=1)
        
