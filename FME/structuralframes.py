import numpy as np

from .dsi_interpolator import DSI
from .geological_points import IPoint,GPoint,TPoint,IePoint
class StructuralFrame:
    """
    Class for representing a slip event of a fault
    """
    def __init__(self,**kwargs):#mesh,fault_event,data,name,region):
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
        self.interpolators['gx'] = DSI(self.mesh,itype='gx',propertyname=self.name+'_'+'gx',**kwargs)
        self.interpolators['gy'] = DSI(self.mesh,itype='gy',propertyname=self.name+'_'+'gy',**kwargs)
        self.interpolators['gz'] = DSI(self.mesh,itype='gz',propertyname=self.name+'_'+'gz',**kwargs)
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
        
