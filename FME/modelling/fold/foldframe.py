from FME.modelling.structuralframes import StructuralFrame
from FME.interpolators.dsi_helper import fold_cg, cg_cstr_to_coo_sq
import numpy as np
class FoldFrame(StructuralFrame):
    """
    Class for representing a slip event of a fault
    """
    def __init__(self,**kwargs):#mesh,fault_event,data,name,region):
        """
        mesh:  support for interpolation
        fault_event: the parent fault that this segment belongs to
        data: the data that are associated with this segment
        name: a name for this segment
        region: region where the interpolation occurs for this segment
        """
        super().__init__(**kwargs)
        self.fold_event = None

        if 'fold_event' in kwargs:
            self.fold_event = kwargs['fault_event']
        if 'name' in kwargs:
            self.name = kwargs['name']
    def buildFrame(self,solver='lsqr',**kwargs):
        """
        Build the fault frame for this segment using the solver specified, default is scipy lsqr
        
        """
        #determine region
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
        gz = False
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
        if 'gz' in kwargs:
            gx = True
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
        self.interpolators['gx'].setup_interpolator(cgw=gxcg,cpw=gxcp,gpw=gxgcp)           
        self.interpolators['gx'].solve_system(solver=solver)
        self.mesh.update_property(self.name+'_'+'gx', self.interpolators['gx'].c)
        self.interpolators['gy'].add_elements_gradient_orthogonal_constraint(np.arange(0,self.mesh.n_elements),self.mesh.property_gradients[self.name+'_'+'gx'],w=gxxgy)
        ##project constant gradient constraint of gy onto gx ----------------------------------------

        self.interpolators['gy'].setup_interpolator(cgw=0.0,cpw=gycp,gpw=gygcp)#cgw=0.1)#
        eg = self.mesh.get_elements_gradients(np.arange(self.mesh.n_elements))
        dgx = self.mesh.property_gradients[self.interpolators['gx'].propertyname]
        idc,c,ncons = fold_cg(eg,dgx,self.mesh.neighbours,self.mesh.elements,self.mesh.nodes)
        a,r,c = cg_cstr_to_coo_sq(c,idc,ncons)
        a = np.array(a)
        a*=gycg
        A = []
        row = []
        col = []
        A.extend(a.tolist())
        row.extend(np.array(r).tolist())
        col.extend(np.array(c).tolist())
        self.interpolators['gy'].add_fold_constraints(A,np.zeros(self.mesh.n_nodes),col,row)
        ####-----------------------------------------------------------------------------------------
        self.interpolators['gy'].solve_system(solver=solver)
        self.mesh.update_property(self.name+'_'+'gy', self.interpolators['gy'].c)

        if gz:
            self.interpolators['gz'].add_elements_gradient_orthogonal_constraint(np.arange(0,self.mesh.n_elements),self.mesh.property_gradients[self.name+'_'+'gx'],w=gxxgz)
            self.interpolators['gz'].add_elements_gradient_orthogonal_constraint(np.arange(0,self.mesh.n_elements),self.mesh.property_gradients[self.name+'_'+'gy'],w=gyxgz)

            self.interpolators['gz'].setup_interpolator(cgw=gzcg,cpw=gzcp,gpw=gzgcp)#cgw=0.1)
            self.interpolators['gz'].solve_system(solver=solver)
            self.mesh.update_property(self.name+'_'+'gz', self.interpolators['gz'].c)  
    def calculate_fold_axis_rotation(self,points):
        s1g = self.get_gx(points[:,:3],grad=True)
        s1g /= np.linalg.norm(s1g,axis=1)[:,None]
        s1 = self.get_gx(points[:,:3],grad=False)
        s1gyg = self.get_gy(points[:,:3],grad=True)
        s1gyg /= np.linalg.norm(s1gyg,axis=1)[:,None]
        s1gy = self.get_gy(points[:,:3],grad=False)
        l1 = points[:,3:]#np.cross(s1g,s0g,axisa=1,axisb=1)
        l1 /= np.linalg.norm(l1,axis=1)[:,None]
        #einsum dot product
        far = np.einsum('ij,ij->i',s1gyg,l1)
        far = np.rad2deg(np.arccos(far))
        #scalar triple product
        stp = np.einsum('ij,ij->i',np.cross(l1,s1gyg,axisa=1,axisb=1),s1g)
        #check bounds
        far[stp<0] = 360.-far[stp<0]
        far[far>90] = far[far>90]+-180
        far[far<-90] = far[far<-90]+180
        return far
    def calculate_fold_limb_rotation(self,points,axis):
        s0g = points[:,3:]
        s1g = self.get_gx(points[:,:3],grad=True)
        s1g /= np.linalg.norm(s1g,axis=1)[:,None]
        s1 = self.get_gx(points[:,:3],grad=False)
        
        projected_s0 = s0g - np.einsum('ij,ij->i',axis,s0g)[:,None]*s0g
        projected_s1 = s1g - np.einsum('ij,ij->i',axis,s1g)[:,None]*s1g
        projected_s0/=np.linalg.norm(projected_s0,axis=1)[:,None]
        projected_s1/=np.linalg.norm(projected_s1,axis=1)[:,None]
        r2 = np.einsum('ij,ij->i',projected_s1,projected_s0)#s1g,s0g)#

        vv = np.cross(s1g,s0g,axisa=1,axisb=1)
        ds = np.einsum('ik,ij->i',axis,vv)
        flr = np.where(ds>0, np.rad2deg(np.arcsin(r2)), (- np.rad2deg(np.arcsin(r2))))
        flr = np.where(flr<-90, (180.+flr),flr)
        flr = np.where(flr>90, -(180.-flr),flr)
        return flr
    def calculate_intersection_lineation(self,points):
        s1g = self.get_gx(points[:,:3],grad=True)
        s1g /= np.linalg.norm(points[:,:3],axis=1)[:,None]
        s0g = points[:,3:]
        s0g /= np.linalg.norm(s0g,axis=1)[:,None]
        l1 = np.cross(s1g,s0g,axisa=1,axisb=1)
        l1 /= np.linalg.norm(l1,axis=1)[:,None]
        return l1
