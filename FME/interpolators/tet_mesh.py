import meshpy.tet
from meshpy import *
#from sympy import Symbol, Matrix
#from sympy.utilities.autowrap import autowrap
import numpy as np
from numpy import linalg as la
from cdot import dot
from dsi_helper import get_element, compute_cg_regularisation_constraint, pointintetra, cg, cg_cstr_to_coo_sq
from mesh_helper import calculate_tetra
from scipy.spatial import cKDTree
import os
class TetMesh:
    def __init__(self,name='TetMesh',path=os.environ['FaultDATA']+'/'):
        self.path=path
        self.name=name
        self.mesh = None
        self.shared_idxs = np.zeros(3,dtype=np.int)
        self.dinfo = {}
        self.cg_calculated = {}
        self.m_inv = autofunc_c#autowrap(minv,backend="cython",args=(x0,x1,x2,x3,y0,y1,y2,y3,z0,z1,z2,z3))
        self.properties = {}
        self.property_gradients = {}
        self.property_gradients_nodes = {}
        self.element_gradients = {}
        self.element_nodes_to_id = {}
        self.regions = {}
        
    def setup_mesh(self,boundary_points,**kwargs):#maxvol=0.001,nstep=100):
       
        minx = boundary_points[0,0]
        miny = boundary_points[0,1]
        minz = boundary_points[0,2]
        maxx = boundary_points[1,0]
        maxy = boundary_points[1,1]
        maxz = boundary_points[1,2]
        points = [
        (minx, miny, minz),
        ( maxx, miny, minz),
        ( maxx,  maxy, minz),
        (minx,  maxy, minz),
        (minx, miny,  maxz),
        ( maxx, miny,  maxz),
        ( maxx,  maxy,  maxz),
        (minx,  maxy,  maxz)
        ]
        n_tetra = 4000#maxvol=0.001
        if 'n_tetra' in kwargs:
            n_tetra = kwargs['n_tetra']
        #calculate the volume of the bounding box
        maxvol = 0.001
        if 'maxvol' not in kwargs:
            lengthU = maxx-minx
            lengthV = maxy-miny
            lengthW = maxz-minz
            boxVol = lengthU*lengthW*lengthV
            correction_factor = 1.91
            maxvol = correction_factor*boxVol / n_tetra 
        if 'maxvol' in kwargs:
            maxvol = kwargs['maxvol']
        
        facets = [
                [0, 1, 2, 3],
                [4, 5, 6, 7],
                [0, 4, 5, 1],
                [1, 5, 6, 2],
                [2, 6, 7, 3],
                [3, 7, 4, 0]
                ]

        # create the mesh
        info = meshpy.tet.MeshInfo()
        info.set_points(points)
        info.set_facets(facets)
        meshpy_mesh = meshpy.tet.build(info, max_volume=maxvol,options=meshpy.tet.Options('pqn'))
        self.nodes = np.array(meshpy_mesh.points)
        self.elements = np.array(meshpy_mesh.elements)
        self.neighbours = np.array(meshpy_mesh.neighbors)
        self.n_nodes = len(self.nodes)
        self.n_elements = len(self.elements)
        self.barycentre = np.sum(self.nodes[self.elements][:,:,:],axis=1) / 4.
        self.tree = cKDTree(self.barycentre)
        self.minx = minx
        self.miny = miny
        self.minz = minz
        self.maxx = maxx
        self.maxy = maxy
        self.maxz = maxz
        self.regions['everywhere'] = np.ones(self.n_nodes).astype(bool)
    def add_region(self,region,name):
        self.regions[name] = region
        self.properties['REGION_'+name] = region.astype(float)
    def update_property(self,name,value,save=True):
        self.properties[name] = value
        self
        #grad = np.zeros((self.elements.shape[0],3))
        grads = self.get_elements_gradients(np.arange(self.n_elements))
        props = self.properties[name][self.elements[np.arange(self.n_elements)]]
        grad = np.einsum('ikj,ij->ik',grads,props)
        self.property_gradients[name] = grad
        if save:
            self.save()
        #//self.property_v = p
    def transfer_gradient_to_nodes(self,propertyname):
        grad = np.zeros((self.nodes.shape[0],3))
        for i, e in enumerate(self.elements):
            for n in e:
                grad[n,:]+=self.property_gradients[propertyname][i,:]
        grad/=4
        self.property_gradients_nodes[propertyname] = grad
            
    def get_neighbours(self,t):
        return self.neighbours[t]
    def get_neighbors(self,t):
        return self.get_neighbours(t) #spelling >.<
    def get_shared_face(self,t1,t2):
        self.shared_idxs[:] = 0
        i = 0
        for p in t1: #for nodes in tetra1
            if p in t2: #if node in tetra2
                self.shared_idxs[i] = p
                i+=1
        shared_points = self.nodes[self.shared_idxs[:i]] #find the locations of shared points
        v1 = shared_points[0,:]-shared_points[1,:]
        v2 = shared_points[2,:]-shared_points[1,:]
        n = np.cross(v1,v2)#v1 * v2
        return n,shared_points,self.shared_idxs
    def calculate_constant_gradient(self,region,shape='rectangular'):
        self.dinfo = np.zeros(self.n_elements).astype(bool)
        A = []
        B = []
        row = []
        col = []
        c_ = 0
        
        #add cg constraint for all of the 
        EG = self.get_elements_gradients(np.arange(self.n_elements))
        idc, c,ncons = cg(EG,self.neighbours,self.elements,self.nodes)
        idc=np.array(idc)
        c = np.array(c)
        if shape == 'rectangular':
            for i in range(ncons):
                for j in range(5):
                    A.append(c[i,j])
                    row.append(i)
                    col.append(idc[i,j])
                B.append(0.)
        c_=ncons            
        if shape == 'square':
            A, row, col = cg_cstr_to_coo_sq(c,idc,ncons)
            A = np.array(A).tolist()
            row = np.array(row).tolist()
            col = np.array(col).tolist()
        self.cg_calculated[shape] = True
        self.cg = [A,B,row,col,c_]
        return True
    def get_constant_gradient(self,w=0.1,region='everywhere',**kwargs):
        shape = 'rectangular'
        if 'shape' in kwargs:
            shape = 'square'
        self.calculate_constant_gradient(region,**kwargs)
        return self.cg[0],self.cg[1],self.cg[2],self.cg[3],self.cg[4]
    def get_element(self,p):
        d,e = self.tree.query(p)
        return self.elements[e]
    def get_element_gradient(self,t):
        tu = tuple(t)
        #if tu in self.element_gradients:
        #    return self.element_gradients[tu]
        ps = self.nodes[t]
        ps -= ps[0,:]
        m= np.array(
                    [[(ps[1,0]-ps[0,0]),(ps[1,1]-ps[0,1]),(ps[1,2]-ps[0,2])],
                    [(ps[2,0]-ps[0,0]),(ps[2,1]-ps[0,1]),(ps[2,2]-ps[0,2])],
                    [(ps[3,0]-ps[0,0]),(ps[3,1]-ps[0,1]),(ps[3,2]-ps[0,2])]])
        I = np.array(
                  [[-1.,1.,0.,0.],
                   [-1.,0.,1.,0.],
                   [-1.,0.,0.,1.]]) 
        return la.inv(m) @ I
        tu = tuple(t)
        if tu in self.element_gradients:
            return self.element_gradients[tu]
        ps = self.nodes[t]
        ps -= ps[0,:]
        #equation 14 mallet2003 gocad
        J = np.ones((3,1))
        O = np.zeros((1,3))
        m = np.zeros((3,3))
        for j in range(3):
            for i in range(3):
                m[i,j] = ps[i+1,j]
        Minv = np.zeros((4,4))
        Minv[0,0] = 1
        Minv[0,1:] = 0
        minv = la.inv(m)
        Minv[1:,0] = (-minv@J)[:,0]
        Minv[1:,1:] = minv
        I = np.zeros((3,4))
        for i in range(3):
            I[i,i+1] = 1
        return I@Minv
    def get_elements_gradients(self,e):
        """
        Returns the gradient of the elements using their global index.
        Speed up using numpy
        """
        ps = self.nodes[self.elements[e]]
        ps -= ps[:,0,:][:,None]
        J = np.ones((len(e),3,1))
        O = np.ones((len(e),1,3))
        m = ps[:,1:,:]
        Minv = np.zeros((len(e),4,4))
        Minv[:,0,0] = 1
        Minv[:,0,1:] = 0
        minv = np.linalg.inv(m)
        Minv[:,1:,0] = (-minv@J)[:,:,0]
        Minv[:,1:,1:] = minv
        I = np.zeros((3,4))
        I[np.arange(3),np.arange(3)+1] = 1
        return I@Minv
    def get_property_value(self,pos=None,element=None,prop=None):
        if element is None:
            if pos is None: 
                print('element and pos None')
                return False
            element = self.get_element(pos)
            if element is None:
                print('Could not find triangle for x:%f y:%f z:%f'%(p.pos[0],p.pos[1],p.pos[2]))
                return False
        points = self.nodes[element]
        c = np.zeros(4)
        #if no position is specified just use the value at the centroid of the element
        if pos is None:
            c[:] = 0.25
            return np.sum(c*self.properties[prop][element])
            #pos = np.mean(points,axis=0)
        M = np.ones((points.shape[0],points.shape[0]))

        for i in range(4):              
            M[i,1] = points[i,0]
            M[i,2] = points[i,1]
            M[i,3] = points[i,2]
        cp = np.ones(points.shape[0])
        cp[1:] = pos[:]
        c = np.dot(cp,la.inv(M))
        #print(c)
        c[:] = 0.25
        return np.sum(c*self.properties[prop][element])
    def get_element_property_gradient(self,element=None,pos=None,prop=None):
        if element is None:
            if pos is None: 
                return False
            element = self.get_element(pos)
            if element is None:
                print('Could not find triangle for x:%f y:%f z:%f'%(p.pos[0],p.pos[1],p.pos[2]))
                return False
        if tuple(element) in self.element_nodes_to_id:
            if prop in self.property_gradients:
                return self.property_gradients[prop][self.element_nodes_to_id[tuple(element)]]
        d_n = self.get_element_gradient(element)#calculate_d(n_points)
        grad = np.sum(d_n*self.properties[prop][element],axis=1)
        grad/=4.
        return grad
    def calc_bary_c(self,e,p):
        """
        Calculate the barycentre coordinates for an array of n elements 
        and n points
        """ 
        points = self.nodes[self.elements[e]]
        npts = len(e)
        vap = np.zeros((3,npts))
        vbp = np.zeros((3,npts))
    
        vcp = np.zeros((3,npts))
        vdp = np.zeros((3,npts))
        vab = np.zeros((3,npts))
        vac = np.zeros((3,npts))
        vad = np.zeros((3,npts))
        vbc = np.zeros((3,npts))
        vbd = np.zeros((3,npts))
        ##bp = 
        vap= p - points[:,0,:]
        vbp= p - points[:,1,:]
        vcp= p - points[:,2,:]
        vdp= p - points[:,3,:]
        vab= points[:,1,:] - points[:,0,:]
        vac= points[:,2,:] - points[:,0,:] 
        vad= points[:,3,:] - points[:,0,:] 
        vbc= points[:,2,:] - points[:,1,:] 
        vbd= points[:,3,:] - points[:,1,:] 
        
        va = np.sum(vbp * np.cross(vbd,vbc,axisa=1,axisb=1), axis=1) /6.
        vb = np.sum(vap * np.cross(vac,vad,axisa=1,axisb=1), axis=1) /6.
        vc = np.sum(vap * np.cross(vad,vab,axisa=1,axisb=1), axis=1) /6.
        vd = np.sum(vap * np.cross(vab,vac,axisa=1,axisb=1), axis=1) /6.
        v = np.sum(vab * np.cross(vac,vad,axisa=1,axisb=1), axis=1) /6.
        c = np.zeros((4,npts))
        c[0,:] = va / v
        c[1,:] = vb / v
        c[2,:] = vc / v
        c[3,:] = vd / v
        return c
    def elements_for_array(self,array,k=10):
        #\TODO could add check to see which points are inside mesh bounding box
        #find the nearest k elements for the arrays
        #reducing k could speed up calculation but migh add errors
        
        d,ee = self.tree.query(array) 
        #d,e = mesh.tree.query(mesh2.nodes,k=k)
        #make a container to find out if the nearest element barycenter actually contains
        #the point
        inside = np.zeros(array.shape[0]).astype(bool)
        inside[:] = True
        inside = array[:,0] < self.maxx[None] 
        inside+= array[:,0] > self.minx[None]
        inside+= array[:,1] > self.miny[None]
        inside+= array[:,1] < self.maxy[None]
        inside+= array[:,2] < self.maxz[None]
        inside+= array[:,2] > self.minz[None]
        #pind = np.array([(array[:,0]-self.origin[0]) // self.step,\
        #                 (array[:,1]-self.origin[1]) // self.step,\
        #                 (array[:,2]-self.origin[2]) // self.step \
        #                ]).astype(int).T
        #    #ix,iy = np.where(np.any(indx[:,None]==pind,axis=2)==True)
        #elements = self.element_map[pind[inside,0],pind[inside,1],pind[inside,2],:]
        #nodes = self.nodes[self.elements[np.array(range(0,self.n_elements))]]
        #elements = np.array(calculate_tetra(elements,array,nodes))
        return ee, inside
        #for i in range(k):
        #    #now loop over the possibilities, if all bc coords > 0 then inside
        #    inside[ee[:,i]<self.n_elements,i] = np.all(self.calc_bary_c(ee[ee[:,i]\
        #          <self.n_elements,i],array[ee[:,i]<len(self.elements),:])>0,axis=0)
        ##find indices for e that are true
        #inside[ee==len(self.elements)] == False
        #ix,iy = np.where(inside==True)
        #e = np.zeros(ee.shape[0]).astype(int)
        #e[ix] = ee[ix,iy]
        ##
        ###just update to check if none of the 5 nearest tetra contian point
        #inside=np.any(inside,axis=1)#inside[np.any
        #return ee,inside
    def elements_for_array_old(self,array,k):
        #\TODO could add check to see which points are inside mesh bounding box
        #find the nearest k elements for the arrays
        #reducing k could speed up calculation but migh add errors
       
        d,ee = self.tree.query(array,k) 
        #d,e = mesh.tree.query(mesh2.nodes,k=k)
        #make a container to find out if the nearest element barycenter actually contains
        #the point
        inside = np.zeros(ee.shape).astype(bool)
        for i in range(k):
            #now loop over the possibilities, if all bc coords > 0 then inside
            inside[ee[:,i]<self.n_elements,i] = np.all(self.calc_bary_c(ee[ee[:,i]\
                  <self.n_elements,i],array[ee[:,i]<len(self.elements),:])>0,axis=0)
        #find indices for e that are true
        inside[ee==len(self.elements)] == False
        ix,iy = np.where(inside==True)
        e = np.zeros(ee.shape[0]).astype(int)
        e[ix] = ee[ix,iy]
        ##
        ###just update to check if none of the 5 nearest tetra contian point
        inside=np.any(inside,axis=1)#inside[np.any
        return e,inside   
    def eval_interpolant(self,array,prop,k=5,e=None,region='everywhere'):
        """Evaluate an interpolant from property on an array of points.
        Uses numpy to speed up calculations but could be expensive for large mesh/points
        """
        if e == None:
            e, inside = self.elements_for_array(array,k)
        else:
            inside = np.array(e.shape).astype(bool)
            inside[:] = True

        bc = self.calc_bary_c(e[inside],array[inside])
        prop_int = np.zeros(e.shape)

        props = self.properties[prop][self.elements[e[inside]]]
        prop_int[inside]=np.sum((bc.T*props),axis=1)
        prop_int[~inside]=np.nan
        return prop_int
    def element_property_value(self,prop):
        bc = np.zeros(4)
        bc[:] = 0.25
        e = np.arange(self.n_elements)
        prop_int = np.zeros(e.shape)
        
        props = self.properties[prop][self.elements[e]]
        prop_int=np.sum((bc.T*props),axis=1)
        return prop_int
    def element_property_gradient(self,prop):
        e = np.arange(self.n_elements)
        
        ps = self.nodes[self.elements[e]]
        m= np.array(
        [[(ps[:,1,0]-ps[:,0,0]),(ps[:,1,1]-ps[:,0,1]),(ps[:,1,2]-ps[:,0,2])],
            [(ps[:,2,0]-ps[:,0,0]),(ps[:,2,1]-ps[:,0,1]),(ps[:,2,2]-ps[:,0,2])],
            [(ps[:,3,0]-ps[:,0,0]),(ps[:,3,1]-ps[:,0,1]),(ps[:,3,2]-ps[:,0,2])]])
        I = np.array(
                  [[-1.,1.,0.,0.],
                   [-1.,0.,1.,0.],
                   [-1.,0.,0.,1.]])
        m = np.swapaxes(m,0,2)
        grads = la.inv(m)

        grads = grads.swapaxes(1,2)
        grads = grads @ I
        vals = self.properties[prop][self.elements[e]]
        #grads = np.swapaxes(grads,1,2)
        a = np.zeros((self.n_elements,3))#array.shape)
        a = (grads*vals[:,None,:]).sum(2)/4.# @ vals.T/
        a /= np.sum(a,axis=1)[:,None]
        return a
        

    def eval_gradient(self,array,prop,k=5,region='everywhere'):
        """Evaluate an interpolant from property on an array of points.
        Uses numpy to speed up calculations but could be expensive for large mesh/points
        """    
        e, inside = self.elements_for_array(array,k)
        ps = self.nodes[self.elements[e[inside]]]
        m= np.array(
        [[(ps[:,1,0]-ps[:,0,0]),(ps[:,1,1]-ps[:,0,1]),(ps[:,1,2]-ps[:,0,2])],
            [(ps[:,2,0]-ps[:,0,0]),(ps[:,2,1]-ps[:,0,1]),(ps[:,2,2]-ps[:,0,2])],
            [(ps[:,3,0]-ps[:,0,0]),(ps[:,3,1]-ps[:,0,1]),(ps[:,3,2]-ps[:,0,2])]])
        I = np.array(
                  [[-1.,1.,0.,0.],
                   [-1.,0.,1.,0.],
                   [-1.,0.,0.,1.]])
        m = np.swapaxes(m,0,2)
        grads = la.inv(m)

        grads = grads.swapaxes(1,2)
        grads = grads @ I
        vals = self.properties[prop][self.elements[e[inside]]]
        #grads = np.swapaxes(grads,1,2)
        a = np.zeros(array.shape)
        a[inside] = (grads*vals[:,None,:]).sum(2)/4.# @ vals.T/
        a[~inside,:]=np.nan
        a /= np.sum(a,axis=1)[:,None]

        return a
        
    def export_to_vtk(self,name='mesh.vtk'):
        import meshio
        meshio.write_points_cells(name,
                          self.nodes,{"tetra": self.elements},
                          point_data=self.properties,
                         cell_data={'tetra':self.property_gradients} 
                         )
        
        meshio.write_points_cells('test.vtk',
                          self.nodes,{"tetra": self.elements},
                          point_data=self.property_gradients_nodes
                         )
    def save(self):
        self.export_to_vtk(self.path+self.name+'.vtk')
