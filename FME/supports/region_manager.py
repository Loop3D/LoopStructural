import numpy as np
from sklearn.decomposition import PCA
class RegionManager:
    def __init__(self,mesh):
        self.mesh = mesh
    def create_region_from_cuboid(self,corners,name):
        """
        corners[0,:] should be minx,miny,minz
        corners[1,:] should be maxx,miny,minz
        corners[2,:] should be minx,maxy,minz
        corners[3,:] should be minx,miny,maxz
        
        """
        u = corners[0,:] - corners[1,:]
        v = corners[0,:] - corners[2,:]
        w = corners[0,:] - corners[3,:]
        
        ux = np.einsum('ij,j->i',self.mesh.nodes,u)
        vx = np.einsum('ij,j->i',self.mesh.nodes,v)
        wx = np.einsum('ij,j->i',self.mesh.nodes,w)
        vp1 = np.dot(v,corners[0,:])
        up1 = np.dot(u,corners[0,:])

        up2 = np.dot(u,corners[1,:])
        vp4 = np.dot(v,corners[2,:])
        wp1 = np.dot(w,corners[0,:])
        wp5 = np.dot(w,corners[3,:])
        logic = np.logical_and(ux<up1,ux>up2)#(condition = []
        logic = np.logical_and(logic,vx<vp1)
        logic = np.logical_and(logic,vx>vp4)
        logic = np.logical_and(logic,wx<wp1)
        logic = np.logical_and(logic,wx>wp5)
        region = np.zeros(self.mesh.n_nodes).astype(bool)
        region[logic] = 1
        self.mesh.regions[name] = region
        
    def create_region_from_object_aligned_box(self,points,name,wscale=1.,hscale=1.,lscale=1.):
        """
        calculates PCA of the points, then calculates bounding box, applies w scale and h scale
        then reporjects back into cartesian space and creates a region
        """
        if points.shape[0] < 3:
            return False
        pca = PCA(n_components=3)
        pca.fit(points)
        
        newp = pca.transform(points)
        
        corners=np.zeros((4,3))
        corners[0,0] = np.min(newp[:,0])
        corners[0,1] = np.min(newp[:,1])
        corners[0,2] = np.min(newp[:,2])

        corners[1,0] = np.max(newp[:,0])
        corners[1,1] = np.min(newp[:,1])
        corners[1,2] = np.min(newp[:,2])

        corners[2,0] = np.min(newp[:,0])
        corners[2,1] = np.max(newp[:,1])
        corners[2,2] = np.min(newp[:,2])

        corners[3,0] = np.min(newp[:,0])
        corners[3,1] = np.min(newp[:,1])
        corners[3,1] = np.max(newp[:,2])
        
        #transform back
        cornerst = pca.inverse_transform(corners)
        self.create_region_from_cuboid(cornerst,name)
    def create_region_from_map_object_aligned_bounding_box(self,points,minz,maxz,name,pc1buffer=0.,pc2buffer=0.):
        pca = PCA(n_components=2)
        pca.fit(points)
        
        newp = pca.transform(points)
        
        corners=np.zeros((4,3))
        corners[0,0] = np.min(newp[:,0])-pc1buffer
        corners[0,1] = np.min(newp[:,1])-pc2buffer
        corners[0,2] = minz

        corners[1,0] = np.max(newp[:,0])+pc1buffer
        corners[1,1] = np.min(newp[:,1])-pc2buffer
        corners[1,2] = minz
        
        corners[2,0] = np.min(newp[:,0])-pc1buffer
        corners[2,1] = np.max(newp[:,1])+pc2buffer
        corners[2,2] = minz

        corners[3,0] = np.min(newp[:,0])-pc1buffer
        corners[3,1] = np.min(newp[:,1])-pc2buffer
        corners[3,2] = maxz
        #transform back
        cornerst=np.zeros((4,3))
        cornerst[:,:2] = pca.inverse_transform(corners[:,:2])
        cornerst[:,2] = corners[:,2]
        
        self.create_region_from_cuboid(cornerst,name)
    def create_region_from_boundary_box(self,boundary_points,name):
        corners = np.zeros((4,3))
        corners[0,:] = np.array([boundary_points[0,0],boundary_points[0,1],boundary_points[0,2]])
        corners[1,:] = np.array([boundary_points[0,0],boundary_points[1,1],boundary_points[0,2]])
        corners[2,:] = np.array([boundary_points[1,0],boundary_points[0,1],boundary_points[0,2]])
        corners[3,:] = np.array([boundary_points[0,0],boundary_points[0,1],boundary_points[1,2]])
        self.create_region_from_cuboid(corners,name)
    
    def create_region_from_property_value(self,propertyname,propertyvalue,region_name,sign=1):
        region = np.zeros(mesh.n_nodes).astype(bool)
        region[mesh.properties[propertyname]*sign>propertyvalue*sign] = 1
        self.mesh.regions[region_name] = region
    def create_properties_for_regions(self):
        for region in self.mesh.regions.keys():
            self.mesh.properties['REGION_%s'%region] = self.mesh.regions[region].astype(float)
