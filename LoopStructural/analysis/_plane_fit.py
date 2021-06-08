import numpy as np

from LoopStructural.utils import getLogger
logger = getLogger(__name__)
def plane_fit(points):
    
    soln = []
#     points = points.loc[~np.isnan(points['val']),['X','Y','Z','val']]
    A=np.zeros((len(points)*4,4))
    B=np.zeros((len(points)*4))
    val_mask = ~np.isnan(points['val'])
    vec_mask = ~np.isnan(points['nx'])
    c = 0
    xyz = points.loc[val_mask,['X','Y','Z']].to_numpy()
    A[:len(points[val_mask]),:3] = xyz[:,:3]
    A[:len(points[val_mask]),3]=1.
    B[:len(points[val_mask])] = points.loc[val_mask,'val'].to_numpy()#v#np.zeros(A.shape[0])
    c = len(points[val_mask])
    A[c:len(points[vec_mask])+c,0] = 1
    B[c:len(points[vec_mask])+c] = points.loc[vec_mask,'nx'].to_numpy()
    c+=len(points[vec_mask])
    A[c:len(points[vec_mask])+c,1] = 1
    B[c:len(points[vec_mask])+c] = points.loc[vec_mask,'ny'].to_numpy()
    c+=len(points[vec_mask])
    A[c:len(points[vec_mask])+c,2] = 1
    B[c:len(points[vec_mask])+c] = points.loc[vec_mask,'nz'].to_numpy()
    c+=len(points[vec_mask])
    x = np.linalg.lstsq(A[:c,:],B[:c],rcond=None)[0]
    soln.append(x)
#     print(B)
    return soln

class PlaneFitFeature:
    def __init__(self,data,name='Plane',model=None):
        self.soln = np.array(plane_fit(data)).T
        self.name = name
        self.model = model
    def evaluate_value(self,xyz):
        if self.model:
            xyz = self.model.rescale(xyz,inplace=False)
        return xyz[:,0]*self.soln[None,0]+xyz[:,1]*self.soln[None,1]+xyz[:,2]*self.soln[None,2]+self.soln[None,3]
    def min(self):

        return -np.inf
    def max(self):
        return np.inf
    def evaluate_gradient(self,xyz):
        grad = np.zeros_like(xyz)
        grad[:,0] = self.soln[None,0]
        grad[:,1] = self.soln[None,1]
        grad[:,2] = self.soln[None,2]
        return grad
