import numpy as np
import pandas as pd
from LoopStructural.interpolators import FiniteDifferenceInterpolator as FDI, PiecewiseLinearInterpolator as PLI
from LoopStructural.interpolators import StructuredGrid,TetMesh
def generate_data(value=True,gradient=False):
    data_list = []
    if value:
        xy = np.array(np.meshgrid(np.linspace(0,1,50),np.linspace(0,1,50))).T.reshape(-1,2)
        xyz = np.hstack([xy,np.zeros((xy.shape[0],1))])
        data = pd.DataFrame(xyz,columns=['X','Y','Z'])
        data['val'] = np.sin(data['X'])
        data['w'] = 1
        data['feature_name'] = 'strati'
        data_list.append(data)
    if gradient:
        data = pd.DataFrame([[0.5,0.5,0.5,0,0,1]],columns=['X','Y','Z','nx','ny','nz'])
        data['w'] = 1
        data['feature_name'] = 'strati'
        data_list.append(data)
    return pd.concat(data_list,ignore_index=True)

def generate_interpolator(interpolator='FDI'):
    origin = np.array([-0.1,-0.1,-0.1])
    maximum = np.array([1.1,1.1,1.1])
    nsteps = np.array([20,20,20])
    step_vector = (maximum-origin)/nsteps
    if type == 'FDI':
        grid = StructuredGrid(origin=origin,nsteps=nsteps,step_vector=step_vector)
        interpolator = FDI(grid)
        return interpolator
    if type == 'PLI':
        grid = TetMesh(origin=origin,nsteps=nsteps,step_vector=step_vector)
        interpolator = PLI(grid)
        return interpolator