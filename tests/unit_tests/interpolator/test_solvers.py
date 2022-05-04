import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from LoopStructural.interpolators import FiniteDifferenceInterpolator as FDI
from LoopStructural.interpolators import StructuredGrid
import matplotlib.pyplot as plt
# def test_FDI():
#     xy = np.array(np.meshgrid(np.linspace(0,1,50),np.linspace(0,1,50))).T.reshape(-1,2)
#     xyz = np.hstack([xy,np.zeros((xy.shape[0],1))])
#     data = pd.DataFrame(xyz,columns=['X','Y','Z'])
#     data['val'] = np.sin(data['X'])
#     data['w'] = 1
#     data['feature_name'] = 'strati'
#     randind = np.arange(0,len(data))
#     np.random.shuffle(randind)
#     data.loc[randind[:int(50*50*.5)],'val'] = np.nan
#     np.random.shuffle(randind)
#     data.loc[randind[:int(50*50*.1)],'val'] = 0
#     origin = np.array([-0.1,-0.1,-0.1])
#     maximum = np.array([1.1,1.1,1.1])
#     nsteps = np.array([10,10,10])
#     step_vector = (maximum-origin)/nsteps
#     grid = StructuredGrid(origin=origin,nsteps=nsteps,step_vector=step_vector)
#     interpolator = FDI(grid)
#     interpolator.set_value_constraints(data[['X','Y','Z','val','w']].to_numpy())
#     interpolator._setup_interpolator()
#     interpolator.solve_system()
#     print(interpolator.constraints['value'])

def test_FDI():
    xy = np.array(np.meshgrid(np.linspace(0,1,50),np.linspace(0,1,50))).T.reshape(-1,2)
    xyz = np.hstack([xy,np.zeros((xy.shape[0],1))])
    data = pd.DataFrame(xyz,columns=['X','Y','Z'])
    data['val'] = np.sin(data['X'])
    data['w'] = 1
    data['feature_name'] = 'strati'
    randind = np.arange(0,len(data))
    # np.random.shuffle(randind)
    # # data.loc[randind[:int(50*50*.5)],'val'] = np.nan
    # np.random.shuffle(randind)
    # data.loc[randind[:int(50*50*.1)],'val'] = 0
    origin = np.array([-0.1,-0.1,-0.1])
    maximum = np.array([1.1,1.1,1.1])
    nsteps = np.array([20,20,20])
    step_vector = (maximum-origin)/nsteps
    grid = StructuredGrid(origin=origin,nsteps=nsteps,step_vector=step_vector)
    interpolator = FDI(grid)
    interpolator.set_value_constraints(data[['X','Y','Z','val','w']].to_numpy())
    interpolator._setup_interpolator()
    interpolator.solve_system()
    assert np.sum(interpolator.evaluate_value(data[['X','Y','Z']].to_numpy())-data[['val']].to_numpy())/len(data) < 0.5
def test_inequality_FDI():
    try:
        import osqp
    except ImportError:
        print('osqp not installed')
        return
    xy = np.array(np.meshgrid(np.linspace(0,1,50),np.linspace(0,1,50))).T.reshape(-1,2)
    xyz = np.hstack([xy,np.zeros((xy.shape[0],1))])
    data = pd.DataFrame(xyz,columns=['X','Y','Z'])
    data['val'] = np.sin(data['X'])
    data['w'] = 1
    data['feature_name'] = 'strati'
    data['l'] = -3
    data['u'] = 10
    randind = np.arange(0,len(data))
    # np.random.shuffle(randind)
    # # data.loc[randind[:int(50*50*.5)],'val'] = np.nan
    # np.random.shuffle(randind)
    # data.loc[randind[:int(50*50*.1)],'val'] = 0
    origin = np.array([-0.1,-0.1,-0.1])
    maximum = np.array([1.1,1.1,1.1])
    nsteps = np.array([20,20,20])
    step_vector = (maximum-origin)/nsteps
    grid = StructuredGrid(origin=origin,nsteps=nsteps,step_vector=step_vector)
    interpolator = FDI(grid)
    interpolator.set_value_constraints(data[['X','Y','Z','val','w']].to_numpy())
    interpolator.set_inequality_constraints(data[['X','Y','Z','l','u']].to_numpy())
    interpolator._setup_interpolator()
    # col = np.arange(0,interpolator.nx,dtype=int)
    # col = np.tile(col, (interpolator.nx, 1)).T
    # interpolator.add_inequality_constraints_to_matrix(np.eye(interpolator.nx),
    #                                                 np.zeros(interpolator.nx)-4,
    #                                                 np.zeros(interpolator.nx)+np.inf,
    #                                                 col
    #                                                 )
    interpolator.solve_system(solver='osqp')

    # print(np.sum(interpolator.evaluate_value(data[['X','Y','Z']].to_numpy())-data[['val']].to_numpy())/len(data))
    # assert np.sum(interpolator.evaluate_value(data[['X','Y','Z']].to_numpy())-data[['val']].to_numpy())/len(data) < 0.5
def test_inequality_FDI_nodes():
    try:
        import osqp
    except ImportError:
        print('osqp not installed')
        return
    xy = np.array(np.meshgrid(np.linspace(0,1,50),np.linspace(0,1,50))).T.reshape(-1,2)
    xyz = np.hstack([xy,np.zeros((xy.shape[0],1))])
    data = pd.DataFrame(xyz,columns=['X','Y','Z'])
    data['val'] = np.sin(data['X'])
    data['w'] = 1
    data['feature_name'] = 'strati'
    data['l'] = -3
    data['u'] = 10
    randind = np.arange(0,len(data))

    origin = np.array([-0.1,-0.1,-0.1])
    maximum = np.array([1.1,1.1,1.1])
    nsteps = np.array([20,20,20])
    step_vector = (maximum-origin)/nsteps
    grid = StructuredGrid(origin=origin,nsteps=nsteps,step_vector=step_vector)
    interpolator = FDI(grid)
    interpolator.set_value_constraints(data[['X','Y','Z','val','w']].to_numpy())
    interpolator.set_inequality_constraints(data[['X','Y','Z','l','u']].to_numpy())
    interpolator._setup_interpolator()
    interpolator.solve_system(solver='osqp')

def test_equality_FDI_nodes():
    xy = np.array(np.meshgrid(np.linspace(0,1,50),np.linspace(0,1,50))).T.reshape(-1,2)
    xyz = np.hstack([xy,np.zeros((xy.shape[0],1))])
    data = pd.DataFrame(xyz,columns=['X','Y','Z'])
    data['val'] = np.sin(data['X'])
    data['w'] = 1
    data['feature_name'] = 'strati'
    origin = np.array([-0.1,-0.1,-0.1])
    maximum = np.array([1.1,1.1,1.1])
    nsteps = np.array([20,20,20])
    step_vector = (maximum-origin)/nsteps
    grid = StructuredGrid(origin=origin,nsteps=nsteps,step_vector=step_vector)
    interpolator = FDI(grid)
    interpolator.set_value_constraints(data[['X','Y','Z','val','w']].to_numpy())

    node_idx = np.arange(0,interpolator.nx)[interpolator.support.nodes[:,2]>.9]
    interpolator.add_equality_constraints(node_idx,np.ones(node_idx.shape[0]),name='top')
    interpolator._setup_interpolator()
    interpolator.solve_system(solver='cg')

if __name__ == '__main__':
    test_inequality_FDI()
    test_inequality_FDI_nodes()
    test_FDI()
    test_equality_FDI_nodes()