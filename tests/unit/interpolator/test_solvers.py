import numpy as np
import pandas as pdB
from LoopStructural.interpolators import FiniteDifferenceInterpolator as FDI
from LoopStructural.interpolators import StructuredGrid
import pytest
def test_FDI(interpolator,data):
    interpolator.set_value_constraints(data.loc[~data['val'].isna(),['X','Y','Z','val','w']].to_numpy())
    interpolator.set_gradient_constraints(data.loc[~data['nx'].isna(),['X','Y','Z','nx','ny','nz','w']].to_numpy())

    interpolator._setup_interpolator()
    interpolator.solve_system()
    if np.sum(~data['val'].isna())>0:
        assert np.sum(interpolator.evaluate_value(data.loc[~data['val'].isna(),['X','Y','Z']].to_numpy())-data.loc[~data['val'].isna(),['val']].to_numpy())/len(data) < 0.5
    else:
        pass #assert np.mean(np.einsum('ij,ij->i',interpolator.evaluate_gradient(data.loc[~data['nx'].isna(),['X','Y','Z']].to_numpy())-data.loc[~data['nx'].isna(),['nx','ny','nz']].to_numpy())) < 0.2

@pytest.mark.skip('not currently working 29-09-2022')
def test_inequality_FDI(interpolator,data):
    try:
        import osqp
    except ImportError:
        print('osqp not installed')
        return
    data['l'] = -3
    data['u'] = 10
    interpolator.set_value_constraints(data.loc[~data['val'].isna(),['X','Y','Z','val','w']].to_numpy())
    interpolator.set_gradient_constraints(data.loc[~data['nx'].isna(),['X','Y','Z','nx','ny','nz','w']].to_numpy())    
    interpolator.set_inequality_constraints(data[['X','Y','Z','l','u']].to_numpy())
    interpolator._setup_interpolator()

    interpolator.solve_system(solver='osqp')

    # print(np.sum(interpolator.evaluate_value(data[['X','Y','Z']].to_numpy())-data[['val']].to_numpy())/len(data))
    # assert np.sum(interpolator.evaluate_value(data[['X','Y','Z']].to_numpy())-data[['val']].to_numpy())/len(data) < 0.5
@pytest.mark.skip('not currently working 29-09-2022')
def test_inequality_FDI_nodes(interpolator,data):
    try:
        import osqp
    except ImportError:
        print('osqp not installed')
        return
    data['l'] = -3
    data['u'] = 10
    randind = np.arange(0,len(data))

    interpolator.set_value_constraints(data.loc[~data['val'].isna(),['X','Y','Z','val','w']].to_numpy())
    interpolator.set_gradient_constraints(data.loc[~data['nx'].isna(),['X','Y','Z','nx','ny','nz','w']].to_numpy())
    interpolator.set_inequality_constraints(data[['X','Y','Z','l','u']].to_numpy())
    interpolator._setup_interpolator()
    interpolator.solve_system(solver='osqp')
    
def test_equality_FDI_nodes(interpolator,data):
    interpolator.set_value_constraints(data.loc[~data['val'].isna(),['X','Y','Z','val','w']].to_numpy())
    interpolator.set_gradient_constraints(data.loc[~data['nx'].isna(),['X','Y','Z','nx','ny','nz','w']].to_numpy())

    node_idx = np.arange(0,interpolator.nx)[interpolator.support.nodes[:,2]>.9]
    interpolator.add_equality_constraints(node_idx,np.ones(node_idx.shape[0]),name='top')
    interpolator._setup_interpolator()
    interpolator.solve_system(solver='cg')

if __name__ == '__main__':
    test_inequality_FDI()
    test_inequality_FDI_nodes()
    test_FDI()
    test_equality_FDI_nodes()