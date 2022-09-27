import numpy as np
import pandas as pdB
from LoopStructural.interpolators import FiniteDifferenceInterpolator as FDI
from LoopStructural.interpolators import StructuredGrid
from generate_interpolator import generate_interpolator, generate_data
import pytest
@pytest.mark.parametrize("interpolator",[('FDI'),('PLI')])
def test_FDI(interpolator):
    data = generate_data()
    interpolator = generate_interpolator(interpolator)
    interpolator.set_value_constraints(data[['X','Y','Z','val','w']].to_numpy())
    interpolator._setup_interpolator()
    interpolator.solve_system()
    assert np.sum(interpolator.evaluate_value(data[['X','Y','Z']].to_numpy())-data[['val']].to_numpy())/len(data) < 0.5

@pytest.mark.parametrize("interpolator",[('FDI'),('PLI')])
def test_inequality_FDI(interpolator):
    try:
        import osqp
    except ImportError:
        print('osqp not installed')
        return
    data = generate_data()
    data['l'] = -3
    data['u'] = 10
    randind = np.arange(0,len(data))
    interpolator = generate_interpolator(interpolator)
    interpolator.set_value_constraints(data[['X','Y','Z','val','w']].to_numpy())
    interpolator.set_inequality_constraints(data[['X','Y','Z','l','u']].to_numpy())
    interpolator._setup_interpolator()

    interpolator.solve_system(solver='osqp')

    # print(np.sum(interpolator.evaluate_value(data[['X','Y','Z']].to_numpy())-data[['val']].to_numpy())/len(data))
    # assert np.sum(interpolator.evaluate_value(data[['X','Y','Z']].to_numpy())-data[['val']].to_numpy())/len(data) < 0.5
@pytest.mark.parametrize("interpolator",[('FDI'),('PLI')])
def test_inequality_FDI_nodes(interpolator):
    try:
        import osqp
    except ImportError:
        print('osqp not installed')
        return
    data = generate_data()
    data['l'] = -3
    data['u'] = 10
    randind = np.arange(0,len(data))

    interpolator = generate_interpolator(interpolator)
    interpolator.set_value_constraints(data[['X','Y','Z','val','w']].to_numpy())
    interpolator.set_inequality_constraints(data[['X','Y','Z','l','u']].to_numpy())
    interpolator._setup_interpolator()
    interpolator.solve_system(solver='osqp')
    
@pytest.mark.parametrize("interpolator",[('FDI'),('PLI')])
def test_equality_FDI_nodes(interpolator):
    data = generate_data()
    interpolator = generate_interpolator(interpolator)
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