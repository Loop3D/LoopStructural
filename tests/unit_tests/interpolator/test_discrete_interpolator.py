from generate_interpolator import generate_data, generate_finite_difference_interpolator
import pytest
def test_nx():
    pass
@pytest.mark.parameterize("interpolator",[('FDI','PLI')])
def test_region(interpolator):
    """Test to see whether restricting the interpolator to a region works
    """
    data  = generate_data()
    interpolator = generate_finite_difference_interpolator(interpolator)
    interpolator.set_value_constraints(data[['X','Y','Z','val','w']].to_numpy())
    interpolator._setup_interpolator()
    region_func = lambda xyz: xyz[:,0] > 0.5
    interpolator.set_region(region_func) 
    interpolator.solve_system()
@pytest.mark.parameterize("interpolator",[('FDI','PLI')])
def test_add_constraint_to_least_squares(interpolator):
    """make sure that when incorrect sized arrays are passed it doesn't get added
    """
    interpolator = generate_finite_difference_interpolator(interpolator)
    interpolator.add_constraints_to_least_squares(A,B,idc,name='test')

def test_update_interpolator():
    pass