import pytest
from tests.generators.generate_interpolator import generate_data, generate_interpolator
@pytest.mark.parameterize("interpolator,value,grad",[('FDI',True,False),('FDI',False,True),('FDI',True,True),('FDI',True,False)])
def test_get_data_locations(interpolator,grad,value):
    data = generate_data(value,grad)
    interpolator = generate_interpolator(interpolator)
    interpolator.set_value_constraints(data.loc[~data['val'].isna(),['X','Y','Z','val','w']].to_numpy())
    interpolator.set_normal_constraints(data.loc[~data['nx'].isna(),['X','Y','Z','nx','ny','nz','w']].to_numpy())
    interpolator.get_data_locations()

