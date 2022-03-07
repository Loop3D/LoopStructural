import numpy as np

# Loop library
from LoopStructural import GeologicalModel
from LoopStructural.modelling.intrusions import IntrusionNetwork
from LoopStructural.modelling.intrusions import IntrusionBuilder
from LoopStructural.modelling.intrusions import IntrusionBody
from LoopStructural.modelling.intrusions import IntrusionFeature
from LoopStructural.modelling.features import StructuralFrame
from LoopStructural.modelling.intrusions import rectangle_function, parallelepiped_function

from LoopStructural.datasets import load_tabular_intrusion
data, boundary_points = load_tabular_intrusion()

def test_intrusion_network():
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.data = data
    model.nsteps = [10,10,10]
    
    intrusion_data = data[data['feature_name'] == 'tabular_intrusion']
    conformable_feature = model.create_and_add_foliation('stratigraphy')

    INet = IntrusionNetwork(feature_data=intrusion_data,
                            intrusion_network_contact='roof',
                            intrusion_network_type='shortest path',
                            model=model,
                            )
    delta_c = 2
    INet.set_data()
    INet.set_contact_anisotropies([conformable_feature])
    INet.set_sequence_of_exploited_anisotropies([conformable_feature])
    INet.set_velocity_parameters()
    INet.set_sections_axis('X')
    intrusion_network_points = INet.build(delta_c=[delta_c])[:,:3]
    
    #test if points lie in the contact of interest

    mean = INet.anisotropies_series_parameters['stratigraphy_0'][1]
#     mean = -10
    stdv = INet.anisotropies_series_parameters['stratigraphy_0'][2]
    evaluated_inet_points =  conformable_feature['feature'].evaluate_value(model.scale(intrusion_network_points))

    assert np.all(np.logical_and((mean - stdv*delta_c)<= evaluated_inet_points,(mean + stdv*delta_c)>= evaluated_inet_points))


def test_intrusion_body(lateral_conceptual_model, vertical_conceptual_model):
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.data = data
    model.nsteps = [10,10,10]
    
    intrusion_data = data[data['feature_name'] == 'tabular_intrusion']
    conformable_feature = model.create_and_add_foliation('stratigraphy')
    
    # create intrusion network 
    INet = IntrusionNetwork(feature_data=intrusion_data,
                        intrusion_network_contact='roof',
                        intrusion_network_type='interpolated',
                        model=model,
                        )

    INet.set_data()
    intrusion_network_points = INet.build(delta_c=[1e15])[:,:3]

    weights = [0,0,0]
    frame_data = model.data[model.data["feature_name"] == 'tabular_intrusion_frame'].copy()
    interpolator = model.get_interpolator(interpolatortype='FDI')
    IFrame_builder = IntrusionBuilder(interpolator, 
                                      model=model, 
                                      feature_name='tabular_intrusion_frame')

    IFrame_builder.set_data(frame_data, INet.intrusion_network_outcome)
    IFrame_builder.setup(
        nelements = 1e2,
        w2=weights[0],
        w1=weights[1],
        w3=weights[2],
    )

    IFrame = IFrame_builder.frame
    
    assert isinstance(IFrame, StructuralFrame)
    
    IBody = IntrusionBody(
        intrusion_data,
        name='tabular_intrusion',
        intrusion_network=INet,
        intrusion_frame=IFrame,
        model=model,
    )


    IBody.set_data_for_s_simulation()
    IBody.set_lateral_extent_conceptual_model(lateral_conceptual_model)
    IBody.set_s_simulation_GSLIBparameters({})
    IBody.make_s_simulation_variogram({})
    IBody.create_grid_for_simulation()
    IBody.simulate_s_thresholds()
    
    assert len(IBody.simulated_s_thresholds) > 0

    IBody.set_data_for_g_simulation()
    IBody.set_vertical_extent_conceptual_model(vertical_conceptual_model)
    IBody.set_g_simulation_GSLIBparameters({})
    IBody.make_g_simulation_variogram({})
    IBody.simulate_g_thresholds()
    
    assert len(IBody.simulated_g_thresholds) > 0
    