import numpy as np

# Loop library
from LoopStructural import GeologicalModel
from LoopStructural.modelling.intrusions import IntrusionFrameBuilder
from LoopStructural.modelling.intrusions import IntrusionBuilder
from LoopStructural.modelling.intrusions import IntrusionFeature
from LoopStructural.modelling.features import StructuralFrame
from LoopStructural.modelling.intrusions import rectangle_function, parallelepiped_function

from LoopStructural.datasets import load_tabular_intrusion
data, boundary_points = load_tabular_intrusion()

def test_intrusion_freame_builder():
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.data = data
    model.nsteps = [10,10,10]

    intrusion_data = data[data['feature_name'] == 'tabular_intrusion']
    frame_data = model.data[model.data["feature_name"] == 'tabular_intrusion_frame'].copy()

    conformable_feature = model.create_and_add_foliation('stratigraphy')

    intrusion_network_parameters = {'type': 'shortest path',
                                    'contact' :'roof',
                                    'delta_c':[2],
                                    'contact_anisotropies' : [conformable_feature],
                                    'shortest_path_sequence':[conformable_feature],
                                    'shortest_path_axis':'X'}
    delta_c = intrusion_network_parameters.get('delta_c')[0]
    # -- get variables for intrusion frame interpolation
    interpolatortype = "FDI"
    nelements = 1e2
    weights = [0,0,0]
    interpolator = model.get_interpolator(interpolatortype=interpolatortype)

    intrusion_frame_builder = IntrusionFrameBuilder(
        interpolator, name='tabular_intrusion_frame', model=model)

    # -- create intrusion network

    intrusion_frame_builder.set_intrusion_network_parameters(intrusion_data, intrusion_network_parameters)
    intrusion_network_geometry = intrusion_frame_builder.create_intrusion_network()

    keys = list(intrusion_frame_builder.anisotropies_series_parameters.keys())
    # #test if points lie in the contact of interest
    mean = intrusion_frame_builder.anisotropies_series_parameters[keys[0]][1]
    # mean = -10
    stdv = intrusion_frame_builder.anisotropies_series_parameters[keys[0]][2]
    evaluated_inet_points =  conformable_feature['feature'].evaluate_value(model.scale(intrusion_network_geometry[:,:3]))

    assert np.all(np.logical_and((mean - stdv*delta_c)<= evaluated_inet_points,(mean + stdv*delta_c)>= evaluated_inet_points))

    # -- create intrusion frame using intrusion network points and flow/inflation measurements
    intrusion_frame_builder.set_intrusion_frame_data(frame_data, intrusion_network_geometry)

    ## -- create intrusion frame
    intrusion_frame_builder.setup(
        nelements=nelements,
        w2=weights[0],
        w1=weights[1],
        gxygz=weights[2],
    )

    intrusion_frame = intrusion_frame_builder.frame

    assert isinstance(intrusion_frame, StructuralFrame)


def test_intrusion_builder():

    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.data = data
    model.nsteps = [10,10,10]

    intrusion_data = data[data['feature_name'] == 'tabular_intrusion']
    frame_data = model.data[model.data["feature_name"] == 'tabular_intrusion_frame'].copy()

    conformable_feature = model.create_and_add_foliation('stratigraphy')

    intrusion_network_parameters = {'type': 'interpolated',
                                    'contact' :'roof'}
    interpolator = model.get_interpolator(interpolatortype='FDI')

    intrusion_frame_builder = IntrusionFrameBuilder(
        interpolator, name='tabular_intrusion_frame', model=model)

    # -- create intrusion network

    intrusion_frame_builder.set_intrusion_network_parameters(intrusion_data, intrusion_network_parameters)
    intrusion_network_geometry = intrusion_frame_builder.create_intrusion_network()

    # -- create intrusion frame using intrusion network points and flow/inflation measurements
    intrusion_frame_builder.set_intrusion_frame_data(frame_data, intrusion_network_geometry)

    ## -- create intrusion frame
    intrusion_frame_builder.setup()
    intrusion_frame = intrusion_frame_builder.frame

    # -- create intrusion builder to simulate distance thresholds along frame coordinates
    intrusion_builder = IntrusionBuilder(
        intrusion_frame, model=model, name="tabular intrusion"
    )
    intrusion_builder.lateral_extent_model =  rectangle_function #intrusion_lateral_extent_model
    intrusion_builder.vertical_extent_model = parallelepiped_function #intrusion_vertical_extent_model

    intrusion_builder.set_data_for_extent_simulation(intrusion_data)
    intrusion_builder.build_arguments = {
        "lateral_extent_sgs_parameters": {},
        "vertical_extent_sgs_parameters": {}}

    intrusion_feature = intrusion_builder.feature
    intrusion_builder.update()

    assert len(intrusion_feature._lateral_simulated_thresholds) > 0
    assert len(intrusion_feature._growth_simulated_thresholds) > 0
    