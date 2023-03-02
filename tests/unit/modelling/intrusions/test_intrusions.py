import pytest
import numpy as np

# Loop library
from LoopStructural import GeologicalModel
from LoopStructural.modelling.intrusions import IntrusionFrameBuilder
from LoopStructural.modelling.intrusions import IntrusionBuilder
from LoopStructural.modelling.intrusions import IntrusionFeature
from LoopStructural.modelling.features import StructuralFrame
from LoopStructural.modelling.intrusions import (
    ellipse_function,
    constant_function,
)

from LoopStructural.datasets import load_tabular_intrusion

data, boundary_points = load_tabular_intrusion()


# @pytest.mark.skip("not currently working 4-10-2022")
def test_intrusion_frame_builder():
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.data = data
    model.nsteps = [10, 10, 10]

    intrusion_data = data[data["feature_name"] == "tabular_intrusion"]
    intrusion_frame_data = model.data[
        model.data["feature_name"] == "tabular_intrusion_frame"
    ]

    conformable_feature = model.create_and_add_foliation("stratigraphy")

    intrusion_frame_parameters = {
        "contact": "roof",
        "contact_anisotropies": [conformable_feature],
    }

    weights = [0, 0, 0]
    interpolator = model.get_interpolator(interpolatortype="PLI")

    intrusion_frame_builder = IntrusionFrameBuilder(
        interpolator, name="tabular_intrusion_frame", model=model
    )

    intrusion_frame_builder.set_intrusion_frame_parameters(
        intrusion_data, intrusion_frame_parameters
    )

    intrusion_frame_builder.create_constraints_for_c0()

    intrusion_frame_builder.set_intrusion_frame_data(intrusion_frame_data)

    ## -- create intrusion frame
    intrusion_frame_builder.setup(
        w2=weights[0],
        w1=weights[1],
        gxygz=weights[2],
    )

    intrusion_frame = intrusion_frame_builder.frame

    assert isinstance(intrusion_frame, StructuralFrame)


def test_intrusion_builder():
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.data = data
    model.nsteps = [10, 10, 10]

    intrusion_data = data[data["feature_name"] == "tabular_intrusion"]
    intrusion_frame_data = model.data[
        model.data["feature_name"] == "tabular_intrusion_frame"
    ]

    conformable_feature = model.create_and_add_foliation("stratigraphy")

    intrusion_frame_parameters = {
        "contact": "roof",
        "contact_anisotropies": [conformable_feature],
    }

    weights = [0, 0, 0]
    interpolator = model.get_interpolator(interpolatortype="PLI")

    intrusion_frame_builder = IntrusionFrameBuilder(
        interpolator, name="tabular_intrusion_frame", model=model
    )

    intrusion_frame_builder.set_intrusion_frame_parameters(
        intrusion_data, intrusion_frame_parameters
    )

    intrusion_frame_builder.create_constraints_for_c0()

    intrusion_frame_builder.set_intrusion_frame_data(intrusion_frame_data)

    ## -- create intrusion frame
    intrusion_frame_builder.setup(
        w2=weights[0],
        w1=weights[1],
        gxygz=weights[2],
    )

    intrusion_frame = intrusion_frame_builder.frame

    assert isinstance(intrusion_frame, StructuralFrame)

    # -- create intrusion builder to simulate distance thresholds along frame coordinates
    intrusion_builder = IntrusionBuilder(
        intrusion_frame,
        model=model,
        interpolator=interpolator,
        lateral_extent_model=ellipse_function,
        vertical_extent_model=constant_function,
        name="tabular_intrusion",
    )

    intrusion_builder.set_data_for_extent_calculation(intrusion_data)

    intrusion_builder.build_arguments = {
        "geometric_scaling_parameters": {},
    }

    intrusion_feature = intrusion_builder.feature
    intrusion_builder.update()

    assert len(intrusion_builder.data_for_lateral_extent_calculation[0]) > 0
    assert len(intrusion_builder.data_for_lateral_extent_calculation[1]) > 0
    assert len(intrusion_builder.data_for_vertical_extent_calculation[0]) > 0
    assert len(intrusion_builder.data_for_vertical_extent_calculation[1]) > 0


# if __name__ == "__main__":
#     test_intrusion_freame_builder()
#     test_intrusion_builder()
# if __name__ == "__main__":
#     test_intrusion_freame_builder()
#     test_intrusion_builder()
