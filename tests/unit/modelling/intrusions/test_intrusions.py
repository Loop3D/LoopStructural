import numpy as np
import pandas as pd
import pytest

# Loop library
from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_tabular_intrusion
from LoopStructural.modelling.features import StructuralFrame
from LoopStructural.modelling.intrusions import (
    IntrusionBuilder,
    IntrusionFrameBuilder,
    constant_function,
    ellipse_function,
)

data, boundary_points = load_tabular_intrusion()


# @pytest.mark.skip("not currently working 4-10-2022")
def test_intrusion_frame_builder():
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.data = data
    model.nsteps = [10, 10, 10]

    intrusion_data = data[data["feature_name"] == "tabular_intrusion"]
    intrusion_frame_data = model.data[model.data["feature_name"] == "tabular_intrusion_frame"]

    conformable_feature = model.create_and_add_foliation("stratigraphy")

    intrusion_frame_parameters = {
        "contact": "roof",
        "contact_anisotropies": [conformable_feature],
    }

    weights = [0, 0, 0]
    # interpolator = model.get_interpolator(interpolatortype="PLI")

    intrusion_frame_builder = IntrusionFrameBuilder(
        interpolatortype="PLI",
        bounding_box=model.bounding_box,
        nelements=1000,
        name="tabular_intrusion_frame",
        model=model,
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

    intrusion_data = model.prepare_data(data[data["feature_name"] == "tabular_intrusion"])
    intrusion_frame_data = model.prepare_data(model.data[model.data["feature_name"] == "tabular_intrusion_frame"])
    
    conformable_feature = model.create_and_add_foliation("stratigraphy")

    intrusion_frame_parameters = {
        "contact": "roof",
        "contact_anisotropies": [conformable_feature],
    }

    weights = [0, 0, 0]

    intrusion_frame_builder = IntrusionFrameBuilder(
        interpolatortype="PLI",
        bounding_box=model.bounding_box,
        nelements=1000,
        name="tabular_intrusion_frame",
        model=model,
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
        lateral_extent_model=ellipse_function,
        vertical_extent_model=constant_function,
        name="tabular_intrusion",
    )

    intrusion_builder.set_data_for_extent_calculation(intrusion_data)

    intrusion_builder.update_build_arguments(
        {
            "geometric_scaling_parameters": {},
        }
    )

    intrusion_builder.update()

    assert len(intrusion_builder.data_for_lateral_extent_calculation[0]) > 0
    assert len(intrusion_builder.data_for_lateral_extent_calculation[1]) > 0
    assert len(intrusion_builder.data_for_vertical_extent_calculation[0]) > 0
    assert len(intrusion_builder.data_for_vertical_extent_calculation[1]) > 0

    # regression test: up_to_date() must not rebuild the intrusion geometry
    # once it has already been built and nothing has changed (previously
    # IntrusionBuilder never set _up_to_date=True, so this rebuilt on every call)
    call_count = {"n": 0}
    original_prepare_data = intrusion_builder.prepare_data

    def counting_prepare_data(*args, **kwargs):
        call_count["n"] += 1
        return original_prepare_data(*args, **kwargs)

    intrusion_builder.prepare_data = counting_prepare_data

    assert intrusion_builder._up_to_date is True
    intrusion_builder.up_to_date()
    intrusion_builder.up_to_date()
    assert call_count["n"] == 0

    intrusion_builder._up_to_date = False
    intrusion_builder.up_to_date()
    assert call_count["n"] == 1
    assert intrusion_builder._up_to_date is True


def test_intrusion_gyxgz_weight_reaches_frame_build():
    """Regression test for a `gxygz`/`gyxgz` typo in
    `GeologicalModel._build_intrusion` that silently dropped the caller's
    coordinate-2 orthogonality weight: it was passed to
    `IntrusionFrameBuilder.build` under the wrong keyword (`gxygz`), which
    `StructuralFrameBuilder.build` doesn't recognise, so it fell into
    `**kwargs` and `w3` silently stayed at its default. See INTRUSIONS.md
    finding 3.
    """
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.data = data
    model.nsteps = [10, 10, 10]

    conformable_feature = model.create_and_add_foliation("stratigraphy")

    captured_kwargs = {}
    original_build = IntrusionFrameBuilder.build

    def capturing_build(self, *args, **kwargs):
        captured_kwargs.update(kwargs)
        return original_build(self, *args, **kwargs)

    IntrusionFrameBuilder.build = capturing_build
    try:
        model.create_and_add_intrusion(
            "tabular_intrusion",
            "tabular_intrusion_frame",
            intrusion_frame_parameters={
                "contact": "roof",
                "contact_anisotropies": [conformable_feature],
            },
            intrusion_lateral_extent_model=ellipse_function,
            intrusion_vertical_extent_model=constant_function,
            gyxgz=5,
        )
    finally:
        IntrusionFrameBuilder.build = original_build

    assert captured_kwargs.get("gyxgz") == 5


def test_intrusion_geometric_scaling_not_implemented():
    """`geometric_scaling_parameters` is only reached when one of the two
    contacts (roof/floor) has no data. Before this test, that path always
    raised `NotImplementedError` regardless of what was passed in, several
    calls deep inside a helper that looked like it partially worked (see
    INTRUSIONS.md finding 2). This pins the simplified, immediate failure
    so the behaviour stays an explicit "not supported" rather than silently
    regressing to something that looks like it works.
    """
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    single_contact_data = data[
        ~(
            (data["feature_name"] == "tabular_intrusion")
            & (data["intrusion_contact_type"] == "floor")
        )
    ].copy()
    # `tabular_intrusion.csv` doesn't have nx/ny/nz/tx/ty/tz columns, so it
    # must be routed through `prepare_data` before evaluating any feature
    # end-to-end (`model.data = ...` alone does not normalise columns) -- see
    # INTRUSIONS.md finding 1b.
    model.data = model.prepare_data(single_contact_data)
    model.nsteps = [10, 10, 10]

    conformable_feature = model.create_and_add_foliation("stratigraphy")

    intrusion_feature = model.create_and_add_intrusion(
        "tabular_intrusion",
        "tabular_intrusion_frame",
        intrusion_frame_parameters={
            "contact": "roof",
            "contact_anisotropies": [conformable_feature],
        },
        intrusion_lateral_extent_model=ellipse_function,
        intrusion_vertical_extent_model=constant_function,
        geometric_scaling_parameters={"thickness": 5.0},
    )

    with pytest.raises(
        NotImplementedError, match="geometric_scaling_parameters is not currently supported"
    ):
        intrusion_feature.evaluate_value(np.array([[2.5, 2.5, 1.5]]))


def test_intrusion_missing_data_raises_clear_error():
    """Before this test, an intrusion_name/intrusion_frame_name typo (or a
    feature_name with no matching rows) would fail deep inside
    `IntrusionFrameBuilder`/`IntrusionBuilder` with a bare `KeyError`, not at
    the `create_and_add_intrusion` boundary. See INTRUSIONS.md finding 4.
    """
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.data = model.prepare_data(data)
    model.nsteps = [10, 10, 10]

    conformable_feature = model.create_and_add_foliation("stratigraphy")

    with pytest.raises(ValueError, match="No data found for intrusion 'not_a_real_feature'"):
        model.create_and_add_intrusion(
            "not_a_real_feature",
            "tabular_intrusion_frame",
            intrusion_frame_parameters={
                "contact": "roof",
                "contact_anisotropies": [conformable_feature],
            },
            intrusion_lateral_extent_model=ellipse_function,
            intrusion_vertical_extent_model=constant_function,
        )

    with pytest.raises(ValueError, match="No data found for intrusion frame 'not_a_real_frame'"):
        model.create_and_add_intrusion(
            "tabular_intrusion",
            "not_a_real_frame",
            intrusion_frame_parameters={
                "contact": "roof",
                "contact_anisotropies": [conformable_feature],
            },
            intrusion_lateral_extent_model=ellipse_function,
            intrusion_vertical_extent_model=constant_function,
        )


def test_intrusion_missing_contact_type_column_raises_clear_error():
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    incomplete_data = model.prepare_data(data).drop(columns=["intrusion_contact_type"])
    model.data = incomplete_data
    model.nsteps = [10, 10, 10]

    conformable_feature = model.create_and_add_foliation("stratigraphy")

    with pytest.raises(ValueError, match="missing required column"):
        model.create_and_add_intrusion(
            "tabular_intrusion",
            "tabular_intrusion_frame",
            intrusion_frame_parameters={
                "contact": "roof",
                "contact_anisotropies": [conformable_feature],
            },
            intrusion_lateral_extent_model=ellipse_function,
            intrusion_vertical_extent_model=constant_function,
        )


def test_intrusion_missing_contact_anisotropies_raises_clear_error():
    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.data = model.prepare_data(data)
    model.nsteps = [10, 10, 10]

    with pytest.raises(ValueError, match="contact_anisotropies"):
        model.create_and_add_intrusion(
            "tabular_intrusion",
            "tabular_intrusion_frame",
            intrusion_frame_parameters={"contact": "roof"},
            intrusion_lateral_extent_model=ellipse_function,
            intrusion_vertical_extent_model=constant_function,
        )


def _build_marginal_fault_model():
    """A minimal sill offset by a single bounding fault, to exercise the
    `marginal_faults` code path in `IntrusionFrameBuilder`
    (`set_marginal_faults_parameters`/`create_constraints_for_c0`), which
    -- unlike the plain tabular-intrusion path -- has no test or example
    anywhere in the codebase (see INTRUSIONS.md finding 1). Reuses the
    `stratigraphy` feature and the `tabular_intrusion`/`tabular_intrusion_frame`
    geometry from `load_tabular_intrusion()` (proven to build end-to-end),
    just renamed, plus one fault.
    """
    stratigraphy_rows = data[data["feature_name"] == "stratigraphy"].copy()

    fault_rows = pd.DataFrame(
        [
            [2.5, 2.5, 2.5, 0, 1, 0, 0, "marginal_fault", 0],
            [2.5, 2.5, 2.5, 1, 0, 0, 1, "marginal_fault", 0],
            [2.5, 2.5, 2.5, 0, 0, 1, 2, "marginal_fault", 0],
        ],
        columns=["X", "Y", "Z", "nx", "ny", "nz", "coord", "feature_name", "val"],
    )

    sill_frame_rows = pd.DataFrame(
        [
            [2.00, 2.00, 2.00, 0, np.nan, 0, 0, -1, "sill_frame"],
            [3.00, 1.00, 2.00, 0, np.nan, 0, 0, -1, "sill_frame"],
            [1.00, 3.00, 2.00, 0, np.nan, 0, 0, -1, "sill_frame"],
            [3.00, 2.00, 1.00, 1, 0, np.nan, np.nan, np.nan, "sill_frame"],
            [3.00, 2.00, 1.00, 1, np.nan, 0, 1, 0, "sill_frame"],
            [2.50, 1.00, 1.00, 2, 0, np.nan, np.nan, np.nan, "sill_frame"],
            [2.50, 2.00, 1.00, 2, 0, np.nan, np.nan, np.nan, "sill_frame"],
            [2.50, 2.00, 1.00, 2, np.nan, 1, 0, 0, "sill_frame"],
        ],
        columns=["X", "Y", "Z", "coord", "val", "gx", "gy", "gz", "feature_name"],
    )

    # Roof/floor contact points spanning both sides of the Y=2.5 fault, and
    # a handful flagged for the lateral (side) extent -- same shape as
    # `tabular_intrusion`'s own contact data, just relabelled.
    sill_contact_rows = pd.DataFrame(
        [
            [3.04, 2.12, 2.18, "roof", False],
            [4.02, 3.85, 2.14, "roof", True],
            [2.02, 3.16, 2.02, "roof", False],
            [2.12, 2.16, 2.02, "roof", False],
            [1.14, 3.50, 2.04, "roof", True],
            [4.08, 3.02, 1.18, "floor", True],
            [1.14, 1.06, 1.16, "floor", True],
            [4.20, 2.18, 1.12, "floor", True],
            [1.02, 3.02, 1.12, "floor", True],
        ],
        columns=["X", "Y", "Z", "intrusion_contact_type", "intrusion_side"],
    )
    sill_contact_rows["feature_name"] = "sill"

    model_data = pd.concat(
        [stratigraphy_rows, fault_rows, sill_frame_rows, sill_contact_rows],
        ignore_index=True,
    )

    model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
    model.data = model.prepare_data(model_data)
    model.nsteps = [10, 10, 10]

    stratigraphy = model.create_and_add_foliation("stratigraphy")
    fault = model.create_and_add_fault("marginal_fault", 0.5, nelements=1e4)

    return model, stratigraphy, fault


def test_intrusion_marginal_faults():
    model, stratigraphy, fault = _build_marginal_fault_model()

    intrusion_feature = model.create_and_add_intrusion(
        "sill",
        "sill_frame",
        intrusion_frame_parameters={
            "contact": "roof",
            "contact_anisotropies": [stratigraphy],
            "marginal_faults": {
                0: {
                    "structure": fault,
                    "block": "hanging wall",
                    "series": stratigraphy,
                }
            },
        },
        intrusion_lateral_extent_model=ellipse_function,
        intrusion_vertical_extent_model=constant_function,
    )

    value = intrusion_feature.evaluate_value(np.array([[2.5, 2.5, 1.5]]))
    assert np.isfinite(value).all()


def test_intrusion_steps_broken_with_current_stratigraphic_column():
    """`intrusion_steps` is currently unusable, not just untested.

    `IntrusionFrameBuilder.set_intrusion_steps_parameters`
    (`intrusion_frame_builder.py`) reads
    `self.model.stratigraphic_column[series_from_name.name][unit_from_name]`,
    i.e. it expects `model.stratigraphic_column` to be the old nested
    `{group_name: {unit_name: {...}}}` dict. `GeologicalModel.stratigraphic_column`
    is now a `StratigraphicColumn` object whose `__getitem__` looks up a
    single element by `uuid` (not by group/series name), and
    `GeologicalModel.set_stratigraphic_column` (the old dict-based setter)
    unconditionally raises `DeprecationWarning` -- so there is no longer any
    way to put the model into the shape `intrusion_steps` expects. This is a
    real regression from an unrelated stratigraphic-column refactor that was
    never propagated to the intrusions module (see INTRUSIONS.md finding 1,
    revised). This test pins the current failure so a future fix (Stage 6)
    has a clear target and so this doesn't regress further/silently.
    """
    model, stratigraphy, fault = _build_marginal_fault_model()
    model.stratigraphic_column.add_unit("unitA", thickness=1)
    model.stratigraphic_column.add_unit("unitB", thickness=1)

    with pytest.raises(KeyError, match="No element found with uuid"):
        model.create_and_add_intrusion(
            "sill",
            "sill_frame",
            intrusion_frame_parameters={
                "contact": "roof",
                "contact_anisotropies": [stratigraphy],
                "intrusion_steps": {
                    0: {
                        "structure": fault,
                        "unit_from": "unitA",
                        "series_from": stratigraphy,
                        "unit_to": "unitB",
                        "series_to": stratigraphy,
                    }
                },
            },
            intrusion_lateral_extent_model=ellipse_function,
            intrusion_vertical_extent_model=constant_function,
        )


# if __name__ == "__main__":
#     test_intrusion_freame_builder()
#     test_intrusion_builder()
# if __name__ == "__main__":
#     test_intrusion_freame_builder()
#     test_intrusion_builder()
