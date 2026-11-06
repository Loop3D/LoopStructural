"""Backend matrix: cells, optional dependencies and pure checker functions (FR-030, SC-018).

A cell is a group of tests for one backend or capability. A mandatory cell must pass.
An optional cell needs a named dependency. It is "excluded" only when the dependency is
absent and the cell was skipped.

Run ``python -m tests.fixtures.backend_matrix <junit>...`` from the repository root. The
command prints the table and exits with 1 on a violation.
"""

from __future__ import annotations

import importlib.util
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

try:  # Literal exists in typing from Python 3.8.
    from typing import Literal
except ImportError:  # pragma: no cover
    from typing_extensions import Literal  # type: ignore

Outcome = Literal["passed", "failed", "skipped", "absent"]

# The 7 rows of "Required backend coverage" in spec.md. The value "extra" is not a row.
SPEC_ROWS: Tuple[str, ...] = (
    "Python FDI and PLI/P1",
    "Python P2",
    "Python FDI-fold, PLI discrete-fold, and P2-fold",
    "Python FDI and PLI constant-norm variants",
    "Common 2D and 3D geometry",
    "Existing native FDI or fold integration",
    "Optional Surfe integration",
)

_SEVERITY = {"passed": 0, "skipped": 1, "failed": 2}


@dataclass(frozen=True)
class OptionalDependency:
    """One optional dependency.

    Parameters
    ----------
    key : str
        Key in ``OPTIONAL_DEPENDENCIES``.
    module : str
        Exact import name that ``find_spec`` looks for.
    distribution : str
        Name to install (pip name or extra).
    capability : str
        What the dependency gives.
    """

    key: str
    module: str
    distribution: str
    capability: str


OPTIONAL_DEPENDENCIES: Dict[str, OptionalDependency] = {
    d.key: d
    for d in (
        OptionalDependency(
            "native",
            "loop_interpolation_core",
            "loop-interpolation[native]",
            "native FDI and fold extension",
        ),
        # The pip name of Surfe is unverified (see T092).
        OptionalDependency("surfe", "surfepy", "surfepy", "Surfe RBF extension"),
        OptionalDependency(
            "viewer", "loopstructuralvisualisation", "loopstructuralvisualisation", "2D viewer"
        ),
        OptionalDependency("matplotlib", "matplotlib", "matplotlib", "plot output"),
        OptionalDependency("pyevtk", "pyevtk", "pyevtk", "VTK export"),
    )
}


@dataclass(frozen=True)
class OptionalStatus:
    """Status of one optional dependency.

    Parameters
    ----------
    key : str
        Key in ``OPTIONAL_DEPENDENCIES``.
    module : str
        Exact module name.
    available : bool
        True when the module can be imported.
    reason : str
        Empty when available. Else the text that names the module and the distribution.
    """

    key: str
    module: str
    available: bool
    reason: str


@dataclass(frozen=True)
class Cell:
    """One cell of the backend matrix.

    Parameters
    ----------
    cell_id : str
        Unique id of the cell.
    spec_row : str
        A value of ``SPEC_ROWS`` or ``"extra"``.
    mandatory : bool
        True when the cell must pass with no skip.
    stub : bool
        True when the cell uses a recording stub. It is not a result for a spec row.
    refs : tuple of (str, str, str or None)
        Test refs: path with "/", test function, id selector or None.
    min_cases : int
        Least number of passed cases for the outcome "passed".
    dependency : str or tuple of str or None
        Key or keys of ``OPTIONAL_DEPENDENCIES``. The cell is "dependency absent" when
        ANY key is not available.
    """

    cell_id: str
    spec_row: str
    mandatory: bool
    stub: bool
    refs: Tuple[Tuple[str, str, Optional[str]], ...]
    min_cases: int
    dependency: Union[str, Tuple[str, ...], None]


# Refs, selectors and min_cases come from `pytest --collect-only -q` (2026-11-06).
# A selector is a literal in the collected ids. min_cases is the collected number of cases.
MATRIX: Tuple[Cell, ...] = (
    Cell(
        "py-fdi",
        "Python FDI and PLI/P1",
        True,
        False,
        (
            (
                "tests/unit/modelling/test_frame_persistence.py",
                "test_round_trip_world_values",
                "[fdi-",
            ),
            (
                "tests/unit/modelling/test_frame_persistence.py",
                "test_payload_declares_frame",
                "[fdi-",
            ),
            (
                "tests/unit/modelling/test_frame_copy_and_lifecycle.py",
                "test_copy_keeps_frame",
                "[fdi-",
            ),
            (
                "tests/unit/modelling/test_frame_copy_and_lifecycle.py",
                "test_stale_frame_rejected_without_state_change",
                "[fdi]",
            ),
            (
                "tests/unit/modelling/test_frame_copy_and_lifecycle.py",
                "test_rebuild_in_new_frame_preserves_world_evaluation",
                "[fdi]",
            ),
            (
                "tests/unit/modelling/test_fault_frame_analytic.py",
                "test_planar_fault_frame_matches_analytic",
                "FDI",
            ),
            (
                "tests/unit/modelling/test_fault_frame_validation.py",
                "test_scalar_displacement_matches_analytic_frame",
                "FDI",
            ),
            (
                "tests/unit/modelling/test_fault_frame_validation.py",
                "test_vector_displacement_matches_analytic_frame",
                "FDI",
            ),
            (
                "tests/unit/modelling/test_reference_constraint_frames.py",
                "test_orthogonal_rows_affine_residual",
                "FDI",
            ),
            (
                "tests/unit/modelling/test_reference_constraint_frames.py",
                "test_world_local_conversion_methods",
                "FDI",
            ),
            (
                "tests/unit/modelling/test_fault_buffer_bounds.py",
                "test_buffer_bounds_match_request",
                "FDI",
            ),
        ),
        152,
        None,
    ),
    Cell(
        "py-pli",
        "Python FDI and PLI/P1",
        True,
        False,
        (
            (
                "tests/unit/modelling/test_frame_persistence.py",
                "test_round_trip_world_values",
                "[pli-",
            ),
            (
                "tests/unit/modelling/test_frame_persistence.py",
                "test_payload_declares_frame",
                "[pli-",
            ),
            (
                "tests/unit/modelling/test_frame_copy_and_lifecycle.py",
                "test_stale_frame_rejected_without_state_change",
                "[pli]",
            ),
            (
                "tests/unit/modelling/test_frame_copy_and_lifecycle.py",
                "test_rebuild_in_new_frame_preserves_world_evaluation",
                "[pli]",
            ),
            (
                "tests/unit/modelling/test_fault_frame_analytic.py",
                "test_planar_fault_frame_matches_analytic",
                "PLI",
            ),
            (
                "tests/unit/modelling/test_fault_frame_validation.py",
                "test_scalar_displacement_matches_analytic_frame",
                "PLI",
            ),
            (
                "tests/unit/modelling/test_fault_frame_validation.py",
                "test_vector_displacement_matches_analytic_frame",
                "PLI",
            ),
            (
                "tests/unit/modelling/test_reference_constraint_frames.py",
                "test_orthogonal_rows_affine_residual",
                "PLI",
            ),
            (
                "tests/unit/modelling/test_reference_constraint_frames.py",
                "test_world_local_conversion_methods",
                "PLI",
            ),
            (
                "tests/unit/modelling/test_fault_buffer_bounds.py",
                "test_buffer_bounds_match_request",
                "PLI",
            ),
        ),
        140,
        None,
    ),
    Cell(
        "py-p2",
        "Python P2",
        True,
        False,
        (
            (
                "tests/unit/modelling/test_frame_persistence.py",
                "test_round_trip_world_values",
                "[p2-",
            ),
            (
                "tests/unit/modelling/test_frame_persistence.py",
                "test_payload_declares_frame",
                "[p2-",
            ),
            (
                "tests/unit/modelling/test_frame_copy_and_lifecycle.py",
                "test_copy_keeps_frame",
                "[p2-",
            ),
            (
                "tests/unit/modelling/test_frame_copy_and_lifecycle.py",
                "test_stale_frame_rejected_without_state_change",
                "[p2]",
            ),
            (
                "tests/unit/modelling/test_frame_copy_and_lifecycle.py",
                "test_rebuild_in_new_frame_preserves_world_evaluation",
                "[p2]",
            ),
            (
                "tests/unit/modelling/test_reference_constraint_frames.py",
                "test_orthogonal_rows_affine_residual",
                "P2",
            ),
            (
                "tests/unit/modelling/test_reference_constraint_frames.py",
                "test_world_local_conversion_methods",
                "P2",
            ),
            (
                "tests/unit/modelling/test_diagnostics_backends.py",
                "test_p2_misfit_leaves_state_unchanged",
                None,
            ),
            (
                "tests/unit/modelling/test_diagnostics_backends.py",
                "test_p2_report_blocks_are_separate",
                None,
            ),
        ),
        110,
        None,
    ),
    Cell(
        "py-fold-fdi",
        "Python FDI-fold, PLI discrete-fold, and P2-fold",
        True,
        False,
        (
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_rows_affine_residual",
                "fdi_fold",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_geometry_sampled_at_world_points",
                "fdi_fold",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_mask_receives_world_points",
                "fdi_fold",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_norm_alignment_uses_world_normals",
                "fdi_fold",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_directional_regularisation_world_contract",
                "-fdi-",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_rows_follow_curved_directions_fdi",
                "fdi_fold",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_regularisation_world_contract",
                "fdi_fold",
            ),
        ),
        108,
        None,
    ),
    Cell(
        "py-fold-pli",
        "Python FDI-fold, PLI discrete-fold, and P2-fold",
        True,
        False,
        (
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_rows_affine_residual",
                "pli_fold",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_geometry_sampled_at_world_points",
                "pli_fold",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_mask_receives_world_points",
                "pli_fold",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_norm_alignment_uses_world_normals",
                "pli_fold",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_directional_regularisation_world_contract",
                "-pli-",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_rows_follow_curved_directions_pli_p2",
                "pli_fold",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_regularisation_world_contract",
                "pli_fold",
            ),
        ),
        108,
        None,
    ),
    Cell(
        "py-fold-p2",
        "Python FDI-fold, PLI discrete-fold, and P2-fold",
        True,
        False,
        (
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_rows_affine_residual",
                "p2_fold",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_geometry_sampled_at_world_points",
                "p2_fold",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_mask_receives_world_points",
                "p2_fold",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_norm_alignment_uses_world_normals",
                "p2_fold",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_directional_regularisation_world_contract",
                "-p2-",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_fold_rows_follow_curved_directions_pli_p2",
                "p2_fold",
            ),
            (
                "tests/unit/modelling/test_fold_constraint_frames.py",
                "test_p2_steepness_wtfunc_receives_world_barycentres",
                None,
            ),
        ),
        156,
        None,
    ),
    Cell(
        "py-cn-fdi",
        "Python FDI and PLI constant-norm variants",
        True,
        False,
        (
            (
                "tests/unit/modelling/test_constant_norm_frames.py",
                "test_constant_norm_fit_recovers_r1_in_rotated_frame",
                "cn_fdi",
            ),
            (
                "tests/unit/modelling/test_frame_copy_and_lifecycle.py",
                "test_stale_frame_rejected_without_state_change",
                "[cn_fdi]",
            ),
            (
                "tests/unit/modelling/test_frame_copy_and_lifecycle.py",
                "test_rebuild_in_new_frame_preserves_world_evaluation",
                "[cn_fdi]",
            ),
        ),
        14,
        None,
    ),
    Cell(
        "py-cn-pli",
        "Python FDI and PLI constant-norm variants",
        True,
        False,
        (
            (
                "tests/unit/modelling/test_constant_norm_frames.py",
                "test_constant_norm_fit_recovers_r1_in_rotated_frame",
                "cn_pli",
            ),
            (
                "tests/unit/modelling/test_frame_copy_and_lifecycle.py",
                "test_stale_frame_rejected_without_state_change",
                "[cn_pli]",
            ),
            (
                "tests/unit/modelling/test_frame_copy_and_lifecycle.py",
                "test_rebuild_in_new_frame_preserves_world_evaluation",
                "[cn_pli]",
            ),
        ),
        14,
        None,
    ),
    Cell(
        "common-geometry",
        "Common 2D and 3D geometry",
        True,
        False,
        (
            (
                "packages/loop_common/tests/test_rigid_transform_contract.py",
                "test_invalid_rotation_keeps_prior_state",
                None,
            ),
            (
                "packages/loop_common/tests/test_rigid_transform_contract.py",
                "test_local_enclosure_contains_projected_corners",
                None,
            ),
            (
                "packages/loop_common/tests/test_rigid_transform_contract.py",
                "test_point_and_direction_round_trips",
                None,
            ),
            (
                "packages/loop_common/tests/test_rigid_transform_contract.py",
                "test_reapply_unchanged_bounds_is_stable",
                None,
            ),
            (
                "packages/loop_common/tests/test_rigid_transform_contract.py",
                "test_validate_rigid_rotation_accepts_proper_rotations",
                None,
            ),
            (
                "packages/loop_common/tests/test_rigid_transform_contract.py",
                "test_validate_rigid_rotation_rejects_invalid",
                None,
            ),
            ("packages/loop_common/tests/test_resize_paths.py", "test_2d_counts_are_derived", None),
            (
                "packages/loop_common/tests/test_resize_paths.py",
                "test_nsteps_below_two_is_rejected_without_change",
                None,
            ),
            (
                "packages/loop_common/tests/test_resize_paths.py",
                "test_nsteps_not_whole_is_rejected_without_change",
                None,
            ),
            (
                "packages/loop_common/tests/test_resize_paths.py",
                "test_nsteps_setter_keeps_extent",
                None,
            ),
            (
                "packages/loop_common/tests/test_resize_paths.py",
                "test_reapplying_same_value_10_times_changes_nothing_3d",
                None,
            ),
            (
                "packages/loop_common/tests/test_resize_paths.py",
                "test_rectilinear_resize_is_rejected",
                None,
            ),
            (
                "packages/loop_common/tests/test_resize_paths.py",
                "test_set_nelements_keeps_extent",
                None,
            ),
            (
                "packages/loop_common/tests/test_resize_paths.py",
                "test_set_nelements_small_request_is_clamped",
                None,
            ),
            (
                "packages/loop_common/tests/test_resize_paths.py",
                "test_step_vector_a_b_a_extent_drift",
                None,
            ),
            (
                "packages/loop_common/tests/test_resize_paths.py",
                "test_step_vector_setter_encloses_extent",
                None,
            ),
            (
                "packages/loop_common/tests/test_cell_centre_sampling.py",
                "test_centres_are_cell_midpoints_in_world",
                None,
            ),
            (
                "packages/loop_common/tests/test_cell_centre_sampling.py",
                "test_centres_inside_box",
                None,
            ),
            (
                "packages/loop_common/tests/test_cell_centre_sampling.py",
                "test_centres_match_exported_cells",
                None,
            ),
            (
                "packages/loop_common/tests/test_cell_centre_sampling.py",
                "test_centres_round_trip_and_index_own_cells",
                None,
            ),
            (
                "packages/loop_common/tests/test_cell_centre_sampling.py",
                "test_local_coordinates_grid_associates_own_centres",
                None,
            ),
            (
                "packages/loop_common/tests/test_euclidean_transformation.py",
                "test_export_from_loop_common_math",
                None,
            ),
            (
                "packages/loop_common/tests/test_euclidean_transformation.py",
                "test_fit_without_rotation_needs_no_scikit_learn",
                None,
            ),
            (
                "packages/loop_common/tests/test_euclidean_transformation.py",
                "test_fit_without_scikit_learn_raises_import_error",
                None,
            ),
            (
                "packages/loop_common/tests/test_euclidean_transformation.py",
                "test_module_imports",
                None,
            ),
            (
                "packages/loop_common/tests/test_euclidean_transformation.py",
                "test_round_trip_without_fit",
                None,
            ),
            (
                "packages/loop_common/tests/test_euclidean_transformation.py",
                "test_three_dimensional_rotation_keeps_the_z_direction",
                None,
            ),
            (
                "packages/loop_common/tests/test_euclidean_transformation.py",
                "test_two_dimensional_rotation_of_known_angle",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_call_is_equivalent_to_transform",
                None,
            ),
            ("tests/unit/utils/test_transformation.py", "test_default_construction", None),
            (
                "tests/unit/utils/test_transformation.py",
                "test_entry_points_are_one_class",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_fit_aligns_a_line_with_negative_slope",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_fit_finds_rotation_that_aligns_main_axis_with_x",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_fit_raises_if_points_have_too_few_columns",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_fit_rotation_false_keeps_angle_zero",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_fit_three_dimensions",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_integer_input_is_not_truncated",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_new_pickle_names_the_loop_common_path",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_one_dimensional_input_raises_value_error",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_pickle_with_the_old_module_path_loads",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_repr_html_contains_translation_and_angle",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_rotation_and_inverse_rotation_are_transposes_in_plane",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_round_trip_after_fit",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_round_trip_explicit_parameters",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_too_few_columns_raises_value_error",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_transform_aligns_variance_with_x_axis",
                None,
            ),
            (
                "tests/unit/utils/test_transformation.py",
                "test_transform_raises_if_points_have_too_few_columns",
                None,
            ),
        ),
        296,
        None,
    ),
    Cell(
        "native-boundary",
        "Existing native FDI or fold integration",
        True,
        True,
        (
            (
                "tests/unit/interpolator/test_backend_world_contract.py",
                "test_native_adapter_sends_local_rows_and_grid",
                None,
            ),
            (
                "tests/unit/interpolator/test_backend_world_contract.py",
                "test_native_solved_field_evaluates_in_world",
                None,
            ),
            (
                "tests/unit/interpolator/test_backend_world_contract.py",
                "test_native_unsupported_settings_keep_rejection_and_state",
                None,
            ),
        ),
        27,
        None,
    ),
    Cell(
        "surfe-boundary",
        "Optional Surfe integration",
        True,
        True,
        (
            (
                "tests/unit/interpolator/test_backend_world_contract.py",
                "test_surfe_adapter_returns_world_values_and_gradients",
                None,
            ),
            (
                "tests/unit/interpolator/test_backend_world_contract.py",
                "test_surfe_adapter_sends_local_points",
                None,
            ),
        ),
        24,
        None,
    ),
    Cell(
        "native-extension",
        "Existing native FDI or fold integration",
        False,
        False,
        (
            (
                "tests/unit/interpolator/test_backend_optional_native_surfe.py",
                "test_native_extension_r1_fit_in_world",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_native_parity.py",
                "test_native_matches_cg_for_constraint_family",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_native_inequality.py",
                "test_native_matches_admm_with_value_bounds",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_native_orthogonal_and_limits.py",
                "test_orthogonal_constraints_give_the_same_field_as_cg",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_native_warm_start_ridge.py",
                "test_native_warm_start_needs_fewer_iterations",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_native_no_rows.py",
                "test_native_before_setup_builds_no_rows",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_native_errors.py",
                "test_native_accepts_normal_data_rows_after_native_setup",
                None,
            ),
            (
                "tests/unit/modelling/test_native_fault_frame.py",
                "test_all_coordinates_of_a_fault_use_the_native_solver",
                None,
            ),
            (
                "tests/unit/modelling/test_builder_solver.py",
                "test_native_solver_set_in_build_arguments_builds_no_rows",
                None,
            ),
        ),
        15,
        "native",
    ),
    Cell(
        "surfe-extension",
        "Optional Surfe integration",
        False,
        False,
        (
            (
                "tests/unit/interpolator/test_backend_optional_native_surfe.py",
                "test_surfe_extension_r1_values_and_gradients_in_world",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_evaluate_gradient_empty_input",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_evaluate_gradient_returns_array",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_evaluate_gradient_shape",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_evaluate_gradient_single_point",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_evaluate_gradient_with_nan_coordinates",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_evaluate_value_all_nan_input",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_evaluate_value_empty_input",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_evaluate_value_single_point",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_evaluate_value_with_nan_coordinates",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_evaluate_value_with_valid_points",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_gradient_and_value_handle_nans_same_way",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_interpolator_is_rbf",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_interpolator_type",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_mixed_valid_and_invalid_coordinates",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_nan_gradient_masking_uses_isnan",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_nan_masking_uses_isnan_not_comparison",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_partial_nan_rows",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_surfe_rbf_creation",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_surfe_rbf_has_surfe_library",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_surfe_rbf_interpolator.py",
                "test_very_large_coordinate_values",
                None,
            ),
            (
                "packages/loop_interpolation/tests/test_p0_surfe_nans.py",
                "test_surfe_evaluate_value_masks_nan_coordinates",
                None,
            ),
        ),
        32,
        "surfe",
    ),
    Cell(
        "display-viewer",
        "extra",
        False,
        False,
        (
            (
                "tests/unit/modelling/test_world_location_consumers.py",
                "test_2d_viewer_fault_ellipse_centre_is_world",
                None,
            ),
        ),
        4,
        ("viewer", "matplotlib"),
    ),
    Cell(
        "export-pyevtk",
        "extra",
        False,
        False,
        (
            (
                "tests/unit/modelling/test_world_location_consumers.py",
                "test_internal_callers_do_not_warn",
                None,
            ),
        ),
        1,
        "pyevtk",
    ),
)

StatusArg = Union[OptionalStatus, Sequence[OptionalStatus], None]


def optional_status(
    key: str, find_spec: Callable[[str], object] = importlib.util.find_spec
) -> OptionalStatus:
    """Find out if an optional dependency is installed.

    Parameters
    ----------
    key : str
        Key in ``OPTIONAL_DEPENDENCIES``.
    find_spec : callable
        Function with the signature of ``importlib.util.find_spec``.

    Returns
    -------
    OptionalStatus
        The status. ``reason`` is empty when the module is available.
    """
    dep = OPTIONAL_DEPENDENCIES[key]
    try:
        found = find_spec(dep.module)
    except (ImportError, ValueError):
        found = None
    if found is not None:
        return OptionalStatus(key, dep.module, True, "")
    reason = "optional dependency '%s' is not installed (%s)" % (dep.module, dep.distribution)
    return OptionalStatus(key, dep.module, False, reason)


def require_optional(key: str) -> None:
    """Skip the running test when an optional dependency is not installed.

    Parameters
    ----------
    key : str
        Key in ``OPTIONAL_DEPENDENCIES``.
    """
    import pytest

    status = optional_status(key)
    if not status.available:
        pytest.skip(status.reason)


def cell_dependency_keys(cell: Cell) -> Tuple[str, ...]:
    """Return the dependency keys of a cell as a tuple.

    Parameters
    ----------
    cell : Cell
        The cell.

    Returns
    -------
    tuple of str
        Empty when the cell has no dependency.
    """
    dep = cell.dependency
    if dep is None:
        return ()
    if isinstance(dep, str):
        return (dep,)
    return tuple(dep)


def classify_cell(cell: Cell, outcome: str, status: StatusArg) -> str:
    """Classify the outcome of a cell.

    Parameters
    ----------
    cell : Cell
        The cell.
    outcome : str
        One of "passed", "failed", "skipped", "absent".
    status : OptionalStatus or sequence of OptionalStatus or None
        Status of the dependencies of the cell.

    Returns
    -------
    str
        "passed", "failed", "excluded: <reason>" or "VIOLATION: <text>".
    """
    if outcome == "failed":
        return "failed"
    if status is None:
        statuses: List[OptionalStatus] = []
    elif isinstance(status, OptionalStatus):
        statuses = [status]
    else:
        statuses = list(status)
    if cell.mandatory:
        if outcome == "passed":
            return "passed"
        return "VIOLATION: mandatory cell '%s' is %s" % (cell.cell_id, outcome)
    absent = [s for s in statuses if not s.available]
    if absent:
        if outcome == "skipped":
            return "excluded: " + "; ".join(s.reason for s in absent)
        return "VIOLATION: optional cell '%s' is %s while a dependency is absent (%s)" % (
            cell.cell_id,
            outcome,
            "; ".join(s.module for s in absent),
        )
    if outcome == "passed":
        return "passed"
    return "VIOLATION: optional cell '%s' is %s while its dependency is installed" % (
        cell.cell_id,
        outcome,
    )


def outcomes_from_junit(paths: Iterable[object]) -> Dict[str, Outcome]:
    """Read junit files and merge the cases.

    Parameters
    ----------
    paths : iterable of path-like
        Junit XML files.

    Returns
    -------
    dict
        Key ``"<classname>::<name>"``, value the outcome. For one key in several files,
        the worse outcome wins (failed, skipped, passed).
    """
    merged: Dict[str, Outcome] = {}
    for path in paths:
        root = ET.parse(str(path)).getroot()
        for case in root.iter("testcase"):
            if case.find("failure") is not None or case.find("error") is not None:
                outcome: Outcome = "failed"
            elif case.find("skipped") is not None:
                outcome = "skipped"
            else:
                outcome = "passed"
            key = "%s::%s" % (case.get("classname", ""), case.get("name", ""))
            old = merged.get(key)
            if old is None or _SEVERITY[outcome] > _SEVERITY[old]:
                merged[key] = outcome
    return merged


def _ref_matches(ref: Tuple[str, str, Optional[str]], classname: str, name: str) -> bool:
    path, func, selector = ref
    module = (path[:-3] if path.endswith(".py") else path).replace("/", ".")
    if classname == "":
        return name == module  # module-level skip
    if classname != module and not classname.startswith(module + "."):
        return False
    if name.split("[", 1)[0] != func:
        return False
    return selector is None or selector in name


def cell_outcome(cell: Cell, cases: Mapping[str, str]) -> Outcome:
    """Compute the outcome of a cell from merged junit cases.

    Parameters
    ----------
    cell : Cell
        The cell.
    cases : mapping
        Result of ``outcomes_from_junit``.

    Returns
    -------
    str
        "failed" if one matching case failed, else "skipped" if one was skipped, else
        "passed" if at least ``min_cases`` passed, else "absent".
    """
    matched: Dict[str, str] = {}
    for key, outcome in cases.items():
        classname, _, name = key.partition("::")
        if any(_ref_matches(ref, classname, name) for ref in cell.refs):
            matched[key] = outcome
    values = list(matched.values())
    if "failed" in values:
        return "failed"
    if "skipped" in values:
        return "skipped"
    if values.count("passed") >= max(cell.min_cases, 1):
        return "passed"
    return "absent"


def rows_without_real_cell(matrix: Iterable[Cell]) -> List[str]:
    """List the spec rows that have no real (non-stub) cell.

    Parameters
    ----------
    matrix : iterable of Cell
        The cells.

    Returns
    -------
    list of str
        Rows of ``SPEC_ROWS`` in order.
    """
    covered = {c.spec_row for c in matrix if not c.stub}
    return [row for row in SPEC_ROWS if row not in covered]


def _cell_statuses(cell: Cell, statuses: Mapping[str, OptionalStatus]) -> List[OptionalStatus]:
    return [statuses[k] for k in cell_dependency_keys(cell)]


def render_table(
    matrix: Iterable[Cell],
    statuses: Mapping[str, OptionalStatus],
    outcomes: Optional[Mapping[str, str]] = None,
) -> str:
    """Render the matrix as text for ``verification.md``.

    Parameters
    ----------
    matrix : iterable of Cell
        The cells.
    statuses : dict
        Status of each optional dependency, by key.
    outcomes : mapping, optional
        Outcome of each cell by ``cell_id``. A missing cell counts as "absent".

    Returns
    -------
    str
        The table. Stub cells are in a separate block.
    """
    cells = list(matrix)
    outcomes = outcomes or {}

    def result(cell: Cell) -> str:
        return classify_cell(
            cell, outcomes.get(cell.cell_id, "absent"), _cell_statuses(cell, statuses)
        )

    lines = ["| Spec row | Cell | Level | Result |", "|---|---|---|---|"]
    for row in SPEC_ROWS:
        real = [c for c in cells if c.spec_row == row and not c.stub]
        if not real:
            lines.append("| %s | (none) | - | VIOLATION: no real cell |" % row)
        for c in real:
            level = "mandatory" if c.mandatory else "optional"
            lines.append("| %s | %s | %s | %s |" % (row, c.cell_id, level, result(c)))
    for c in cells:
        if c.spec_row not in SPEC_ROWS and not c.stub:
            level = "mandatory" if c.mandatory else "optional"
            lines.append("| %s | %s | %s | %s |" % (c.spec_row, c.cell_id, level, result(c)))
    lines.append("")
    lines.append("Boundary stub cells (Python side; not a native or Surfe result)")
    stubs = [c for c in cells if c.stub]
    if not stubs:
        lines.append("(none)")
    for c in stubs:
        lines.append("- %s (%s): %s" % (c.cell_id, c.spec_row, result(c)))
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Print the matrix table for the given junit files.

    Parameters
    ----------
    argv : sequence of str, optional
        Paths of junit files. Default is ``sys.argv[1:]``.

    Returns
    -------
    int
        1 when a cell is a violation, else 0.
    """
    args = list(sys.argv[1:] if argv is None else argv)
    cases = outcomes_from_junit(args)
    statuses = {key: optional_status(key) for key in OPTIONAL_DEPENDENCIES}
    outcomes = {c.cell_id: cell_outcome(c, cases) for c in MATRIX}
    print(render_table(MATRIX, statuses, outcomes))
    print("Stub cells test the Python side of the boundary only. They are not a native or "
          "Surfe pass.")
    violation = any(
        classify_cell(c, outcomes[c.cell_id], _cell_statuses(c, statuses)).startswith(
            "VIOLATION"
        )
        for c in MATRIX
    )
    return 1 if violation else 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
