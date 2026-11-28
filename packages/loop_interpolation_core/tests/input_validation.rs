//! Input validation tests for `Problem` (FR-018, story 1 scenario 10).
//!
//! Each test builds a `Problem` and calls one `add_*` method with bad input.
//! The test checks the error variant and that `constraint` names the family
//! and the row, for example `value[3]`.
//!
//! A rotation cannot be given. The API has no rotation argument. The grid
//! has an origin, node counts and a step only. So no test for rotation
//! exists.

use loop_interpolation_core::{CoreError, Problem, Settings};

/// Make a valid 5 x 5 x 5 grid with step 1 and origin 0.
fn problem() -> Problem {
    Problem::new([0.0; 3], [5, 5, 5], [1.0; 3]).expect("valid grid")
}

/// Make `n` points inside the grid. Each point has 3 values.
fn pts(n: usize) -> Vec<f64> {
    (0..n)
        .flat_map(|i| [1.0 + (i % 3) as f64, 1.5, 1.5])
        .collect()
}

/// Check that the error is `InvalidInput` and that `constraint` has the
/// family name at the start and the row in brackets.
fn assert_invalid(r: Result<(), CoreError>, family: &str, row: usize) {
    match r {
        Err(CoreError::InvalidInput { constraint, .. }) => {
            let tag = format!("[{row}]");
            assert!(
                constraint.starts_with(family) && constraint.contains(&tag),
                "constraint was {constraint:?}, expected {family}{tag}"
            );
        }
        other => panic!("expected InvalidInput, got {other:?}"),
    }
}

/// True if one warning text contains the text `needle`.
fn has_warning(p: &Problem, needle: &str) -> bool {
    p.warnings().iter().any(|w| w.contains(needle))
}

#[test]
fn nan_in_point() {
    let mut p = problem();
    let mut x = pts(5);
    x[3 * 3 + 1] = f64::NAN; // row 3
    let r = p.add_value(&x, &[0.0; 5], &[1.0; 5]);
    assert_invalid(r, "value", 3);
}

#[test]
fn infinite_value() {
    let mut p = problem();
    let mut v = vec![0.0; 5];
    v[3] = f64::INFINITY;
    let r = p.add_value(&pts(5), &v, &[1.0; 5]);
    assert_invalid(r, "value", 3);
}

#[test]
fn array_of_wrong_length() {
    let mut p = problem();
    // 5 values for the points: not a multiple of 3.
    let r = p.add_value(&[1.0; 5], &[0.0; 2], &[1.0; 2]);
    assert!(
        matches!(r, Err(CoreError::InvalidInput { .. })),
        "got {r:?}"
    );
    // Values array has a different length than the points.
    let r = p.add_value(&pts(4), &[0.0; 3], &[1.0; 4]);
    assert!(
        matches!(r, Err(CoreError::InvalidInput { .. })),
        "got {r:?}"
    );
}

fn zero_vector_input() -> Vec<f64> {
    let mut v: Vec<f64> = (0..5).flat_map(|_| [0.0, 0.0, 1.0]).collect();
    v[3 * 3..3 * 3 + 3].copy_from_slice(&[0.0, 0.0, 0.0]); // row 3
    v
}

#[test]
fn zero_vector_gradient() {
    let mut p = problem();
    let r = p.add_gradient(&pts(5), &zero_vector_input(), &[1.0; 5]);
    assert_invalid(r, "gradient", 3);
}

#[test]
fn zero_vector_normal() {
    let mut p = problem();
    let r = p.add_normal(&pts(5), &zero_vector_input(), &[1.0; 5]);
    assert_invalid(r, "normal", 3);
}

#[test]
fn zero_vector_tangent() {
    let mut p = problem();
    let r = p.add_tangent(&pts(5), &zero_vector_input(), &[1.0; 5]);
    assert_invalid(r, "tangent", 3);
}

#[test]
fn weight_zero() {
    let mut p = problem();
    let mut w = vec![1.0; 5];
    w[3] = 0.0;
    let r = p.add_value(&pts(5), &[0.0; 5], &w);
    assert_invalid(r, "value", 3);
}

#[test]
fn weight_negative() {
    let mut p = problem();
    let mut w = vec![1.0; 5];
    w[3] = -1.0;
    let r = p.add_value(&pts(5), &[0.0; 5], &w);
    assert_invalid(r, "value", 3);
}

#[test]
fn value_bound_lower_above_upper() {
    let mut p = problem();
    let lower = vec![0.0; 5];
    let mut upper = vec![1.0; 5];
    upper[3] = -1.0; // lower > upper at row 3
    let r = p.add_value_bound(&pts(5), &lower, &upper, &[1.0; 5]);
    assert_invalid(r, "value_bound", 3);
}

#[test]
fn pair_bound_lower_above_upper() {
    let mut p = problem();
    let lower = vec![0.0; 5];
    let mut upper = vec![1.0; 5];
    upper[3] = -1.0; // lower > upper at row 3
    let r = p.add_pair_bound(&pts(5), &pts(5), &lower, &upper, &[1.0; 5]);
    assert_invalid(r, "pair_bound", 3);
}

#[test]
fn grid_with_two_nodes_on_one_axis() {
    let r = Problem::new([0.0; 3], [5, 2, 5], [1.0; 3]);
    assert!(
        matches!(r, Err(CoreError::GridTooSmall { axis: 1, nodes: 2 })),
        "got {:?}",
        r.map(|_| ())
    );
}

#[test]
fn grid_above_max_nodes() {
    // 2^21 * 2^21 * 2^21 = 2^63. This is more than isize::MAX / 8.
    // The check must use u128 maths and must not allocate.
    let n = 1usize << 21;
    let r = Problem::new([0.0; 3], [n, n, n], [1.0; 3]);
    assert!(
        matches!(r, Err(CoreError::GridTooLarge { .. })),
        "got {:?}",
        r.map(|_| ())
    );
}

#[test]
fn start_with_zero_constraints() {
    let p = problem();
    let r = p.start(Settings::default());
    assert!(
        matches!(r, Err(CoreError::NoConstraints)),
        "got {:?}",
        r.map(|_| ())
    );
}

// Warnings: these are not errors.

#[test]
fn point_outside_grid_is_dropped_with_warning() {
    let mut p = problem();
    let mut x = pts(4);
    x[0..3].copy_from_slice(&[100.0, 100.0, 100.0]); // one point outside
    p.add_value(&x, &[0.0; 4], &[1.0; 4])
        .expect("outside point is not an error");
    assert_eq!(p.warnings().len(), 1, "warnings: {:?}", p.warnings());
    assert!(has_warning(&p, "1"), "warnings: {:?}", p.warnings());
}

#[test]
fn conflicting_values_at_one_point_give_warning() {
    let mut p = problem();
    let x = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0]; // two rows at one point
    p.add_value(&x, &[1.0, 5.0], &[1.0, 1.0])
        .expect("conflict is not an error");
    assert!(has_warning(&p, "1"), "warnings: {:?}", p.warnings());
}

#[test]
fn step_ratio_above_100_gives_warning() {
    let p = Problem::new([0.0; 3], [5, 5, 5], [1.0, 1.0, 150.0]).expect("valid grid");
    assert!(!p.warnings().is_empty(), "expected a step ratio warning");
}
