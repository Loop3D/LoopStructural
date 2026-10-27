//! Tests of the ADMM solve session (T021, story 2).
//! Inequality bounds must hold, a bound that is not active must change
//! nothing, and bounds that cannot hold together must give `Infeasible`.
//! These tests fail until ADMM is implemented (T022).

use loop_interpolation_core::{
    CoreError, Field, Precision, Problem, Progress, Settings, SolveReport, Step, ViolatedConstraint,
};

const N: usize = 20;
/// A bound that never acts.
const BIG: f64 = 1.0e6;

/// Index of the node (i, j, k) in the field values.
fn idx(i: usize, j: usize, k: usize) -> usize {
    i + N * (j + N * k)
}

/// Node (i, j, k) as a point. The grid origin is 0 and the step is 1.
fn pt(i: usize, j: usize, k: usize) -> [f64; 3] {
    [i as f64, j as f64, k as f64]
}

/// Problem with data that rise along x: value = x on a set of nodes.
fn base_problem() -> Problem {
    let mut p = Problem::new([0.0; 3], [N, N, N], [1.0; 3]).unwrap();
    let mut points = Vec::new();
    let mut values = Vec::new();
    let mut weights = Vec::new();
    for i in 2..=17 {
        for &j in &[3usize, 9, 16] {
            for &k in &[3usize, 9, 16] {
                points.extend_from_slice(&pt(i, j, k));
                values.push(i as f64);
                weights.push(1.0);
            }
        }
    }
    p.add_value(&points, &values, &weights).unwrap();
    p
}

fn settings() -> Settings {
    Settings {
        tol: 1e-10,
        max_iter: 20_000,
        ..Settings::default()
    }
}

/// Run the session with small steps. Return the field, the report and
/// every progress that was reported.
fn run(p: Problem) -> (Field, SolveReport, Vec<Progress>) {
    let mut s = p.start(settings()).unwrap();
    let mut seen = Vec::new();
    for _ in 0..1_000_000 {
        match s.step(25).unwrap() {
            Step::Running(pr) => seen.push(pr),
            Step::Done(_) => {
                let (f, r) = s.finish().unwrap();
                return (f, r, seen);
            }
        }
    }
    panic!("the session does not finish");
}

fn range_of(v: &[f64]) -> f64 {
    let lo = v.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    hi - lo
}

/// Violation of a value bound: how far the value is outside [lo, hi].
fn value_violation(v: f64, lo: f64, hi: f64) -> f64 {
    (lo - v).max(v - hi).max(0.0)
}

/// Violation of a pair bound: how far a - b is outside [lo, hi].
fn pair_violation(a: f64, b: f64, lo: f64, hi: f64) -> f64 {
    value_violation(a - b, lo, hi)
}

/// Field without any bound. It shows that the data break the bounds.
fn free_values() -> Vec<f64> {
    let (f, r, _) = run(base_problem());
    assert!(r.converged);
    f.values_f64()
}

/// (1) Value bound: upper bound 10 on a node where the free field is about 15.
#[test]
fn value_bound_holds() {
    let free = free_values();
    let range = range_of(&free);
    let at = idx(15, 9, 9);
    assert!(free[at] > 12.0, "the data must break the bound");
    let mut p = base_problem();
    p.add_value_bound(&pt(15, 9, 9), &[-BIG], &[10.0], &[1.0])
        .unwrap();
    let (f, r, _) = run(p);
    assert!(r.converged, "ADMM did not converge");
    let v = f.values_f64()[at];
    let viol = value_violation(v, -BIG, 10.0);
    assert!(viol <= 1e-4 * range, "violation {viol}, range {range}");
}

/// (2) Pair bound A - B >= 0 against data that rise from A to B.
#[test]
fn pair_bound_holds() {
    let free = free_values();
    let range = range_of(&free);
    let (a, b) = (idx(5, 9, 9), idx(15, 9, 9));
    assert!(free[a] - free[b] < -5.0, "the data must break the bound");
    let mut p = base_problem();
    p.add_pair_bound(&pt(5, 9, 9), &pt(15, 9, 9), &[0.0], &[BIG], &[1.0])
        .unwrap();
    let (f, r, _) = run(p);
    assert!(r.converged, "ADMM did not converge");
    let v = f.values_f64();
    let viol = pair_violation(v[a], v[b], 0.0, BIG);
    assert!(viol <= 1e-4 * range, "violation {viol}, range {range}");
}

/// (5) Lower and upper bound on one row: both hold. The lower bound is
/// active on the first row, the upper bound on the second row.
#[test]
fn lower_and_upper_bound_hold() {
    let free = free_values();
    let range = range_of(&free);
    let (lo_node, hi_node) = (idx(3, 9, 9), idx(15, 9, 9));
    assert!(free[lo_node] < 5.0 && free[hi_node] > 6.0);
    let mut p = base_problem();
    let mut points = Vec::new();
    points.extend_from_slice(&pt(3, 9, 9));
    points.extend_from_slice(&pt(15, 9, 9));
    p.add_value_bound(&points, &[5.0, 5.0], &[6.0, 6.0], &[1.0, 1.0])
        .unwrap();
    let (f, r, _) = run(p);
    assert!(r.converged, "ADMM did not converge");
    let v = f.values_f64();
    for &n in &[lo_node, hi_node] {
        let viol = value_violation(v[n], 5.0, 6.0);
        assert!(viol <= 1e-4 * range, "violation {viol}, range {range}");
    }
}

/// (4) A bound that is not active changes the field by at most 1e-8 of the range.
#[test]
fn inactive_bound_changes_nothing() {
    let free = free_values();
    let range = range_of(&free);
    let mut p = base_problem();
    p.add_value_bound(&pt(15, 9, 9), &[-BIG], &[BIG], &[1.0])
        .unwrap();
    p.add_pair_bound(&pt(5, 9, 9), &pt(15, 9, 9), &[-BIG], &[BIG], &[1.0])
        .unwrap();
    let (f, r, _) = run(p);
    assert!(r.converged, "ADMM did not converge");
    let got = f.values_f64();
    assert_eq!(got.len(), free.len());
    let max_diff = got
        .iter()
        .zip(free.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(
        max_diff <= 1e-8 * range,
        "max diff {max_diff} is over 1e-8 of range {range}"
    );
}

/// (6) Two bounds that cannot hold together: `Infeasible`, at most 20 rows,
/// largest violation first. 30 nodes have the two rows, so at least 30 rows
/// are violated and the list must be cut to 20.
#[test]
fn infeasible_bounds_are_reported() {
    let mut p = base_problem();
    let mut points = Vec::new();
    for i in 2..=16 {
        for &j in &[6usize, 12] {
            points.extend_from_slice(&pt(i, j, 9));
        }
    }
    let m = points.len() / 3;
    assert_eq!(m, 30);
    // Row set 1: value at least 12. Row set 2: value at most 5.
    p.add_value_bound(&points, &vec![12.0; m], &vec![BIG; m], &vec![1.0; m])
        .unwrap();
    p.add_value_bound(&points, &vec![-BIG; m], &vec![5.0; m], &vec![1.0; m])
        .unwrap();
    match p.solve(settings()) {
        Err(CoreError::Infeasible { worst }) => {
            assert!(
                !worst.is_empty() && worst.len() <= 20,
                "len {}",
                worst.len()
            );
            assert_eq!(worst.len(), 20);
            let sorted = worst
                .windows(2)
                .all(|w: &[ViolatedConstraint]| w[0].violation >= w[1].violation);
            assert!(sorted, "not sorted by falling violation: {worst:?}");
            assert!(worst[0].violation > 0.0);
        }
        Err(e) => panic!("wrong error: {e}"),
        Ok(_) => panic!("an infeasible problem must not solve"),
    }
}

/// Every `Running` progress has a finite `max_violation`. The final report has it.
#[test]
fn progress_and_report_have_max_violation() {
    let free = free_values();
    let range = range_of(&free);
    let mut p = base_problem();
    p.add_value_bound(&pt(15, 9, 9), &[-BIG], &[10.0], &[1.0])
        .unwrap();
    let (_f, r, seen) = run(p);
    assert!(!seen.is_empty(), "no progress was reported");
    for pr in &seen {
        assert!(
            pr.max_violation.is_finite() && pr.max_violation >= 0.0,
            "bad max_violation {}",
            pr.max_violation
        );
    }
    assert!(r.max_violation.is_finite() && r.max_violation >= 0.0);
    assert!(
        r.max_violation <= 1e-4 * range,
        "final max_violation {} over 1e-4 of range {range}",
        r.max_violation
    );
}

/// `Settings` has no `hybrid` or `active_set` field. This literal names every
/// field of D8 and has no `..Default::default()`. A new field would not compile.
#[test]
fn settings_has_only_the_d8_fields() {
    let s = Settings {
        precision: Precision::F64,
        tol: 1e-8,
        max_iter: 5000,
        budget_bytes: None,
        rho: 1.0,
        inner_max_iter: 3000,
        inner_tol: 1e-6,
        admm_abs_tol: 1e-4,
        admm_rel_tol: 1e-3,
    };
    assert_eq!(s.max_iter, 5000);
}
