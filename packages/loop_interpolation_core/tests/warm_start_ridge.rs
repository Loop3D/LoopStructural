//! Tests of the start value (`set_initial_guess`) and the ridge factor
//! (`set_ridge_factor`). A warm start must give the same field as a cold
//! start, and a start at the solution must need almost no iterations.

use loop_interpolation_core::{CoreError, Field, Precision, Problem, Settings, SolveReport, Step};

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

/// Problem with an active upper bound, so that the session uses ADMM.
fn bound_problem() -> Problem {
    let mut p = base_problem();
    p.add_value_bound(&pt(15, 9, 9), &[-BIG], &[10.0], &[1.0])
        .unwrap();
    p
}

fn settings() -> Settings {
    Settings {
        tol: 1e-10,
        max_iter: 20_000,
        ..Settings::default()
    }
}

fn run(p: Problem, s: Settings) -> (Field, SolveReport) {
    let mut session = p.start(s).unwrap();
    for _ in 0..1_000_000 {
        if let Step::Done(_) = session.step(50).unwrap() {
            return session.finish().unwrap();
        }
    }
    panic!("the session does not finish");
}

fn range_of(v: &[f64]) -> f64 {
    let lo = v.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    hi - lo
}

fn max_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

/// A start value that is far from the solution, but not zero.
fn rough_start() -> Vec<f64> {
    (0..N * N * N)
        .map(|n| 3.0 * ((n % 7) as f64) - 5.0)
        .collect()
}

#[test]
fn cg_warm_start_at_the_solution_needs_few_iterations() {
    let (cold, cold_report) = run(base_problem(), settings());
    assert!(cold_report.converged);
    let mut p = base_problem();
    p.set_initial_guess(&cold.values_f64()).unwrap();
    let (warm, warm_report) = run(p, settings());
    assert!(warm_report.converged);
    assert!(
        warm_report.iterations * 10 < cold_report.iterations,
        "warm {} cold {}",
        warm_report.iterations,
        cold_report.iterations
    );
    let cold = cold.values_f64();
    assert!(max_diff(&warm.values_f64(), &cold) <= 1e-6 * range_of(&cold));
}

#[test]
fn cg_warm_start_from_a_rough_start_gives_the_same_field() {
    let (cold, _) = run(base_problem(), settings());
    let mut p = base_problem();
    p.set_initial_guess(&rough_start()).unwrap();
    let (warm, report) = run(p, settings());
    assert!(report.converged);
    let cold = cold.values_f64();
    let diff = max_diff(&warm.values_f64(), &cold);
    assert!(diff <= 1e-6 * range_of(&cold), "difference {diff}");
}

#[test]
fn f32_warm_start_at_the_solution_converges() {
    let s = Settings {
        precision: Precision::F32,
        // The f32 residual of this problem stops near 2e-5.
        tol: 1e-4,
        ..settings()
    };
    let (cold, cold_report) = run(base_problem(), s.clone());
    assert!(cold_report.converged, "{cold_report:?}");
    let mut p = base_problem();
    p.set_initial_guess(&cold.values_f64()).unwrap();
    let (_, warm_report) = run(p, s);
    assert!(warm_report.converged);
    assert!(warm_report.iterations <= cold_report.iterations);
}

#[test]
fn admm_warm_start_gives_the_same_field() {
    let (cold, cold_report) = run(bound_problem(), settings());
    assert!(cold_report.converged);
    let cold = cold.values_f64();
    let range = range_of(&cold);
    for start in [cold.clone(), rough_start()] {
        let mut p = bound_problem();
        p.set_initial_guess(&start).unwrap();
        let (warm, report) = run(p, settings());
        assert!(report.converged);
        let warm = warm.values_f64();
        assert!(warm[idx(15, 9, 9)] <= 10.0 + 1e-4 * range);
        // ADMM stops at a tolerance, so the two fields agree less closely.
        let diff = max_diff(&warm, &cold);
        assert!(diff <= 1e-3 * range, "difference {diff}");
    }
}

#[test]
fn initial_guess_is_checked() {
    let mut p = base_problem();
    let short = vec![0.0; N * N * N - 1];
    assert!(matches!(
        p.set_initial_guess(&short),
        Err(CoreError::InvalidInput { .. })
    ));
    let mut bad = vec![0.0; N * N * N];
    bad[7] = f64::NAN;
    match p.set_initial_guess(&bad) {
        Err(CoreError::InvalidInput { constraint, .. }) => {
            assert_eq!(constraint, "initial_guess[7]")
        }
        other => panic!("expected InvalidInput, got {other:?}"),
    }
}

#[test]
fn initial_guess_is_in_the_memory_estimate() {
    let s = settings();
    let mut p = base_problem();
    let without = p.estimate_memory(&s).unwrap().total;
    p.set_initial_guess(&vec![0.0; N * N * N]).unwrap();
    let with = p.estimate_memory(&s).unwrap().total;
    assert!(with >= without + 8 * (N * N * N) as u64);
}

#[test]
fn large_ridge_factor_pulls_the_field_to_zero() {
    let norm = |ridge: f64| {
        let mut p = base_problem();
        p.set_ridge_factor(ridge).unwrap();
        let (f, r) = run(p, settings());
        assert!(r.converged);
        f.values_f64().iter().map(|v| v * v).sum::<f64>().sqrt()
    };
    let none = norm(0.0);
    let strong = norm(1.0);
    assert!(none.is_finite());
    assert!(strong < 0.9 * none, "ridge 0: {none}, ridge 1: {strong}");
}

#[test]
fn ridge_factor_is_checked() {
    let mut p = base_problem();
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            p.set_ridge_factor(bad),
            Err(CoreError::InvalidInput { .. })
        ));
    }
    assert!(p.set_ridge_factor(0.0).is_ok());
}
