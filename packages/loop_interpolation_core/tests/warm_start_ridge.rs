//! Tests of the start value (`set_start`) and the ridge factor
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
    p.set_start(&cold.values_f64()).unwrap();
    let (warm, warm_report) = run(p, settings());
    assert!(warm_report.converged);
    println!(
        "start at the solution: {} iterations, cold: {}",
        warm_report.iterations, cold_report.iterations
    );
    // The recursive residual of the cold solve is below the tolerance. The
    // true residual of the same field can be a little above it, so a few
    // iterations are allowed.
    assert!(
        warm_report.iterations <= 3,
        "warm {} cold {}",
        warm_report.iterations,
        cold_report.iterations
    );
    let cold = cold.values_f64();
    assert!(max_diff(&warm.values_f64(), &cold) <= 1e-6 * range_of(&cold));
}

/// Data for the nested test: value = a smooth function, at points inside
/// [1, 18]^3 (physical units, the same for both grids).
fn add_nested_data(p: &mut Problem) {
    let mut points = Vec::new();
    let mut values = Vec::new();
    for a in 0..6 {
        for b in 0..6 {
            for c in 0..4 {
                let q = [
                    1.0 + 3.2 * a as f64,
                    1.5 + 3.1 * b as f64,
                    2.0 + 4.7 * c as f64,
                ];
                values.push(0.1 * q[0] + 0.05 * q[1] * q[1] / 18.0 + 0.2 * (q[2] * 0.3).sin());
                points.extend_from_slice(&q);
            }
        }
    }
    let weights = vec![1.0; values.len()];
    p.add_value(&points, &values, &weights).unwrap();
}

/// Nested iteration: solve on 10^3, evaluate that field at the nodes of a
/// 20^3 grid with the same extent, and use it as the start of the 20^3
/// solve. The result equals a cold 20^3 solve within SC-003 (1e-6 of the
/// field range).
#[test]
fn nested_10_to_20_equals_a_cold_solve() {
    const EXTENT: f64 = 19.0;
    let mut coarse = Problem::new([0.0; 3], [10; 3], [EXTENT / 9.0; 3]).unwrap();
    add_nested_data(&mut coarse);
    let (coarse_field, coarse_report) = run(coarse, settings());
    assert!(coarse_report.converged);

    // Fine nodes in the node order of the grid (x fastest).
    let mut nodes = Vec::with_capacity(3 * 20 * 20 * 20);
    for k in 0..20 {
        for j in 0..20 {
            for i in 0..20 {
                nodes.extend_from_slice(&[i as f64, j as f64, k as f64]);
            }
        }
    }
    let start = coarse_field.evaluate(&nodes).unwrap();

    let mut cold_problem = Problem::new([0.0; 3], [20; 3], [1.0; 3]).unwrap();
    add_nested_data(&mut cold_problem);
    let (cold, cold_report) = run(cold_problem, settings());
    assert!(cold_report.converged);

    let mut nested_problem = Problem::new([0.0; 3], [20; 3], [1.0; 3]).unwrap();
    add_nested_data(&mut nested_problem);
    nested_problem.set_start(&start).unwrap();
    let (nested, nested_report) = run(nested_problem, settings());
    assert!(nested_report.converged);
    println!(
        "nested: {} iterations (+{} on 10^3), cold: {}",
        nested_report.iterations, coarse_report.iterations, cold_report.iterations
    );

    let cold = cold.values_f64();
    let diff = max_diff(&nested.values_f64(), &cold);
    assert!(diff <= 1e-6 * range_of(&cold), "difference {diff}");
    assert!(nested_report.iterations <= cold_report.iterations);
}

#[test]
fn cg_warm_start_from_a_rough_start_gives_the_same_field() {
    let (cold, _) = run(base_problem(), settings());
    let mut p = base_problem();
    p.set_start(&rough_start()).unwrap();
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
    p.set_start(&cold.values_f64()).unwrap();
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
        p.set_start(&start).unwrap();
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
fn start_is_checked() {
    let mut p = base_problem();
    let short = vec![0.0; N * N * N - 1];
    assert!(matches!(
        p.set_start(&short),
        Err(CoreError::InvalidInput { .. })
    ));
    let mut bad = vec![0.0; N * N * N];
    bad[7] = f64::NAN;
    match p.set_start(&bad) {
        Err(CoreError::InvalidInput { constraint, .. }) => {
            assert_eq!(constraint, "start[7]")
        }
        other => panic!("expected InvalidInput, got {other:?}"),
    }
}

#[test]
fn start_is_in_the_memory_estimate() {
    let s = settings();
    let mut p = base_problem();
    let without = p.estimate_memory(&s).unwrap().total;
    p.set_start(&vec![0.0; N * N * N]).unwrap();
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
