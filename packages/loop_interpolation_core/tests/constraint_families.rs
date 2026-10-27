//! Tests of the constraint families on a solved field (T024).
//! They check: tangent and interface properties, point evaluation, a large
//! grid origin, and the thread count.
//! Some tests use `Field::evaluate` and `Field::evaluate_gradient`.
//! These tests fail until those functions are implemented (T025).

use loop_interpolation_core::{CoreError, Field, Problem, Settings};

const N: usize = 12;

/// Settings with a tight tolerance, so that the solve error is small.
fn settings() -> Settings {
    Settings {
        tol: 1e-12,
        max_iter: 5000,
        ..Settings::default()
    }
}

/// Add `shift` to each coordinate of a flat array of points.
fn shifted(points: &[f64], shift: f64) -> Vec<f64> {
    points.iter().map(|v| v + shift).collect()
}

/// Points for the value constraints: two planes, z = 2 and z = 9.
fn value_data() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let xy = [(2.0, 2.0), (8.0, 2.0), (2.0, 8.0), (8.0, 8.0), (5.0, 5.0)];
    let mut p = vec![];
    let mut v = vec![];
    for &(x, y) in &xy {
        p.extend_from_slice(&[x, y, 2.0]);
        v.push(0.0);
        p.extend_from_slice(&[x, y, 9.0]);
        v.push(7.0);
    }
    let w = vec![1.0; v.len()];
    (p, v, w)
}

/// Points of one interface (same id). They are at z = 5.5.
fn interface_data() -> (Vec<f64>, Vec<f64>) {
    let p = vec![
        3.0, 3.0, 5.5, //
        7.0, 3.0, 5.5, //
        3.0, 7.0, 5.5, //
        7.0, 7.0, 5.5,
    ];
    let ids = vec![1.0; 4];
    (p, ids)
}

/// Points and unit tangent vectors (direction x).
fn tangent_data() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let p = vec![4.0, 4.0, 4.0, 6.0, 6.0, 7.0, 5.0, 3.0, 6.0];
    let t = vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let w = vec![1.0; 3];
    (p, t, w)
}

/// Build and solve a problem. All coordinates move by `shift`.
fn solve_shifted(shift: f64) -> Field {
    let mut p = Problem::new([shift; 3], [N, N, N], [1.0; 3]).unwrap();
    let (vp, vv, vw) = value_data();
    p.add_value(&shifted(&vp, shift), &vv, &vw).unwrap();
    let (ip, ids) = interface_data();
    p.add_interface(&shifted(&ip, shift), &ids, 1.0).unwrap();
    let (tp, tv, tw) = tangent_data();
    p.add_tangent(&shifted(&tp, shift), &tv, &tw).unwrap();
    let (field, _report) = p.solve(settings()).unwrap();
    field
}

fn range_of(values: &[f64]) -> f64 {
    let lo = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    hi - lo
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

/// Scenario 5: tangent vectors and interface pairs.
#[test]
fn tangent_and_interface_properties() {
    let field = solve_shifted(0.0);
    let range = range_of(&field.values_f64());

    // The gradient is normal to each tangent vector.
    let (tp, tv, _) = tangent_data();
    let grad = field.evaluate_gradient(&tp).unwrap();
    for i in 0..tp.len() / 3 {
        let g = &grad[3 * i..3 * i + 3];
        let t = &tv[3 * i..3 * i + 3];
        let len = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        let dot = g[0] * t[0] + g[1] * t[1] + g[2] * t[2];
        assert!(len > 0.0, "gradient length is zero at tangent point {i}");
        assert!(
            (dot / len).abs() < 1e-3,
            "tangent {i}: dot / |g| = {}",
            dot / len
        );
    }

    // The interface points have the same field value.
    let (ip, _) = interface_data();
    let vals = field.evaluate(&ip).unwrap();
    let spread = range_of(&vals);
    assert!(
        spread <= 1e-4 * range,
        "interface spread {spread} is more than 1e-4 of the range {range}"
    );
}

/// Scenario 8: point evaluation.
#[test]
fn evaluate_matches_nodes_and_gradient() {
    let field = solve_shifted(0.0);
    let values = field.values_f64();
    let g = &field.grid;

    // Evaluation at the nodes gives the node values.
    let mut pts = Vec::with_capacity(3 * values.len());
    for flat in 0..g.n_nodes() {
        let (i, j, k) = g.node_ijk(flat);
        pts.push(g.origin[0] + i as f64 * g.step_vector[0]);
        pts.push(g.origin[1] + j as f64 * g.step_vector[1]);
        pts.push(g.origin[2] + k as f64 * g.step_vector[2]);
    }
    let at_nodes = field.evaluate(&pts).unwrap();
    let range = range_of(&values);
    assert!(max_abs_diff(&at_nodes, &values) <= 1e-9 * range.max(1.0));

    // The gradient equals central differences of the field.
    let h = 1e-3;
    let probes = [
        [3.37, 4.21, 5.43],
        [6.61, 2.77, 8.19],
        [9.13, 7.49, 3.31],
        [2.29, 9.83, 6.57],
    ];
    let mut probe_flat = vec![];
    for p in &probes {
        probe_flat.extend_from_slice(p);
    }
    let grad = field.evaluate_gradient(&probe_flat).unwrap();
    let mut fd = vec![0.0; grad.len()];
    for (n, p) in probes.iter().enumerate() {
        for axis in 0..3 {
            let mut hi = *p;
            let mut lo = *p;
            hi[axis] += h;
            lo[axis] -= h;
            let fh = field.evaluate(&hi).unwrap()[0];
            let fl = field.evaluate(&lo).unwrap()[0];
            fd[3 * n + axis] = (fh - fl) / (2.0 * h);
        }
    }
    let largest = grad.iter().cloned().map(f64::abs).fold(0.0, f64::max);
    assert!(
        max_abs_diff(&grad, &fd) <= 1e-3 * largest,
        "gradient and central differences differ"
    );

    // A point outside the grid gives an error with the point index.
    let outside = [5.0, 5.0, 5.0, -1.0, 5.0, 5.0];
    assert_eq!(
        field.evaluate(&outside).unwrap_err(),
        CoreError::PointOutsideGrid { index: 1 }
    );
    assert_eq!(
        field.evaluate_gradient(&outside).unwrap_err(),
        CoreError::PointOutsideGrid { index: 1 }
    );
}

/// Scenario 9: a large grid origin gives the same field.
#[test]
fn large_origin_gives_same_field() {
    let base = solve_shifted(0.0).values_f64();
    let moved = solve_shifted(5_000_000.0).values_f64();
    let range = range_of(&base);
    assert!(
        max_abs_diff(&base, &moved) <= 1e-6 * range,
        "the shifted field differs from the base field"
    );
}

/// Scenario 7: the thread count does not change the field.
#[cfg(feature = "parallel")]
#[test]
fn thread_count_does_not_change_field() {
    use rayon::ThreadPoolBuilder;

    let run = |threads: usize| -> Vec<f64> {
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| solve_shifted(0.0).values_f64())
    };
    let one = run(1);
    let four = run(4);
    let range = range_of(&one);
    assert!(
        max_abs_diff(&one, &four) <= 1e-6 * range,
        "fields from 1 and 4 threads differ"
    );
}
