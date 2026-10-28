//! Memory tests with many active pair bounds (T023, story 2 scenario 3).
//! The pair bounds are active at the free answer: each bound is the free
//! difference of its two nodes. The peak heap use with the bounds must stay
//! near the peak without the bounds.

mod common;

use loop_interpolation_core::{Field, Problem, Settings, SolveReport};
use std::sync::Mutex;
use std::time::Instant;

/// Only one test at a time may use the counting allocator.
static LOCK: Mutex<()> = Mutex::new(());

/// Small deterministic generator (LCG). It returns a number in [0, 1).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }

    /// A number in 0..n.
    fn below(&mut self, n: usize) -> usize {
        ((self.next() * n as f64) as usize).min(n - 1)
    }
}

/// Node pairs. The points are at grid nodes (origin 0, step 1).
struct Pairs {
    a_pts: Vec<f64>,
    b_pts: Vec<f64>,
    a_idx: Vec<usize>,
    b_idx: Vec<usize>,
}

impl Pairs {
    fn len(&self) -> usize {
        self.a_idx.len()
    }
}

/// Make `count` pairs of different nodes. The seed is fixed.
fn make_pairs(n: usize, count: usize) -> Pairs {
    let mut rng = Lcg(777);
    let mut pairs = Pairs {
        a_pts: Vec::with_capacity(3 * count),
        b_pts: Vec::with_capacity(3 * count),
        a_idx: Vec::with_capacity(count),
        b_idx: Vec::with_capacity(count),
    };
    while pairs.len() < count {
        let a = [rng.below(n), rng.below(n), rng.below(n)];
        let b = [rng.below(n), rng.below(n), rng.below(n)];
        if a == b {
            continue;
        }
        pairs.a_pts.extend(a.iter().map(|&x| x as f64));
        pairs.b_pts.extend(b.iter().map(|&x| x as f64));
        pairs.a_idx.push(a[0] + n * (a[1] + n * a[2]));
        pairs.b_idx.push(b[0] + n * (b[1] + n * b[2]));
    }
    pairs
}

/// Pair bounds: the pairs, the lower bounds and the upper bounds.
type Bounds<'a> = Option<(&'a Pairs, &'a [f64], &'a [f64])>;

/// Build the problem like R1: value and gradient data from a plane field.
/// Add the pair bounds when `bounds` is given.
fn build(n: usize, n_value: usize, n_normal: usize, bounds: Bounds) -> Problem {
    let mut p = Problem::new([0.0; 3], [n; 3], [1.0; 3]).unwrap();
    let extent = (n - 1) as f64;
    let mut rng = Lcg(12345);
    let mut pts = Vec::with_capacity(3 * n_value);
    let mut vals = Vec::with_capacity(n_value);
    for _ in 0..n_value {
        let q = [
            rng.next() * extent,
            rng.next() * extent,
            rng.next() * extent,
        ];
        vals.push(q[0] + 0.5 * q[1] + 0.25 * q[2]);
        pts.extend_from_slice(&q);
    }
    p.add_value(&pts, &vals, &vec![1.0; n_value]).unwrap();
    let mut pts = Vec::with_capacity(3 * n_normal);
    let mut dirs = Vec::with_capacity(3 * n_normal);
    for _ in 0..n_normal {
        pts.extend_from_slice(&[
            rng.next() * extent,
            rng.next() * extent,
            rng.next() * extent,
        ]);
        dirs.extend_from_slice(&[1.0, 0.5, 0.25]);
    }
    p.add_gradient(&pts, &dirs, &vec![1.0; n_normal]).unwrap();
    if let Some((pairs, lower, upper)) = bounds {
        p.add_pair_bound(
            &pairs.a_pts,
            &pairs.b_pts,
            lower,
            upper,
            &vec![1.0; pairs.len()],
        )
        .unwrap();
    }
    p
}

fn settings() -> Settings {
    Settings {
        tol: 1e-8,
        max_iter: 20_000,
        ..Settings::default()
    }
}

/// Result of one measured solve.
struct Run {
    peak: usize,
    seconds: f64,
    field: Field,
    report: SolveReport,
}

/// Build and solve. The peak is above the level at the start and includes
/// the build of the Problem.
fn measure(n: usize, n_value: usize, n_normal: usize, bounds: Bounds) -> Run {
    common::reset_peak();
    let start = common::current();
    let t = Instant::now();
    let p = build(n, n_value, n_normal, bounds);
    let (field, report) = p.solve(settings()).unwrap();
    let seconds = t.elapsed().as_secs_f64();
    let peak = common::peak().saturating_sub(start);
    Run {
        peak,
        seconds,
        field,
        report,
    }
}

fn range_of(v: &[f64]) -> f64 {
    let lo = v.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    hi - lo
}

/// Largest violation of the lower bound `a - b >= lower`.
fn largest_violation(values: &[f64], pairs: &Pairs, lower: &[f64]) -> f64 {
    (0..pairs.len())
        .map(|i| (lower[i] - (values[pairs.a_idx[i]] - values[pairs.b_idx[i]])).max(0.0))
        .fold(0.0, f64::max)
}

/// Solve without bounds, then with active bounds. Return both runs, the
/// pairs, the lower bounds and the range of the free field.
fn free_then_bounded(
    n: usize,
    n_value: usize,
    n_normal: usize,
    n_pairs: usize,
) -> (Run, Run, Pairs, Vec<f64>, f64) {
    let free = measure(n, n_value, n_normal, None);
    let values = free.field.values_f64();
    let range = range_of(&values);
    let pairs = make_pairs(n, n_pairs);
    // The lower bound is the free difference, so the bound is active.
    let lower: Vec<f64> = (0..pairs.len())
        .map(|i| values[pairs.a_idx[i]] - values[pairs.b_idx[i]])
        .collect();
    let upper = vec![f64::INFINITY; pairs.len()];
    drop(values);
    let bounded = measure(n, n_value, n_normal, Some((&pairs, &lower, &upper)));
    (free, bounded, pairs, lower, range)
}

/// 50^3 grid with 2 500 active pair bounds: the peak stays near the free peak.
#[test]
fn medium_active_bounds_keep_peak_near_free_peak() {
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let (free, bounded, pairs, lower, range) = free_then_bounded(50, 500, 125, 2500);
    println!(
        "free peak={} bounded peak={} iterations={} max_violation={}",
        free.peak, bounded.peak, bounded.report.iterations, bounded.report.max_violation
    );
    assert!(bounded.report.converged, "ADMM did not converge");
    assert!(
        bounded.report.max_violation <= 1e-4 * range,
        "report violation {}, range {range}",
        bounded.report.max_violation
    );
    let own = largest_violation(&bounded.field.values_f64(), &pairs, &lower);
    assert!(own <= 1e-4 * range, "violation {own}, range {range}");
    assert!(
        bounded.peak as f64 <= 1.2 * free.peak as f64,
        "peak with bounds {}, peak without bounds {}",
        bounded.peak,
        free.peak
    );
}

/// Run the 100^3 problem with `n_pairs` active pair bounds. Print the numbers.
fn big_run(n_pairs: usize) {
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let n = 100usize;
    let (free, bounded, pairs, lower, range) = free_then_bounded(n, 2000, 500, n_pairs);
    let dof = n * n * n;
    let own = largest_violation(&bounded.field.values_f64(), &pairs, &lower);
    println!("free: peak_bytes={} seconds={:.2}", free.peak, free.seconds);
    println!(
        "bounded ({n_pairs} pairs): peak_bytes={} bytes_per_dof={} seconds={:.2}",
        bounded.peak,
        bounded.peak as f64 / dof as f64,
        bounded.seconds
    );
    println!(
        "ratio={:.3} report_peak_bytes={} iterations={} converged={}",
        bounded.peak as f64 / free.peak as f64,
        bounded.report.peak_bytes,
        bounded.report.iterations,
        bounded.report.converged
    );
    println!(
        "max_violation={} (report {}) range={range}",
        own, bounded.report.max_violation
    );
}

/// Problem R2: 100^3 grid, 200 active pair bounds. It only prints.
#[test]
#[ignore]
fn r2() {
    big_run(200);
}

/// Problem R3: 100^3 grid, 20 000 active pair bounds. It only prints.
#[test]
#[ignore]
fn r3() {
    big_run(20_000);
}
