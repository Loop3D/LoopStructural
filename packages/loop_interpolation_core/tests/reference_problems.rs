//! Reference problem tests (SC-001, SC-007, SC-008, SC-009, SC-012).
//! The tests use the counting allocator of `common`.
//! The CI tests use small iteration counts. The big tests are ignored.

mod common;

use std::sync::{Mutex, MutexGuard};
use std::time::Instant;

use loop_interpolation_core::{Precision, Problem, Settings};

/// The counting allocator is global for this test binary.
/// Each test takes this lock, so that only one test measures at a time.
static SERIAL: Mutex<()> = Mutex::new(());

fn lock() -> MutexGuard<'static, ()> {
    // A failed test poisons the lock. Ignore the poison.
    SERIAL.lock().unwrap_or_else(|e| e.into_inner())
}

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
}

/// Build a problem with the data generator of `session_r1` in tests/memory.rs:
/// a plane field, `n_value` value rows and `n_normal` gradient rows.
fn build(nodes: usize, n_value: usize, n_normal: usize) -> Problem {
    let mut p = Problem::new([0.0; 3], [nodes; 3], [1.0; 3]).expect("grid must be valid");
    let mut rng = Lcg(12345);
    let extent = (nodes - 1) as f64;
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
    p.add_value(&pts, &vals, &vec![1.0; n_value])
        .expect("value rows must be valid");
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
    p.add_gradient(&pts, &dirs, &vec![1.0; n_normal])
        .expect("gradient rows must be valid");
    p
}

/// Start the thread pool (feature `parallel`) once, before the first measurement.
/// The pool uses memory one time only. The memory does not belong to a solve.
fn warm_up() {
    static ONCE: std::sync::Once = std::sync::Once::new();
    ONCE.call_once(|| {
        let p = build(12, 20, 5);
        let s = Settings {
            max_iter: 3,
            ..Settings::default()
        };
        let _ = p.solve(s);
    });
}

/// Build and solve. Return the peak bytes above the level at the start
/// (the build of the Problem is included), the estimate total, and the
/// solve time in seconds.
fn measure(nodes: usize, n_value: usize, n_normal: usize, settings: Settings) -> (usize, u64, f64) {
    warm_up();
    common::reset_peak();
    let start = common::current();
    let p = build(nodes, n_value, n_normal);
    let estimate = p
        .estimate_memory(&settings)
        .expect("estimate must succeed")
        .total;
    let t0 = Instant::now();
    let (field, _report) = p.solve(settings).expect("solve must succeed");
    let secs = t0.elapsed().as_secs_f64();
    drop(field);
    let peak = common::peak().saturating_sub(start);
    (peak, estimate, secs)
}

/// Settings for a short solve. All vectors exist from the start,
/// so the peak is the same as for a full solve.
fn short_settings() -> Settings {
    Settings {
        max_iter: 20,
        ..Settings::default()
    }
}

/// Check that the estimate is not lower than the peak and at most 1.25 times it.
fn check_estimate(sizes: &[usize]) {
    for &n in sizes {
        let (peak, estimate, _) = measure(n, 500, 125, short_settings());
        let ratio = estimate as f64 / peak as f64;
        println!("n={n} estimate={estimate} peak={peak} ratio={ratio:.3}");
        assert!(estimate as usize >= peak, "n={n}: estimate below peak");
        assert!(ratio <= 1.25, "n={n}: estimate ratio {ratio} above 1.25");
    }
}

/// Least-squares fit y = a + b x. Return the R^2 of the fit.
fn r_squared(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let sxx: f64 = x.iter().map(|v| (v - mx).powi(2)).sum();
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let syy: f64 = y.iter().map(|v| (v - my).powi(2)).sum();
    let b = sxy / sxx;
    let a = my - b * mx;
    let ss_res: f64 = x
        .iter()
        .zip(y)
        .map(|(xi, yi)| (yi - (a + b * xi)).powi(2))
        .sum();
    if syy == 0.0 { 1.0 } else { 1.0 - ss_res / syy }
}

/// Check that the peak is linear in the dof for the given grid sizes.
fn check_linear(sizes: &[usize]) {
    let mut dofs = Vec::new();
    let mut peaks = Vec::new();
    for &n in sizes {
        let (peak, _, _) = measure(n, 500, 125, short_settings());
        println!("n={n} dof={} peak={peak}", n * n * n);
        dofs.push((n * n * n) as f64);
        peaks.push(peak as f64);
    }
    let r2 = r_squared(&dofs, &peaks);
    println!("r_squared={r2}");
    assert!(r2 >= 0.99, "R^2 {r2} below 0.99");
}

/// SC-009: ten grid sizes from 10^3 to 100^3.
#[test]
fn estimate_matches_peak_ci() {
    let _guard = lock();
    check_estimate(&[10, 15, 20, 30, 40, 50, 60, 70, 85, 100]);
}

/// SC-012: peak against dof is a line (CI sizes).
#[test]
fn peak_is_linear_in_dof_ci() {
    let _guard = lock();
    check_linear(&[30, 45, 60]);
}

/// SC-012: the same check for 50^3, 100^3 and 150^3.
#[test]
#[ignore]
fn peak_is_linear_in_dof_large() {
    let _guard = lock();
    check_linear(&[50, 100, 150]);
}

/// SC-001: R1 peak is at most 100 MB. The time is printed.
#[test]
#[ignore]
fn r1_peak_within_100_mb() {
    let _guard = lock();
    let settings = Settings {
        tol: 1e-8,
        max_iter: 5000,
        ..Settings::default()
    };
    let (peak, estimate, secs) = measure(100, 2000, 500, settings);
    println!("R1 peak={peak} estimate={estimate} time_s={secs:.3}");
    assert!(peak <= 100_000_000, "R1 peak {peak} above 100 MB");
}

/// SC-008: 250^3 in F32 with a budget of 2 GB. The peak is at most 1.5 GB.
#[test]
#[ignore]
fn f32_250_cubed_within_budget() {
    let _guard = lock();
    let settings = Settings {
        precision: Precision::F32,
        tol: 1e-5,
        max_iter: 5000,
        budget_bytes: Some(2_000_000_000),
        ..Settings::default()
    };
    let (peak, estimate, secs) = measure(250, 500, 125, settings);
    println!("250^3 F32 peak={peak} estimate={estimate} time_s={secs:.3}");
    assert!(peak <= 1_500_000_000, "peak {peak} above 1.5 GB");
}

/// SC-009 extension: sizes up to 200^3.
#[test]
#[ignore]
fn estimate_matches_peak_to_200() {
    let _guard = lock();
    check_estimate(&[10, 20, 40, 60, 80, 100, 120, 150, 175, 200]);
}

/// SC-007 (report only): solve R1 fully and write the field as
/// little-endian f64 to `target/r1_field.f64` for the browser comparison.
#[test]
#[ignore]
fn r1_write_field() {
    let _guard = lock();
    let p = build(100, 2000, 500);
    let settings = Settings {
        tol: 1e-8,
        max_iter: 5000,
        ..Settings::default()
    };
    let t0 = Instant::now();
    let (field, report) = p.solve(settings).expect("solve must succeed");
    let secs = t0.elapsed().as_secs_f64();
    println!(
        "R1 time_s={secs:.3} iterations={} converged={}",
        report.iterations, report.converged
    );
    let mut bytes = Vec::with_capacity(8 * 1_000_000);
    for v in field.values_f64() {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target");
    std::fs::create_dir_all(&dir).expect("target directory must be created");
    let path = dir.join("r1_field.f64");
    std::fs::write(&path, &bytes).expect("field file must be written");
    println!("written {} bytes to {}", bytes.len(), path.display());
}
