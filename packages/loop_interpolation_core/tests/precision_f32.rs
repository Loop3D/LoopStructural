//! Single precision tests (T046, story 5, FR-013, FR-022, SC-005).
//! F32 must use at most 60 percent of the F64 memory, must agree with F64,
//! and must report its residual and its convergence honestly.

mod common;

use loop_interpolation_core::{Field, FieldData, Precision, Problem, Settings, SolveReport};
use std::sync::Mutex;

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
}

/// Build the R1-type problem: value and gradient data from a plane field.
/// Each call builds the same problem (the seed is fixed).
fn build(n: usize, n_value: usize, n_normal: usize) -> Problem {
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
    p
}

fn settings(precision: Precision, tol: f64) -> Settings {
    Settings {
        precision,
        tol,
        max_iter: 20_000,
        ..Settings::default()
    }
}

/// Result of one measured solve.
struct Run {
    peak: usize,
    field: Field,
    report: SolveReport,
}

/// Build and solve. The peak is above the level at the start. It includes
/// the build of the Problem, for F32 and for F64 in the same way.
fn measure(n: usize, precision: Precision, tol: f64) -> Run {
    common::reset_peak();
    let start = common::current();
    let p = build(n, 500, 125);
    let (field, report) = p.solve(settings(precision, tol)).unwrap();
    let peak = common::peak().saturating_sub(start);
    Run {
        peak,
        field,
        report,
    }
}

fn range_of(v: &[f64]) -> f64 {
    let lo = v.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    hi - lo
}

/// Relative difference of two numbers.
fn rel_diff(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(f64::MIN_POSITIVE)
}

/// Tolerance for F32. F32 cannot reach the F64 tolerance.
const TOL_F32: f64 = 1e-5;
const TOL_F64: f64 = 1e-8;

/// (1) The peak of F32 is at most 0.6 times the peak of F64.
#[test]
fn f32_peak_is_at_most_60_percent_of_f64_peak() {
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let r64 = measure(50, Precision::F64, TOL_F64);
    let r32 = measure(50, Precision::F32, TOL_F32);
    println!(
        "peak f64={} f32={} report peak f64={} f32={}",
        r64.peak, r32.peak, r64.report.peak_bytes, r32.report.peak_bytes
    );
    assert!(
        r32.peak as f64 <= 0.6 * r64.peak as f64,
        "F32 peak {}, F64 peak {}",
        r32.peak,
        r64.peak
    );
}

/// (2) The F32 field differs from the F64 field by at most 1e-3 of the range.
#[test]
fn f32_field_agrees_with_f64_field() {
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let r64 = measure(50, Precision::F64, TOL_F64);
    let r32 = measure(50, Precision::F32, TOL_F32);
    let v64 = r64.field.values_f64();
    let v32 = r32.field.values_f64();
    assert_eq!(v64.len(), v32.len());
    let range = range_of(&v64);
    let max_diff = v64
        .iter()
        .zip(v32.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(
        max_diff <= 1e-3 * range,
        "largest difference {max_diff}, range {range}"
    );
}

/// (3) The report residual matches the residual that `residual_f64` computes,
/// in both modes.
#[test]
fn report_residual_matches_residual_f64() {
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    for (precision, tol) in [(Precision::F64, TOL_F64), (Precision::F32, TOL_F32)] {
        // `solve` consumes the problem. Build a second copy for the check.
        let kept = build(50, 500, 125);
        let (field, report) = build(50, 500, 125).solve(settings(precision, tol)).unwrap();
        let own = kept.residual_f64(&field);
        let rel = rel_diff(report.residual, own);
        assert!(
            rel <= 1e-6,
            "{precision:?}: report residual {}, residual_f64 {own}, relative difference {rel}",
            report.residual
        );
    }
}

/// (4) F32 with tol = 1e-12 does not converge, and the report says so.
#[test]
fn f32_with_tiny_tolerance_reports_not_converged() {
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let (_field, report) = build(20, 100, 25)
        .solve(Settings {
            precision: Precision::F32,
            tol: 1e-12,
            max_iter: 2000,
            ..Settings::default()
        })
        .unwrap();
    assert!(
        !report.converged,
        "F32 cannot reach 1e-12, residual {}",
        report.residual
    );
    assert!(
        report.residual > 1e-12,
        "residual {} is below the tolerance but converged is false",
        report.residual
    );
}

/// (5) The Field in F32 mode holds a Vec<f32>. The F64 field holds a Vec<f64>.
#[test]
fn f32_field_holds_f32_values() {
    let _guard = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let (field, _) = build(20, 100, 25)
        .solve(settings(Precision::F32, TOL_F32))
        .unwrap();
    assert!(matches!(field.data, FieldData::F32(_)));
    let (field, _) = build(20, 100, 25)
        .solve(settings(Precision::F64, TOL_F64))
        .unwrap();
    assert!(matches!(field.data, FieldData::F64(_)));
}
