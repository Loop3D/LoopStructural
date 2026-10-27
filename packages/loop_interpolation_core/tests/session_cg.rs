//! Tests of the CG solve session (T019).
//! The session must equal the one-shot solve, give progress, be able to
//! stop early, and free its memory when it is dropped.

mod common;

use std::sync::Mutex;

use loop_interpolation_core::{
    FiniteDifferenceInterpolator, Problem, Settings, SolveSession, Stage, Step, StructuredGrid3D,
};

/// The tests that read the allocator counters must not run at the same time.
/// Every test takes this lock. The allocator counters are global, so a test that
/// runs at the same time would change the memory level of another test.
static MEM_LOCK: Mutex<()> = Mutex::new(());

const N: usize = 30;

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

/// Fixed data: flat arrays of points, values and normals.
struct Data {
    vp: Vec<f64>,
    vv: Vec<f64>,
    vw: Vec<f64>,
    np: Vec<f64>,
    nv: Vec<f64>,
    nw: Vec<f64>,
}

fn data(n: usize, n_value: usize, n_normal: usize) -> Data {
    let extent = (n - 1) as f64;
    let mut rng = Lcg(2024);
    let mut d = Data {
        vp: vec![],
        vv: vec![],
        vw: vec![],
        np: vec![],
        nv: vec![],
        nw: vec![],
    };
    for _ in 0..n_value {
        let p = [
            1.0 + rng.next() * (extent - 2.0),
            1.0 + rng.next() * (extent - 2.0),
            1.0 + rng.next() * (extent - 2.0),
        ];
        d.vv.push(p[0] + 0.5 * p[1] + 0.25 * p[2] + 0.1 * (p[0] * 0.3).sin());
        d.vp.extend_from_slice(&p);
        d.vw.push(1.0);
    }
    for _ in 0..n_normal {
        let p = [
            1.0 + rng.next() * (extent - 2.0),
            1.0 + rng.next() * (extent - 2.0),
            1.0 + rng.next() * (extent - 2.0),
        ];
        d.np.extend_from_slice(&p);
        d.nv.extend_from_slice(&[1.0, 0.5, 0.25]);
        d.nw.push(1.0);
    }
    d
}

/// Build the problem with the flat API.
fn problem(n: usize, d: &Data) -> Problem {
    let mut p = Problem::new([0.0; 3], [n, n, n], [1.0; 3]).unwrap();
    p.add_value(&d.vp, &d.vv, &d.vw).unwrap();
    p.add_normal(&d.np, &d.nv, &d.nw).unwrap();
    p
}

fn settings(tol: f64, max_iter: usize) -> Settings {
    Settings {
        tol,
        max_iter,
        ..Settings::default()
    }
}

/// Reference: the ported one-shot solve with the same data.
fn reference(n: usize, d: &Data) -> loop_interpolation_core::CgResult {
    let grid = StructuredGrid3D::new([0.0; 3], [n, n, n], [1.0; 3]);
    let mut interp = FiniteDifferenceInterpolator::new(grid);
    for i in 0..d.vv.len() {
        interp.add_value_constraint(
            [d.vp[3 * i], d.vp[3 * i + 1], d.vp[3 * i + 2]],
            d.vv[i],
            d.vw[i],
        );
    }
    for i in 0..d.nw.len() {
        interp.add_norm_constraint(
            [d.np[3 * i], d.np[3 * i + 1], d.np[3 * i + 2]],
            [d.nv[3 * i], d.nv[3 * i + 1], d.nv[3 * i + 2]],
            d.nw[i],
        );
    }
    interp.solve(5000, 1e-10)
}

/// Run `step(k)` until `Done`. Return the report.
fn run_to_done(s: &mut SolveSession, k: usize) -> loop_interpolation_core::SolveReport {
    for _ in 0..100_000 {
        if let Step::Done(r) = s.step(k).unwrap() {
            return r;
        }
    }
    panic!("the session does not finish");
}

/// (1) The session field equals the one-shot field.
#[test]
fn session_equals_one_shot() {
    let _g = MEM_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let d = data(N, 300, 60);
    let want = reference(N, &d);
    assert!(want.converged);
    let mut s = problem(N, &d).start(settings(1e-10, 5000)).unwrap();
    run_to_done(&mut s, 7);
    let (field, report) = s.finish().unwrap();
    assert!(report.converged);
    let got = field.values_f64();
    assert_eq!(got.len(), want.x.len());
    let lo = want.x.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = want.x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let range = hi - lo;
    let max_diff = got
        .iter()
        .zip(want.x.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(
        max_diff <= 1e-8 * range,
        "max diff {max_diff} is over 1e-8 of range {range}"
    );
}

/// (2) Progress: stage Solve, growing iteration, falling residual.
#[test]
fn progress_is_reported() {
    let _g = MEM_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let d = data(N, 300, 60);
    let mut s = problem(N, &d).start(settings(1e-10, 5000)).unwrap();
    let mut seen = Vec::new();
    for _ in 0..5 {
        match s.step(7).unwrap() {
            Step::Running(p) => {
                assert!(matches!(p.stage, Stage::Solve));
                seen.push((p.iteration, p.residual));
            }
            Step::Done(_) => panic!("the solve ended in 35 iterations"),
        }
    }
    assert_eq!(seen.len(), 5);
    for w in seen.windows(2) {
        assert!(w[1].0 > w[0].0, "iteration does not grow: {seen:?}");
    }
    assert!(
        seen[4].1 < seen[0].1,
        "residual does not fall over 5 steps: {seen:?}"
    );
}

/// Warm up lazy allocations (thread pool) so that they do not count.
fn warm_up() {
    let d = data(8, 20, 5);
    let mut s = problem(8, &d).start(settings(1e-8, 100)).unwrap();
    let _ = s.step(3);
    drop(s);
}

/// (3) Drop after 3 steps: memory returns to the level before `start`.
#[test]
fn drop_after_steps_frees_memory() {
    let _g = MEM_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let d = data(N, 300, 60);
    warm_up();
    let p = problem(N, &d);
    let before = common::current();
    let mut s = p.start(settings(1e-10, 5000)).unwrap();
    for _ in 0..3 {
        let _ = s.step(7).unwrap();
    }
    drop(s);
    let after = common::current();
    // The session owns the problem, so the level may be lower. It must not be higher.
    assert!(
        after <= before + 64 * 1024,
        "leak: before {before}, after {after}"
    );
}

/// (4) Drop right after `start` (cancel during setup): same check.
#[test]
fn drop_after_start_frees_memory() {
    let _g = MEM_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let d = data(N, 300, 60);
    warm_up();
    let p = problem(N, &d);
    let before = common::current();
    let s = p.start(settings(1e-10, 5000)).unwrap();
    drop(s);
    let after = common::current();
    assert!(
        after <= before + 64 * 1024,
        "leak: before {before}, after {after}"
    );
}

/// (5) `finish()` after 3 steps gives a field and `converged == false`.
/// `best_so_far()` gives the same.
#[test]
fn finish_early_is_not_converged() {
    let _g = MEM_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let d = data(N, 300, 60);
    let mut s = problem(N, &d).start(settings(1e-10, 5000)).unwrap();
    for _ in 0..3 {
        assert!(matches!(s.step(7).unwrap(), Step::Running(_)));
    }
    let (best, best_report) = s.best_so_far().unwrap();
    assert!(!best_report.converged);
    let (field, report) = s.finish().unwrap();
    assert!(!report.converged);
    assert_eq!(report.iterations, best_report.iterations);
    let a = field.values_f64();
    let b = best.values_f64();
    assert_eq!(a.len(), N * N * N);
    assert_eq!(a, b);
}

/// (6) `max_iter = 3` gives `Done` with `converged == false` and a residual.
#[test]
fn max_iter_stops_not_converged() {
    let _g = MEM_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let d = data(N, 300, 60);
    let mut s = problem(N, &d).start(settings(1e-10, 3)).unwrap();
    let r = run_to_done(&mut s, 7);
    assert!(!r.converged);
    assert!(r.iterations <= 3);
    assert!(r.residual.is_finite() && r.residual > 0.0);
}

/// (7) The report has `peak_bytes > 0` and the true count of iterations.
#[test]
fn report_has_peak_and_iterations() {
    let _g = MEM_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let d = data(N, 300, 60);
    let want = reference(N, &d);
    let mut s = problem(N, &d).start(settings(1e-10, 5000)).unwrap();
    let r = run_to_done(&mut s, 7);
    assert!(r.peak_bytes > 0);
    assert!(r.iterations > 0);
    // The same algorithm gives about the same count as the one-shot solve.
    let diff = (r.iterations as i64 - want.iterations as i64).abs();
    assert!(
        diff <= 5,
        "iterations {} and one-shot {}",
        r.iterations,
        want.iterations
    );
    let (_f, r2) = s.finish().unwrap();
    assert_eq!(r2.iterations, r.iterations);
}
