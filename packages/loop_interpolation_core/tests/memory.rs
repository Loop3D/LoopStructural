//! Memory tests. They measure the peak heap use of a full solve.

mod common;

use loop_interpolation_core::{FiniteDifferenceInterpolator, StructuredGrid3D};

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

/// Build and solve a problem on an n x n x n grid.
/// Return the peak bytes above the level at the start, and the dof.
fn run(n: usize, n_value: usize, n_normal: usize) -> (usize, usize) {
    common::reset_peak();
    let start = common::current();

    let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [n, n, n], [1.0, 1.0, 1.0]);
    let mut interp = FiniteDifferenceInterpolator::new(grid);
    let dof = interp.dof();
    let extent = (n - 1) as f64;
    let mut rng = Lcg(12345);

    for _ in 0..n_value {
        let p = [
            rng.next() * extent,
            rng.next() * extent,
            rng.next() * extent,
        ];
        // The value is a plane field, so the data are consistent.
        let v = p[0] + 0.5 * p[1] + 0.25 * p[2];
        interp.add_value_constraint(p, v, 1.0);
    }
    for _ in 0..n_normal {
        let p = [
            rng.next() * extent,
            rng.next() * extent,
            rng.next() * extent,
        ];
        interp.add_gradient_constraint(p, [1.0, 0.5, 0.25], 1.0);
    }

    let result = interp.solve(5000, 1e-8);
    println!(
        "iterations={} converged={} residual={}",
        result.iterations, result.converged, result.residual_norm
    );
    drop(result);

    let peak = common::peak().saturating_sub(start);
    (peak, dof)
}

/// Problem R1: 100^3 grid, 2000 value and 500 normal constraints.
#[test]
#[ignore]
fn baseline_r1() {
    let (peak, dof) = run(100, 2000, 500);
    let bytes_per_dof = peak as f64 / dof as f64;
    println!("bytes_per_dof={}", bytes_per_dof);
    assert!(bytes_per_dof > 0.0);
}

/// Small problem on a 30^3 grid. It only prints.
#[test]
fn baseline_small() {
    let (peak, dof) = run(30, 200, 50);
    println!("small bytes_per_dof={}", peak as f64 / dof as f64);
}
