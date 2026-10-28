//! The per-class stencil (the default operator on grids with 7 or more
//! nodes on every axis) must equal the unfused reference operator to about
//! 1e-15, on random vectors and on several grid shapes.

use loop_interpolation_core::{
    BorderRegularisationWeights, FiniteDifferenceInterpolator, RegularisationWeights,
    StructuredGrid3D,
};

/// Small deterministic generator (LCG). It returns a number in [-0.5, 0.5).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
    }
}

fn interpolator(
    shape: [usize; 3],
    weights: RegularisationWeights,
    border: BorderRegularisationWeights,
    ridge: f64,
) -> FiniteDifferenceInterpolator {
    let grid = StructuredGrid3D::new([0.0; 3], shape, [1.0; 3]);
    let mut interp = FiniteDifferenceInterpolator::new(grid)
        .with_regularisation(weights)
        .with_border_regularisation(border);
    interp.ridge_factor = ridge;
    // Data rows at a corner, on a face and inside, so the data term is
    // checked with the sweep.
    interp.add_value_constraint([0.0, 0.0, 0.0], 1.0, 2.0);
    interp.add_value_constraint([2.5, 1.5, 3.0], -0.4, 1.0);
    interp.add_gradient_constraint([1.0, 2.0, 0.5], [1.0, 0.3, -0.2], 1.5);
    interp
}

/// Largest difference, relative to the largest entry of the reference.
fn relative_difference(a: &[f64], b: &[f64]) -> f64 {
    let scale = b.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1e-300);
    let diff = a
        .iter()
        .zip(b.iter())
        .fold(0.0f64, |m, (x, y)| m.max((x - y).abs()));
    diff / scale
}

fn check(shape: [usize; 3], weights: RegularisationWeights, border: BorderRegularisationWeights) {
    for ridge in [0.0, 1e-8, 0.3] {
        let interp = interpolator(shape, weights, border, ridge);
        // The class path needs no scratch. The fused path does.
        let class_path = shape.iter().all(|&n| n >= 7);
        assert_eq!(interp.scratch_len() == 0, class_path, "shape {shape:?}");
        let mut rng = Lcg(shape[0] as u64 * 131 + shape[1] as u64 * 17 + shape[2] as u64);
        for _ in 0..3 {
            let x: Vec<f64> = (0..interp.dof()).map(|_| rng.next()).collect();
            let class = interp.apply_normal_operator(&x);
            let reference = interp.apply_normal_operator_unfused(&x);
            let rel = relative_difference(&class, &reference);
            println!("shape {shape:?} ridge {ridge}: relative difference {rel:e}");
            assert!(
                rel <= 2e-15,
                "shape {shape:?} ridge {ridge} weights {weights:?} border {border:?}: \
                 relative difference {rel:e}"
            );
        }
    }
}

const SHAPES: [[usize; 3]; 9] = [
    [7, 7, 7],
    [7, 8, 9],
    [9, 7, 16],
    [16, 16, 16],
    [12, 10, 7],
    [8, 15, 11],
    [30, 9, 7],
    [61, 61, 31],
    [40, 33, 7],
];

#[test]
fn class_stencil_equals_the_reference_with_default_weights() {
    for shape in SHAPES {
        check(
            shape,
            RegularisationWeights::default(),
            BorderRegularisationWeights::default(),
        );
    }
}

#[test]
fn class_stencil_equals_the_reference_with_other_weights() {
    let weights = RegularisationWeights {
        dxx: 2.0,
        dyy: 0.5,
        dzz: 1.5,
        dxy: 0.0,
        dxz: 3.0,
        dyz: 0.25,
    };
    let border = BorderRegularisationWeights {
        dx: 0.0,
        dy: 2.0,
        dz: 0.7,
    };
    for shape in SHAPES {
        check(shape, weights, border);
    }
}

#[test]
fn class_stencil_equals_the_reference_without_curvature_or_border() {
    let no_curvature = RegularisationWeights {
        dxx: 0.0,
        dyy: 0.0,
        dzz: 0.0,
        dxy: 0.0,
        dxz: 0.0,
        dyz: 0.0,
    };
    let no_border = BorderRegularisationWeights {
        dx: 0.0,
        dy: 0.0,
        dz: 0.0,
    };
    for shape in [[7, 7, 7], [11, 8, 13]] {
        check(shape, no_curvature, BorderRegularisationWeights::default());
        check(shape, RegularisationWeights::default(), no_border);
    }
}

/// Below 7 nodes on an axis the operator keeps the fused path. It must
/// still equal the reference.
#[test]
fn small_grids_keep_the_fused_path() {
    for shape in [[6, 9, 9], [9, 5, 9], [4, 4, 4], [6, 6, 6]] {
        check(
            shape,
            RegularisationWeights::default(),
            BorderRegularisationWeights::default(),
        );
    }
}

/// The class path is symmetric: <N x, y> = <x, N y>.
#[test]
fn class_stencil_is_symmetric() {
    let interp = interpolator(
        [13, 9, 11],
        RegularisationWeights::default(),
        BorderRegularisationWeights::default(),
        1e-8,
    );
    let mut rng = Lcg(7);
    let x: Vec<f64> = (0..interp.dof()).map(|_| rng.next()).collect();
    let y: Vec<f64> = (0..interp.dof()).map(|_| rng.next()).collect();
    let nx = interp.apply_normal_operator(&x);
    let ny = interp.apply_normal_operator(&y);
    let a: f64 = nx.iter().zip(y.iter()).map(|(p, q)| p * q).sum();
    let b: f64 = x.iter().zip(ny.iter()).map(|(p, q)| p * q).sum();
    assert!((a - b).abs() <= 1e-12 * a.abs().max(b.abs()), "{a} vs {b}");
}
