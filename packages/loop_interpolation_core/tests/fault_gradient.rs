//! Tests of the native gradient of a faulted feature on hand-built fault chains.

use loop_interpolation_core::fault::gradient::{feature_gradient, feature_gradient_threads};
use loop_interpolation_core::fault::profile::{Profile, Profile1};
use loop_interpolation_core::fault::{
    BlockSettings, Extension, Fault, FaultChain, GridField, Side,
};

const N: usize = 11;

/// A grid on [-5, 5]^3 with node values `f(x, y, z)`.
fn grid(f: impl Fn(f64, f64, f64) -> f64) -> GridField {
    let mut values = Vec::new();
    for k in 0..N {
        for j in 0..N {
            for i in 0..N {
                values.push(f(i as f64 - 5.0, j as f64 - 5.0, k as f64 - 5.0));
            }
        }
    }
    GridField {
        origin: [-5.0; 3],
        nsteps: [N; 3],
        step_vector: [1.0; 3],
        local_origin: None,
        local_rotation: None,
        values,
    }
}

fn chain() -> FaultChain {
    let gx = Profile1::Composite {
        positive: Box::new(Profile1::Ones),
        negative: Box::new(Profile1::Zeros),
    };
    FaultChain {
        faults: vec![Fault {
            frame: [grid(|x, _, _| x), grid(|_, y, _| y), grid(|_, _, z| z)],
            cutters: vec![],
            regions: [vec![], vec![], vec![]],
            profile: Profile {
                scale: 1.0,
                gx,
                gy: Profile1::Ones,
                gz: Profile1::Ones,
            },
            displacement: 0.7,
            steps: 2,
        }],
        block: BlockSettings {
            side_epsilon: 1e-9,
            extension: Extension::Smooth,
            exact: 0.25,
            limit: 0.5,
            step: 1e-3,
        },
    }
}

fn linear() -> GridField {
    grid(|x, y, z| 2.0 * x - 3.0 * y + 0.5 * z)
}

const S: [Side; 1] = [Side::Free];

#[test]
fn no_fault_equals_field_gradient() {
    let c = chain();
    let f = linear();
    let pts = [[1.0, 0.5, 0.25], [-2.3, 1.1, 3.0]];
    let (g, n) = feature_gradient(&c, &f, &pts, &[], &S, &[false], &[false], 0.5);
    assert_eq!(n, 0);
    for (p, g) in pts.iter().zip(g) {
        let e = f.gradient(*p);
        assert!((0..3).all(|k| (g[k] - e[k]).abs() < 1e-10), "{g:?} {e:?}");
    }
}

#[test]
fn point_on_fault_is_nan_and_far_point_is_not() {
    let c = chain();
    let f = linear();
    let pts = [[0.0, 0.5, 0.25], [3.0, 0.5, 0.25]];
    let (g, n) = feature_gradient(&c, &f, &pts, &[0], &S, &[false], &[true], 0.5);
    assert_eq!(n, 1);
    assert!(g[0].iter().all(|v| v.is_nan()));
    assert!(g[1].iter().all(|v| v.is_finite()));
    // A point near the fault is refined and then fine.
    let (g, n) = feature_gradient(
        &c,
        &f,
        &[[0.1, 0.5, 0.25]],
        &[0],
        &S,
        &[false],
        &[true],
        0.5,
    );
    assert_eq!(n, 1);
    assert!(g[0].iter().all(|v| v.is_finite()));
}

#[test]
fn threads_variant_equals_plain_and_empty() {
    let c = chain();
    let f = linear();
    let pts: Vec<[f64; 3]> = (0..40).map(|i| [i as f64 * 0.1 - 2.0, 0.2, 0.1]).collect();
    let a = feature_gradient(&c, &f, &pts, &[0], &S, &[false], &[true], 0.5);
    for t in [0, 1, 3] {
        let b = feature_gradient_threads(&c, &f, &pts, &[0], &S, &[false], &[true], 0.5, t);
        assert_eq!(a.1, b.1);
        for (x, y) in a.0.iter().zip(&b.0) {
            for k in 0..3 {
                assert!(x[k].to_bits() == y[k].to_bits());
            }
        }
    }
    let e = feature_gradient(&c, &f, &[], &[0], &S, &[false], &[true], 0.5);
    assert!(e.0.is_empty() && e.1 == 0);
}
