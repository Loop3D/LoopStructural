//! Tests of the native fault restore on hand-built fault chains.
//!
//! Not done: a comparison with `FaultSegment.apply_to_points`. The frame grids of a
//! LoopStructural model are too big to copy into a test. T013 and the Python
//! tests cover that comparison.

use loop_model_core::fault::profile::{Profile, Profile1};
use loop_model_core::fault::{
    BlockSettings, Extension, Fault, FaultChain, GridField, Region, Side,
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

fn fault(x0: f64, gx: Profile1, displacement: f64, steps: usize) -> Fault {
    Fault {
        frame: [grid(|x, _, _| x - x0), grid(|_, y, _| y), grid(|_, _, z| z)],
        cutters: vec![],
        regions: [vec![], vec![], vec![]],
        profile: Profile {
            scale: 1.0,
            gx,
            gy: Profile1::Ones,
            gz: Profile1::Ones,
        },
        displacement,
        steps,
    }
}

fn chain(faults: Vec<Fault>) -> FaultChain {
    FaultChain {
        faults,
        block: BlockSettings {
            side_epsilon: 1e-9,
            extension: Extension::Smooth,
            exact: 0.25,
            limit: 0.5,
            step: 1e-3,
        },
    }
}

fn close(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
    (0..3).all(|i| (a[i] - b[i]).abs() <= tol)
}

#[test]
fn constant_displacement_shifts_along_slip() {
    let c = chain(vec![fault(0.0, Profile1::Ones, 0.7, 4)]);
    let s = [Side::Free];
    let r = c.restore(&[0], [1.0, 0.5, 0.25], false, &s, &[false]);
    assert!(close(r, [1.0, 1.2, 0.25], 1e-12), "{r:?}");
}

#[test]
fn reverse_undoes_restore() {
    let c = chain(vec![fault(0.0, Profile1::Ones, 0.7, 3)]);
    let s = [Side::Free];
    let p = [1.0, 0.5, 0.25];
    let r = c.restore(&[0], p, false, &s, &[false]);
    let b = c.restore(&[0], r, true, &s, &[false]);
    assert!(close(b, p, 1e-9), "{b:?}");
}

#[test]
fn nan_in_gives_nan_out() {
    let c = chain(vec![fault(0.0, Profile1::Ones, 0.7, 2)]);
    let r = c.restore(&[0], [f64::NAN, 0.0, 0.0], false, &[Side::Free], &[false]);
    assert!(r[0].is_nan());
    let r = c.restore(&[0], [1.0, f64::NAN, 0.0], false, &[Side::Free], &[false]);
    assert!(r[1].is_nan());
}

#[test]
fn absent_side_moves_nothing() {
    let c = chain(vec![fault(0.0, Profile1::Ones, 0.7, 2)]);
    let p = [1.0, 0.5, 0.25];
    assert_eq!(c.apply(0, p, false, &[Side::Absent], &[false]), p);
}

#[test]
fn footwall_point_does_not_move() {
    let gx = Profile1::Composite {
        positive: Box::new(Profile1::Ones),
        negative: Box::new(Profile1::Zeros),
    };
    let c = chain(vec![fault(0.0, gx, 0.7, 2)]);
    let s = [Side::Free];
    let foot = [-1.0, 0.5, 0.25];
    assert_eq!(c.apply(0, foot, false, &s, &[false]), foot);
    let hang = c.apply(0, [1.0, 0.5, 0.25], false, &s, &[false]);
    assert!(close(hang, [1.0, 1.2, 0.25], 1e-12));
}

#[test]
fn two_faults_equal_two_applies_and_restored_skips() {
    let c = chain(vec![
        fault(0.0, Profile1::Ones, 0.5, 2),
        fault(0.0, Profile1::Ones, 0.25, 1),
    ]);
    let s = [Side::Free, Side::Free];
    let p = [1.0, 0.0, 0.0];
    let a = c.apply(0, p, false, &s, &[false, false]);
    let b = c.apply(1, a, false, &s, &[true, false]);
    let r = c.restore(&[0, 1], p, false, &s, &[false, false]);
    assert!(close(r, b, 1e-15));
    assert!(close(r, [1.0, 0.75, 0.0], 1e-12));
    let skip = c.restore(&[0, 1], p, false, &s, &[true, false]);
    assert!(close(skip, [1.0, 0.25, 0.0], 1e-12));
}

#[test]
fn sign_region_excludes_point() {
    // The region is the positive side of fault 1 (x - 3 > 0); the point has x = 1.
    let mut a = fault(0.0, Profile1::Ones, 0.7, 2);
    let b = fault(3.0, Profile1::Ones, 0.0, 1);
    let region = Region::Sign {
        positive: true,
        fault: 1,
        vector: Some([1.0, 0.0, 0.0]),
        point: Some([3.0, 0.0, 0.0]),
    };
    a.regions[0].push(region);
    let c = chain(vec![a, b]);
    let s = [Side::Free, Side::Free];
    let p = [1.0, 0.5, 0.25];
    assert!(c.frame_value(0, 0, p, &s, &[false; 2]).is_nan());
    assert_eq!(c.restore(&[0], p, false, &s, &[false; 2]), p);
    let q = [4.0, 0.5, 0.25];
    assert!(!c.frame_value(0, 0, q, &s, &[false; 2]).is_nan());
}

#[test]
fn unconformity_region() {
    let mut a = fault(0.0, Profile1::Ones, 0.7, 1);
    a.regions[0].push(Region::Unconformity {
        field: grid(|_, _, z| z),
        faults: vec![],
        value: 1.0,
        sign: true,
    });
    let c = chain(vec![a]);
    let s = [Side::Free];
    assert!(!c.frame_value(0, 0, [1.0, 0.0, 0.5], &s, &[false]).is_nan());
    assert!(c.frame_value(0, 0, [1.0, 0.0, 2.0], &s, &[false]).is_nan());
    assert!(
        c.frame_value(0, 0, [1.0, 0.0, f64::NAN], &s, &[false])
            .is_nan()
    );
}

#[test]
fn batch_equals_single_and_empty() {
    let c = chain(vec![fault(0.0, Profile1::Ones, 0.7, 3)]);
    let s = [Side::Free];
    let pts: Vec<[f64; 3]> = (0..50)
        .map(|i| [i as f64 * 0.1 - 2.5, 0.01 * i as f64, 0.3])
        .chain([[f64::NAN, 0.0, 0.0]])
        .collect();
    for threads in [0, 1, 3] {
        let b = c.restore_batch(&[0], &pts, false, &s, &[false], threads);
        let v = c.frame_value_batch(0, 0, &pts, &s, &[false], threads);
        for (i, p) in pts.iter().enumerate() {
            let r = c.restore(&[0], *p, false, &s, &[false]);
            for k in 0..3 {
                assert!(b[i][k] == r[k] || (b[i][k].is_nan() && r[k].is_nan()));
            }
            let e = c.frame_value(0, 0, *p, &s, &[false]);
            assert!(v[i] == e || (v[i].is_nan() && e.is_nan()));
        }
    }
    assert!(
        c.restore_batch(&[0], &[], false, &s, &[false], 2)
            .is_empty()
    );
    assert!(c.frame_value_batch(0, 0, &[], &s, &[false], 2).is_empty());
}
