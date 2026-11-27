//! Tests of the structured grid lookup of `GridField`.

use loop_interpolation_core::fault::grid_field::GridField;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * b.abs().max(1.0)
}

/// Build a grid with node values `f(x, y, z)` at the local node coordinates (x runs fastest).
fn grid(
    origin: [f64; 3],
    nsteps: [usize; 3],
    step: [f64; 3],
    f: impl Fn([f64; 3]) -> f64,
) -> GridField {
    let mut values = Vec::new();
    for k in 0..nsteps[2] {
        for j in 0..nsteps[1] {
            for i in 0..nsteps[0] {
                let n = [i, j, k];
                values.push(f(std::array::from_fn(|a| {
                    origin[a] + step[a] * n[a] as f64
                })));
            }
        }
    }
    GridField {
        origin,
        nsteps,
        step_vector: step,
        local_origin: None,
        local_rotation: None,
        values,
    }
}

#[test]
fn linear_field_is_exact() {
    let g = grid([1.0, -2.0, 0.5], [4, 5, 3], [0.5, 0.8, 1.0], |p| {
        2.0 * p[0] - 3.0 * p[1] + 0.5 * p[2] + 1.0
    });
    let p = [1.7, 0.3, 1.2];
    assert!(close(g.value(p), 2.0 * 1.7 - 0.9 + 0.6 + 1.0));
    let d = g.gradient(p);
    assert!(close(d[0], 2.0) && close(d[1], -3.0) && close(d[2], 0.5));
}

#[test]
fn rotated_local_frame() {
    let (c, s) = (0.3_f64.cos(), 0.3_f64.sin());
    let r = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
    let o = [10.0, 20.0, 30.0];
    let mut g = grid([0.0; 3], [5, 5, 5], [1.0; 3], |q| q[0] + 2.0 * q[1] - q[2]);
    g.local_origin = Some(o);
    g.local_rotation = Some(r);
    // World point whose local coordinates are (1.5, 2.5, 3.5).
    let q = [1.5, 2.5, 3.5];
    let p: [f64; 3] = std::array::from_fn(|i| o[i] + (0..3).map(|j| r[j][i] * q[j]).sum::<f64>());
    assert!(close(g.value(p), 1.5 + 5.0 - 3.5));
    // The gradient in world axes is R^T (1, 2, -1).
    let a = [1.0, 2.0, -1.0];
    let d = g.gradient(p);
    for i in 0..3 {
        assert!(close(d[i], (0..3).map(|j| r[j][i] * a[j]).sum::<f64>()));
    }
}

#[test]
fn nan_and_outside_points() {
    let g = grid([0.0; 3], [3, 3, 3], [1.0; 3], |p| p[0]);
    for p in [
        [f64::NAN, 1.0, 1.0],
        [1.0, f64::NAN, 1.0],
        [3.0, 1.0, 1.0],
        [1.0, 1.0, -0.1],
    ] {
        assert!(g.value(p).is_nan());
        // Outside finite points give a zero gradient. NaN points give NaN.
        let is_nan_point = p.iter().any(|v| v.is_nan());
        for v in g.gradient(p) {
            assert!(if is_nan_point { v.is_nan() } else { v == 0.0 });
        }
    }
    // The upper boundary is inside.
    assert!(close(g.value([2.0, 2.0, 2.0]), 2.0));
}

/// Expected values come from this Python run (`.venv/bin/python -I ref.py`):
///
/// ```python
/// import numpy as np
/// from loop_common.supports import StructuredGrid
/// g = StructuredGrid(origin=np.array([1.0, -2.0, 0.5]), nsteps=np.array([5, 4, 6]),
///                    step_vector=np.array([0.7, 1.1, 0.4]))  # nsteps here = cells, nodes [6, 5, 7]
/// n = g.nodes
/// c = np.sin(n[:, 0]) + np.cos(1.3 * n[:, 1]) * n[:, 2]
/// pts = np.array([[1.3, -1.5, 0.9], [2.2, 0.4, 1.7], [3.1, 1.0, 2.2],
///                 [3.4, -0.3, 1.1], [4.5, 1.2, 2.9]])
/// lk = g._lookup_points(pts)
/// val = g._values_from_coefficients(lk[0], lk[1], g.trilinear(lk[2][lk[1]]), c)
/// gr = g._evaluate_gradient_from_lookup(lk, c)
/// ```
#[test]
fn matches_python_grid() {
    let g = grid([1.0, -2.0, 0.5], [6, 5, 7], [0.7, 1.1, 0.4], |p| {
        p[0].sin() + (1.3 * p[1]).cos() * p[2]
    });
    let pts = [
        [1.3, -1.5, 0.9],
        [2.2, 0.4, 1.7],
        [3.1, 1.0, 2.2],
        [3.4, -0.3, 1.1],
        [4.5, 1.2, 2.9],
    ];
    let val = [
        0.644792795517428,
        2.0732095333664233,
        0.43114013843941873,
        0.5364453971641094,
        -1.036275078847976,
    ];
    let gr = [
        [
            0.21456260806367436,
            1.0203058126449638,
            -0.29005219078841216,
        ],
        [-0.4517166141447394, -1.6772996642117397, 0.7690606058743095],
        [
            -0.9334836476800142,
            -2.1706230948622522,
            0.17707248909369477,
        ],
        [-0.9334836476800141, 0.5762382938262829, 0.7044634809407484],
        [
            -0.5223888953176828,
            -2.861275897772968,
            -0.02025688316650953,
        ],
    ];
    for i in 0..5 {
        assert!(close(g.value(pts[i]), val[i]), "value {i}");
        let d = g.gradient(pts[i]);
        for a in 0..3 {
            assert!(close(d[a], gr[i][a]), "gradient {i} axis {a}");
        }
    }
}
