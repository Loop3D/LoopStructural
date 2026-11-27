//! Gradient of a faulted field.

use super::nearest::KdTree;
use super::{FaultChain, GridField, Side};

type Tet = [[f64; 3]; 4];

/// Return the gradient of `field` at `points` after restore through `order`.
/// The second result is the number of points refined at least once.
#[allow(clippy::too_many_arguments)]
pub fn feature_gradient(
    chain: &FaultChain,
    field: &GridField,
    points: &[[f64; 3]],
    order: &[usize],
    sides: &[Side],
    restored: &[bool],
    jumping: &[bool],
    scale: f64,
) -> (Vec<[f64; 3]>, usize) {
    feature_gradient_threads(
        chain, field, points, order, sides, restored, jumping, scale, 1,
    )
}

/// Same as `feature_gradient`. Use `threads` threads (0 is all cores) for the
/// work on each corner with the feature `parallel`; otherwise `threads` has no effect.
#[allow(clippy::too_many_arguments)]
pub fn feature_gradient_threads(
    chain: &FaultChain,
    field: &GridField,
    points: &[[f64; 3]],
    order: &[usize],
    sides: &[Side],
    restored: &[bool],
    jumping: &[bool],
    scale: f64,
    threads: usize,
) -> (Vec<[f64; 3]>, usize) {
    let mut tets: Vec<Tet> = points.iter().map(|&p| tetrahedron(p, scale)).collect();
    let mut active: Vec<usize> = (0..points.len()).collect();
    active = straddling(chain, jumping, sides, &tets, &active, threads);
    let count = active.len();
    let mut scale = scale;
    for _ in 0..10 {
        if active.is_empty() {
            break;
        }
        scale *= 0.5;
        for &i in &active {
            tets[i] = tetrahedron(points[i], scale);
        }
        active = straddling(chain, jumping, sides, &tets, &active, threads);
    }
    let mut bad = vec![false; points.len()];
    for &i in &active {
        bad[i] = true;
    }
    let grads = map(&tets, threads, |i, t| {
        if bad[i] {
            return [f64::NAN; 3];
        }
        let v: [f64; 4] = std::array::from_fn(|c| {
            field.value(chain.restore(order, t[c], false, sides, restored))
        });
        gradient_from_tetrahedron(t, &v)
    });
    (grads, count)
}

/// Run `f(index, item)` for each item, in order. Parallel with the feature `parallel`.
#[cfg(feature = "parallel")]
fn map<T: Sync, R: Send>(
    items: &[T],
    threads: usize,
    f: impl Fn(usize, &T) -> R + Sync + Send,
) -> Vec<R> {
    use rayon::prelude::*;
    super::pool::pool(threads)
        .install(|| items.par_iter().enumerate().map(|(i, t)| f(i, t)).collect())
}

#[cfg(not(feature = "parallel"))]
fn map<T, R>(items: &[T], _threads: usize, f: impl Fn(usize, &T) -> R) -> Vec<R> {
    items.iter().enumerate().map(|(i, t)| f(i, t)).collect()
}

/// Return the regular tetrahedron of size `scale` around `p`.
fn tetrahedron(p: [f64; 3], scale: f64) -> Tet {
    let a = (8.0f64 / 9.0).sqrt();
    let b = (2.0f64 / 9.0).sqrt();
    let c = (2.0f64 / 3.0).sqrt();
    let r = [
        [a, 0.0, -1.0 / 3.0],
        [-b, c, -1.0 / 3.0],
        [-b, -c, -1.0 / 3.0],
        [0.0, 0.0, 1.0],
    ];
    // The regular tetrahedron is scaled first, then moved to the centre.
    std::array::from_fn(|i| std::array::from_fn(|k| p[k] + r[i][k] * scale))
}

/// Return the points (from `active`) whose tetrahedron has a fault between its corners.
fn straddling(
    chain: &FaultChain,
    jumping: &[bool],
    sides: &[Side],
    tets: &[Tet],
    active: &[usize],
    threads: usize,
) -> Vec<usize> {
    let corners: Vec<[f64; 3]> = active.iter().flat_map(|&i| tets[i]).collect();
    let all_false = vec![false; chain.faults.len()];
    let mut out = vec![false; active.len()];
    for j in (0..jumping.len()).filter(|&j| jumping[j]) {
        let cutters = &chain.faults[j].cutters;
        let q = map(&corners, threads, |_, &c| {
            chain.restore(cutters, c, false, sides, &all_false)
        });
        let mut v = map(&corners, threads, |_, &c| {
            chain.frame_value(j, 0, c, sides, &all_false)
        });
        fill_nearest(&mut v, &q);
        for (k, o) in out.iter_mut().enumerate() {
            let s = &v[4 * k..4 * k + 4];
            *o |= !(s.iter().all(|&x| x > 0.0) || s.iter().all(|&x| x < 0.0));
        }
    }
    active
        .iter()
        .zip(out)
        .filter(|(_, o)| *o)
        .map(|(&i, _)| i)
        .collect()
}

/// Give each NaN value the value of the nearest known point of the same batch.
fn fill_nearest(v: &mut [f64], q: &[[f64; 3]]) {
    let known: Vec<usize> = (0..v.len()).filter(|&i| !v[i].is_nan()).collect();
    let nan: Vec<usize> = (0..v.len())
        .filter(|&i| v[i].is_nan() && q[i].iter().all(|x| x.is_finite()))
        .collect();
    if known.is_empty() || nan.is_empty() {
        return;
    }
    let tree = KdTree::new(&known.iter().map(|&i| q[i]).collect::<Vec<_>>());
    for i in nan {
        v[i] = v[known[tree.nearest(q[i])]];
    }
}

/// Return the gradient of the linear field with `value` at the corners of `t`.
fn gradient_from_tetrahedron(t: &Tet, value: &[f64; 4]) -> [f64; 3] {
    // The rows of `e` are the edges from corner 0.
    let e: [[f64; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|k| t[i + 1][k] - t[0][k]));
    let inv = invert(&e);
    // G = inv * I, with I the matrix of the corner differences.
    let i_mat = [
        [-1.0, 1.0, 0.0, 0.0],
        [-1.0, 0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0, 1.0],
    ];
    std::array::from_fn(|a| {
        (0..4)
            .map(|c| (0..3).map(|k| inv[a][k] * i_mat[k][c]).sum::<f64>() * value[c])
            .sum()
    })
}

/// Return the inverse of a 3x3 matrix. A singular matrix gives NaN.
fn invert(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let c = |a: usize, b: usize, c: usize, d: usize| m[a][b] * m[c][d] - m[a][d] * m[c][b];
    let adj = [
        [c(1, 1, 2, 2), -c(0, 1, 2, 2), c(0, 1, 1, 2)],
        [-c(1, 0, 2, 2), c(0, 0, 2, 2), -c(0, 0, 1, 2)],
        [c(1, 0, 2, 1), -c(0, 0, 2, 1), c(0, 0, 1, 1)],
    ];
    let det = m[0][0] * adj[0][0] + m[0][1] * adj[1][0] + m[0][2] * adj[2][0];
    if det == 0.0 {
        return [[f64::NAN; 3]; 3];
    }
    adj.map(|r| r.map(|x| x / det))
}
