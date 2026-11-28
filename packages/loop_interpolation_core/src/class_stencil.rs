//! Per-class regularisation kernels.
//!
//! The regularisation part of the normal operator is
//! `sum_family w^2 S^T S + sum_border w^2 B^T B`. Its row at node `p` is a
//! stencil of radius 2. The stencil depends only on which rows of `S` and
//! `B` touch `p`:
//!
//! * A curvature row at `r` is active when `m <= r[a] <= n[a] - 1 - m` on
//!   every axis `a` (the margin `m` is 1). The test is separable by axis.
//! * A border row lives on one layer of one axis, at every position on the
//!   other two axes. Only the coordinate on its own axis matters.
//!
//! So the row at `p` depends only on the class of each coordinate. There
//! are 7 classes per axis: the first, second and third node, the interior,
//! and the mirrored three at the other end. That gives 7^3 = 343 kernels.
//! Each kernel is a sum of products of the family taps and weights, so it
//! is built directly from the stencils (not by probing the operator).
//!
//! The kernels replace the fused-interior pass, the shell correction and
//! the border passes, and their scratch buffers. They need at least
//! [`MIN_AXIS_NODES`] nodes on every axis, so that the interior class is
//! not empty. They do not cover directional regularisation, whose weights
//! change from node to node.

use crate::PAR_MIN_LEN;
use crate::grid::StructuredGrid3D;
use crate::par::*;
use crate::scalar::Scalar;
use crate::stencil::{EdgeStencil, Stencil};

/// Classes on one axis.
const N_CLASS: usize = 7;
/// The interior class.
const INTERIOR: usize = 3;
/// Least number of nodes on an axis for the class kernels.
pub(crate) const MIN_AXIS_NODES: usize = 7;
/// Taps of the dense 5 x 5 x 5 work array (offsets -2..=2 on each axis).
const DENSE: usize = 125;
/// Nodes in one block of an interior run.
const BLOCK: usize = 64;
/// Upper limit of [`ClassKernels::bytes`] (83,704 bytes with all weights
/// on; a zero weight only removes taps). The memory estimate uses it.
pub(crate) const MAX_BYTES: u64 = 86_016;

/// Class of coordinate `p` on an axis with `n` nodes (`n >= MIN_AXIS_NODES`).
#[inline]
fn axis_class(p: usize, n: usize) -> usize {
    if p < INTERIOR {
        p
    } else if p + INTERIOR >= n {
        N_CLASS - (n - p)
    } else {
        INTERIOR
    }
}

/// A coordinate of class `c` on an axis with `n` nodes.
fn class_coordinate(c: usize, n: usize) -> usize {
    if c <= INTERIOR { c } else { n - (N_CLASS - c) }
}

/// Index of an offset in the dense work array.
#[inline]
fn dense_index(e: [i32; 3]) -> usize {
    (e[0] + 2) as usize + 5 * ((e[1] + 2) as usize + 5 * (e[2] + 2) as usize)
}

/// The 343 class kernels of one grid shape and one set of weights.
pub(crate) struct ClassKernels {
    nsteps: [usize; 3],
    /// Kernel of class `(cx, cy, cz)` at `cx + 7 * (cy + 7 * cz)`.
    /// A tap is (flat index offset, coefficient).
    kernels: Vec<Vec<(i64, f64)>>,
}

impl ClassKernels {
    /// Build the kernels. `None` if the grid or the stencils do not fit the
    /// class model: an axis with fewer than `MIN_AXIS_NODES` nodes, or a tap
    /// farther than one step.
    pub(crate) fn build(
        grid: &StructuredGrid3D,
        stencils: &[(Stencil, f64)],
        borders: &[(EdgeStencil, f64)],
    ) -> Option<Self> {
        let n = grid.nsteps;
        if n.iter().any(|&a| a < MIN_AXIS_NODES) {
            return None;
        }
        let near = |d: i32| (-1..=1).contains(&d);
        if stencils.iter().any(|(s, _)| {
            s.offsets
                .iter()
                .any(|d| !(near(d.0) && near(d.1) && near(d.2)))
        }) || borders
            .iter()
            .any(|(b, _)| b.taps.iter().any(|(o, _)| !near(*o)))
        {
            return None;
        }

        let nx = n[0] as i64;
        let nxy = nx * n[1] as i64;
        let mut kernels = Vec::with_capacity(N_CLASS * N_CLASS * N_CLASS);
        for cz in 0..N_CLASS {
            for cy in 0..N_CLASS {
                for cx in 0..N_CLASS {
                    let p = [
                        class_coordinate(cx, n[0]) as i64,
                        class_coordinate(cy, n[1]) as i64,
                        class_coordinate(cz, n[2]) as i64,
                    ];
                    let (sum, abs_sum) = Self::dense_kernel(n, p, stencils, borders);
                    // A tap is kept if it is not the round-off of terms that
                    // cancel. A true zero has no terms at all.
                    let mut taps: Vec<(i64, f64)> = (0..DENSE)
                        .filter(|&t| abs_sum[t] > 0.0 && sum[t].abs() > 1e-13 * abs_sum[t])
                        .map(|t| {
                            let e = [
                                (t % 5) as i64 - 2,
                                ((t / 5) % 5) as i64 - 2,
                                (t / 25) as i64 - 2,
                            ];
                            (e[0] + nx * e[1] + nxy * e[2], sum[t])
                        })
                        .collect();
                    // Keep the table at its exact size (see `MAX_BYTES`).
                    taps.shrink_to_fit();
                    kernels.push(taps);
                }
            }
        }
        Some(Self { nsteps: n, kernels })
    }

    /// The kernel at node `p` as a dense array, with the sum of the
    /// absolute values of the terms of each tap.
    fn dense_kernel(
        n: [usize; 3],
        p: [i64; 3],
        stencils: &[(Stencil, f64)],
        borders: &[(EdgeStencil, f64)],
    ) -> ([f64; DENSE], [f64; DENSE]) {
        let mut sum = [0.0f64; DENSE];
        let mut abs_sum = [0.0f64; DENSE];
        let mut add = |e: [i32; 3], v: f64| {
            let t = dense_index(e);
            sum[t] += v;
            abs_sum[t] += v.abs();
        };
        // Curvature: p is tap d1 of the row at r = p - d1. That row reads
        // x[r + d2] = x[p + (d2 - d1)].
        for (stencil, weight) in stencils {
            if *weight == 0.0 {
                continue;
            }
            let w2 = weight * weight;
            for (d1, c1) in stencil.offsets.iter().zip(stencil.coeffs.iter()) {
                let d1 = [d1.0, d1.1, d1.2];
                let active = (0..3).all(|a| {
                    let r = p[a] - d1[a] as i64;
                    let m = stencil.margin[a] as i64;
                    r >= m && r <= n[a] as i64 - 1 - m
                });
                if !active {
                    continue;
                }
                for (d2, c2) in stencil.offsets.iter().zip(stencil.coeffs.iter()) {
                    let e = [d2.0 - d1[0], d2.1 - d1[1], d2.2 - d1[2]];
                    add(e, w2 * c1 * c2);
                }
            }
        }
        // Border: a row on `layer` of `axis` reads only along that axis.
        for (border, weight) in borders {
            if *weight == 0.0 {
                continue;
            }
            let w2 = weight * weight;
            let a = border.axis;
            for (o1, t1) in &border.taps {
                if p[a] - *o1 as i64 != border.layer as i64 {
                    continue;
                }
                for (o2, t2) in &border.taps {
                    let mut e = [0i32; 3];
                    e[a] = o2 - o1;
                    add(e, w2 * t1 * t2);
                }
            }
        }
        (sum, abs_sum)
    }

    /// Bytes held by the kernels.
    pub(crate) fn bytes(&self) -> usize {
        self.kernels
            .iter()
            .map(|k| k.capacity() * std::mem::size_of::<(i64, f64)>())
            .sum::<usize>()
            + self.kernels.capacity() * std::mem::size_of::<Vec<(i64, f64)>>()
    }

    /// Write `out = (sum_family w^2 S^T S + sum_border w^2 B^T B + r2 I) x`
    /// in one sweep. Each x-line has its edge nodes, which use their class
    /// kernel, and one run of interior nodes that share one kernel. The run
    /// is done as contiguous multiply-add passes, one for each tap. All sums
    /// are in `f64`. `out` is overwritten.
    pub(crate) fn apply_into<S: Scalar>(&self, x: &[S], r2: f64, out: &mut [S]) {
        let [nx, ny, nz] = self.nsteps;
        let min_lines = (PAR_MIN_LEN / nx).max(1);
        out.par_chunks_mut(nx)
            .enumerate()
            .with_min_len(min_lines)
            .for_each(|(line, o)| {
                let cy = axis_class(line % ny, ny);
                let cz = axis_class(line / ny, nz);
                let row = N_CLASS * (cy + N_CLASS * cz);
                let base = line * nx;
                for i in (0..INTERIOR).chain(nx - INTERIOR..nx) {
                    let taps = &self.kernels[axis_class(i, nx) + row];
                    o[i] = node(taps, x, r2, base + i);
                }
                let taps = &self.kernels[INTERIOR + row];
                run(taps, x, r2, base, INTERIOR, nx - INTERIOR, o);
            });
    }
}

/// One node: `r2 x[p] + sum_tap c x[p + offset]`.
#[inline]
fn node<S: Scalar>(taps: &[(i64, f64)], x: &[S], r2: f64, p: usize) -> S {
    let mut acc = r2 * x[p].to_f64();
    for &(off, c) in taps {
        acc += c * x[(p as i64 + off) as usize].to_f64();
    }
    S::from_f64(acc)
}

/// Nodes `start..end` of the line that starts at flat index `base`, which
/// all use `taps`. Blocks of `BLOCK` nodes are summed tap by tap into an
/// `f64` buffer on the stack, so each pass is a contiguous multiply-add.
fn run<S: Scalar>(
    taps: &[(i64, f64)],
    x: &[S],
    r2: f64,
    base: usize,
    start: usize,
    end: usize,
    o: &mut [S],
) {
    let mut buf = [0.0f64; BLOCK];
    let mut s = start;
    while s < end {
        let len = (end - s).min(BLOCK);
        let acc = &mut buf[..len];
        let p0 = base + s;
        for (a, xi) in acc.iter_mut().zip(&x[p0..p0 + len]) {
            *a = r2 * xi.to_f64();
        }
        for &(off, c) in taps {
            let q0 = (p0 as i64 + off) as usize;
            for (a, xi) in acc.iter_mut().zip(&x[q0..q0 + len]) {
                *a += c * xi.to_f64();
            }
        }
        for (oi, a) in o[s..s + len].iter_mut().zip(acc.iter()) {
            *oi = S::from_f64(*a);
        }
        s += len;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kernel_table_size() {
        use crate::stencil::{dxx, dxy, dxz, dyy, dyz, dzz, edge_backward, edge_forward};
        let grid = StructuredGrid3D::new([0.0; 3], [10, 10, 10], [1.0; 3]);
        let stencils = [
            (dxx(), 1.0),
            (dyy(), 1.0),
            (dzz(), 1.0),
            (dxy(), 1.0),
            (dxz(), 1.0),
            (dyz(), 1.0),
        ];
        let mut borders = Vec::new();
        for axis in 0..3 {
            borders.push((edge_forward(axis, &grid.nsteps).unwrap(), 1.0));
            borders.push((edge_backward(axis, &grid.nsteps).unwrap(), 1.0));
        }
        let k = ClassKernels::build(&grid, &stencils, &borders).unwrap();
        let taps: usize = k.kernels.iter().map(|v| v.len()).sum();
        let interior = k.kernels[INTERIOR + N_CLASS * (INTERIOR + N_CLASS * INTERIOR)].len();
        println!("bytes={} taps={taps} interior taps={interior}", k.bytes());
        assert_eq!(interior, 25);
        assert!(k.bytes() as u64 <= MAX_BYTES, "bytes {}", k.bytes());
    }

    #[test]
    fn classes_cover_an_axis() {
        let n = 9;
        let classes: Vec<usize> = (0..n).map(|p| axis_class(p, n)).collect();
        assert_eq!(classes, [0, 1, 2, 3, 3, 3, 4, 5, 6]);
        for c in 0..N_CLASS {
            assert_eq!(axis_class(class_coordinate(c, n), n), c);
        }
        let n = 7;
        let classes: Vec<usize> = (0..n).map(|p| axis_class(p, n)).collect();
        assert_eq!(classes, [0, 1, 2, 3, 4, 5, 6]);
    }
}
