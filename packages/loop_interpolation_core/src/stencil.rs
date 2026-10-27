//! Matrix-free translation-invariant stencil operators.
//!
//! This is the Rust analogue of the Python `FiniteDifferenceInterpolator`'s
//! `_assemble_operator` / `matrix_free_regularisation_blocks` path, but
//! instead of gathering explicit column-index arrays (`neighbour_global_indexes`)
//! and storing them per family, a stencil here is just a short list of
//! `(offset, coefficient)` pairs applied directly against the dense 3D grid.
//! Nothing is ever materialised as a matrix: `forward`/`adjoint` are pure
//! functions of the field array, so the whole regularisation system costs
//! `O(dof)` time and `O(1)` extra memory beyond the field itself (plus
//! whatever scratch the caller supplies for the output).
//!
//! `adjoint` is implemented as a *correlation with the flipped stencil*
//! rather than a scatter-add: `adjoint(y)[c] = sum_offset coeff * y[c -
//! offset]`, treating a source row `c - offset` as contributing zero (rather
//! than reading `y` there) whenever that row is inactive (out of grid, or
//! within `margin` of a boundary) -- exactly mirroring which rows `forward`
//! zeroes out, so this is `A^T y` for any `y`, not just vectors already
//! produced by `forward`. No index collisions to resolve and no atomics
//! required, so both passes parallelise trivially over grid nodes.

use crate::PAR_MIN_LEN;
use crate::grid::StructuredGrid3D;
use rayon::prelude::*;

#[derive(Debug, Clone)]
pub struct Stencil {
    /// (di, dj, dk) neighbour offsets, paired index-for-index with `coeffs`.
    pub offsets: Vec<(i32, i32, i32)>,
    /// L2-normalised coefficients (unit-norm row, matching the Python
    /// `add_constraints_to_least_squares` / `_store_matrix_free_regularisation_block`
    /// row normalisation).
    pub coeffs: Vec<f64>,
    /// Per-axis maximum absolute offset; a row is only evaluated at nodes at
    /// least `margin[axis]` away from that axis's boundary.
    pub margin: [i32; 3],
}

impl Stencil {
    pub fn new(taps: &[((i32, i32, i32), f64)]) -> Self {
        let norm = taps.iter().map(|(_, v)| v * v).sum::<f64>().sqrt();
        let offsets = taps.iter().map(|(o, _)| *o).collect::<Vec<_>>();
        let coeffs = taps
            .iter()
            .map(|(_, v)| if norm > 0.0 { v / norm } else { *v })
            .collect::<Vec<_>>();
        let mut margin = [0i32; 3];
        for (di, dj, dk) in &offsets {
            margin[0] = margin[0].max(di.abs());
            margin[1] = margin[1].max(dj.abs());
            margin[2] = margin[2].max(dk.abs());
        }
        Self {
            offsets,
            coeffs,
            margin,
        }
    }

    /// Whether the stencil's full footprint fits at node `(i, j, k)`, i.e.
    /// whether this is a valid regularisation row there. `pub(crate)` so
    /// `FiniteDifferenceInterpolator`'s fused six-family pass (see
    /// `apply_normal_operator_into`) can reuse it directly instead of going
    /// through `forward`/`adjoint` one family at a time.
    #[inline]
    pub(crate) fn is_row_active(
        &self,
        grid: &StructuredGrid3D,
        i: usize,
        j: usize,
        k: usize,
    ) -> bool {
        let ijk = [i as i64, j as i64, k as i64];
        for (axis, coord) in ijk.into_iter().enumerate() {
            let m = self.margin[axis] as i64;
            let n = grid.nsteps[axis] as i64;
            if coord < m || coord >= n - m {
                return false;
            }
        }
        true
    }

    /// `out = A x` (overwrites `out`, no allocation). Rows too close to a
    /// boundary for the full stencil footprint to fit are zero, matching
    /// the Python "interior-only" regularisation rows.
    pub fn forward_into(&self, grid: &StructuredGrid3D, x: &[f64], out: &mut [f64]) {
        let nx = grid.nsteps[0];
        let ny = grid.nsteps[1];
        out.par_iter_mut()
            .enumerate()
            .with_min_len(PAR_MIN_LEN)
            .for_each(|(flat, o)| {
                let i = flat % nx;
                let j = (flat / nx) % ny;
                let k = flat / (nx * ny);
                *o = if self.is_row_active(grid, i, j, k) {
                    self.offsets
                        .iter()
                        .zip(self.coeffs.iter())
                        .map(|(d, c)| {
                            c * x[grid
                                .neighbour_index(i, j, k, *d)
                                .expect("row within margin")]
                        })
                        .sum()
                } else {
                    0.0
                };
            });
    }

    /// `y = A x`, allocating a fresh output vector. Convenience wrapper
    /// around [`forward_into`]; prefer the `_into` form in a hot loop (see
    /// `FiniteDifferenceInterpolator::solve`).
    pub fn forward(&self, grid: &StructuredGrid3D, x: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; grid.n_nodes()];
        self.forward_into(grid, x, &mut out);
        out
    }

    /// `out += scale * A^T y` (accumulates, no allocation), implemented as
    /// correlation with the flipped stencil (see module docs) so it needs
    /// no scatter-add and no atomics.
    pub fn adjoint_accumulate_into(
        &self,
        grid: &StructuredGrid3D,
        y: &[f64],
        scale: f64,
        out: &mut [f64],
    ) {
        if scale == 0.0 {
            return;
        }
        let nx = grid.nsteps[0];
        let ny = grid.nsteps[1];
        out.par_iter_mut()
            .enumerate()
            .with_min_len(PAR_MIN_LEN)
            .for_each(|(flat, o)| {
                let i = flat % nx;
                let j = (flat / nx) % ny;
                let k = flat / (nx * ny);
                let value: f64 = self
                    .offsets
                    .iter()
                    .zip(self.coeffs.iter())
                    .map(|(d, c)| {
                        let flipped = (-d.0, -d.1, -d.2);
                        match grid.neighbour_index(i, j, k, flipped) {
                            // The source row `src` only actually contributes if
                            // it was an active row in `forward` (A[src, c] is
                            // defined to be 0 for inactive rows regardless of
                            // `y[src]`, since `y` here is an arbitrary covector,
                            // not necessarily `forward`'s output).
                            //
                            // `(si, sj, sk)` is just `(i, j, k) + flipped`,
                            // already validated in-bounds by `neighbour_index`
                            // above -- recover it with that same cheap add
                            // instead of decomposing `src` back via `node_ijk`'s
                            // two divisions + two modulos.
                            Some(src) => {
                                let si = (i as i64 + flipped.0 as i64) as usize;
                                let sj = (j as i64 + flipped.1 as i64) as usize;
                                let sk = (k as i64 + flipped.2 as i64) as usize;
                                if self.is_row_active(grid, si, sj, sk) {
                                    c * y[src]
                                } else {
                                    0.0
                                }
                            }
                            None => 0.0,
                        }
                    })
                    .sum();
                *o += scale * value;
            });
    }

    /// `out = A^T y`, allocating a fresh output vector. Convenience wrapper
    /// around [`adjoint_accumulate_into`].
    pub fn adjoint(&self, grid: &StructuredGrid3D, y: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; grid.n_nodes()];
        self.adjoint_accumulate_into(grid, y, 1.0, &mut out);
        out
    }

    /// `out += A^T (weight2 .* (A x))`, i.e. this stencil's contribution to
    /// the normal-equations operator when the row weight varies *per node*
    /// instead of being one uniform scalar (`weight2[c]` is squared-weight at
    /// row/node `c`) -- the anisotropic/directional-regularisation analogue
    /// of [`accumulate_normal_into`](Self::accumulate_normal_into), which
    /// only supports a uniform scalar weight. `scratch` is reused as the
    /// intermediate `A x` buffer (then overwritten in place with `weight2 .*
    /// (A x)`), so this makes no allocations of its own -- same contract as
    /// `accumulate_normal_into`.
    pub fn accumulate_normal_weighted_into(
        &self,
        grid: &StructuredGrid3D,
        x: &[f64],
        weight2: &[f64],
        scratch: &mut [f64],
        out: &mut [f64],
    ) {
        self.forward_into(grid, x, scratch);
        scratch
            .par_iter_mut()
            .zip(weight2.par_iter())
            .with_min_len(PAR_MIN_LEN)
            .for_each(|(s, w2)| *s *= w2);
        self.adjoint_accumulate_into(grid, scratch, 1.0, out);
    }

    /// `diag(A^T diag(weight2) A)`, the per-node-weighted analogue of
    /// [`diagonal`](Self::diagonal): mirrors its row-activity logic exactly,
    /// but accumulates `weight2[src] * coeff^2` (the *row's* -- i.e. source
    /// node's -- squared weight) instead of a uniform scalar times `coeff^2`.
    pub fn diagonal_weighted(&self, grid: &StructuredGrid3D, weight2: &[f64]) -> Vec<f64> {
        let nx = grid.nsteps[0];
        let ny = grid.nsteps[1];
        (0..grid.n_nodes())
            .into_par_iter()
            .with_min_len(PAR_MIN_LEN)
            .map(|flat| {
                let i = flat % nx;
                let j = (flat / nx) % ny;
                let k = flat / (nx * ny);
                self.offsets
                    .iter()
                    .zip(self.coeffs.iter())
                    .map(|(d, c)| {
                        let flipped = (-d.0, -d.1, -d.2);
                        match grid.neighbour_index(i, j, k, flipped) {
                            Some(src) => {
                                let si = (i as i64 + flipped.0 as i64) as usize;
                                let sj = (j as i64 + flipped.1 as i64) as usize;
                                let sk = (k as i64 + flipped.2 as i64) as usize;
                                if self.is_row_active(grid, si, sj, sk) {
                                    weight2[src] * c * c
                                } else {
                                    0.0
                                }
                            }
                            None => 0.0,
                        }
                    })
                    .sum()
            })
            .collect()
    }

    /// `diag(A^T A)`, used to build a Jacobi preconditioner without ever
    /// forming `A`. Mirrors `adjoint`'s row-activity logic but accumulates
    /// `coeff^2` instead of `coeff * y[src]`. Only called once per solve
    /// (to build the preconditioner), so allocating here is fine.
    pub fn diagonal(&self, grid: &StructuredGrid3D) -> Vec<f64> {
        let nx = grid.nsteps[0];
        let ny = grid.nsteps[1];
        (0..grid.n_nodes())
            .into_par_iter()
            .with_min_len(PAR_MIN_LEN)
            .map(|flat| {
                let i = flat % nx;
                let j = (flat / nx) % ny;
                let k = flat / (nx * ny);
                self.offsets
                    .iter()
                    .zip(self.coeffs.iter())
                    .map(|(d, c)| {
                        let flipped = (-d.0, -d.1, -d.2);
                        match grid.neighbour_index(i, j, k, flipped) {
                            Some(_) => {
                                let si = (i as i64 + flipped.0 as i64) as usize;
                                let sj = (j as i64 + flipped.1 as i64) as usize;
                                let sk = (k as i64 + flipped.2 as i64) as usize;
                                if self.is_row_active(grid, si, sj, sk) {
                                    c * c
                                } else {
                                    0.0
                                }
                            }
                            None => 0.0,
                        }
                    })
                    .sum()
            })
            .collect()
    }

    /// `out += w^2 * A^T A x`, i.e. this stencil's contribution to the
    /// weighted normal-equations operator. `scratch` is reused as the
    /// intermediate `A x` buffer so this makes no allocations of its own,
    /// which matters because it runs on every conjugate-gradient iteration.
    pub fn accumulate_normal_into(
        &self,
        grid: &StructuredGrid3D,
        x: &[f64],
        weight: f64,
        scratch: &mut [f64],
        out: &mut [f64],
    ) {
        if weight == 0.0 {
            return;
        }
        self.forward_into(grid, x, scratch);
        self.adjoint_accumulate_into(grid, scratch, weight * weight, out);
    }

    /// `out += w^2 * A^T A x`, allocating its own scratch buffer. Convenience
    /// wrapper around [`accumulate_normal_into`] for callers outside the CG
    /// hot loop (e.g. tests).
    pub fn accumulate_normal(
        &self,
        grid: &StructuredGrid3D,
        x: &[f64],
        weight: f64,
        out: &mut [f64],
    ) {
        let mut scratch = vec![0.0; grid.n_nodes()];
        self.accumulate_normal_into(grid, x, weight, &mut scratch, out);
    }
}

/// The six canonical second-derivative stencil families used for isotropic
/// curvature regularisation, mirroring `Operator.Dxx_mask` ... `Dyz_mask`.
pub fn dxx() -> Stencil {
    Stencil::new(&[((-1, 0, 0), 1.0), ((0, 0, 0), -2.0), ((1, 0, 0), 1.0)])
}

pub fn dyy() -> Stencil {
    Stencil::new(&[((0, -1, 0), 1.0), ((0, 0, 0), -2.0), ((0, 1, 0), 1.0)])
}

pub fn dzz() -> Stencil {
    Stencil::new(&[((0, 0, -1), 1.0), ((0, 0, 0), -2.0), ((0, 0, 1), 1.0)])
}

pub fn dxy() -> Stencil {
    Stencil::new(&[
        ((-1, -1, 0), 1.0),
        ((1, -1, 0), -1.0),
        ((-1, 1, 0), -1.0),
        ((1, 1, 0), 1.0),
    ])
}

pub fn dxz() -> Stencil {
    Stencil::new(&[
        ((-1, 0, -1), 1.0),
        ((1, 0, -1), -1.0),
        ((-1, 0, 1), -1.0),
        ((1, 0, 1), 1.0),
    ])
}

pub fn dyz() -> Stencil {
    Stencil::new(&[
        ((0, -1, -1), 1.0),
        ((0, 1, -1), -1.0),
        ((0, -1, 1), -1.0),
        ((0, 1, 1), 1.0),
    ])
}

/// A one-sided first-derivative stencil active only on a fixed axis-aligned
/// "shell" layer (every node with `coord[axis] == layer`), rather than
/// everywhere-away-from-a-margin like [`Stencil`].
///
/// This is the Rust analogue of the Python `FiniteDifferenceInterpolator`'s
/// border regularisation terms (`assemble_borders`, built from
/// `Operator.Dx_forward_mask`/`Dx_backward_mask` and their y/z equivalents,
/// applied only at `_boundary_indexes`' single interior layer just inside
/// each face) -- structurally different from the six interior curvature
/// families, which is why it needs its own type instead of reusing
/// [`Stencil`]'s margin-based activation rule.
#[derive(Debug, Clone)]
pub struct EdgeStencil {
    /// Which axis (0=x, 1=y, 2=z) the taps and the fixed layer are along.
    pub axis: usize,
    /// The fixed index along `axis` at which this stencil's rows live (every
    /// node with `coord[axis] == layer`, for any index along the other two
    /// axes).
    pub layer: usize,
    /// (signed step along `axis`, L2-normalised coefficient) taps.
    pub taps: Vec<(i32, f64)>,
}

impl EdgeStencil {
    pub fn new(axis: usize, layer: usize, taps: &[(i32, f64)]) -> Self {
        let norm = taps.iter().map(|(_, v)| v * v).sum::<f64>().sqrt();
        let taps = taps
            .iter()
            .map(|(o, v)| (*o, if norm > 0.0 { v / norm } else { *v }))
            .collect();
        Self { axis, layer, taps }
    }

    #[inline]
    fn offset3(&self, o: i32) -> (i32, i32, i32) {
        match self.axis {
            0 => (o, 0, 0),
            1 => (0, o, 0),
            _ => (0, 0, o),
        }
    }

    #[inline]
    fn is_row_active(&self, ijk: [usize; 3]) -> bool {
        ijk[self.axis] == self.layer
    }

    /// `(i, j, k)` of the `pos`-th node (`pos` in `0..layer_size(grid)`) on
    /// the fixed layer `coord[axis] == layer_value`, in an arbitrary but
    /// fixed row-major enumeration over the two non-`axis` dimensions.
    #[inline]
    fn ijk_at_layer(
        &self,
        nx: usize,
        ny: usize,
        layer_value: usize,
        pos: usize,
    ) -> (usize, usize, usize) {
        match self.axis {
            0 => (layer_value, pos % ny, pos / ny),
            1 => (pos % nx, layer_value, pos / nx),
            _ => (pos % nx, pos / nx, layer_value),
        }
    }

    /// Node count of one fixed-`axis` layer (`grid.n_nodes() /
    /// grid.nsteps[axis]`) -- every position `forward_into` ever writes to,
    /// and the size of one `target_layers` pass in
    /// `adjoint_accumulate_into`/`diagonal`.
    #[inline]
    fn layer_size(&self, grid: &StructuredGrid3D) -> usize {
        grid.n_nodes() / grid.nsteps[self.axis]
    }

    /// The (at most `taps.len()`, deduplicated) axis-coordinates an output
    /// node can ever receive a nonzero adjoint contribution at: `y[src]` is
    /// only read (see `adjoint_accumulate_into`) when the *source* row,
    /// reached by stepping `-off` from the output node, lands exactly on
    /// `layer` -- i.e. `coord[axis] - off == layer`, so the output node's
    /// own `coord[axis]` must be `layer + off` for some tap. Taps are always
    /// single-step (`edge_forward`/`edge_backward` only ever emit `off` in
    /// `{-1, 0, 1}`), so this is at most 2 distinct layers in practice, out
    /// of `nsteps[axis]` possible ones -- the source of this type's whole
    /// speedup over a naive full-`dof` scan.
    fn target_layers(&self, n_axis: usize) -> Vec<usize> {
        let mut layers: Vec<usize> = self
            .taps
            .iter()
            .filter_map(|(off, _)| {
                let t = self.layer as i64 + *off as i64;
                (t >= 0 && (t as usize) < n_axis).then_some(t as usize)
            })
            .collect();
        layers.sort_unstable();
        layers.dedup();
        layers
    }

    /// `out = A x` (overwrites every row this stencil can ever produce;
    /// every other row -- everywhere off this stencil's fixed `layer` -- is
    /// left untouched rather than explicitly zeroed, since this stencil's
    /// rows are zero there by construction and every caller either starts
    /// from an already-zeroed buffer or (`accumulate_normal_into`) only ever
    /// reads `out` back at positions this same call just wrote). Computes
    /// each layer node's value in parallel into a small `layer_size`-length
    /// buffer, then scatter-writes it into `out` serially (`layer_size` is
    /// small enough -- one node layer, not the full grid -- that this beats
    /// touching all `dof` positions to filter out the inactive ones, which
    /// is what a naive full-grid pass would otherwise spend almost all its
    /// time doing).
    pub fn forward_into(&self, grid: &StructuredGrid3D, x: &[f64], out: &mut [f64]) {
        let nx = grid.nsteps[0];
        let ny = grid.nsteps[1];
        let layer_size = self.layer_size(grid);
        let values: Vec<f64> = (0..layer_size)
            .into_par_iter()
            .with_min_len(PAR_MIN_LEN)
            .map(|pos| {
                let (i, j, k) = self.ijk_at_layer(nx, ny, self.layer, pos);
                self.taps
                    .iter()
                    .map(|(off, c)| {
                        c * x[grid
                            .neighbour_index(i, j, k, self.offset3(*off))
                            .expect("edge stencil taps stay within the grid by construction")]
                    })
                    .sum()
            })
            .collect();
        for (pos, value) in values.into_iter().enumerate() {
            let (i, j, k) = self.ijk_at_layer(nx, ny, self.layer, pos);
            out[grid.index(i, j, k)] = value;
        }
    }

    pub fn forward(&self, grid: &StructuredGrid3D, x: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; grid.n_nodes()];
        self.forward_into(grid, x, &mut out);
        out
    }

    /// `out += scale * A^T y`; see [`Stencil::adjoint_accumulate_into`] for
    /// the flipped-stencil-correlation reasoning this mirrors, restricted
    /// (see [`target_layers`](Self::target_layers)) to the handful of axis
    /// layers that can ever receive a nonzero contribution instead of
    /// scanning the full `dof`-length grid.
    pub fn adjoint_accumulate_into(
        &self,
        grid: &StructuredGrid3D,
        y: &[f64],
        scale: f64,
        out: &mut [f64],
    ) {
        if scale == 0.0 {
            return;
        }
        let nx = grid.nsteps[0];
        let ny = grid.nsteps[1];
        let n_axis = grid.nsteps[self.axis];
        let layer_size = self.layer_size(grid);
        for target_layer in self.target_layers(n_axis) {
            let contributions: Vec<f64> = (0..layer_size)
                .into_par_iter()
                .with_min_len(PAR_MIN_LEN)
                .map(|pos| {
                    let (i, j, k) = self.ijk_at_layer(nx, ny, target_layer, pos);
                    let ijk = [i, j, k];
                    self.taps
                        .iter()
                        .map(|(off, c)| {
                            let flipped = self.offset3(-*off);
                            match grid.neighbour_index(i, j, k, flipped) {
                                Some(src) => {
                                    let mut sijk = ijk;
                                    sijk[self.axis] =
                                        (ijk[self.axis] as i64 - *off as i64) as usize;
                                    if self.is_row_active(sijk) {
                                        c * y[src]
                                    } else {
                                        0.0
                                    }
                                }
                                None => 0.0,
                            }
                        })
                        .sum()
                })
                .collect();
            for (pos, value) in contributions.into_iter().enumerate() {
                let (i, j, k) = self.ijk_at_layer(nx, ny, target_layer, pos);
                out[grid.index(i, j, k)] += scale * value;
            }
        }
    }

    pub fn adjoint(&self, grid: &StructuredGrid3D, y: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; grid.n_nodes()];
        self.adjoint_accumulate_into(grid, y, 1.0, &mut out);
        out
    }

    /// `diag(A^T A)`, restricted the same way as
    /// [`adjoint_accumulate_into`](Self::adjoint_accumulate_into) (same
    /// `is_row_active`-gated per-tap sum, just `c * c` instead of `c *
    /// y[src]`).
    pub fn diagonal(&self, grid: &StructuredGrid3D) -> Vec<f64> {
        let nx = grid.nsteps[0];
        let ny = grid.nsteps[1];
        let n_axis = grid.nsteps[self.axis];
        let layer_size = self.layer_size(grid);
        let mut diag = vec![0.0; grid.n_nodes()];
        for target_layer in self.target_layers(n_axis) {
            let contributions: Vec<f64> = (0..layer_size)
                .into_par_iter()
                .with_min_len(PAR_MIN_LEN)
                .map(|pos| {
                    let (i, j, k) = self.ijk_at_layer(nx, ny, target_layer, pos);
                    let ijk = [i, j, k];
                    self.taps
                        .iter()
                        .map(|(off, c)| {
                            let flipped = self.offset3(-*off);
                            match grid.neighbour_index(i, j, k, flipped) {
                                Some(_) => {
                                    let mut sijk = ijk;
                                    sijk[self.axis] =
                                        (ijk[self.axis] as i64 - *off as i64) as usize;
                                    if self.is_row_active(sijk) { c * c } else { 0.0 }
                                }
                                None => 0.0,
                            }
                        })
                        .sum()
                })
                .collect();
            for (pos, value) in contributions.into_iter().enumerate() {
                let (i, j, k) = self.ijk_at_layer(nx, ny, target_layer, pos);
                diag[grid.index(i, j, k)] = value;
            }
        }
        diag
    }

    /// `out += w^2 * A^T A x`, using `scratch` as the intermediate `A x`
    /// buffer (`forward_into`/`adjoint_accumulate_into` do each allocate a
    /// small `layer_size`-length scratch `Vec` of their own -- see their
    /// docs -- but nothing proportional to `dof`).
    pub fn accumulate_normal_into(
        &self,
        grid: &StructuredGrid3D,
        x: &[f64],
        weight: f64,
        scratch: &mut [f64],
        out: &mut [f64],
    ) {
        if weight == 0.0 {
            return;
        }
        self.forward_into(grid, x, scratch);
        self.adjoint_accumulate_into(grid, scratch, weight * weight, out);
    }

    /// `out += w^2 * A^T A x`, allocating its own scratch buffer. Convenience
    /// wrapper around [`accumulate_normal_into`](Self::accumulate_normal_into)
    /// for callers outside the CG hot loop (e.g. tests).
    pub fn accumulate_normal(
        &self,
        grid: &StructuredGrid3D,
        x: &[f64],
        weight: f64,
        out: &mut [f64],
    ) {
        let mut scratch = vec![0.0; grid.n_nodes()];
        self.accumulate_normal_into(grid, x, weight, &mut scratch, out);
    }
}

/// Border regularisation stencil for the layer just inside the *lower* face
/// along `axis` (index 1): taps at offset -1 (the true boundary node, index
/// 0) and 0 (the centre, index 1), matching the Python
/// `Operator.Dx_forward_mask` (see module docs re: sign convention -- these
/// rows always target zero, so the operator's contribution to the
/// least-squares system is invariant to an overall row sign flip). Returns
/// `None` if `nsteps[axis] < 3`, matching the Python
/// `_boundary_indexes`/`assemble_borders`'s skip for axes too thin to have a
/// distinct "one layer in from the boundary" row.
pub fn edge_forward(axis: usize, nsteps: &[usize; 3]) -> Option<EdgeStencil> {
    if nsteps[axis] < 3 {
        return None;
    }
    Some(EdgeStencil::new(axis, 1, &[(-1, -1.0), (0, 1.0)]))
}

/// Border regularisation stencil for the layer just inside the *upper* face
/// along `axis` (index `nsteps[axis] - 2`): taps at offset 0 (centre) and +1
/// (the true boundary node), matching `Operator.Dx_backward_mask`. `None`
/// under the same `nsteps[axis] < 3` condition as [`edge_forward`].
pub fn edge_backward(axis: usize, nsteps: &[usize; 3]) -> Option<EdgeStencil> {
    if nsteps[axis] < 3 {
        return None;
    }
    Some(EdgeStencil::new(
        axis,
        nsteps[axis] - 2,
        &[(0, -1.0), (1, 1.0)],
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid() -> StructuredGrid3D {
        StructuredGrid3D::new([0.0, 0.0, 0.0], [5, 5, 5], [1.0, 1.0, 1.0])
    }

    #[test]
    fn dxx_zero_on_linear_field() {
        let grid = grid();
        // f(i,j,k) = i (a linear function of x), dxx should vanish everywhere.
        let x: Vec<f64> = (0..grid.n_nodes())
            .map(|flat| grid.node_ijk(flat).0 as f64)
            .collect();
        let y = dxx().forward(&grid, &x);
        assert!(y.iter().all(|v| v.abs() < 1e-12));
    }

    #[test]
    fn dxx_nonzero_on_quadratic_field() {
        let grid = grid();
        let x: Vec<f64> = (0..grid.n_nodes())
            .map(|flat| {
                let i = grid.node_ijk(flat).0 as f64;
                i * i
            })
            .collect();
        let y = dxx().forward(&grid, &x);
        // interior nodes (i in [1, nx-2]) should be nonzero
        let (nx, _, _) = (grid.nsteps[0], grid.nsteps[1], grid.nsteps[2]);
        let mid = grid.index(nx / 2, 2, 2);
        assert!(y[mid].abs() > 1e-6);
    }

    #[test]
    fn adjoint_matches_transpose_dot_product() {
        // <A x, y> == <x, A^T y> for a random-ish field.
        let grid = grid();
        let n = grid.n_nodes();
        let x: Vec<f64> = (0..n).map(|i| ((i * 37 + 11) % 97) as f64 / 97.0).collect();
        let y: Vec<f64> = (0..n).map(|i| ((i * 53 + 7) % 89) as f64 / 89.0).collect();

        for stencil in [dxx(), dyy(), dzz(), dxy(), dxz(), dyz()] {
            let ax = stencil.forward(&grid, &x);
            let aty = stencil.adjoint(&grid, &y);
            let lhs: f64 = ax.iter().zip(y.iter()).map(|(a, b)| a * b).sum();
            let rhs: f64 = x.iter().zip(aty.iter()).map(|(a, b)| a * b).sum();
            assert!((lhs - rhs).abs() < 1e-9, "lhs={lhs} rhs={rhs}");
        }
    }

    fn edge_stencils(nsteps: &[usize; 3]) -> Vec<EdgeStencil> {
        let mut v = Vec::new();
        for axis in 0..3 {
            if let Some(s) = edge_forward(axis, nsteps) {
                v.push(s);
            }
            if let Some(s) = edge_backward(axis, nsteps) {
                v.push(s);
            }
        }
        v
    }

    #[test]
    fn edge_stencil_adjoint_matches_transpose_dot_product() {
        let grid = grid();
        let n = grid.n_nodes();
        let x: Vec<f64> = (0..n).map(|i| ((i * 37 + 11) % 97) as f64 / 97.0).collect();
        let y: Vec<f64> = (0..n).map(|i| ((i * 53 + 7) % 89) as f64 / 89.0).collect();

        for stencil in edge_stencils(&grid.nsteps) {
            let ax = stencil.forward(&grid, &x);
            let aty = stencil.adjoint(&grid, &y);
            let lhs: f64 = ax.iter().zip(y.iter()).map(|(a, b)| a * b).sum();
            let rhs: f64 = x.iter().zip(aty.iter()).map(|(a, b)| a * b).sum();
            assert!(
                (lhs - rhs).abs() < 1e-9,
                "axis={} layer={} lhs={lhs} rhs={rhs}",
                stencil.axis,
                stencil.layer
            );
        }
    }

    #[test]
    fn edge_stencil_only_active_on_its_fixed_layer() {
        let grid = grid();
        let n = grid.n_nodes();
        let x: Vec<f64> = (0..n).map(|i| ((i * 13 + 3) % 71) as f64 / 71.0).collect();

        for stencil in edge_stencils(&grid.nsteps) {
            let y = stencil.forward(&grid, &x);
            #[allow(clippy::needless_range_loop)]
            for flat in 0..n {
                let (i, j, k) = grid.node_ijk(flat);
                let ijk = [i, j, k][stencil.axis];
                if ijk != stencil.layer {
                    assert!(y[flat].abs() < 1e-12, "expected zero row away from layer");
                }
            }
            // At least one row on the layer should be nonzero for a
            // non-constant field.
            assert!(y.iter().any(|v| v.abs() > 1e-9));
        }
    }

    #[test]
    fn edge_stencil_zero_on_constant_field() {
        let grid = grid();
        let x = vec![1.0; grid.n_nodes()];
        for stencil in edge_stencils(&grid.nsteps) {
            let y = stencil.forward(&grid, &x);
            assert!(y.iter().all(|v| v.abs() < 1e-12));
        }
    }

    #[test]
    fn edge_stencil_none_for_thin_axis() {
        let nsteps = [2, 5, 5];
        assert!(edge_forward(0, &nsteps).is_none());
        assert!(edge_backward(0, &nsteps).is_none());
        assert!(edge_forward(1, &nsteps).is_some());
    }
}
