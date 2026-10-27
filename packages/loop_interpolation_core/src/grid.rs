//! Regular (axis-aligned, uniform-spacing) 3D grid support.
//!
//! Node layout matches the Python `StructuredGrid`: nodes are flattened in
//! Fortran (column-major) order, `flat = i + nx * (j + ny * k)`, so that a
//! flat field can be reshaped as `arr.reshape((nx, ny, nz), order="F")` on
//! the Python side for cross-checking.

#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct StructuredGrid3D {
    pub origin: [f64; 3],
    pub nsteps: [usize; 3],
    pub step_vector: [f64; 3],
    /// Optional rotation carrying the grid's local (axis-aligned) frame into
    /// world space: `world = origin + R * (local - origin)`. `None` means
    /// identity (the common, unrotated case) and is handled as a fast path
    /// that skips the matrix multiply entirely.
    ///
    /// Mirrors the Python `StructuredGrid`'s `rotation_xy` (a 3x3 matrix, or
    /// a single XY-plane angle in degrees via the setter -- see
    /// [`with_rotation_xy_degrees`](Self::with_rotation_xy_degrees)). Note:
    /// `BaseStructuredSupport.rotate()` is dead code on the Python side (grep
    /// confirms it is never called), so there is no numeric Python reference
    /// behaviour to match bit-for-bit here -- this implements the
    /// geometrically standard convention described above: the grid box is
    /// defined in its own unrotated local frame and then rotated about
    /// `origin` into world space, so any world-space point/direction must be
    /// rotated *back* (via the inverse, i.e. transpose, rotation) before the
    /// axis-aligned trilinear math runs. See [`to_local`](Self::to_local) /
    /// [`direction_to_local`](Self::direction_to_local).
    pub rotation: Option<[[f64; 3]; 3]>,
}

impl StructuredGrid3D {
    pub fn new(origin: [f64; 3], nsteps: [usize; 3], step_vector: [f64; 3]) -> Self {
        Self {
            origin,
            nsteps,
            step_vector,
            rotation: None,
        }
    }

    /// Attach an explicit 3x3 rotation matrix (local -> world).
    pub fn with_rotation_matrix(mut self, matrix: [[f64; 3]; 3]) -> Self {
        self.rotation = Some(matrix);
        self
    }

    /// Attach an XY-plane rotation given as an angle in degrees, matching the
    /// Python `rotation_xy` setter's scalar-angle overload: a standard
    /// counter-clockwise rotation about the Z axis,
    /// `[[cos,-sin,0],[sin,cos,0],[0,0,1]]`.
    pub fn with_rotation_xy_degrees(mut self, degrees: f64) -> Self {
        let radians = degrees.to_radians();
        let (s, c) = radians.sin_cos();
        self.rotation = Some([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]);
        self
    }

    #[inline]
    pub fn n_nodes(&self) -> usize {
        self.nsteps[0] * self.nsteps[1] * self.nsteps[2]
    }

    #[inline]
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.nsteps[0] * (j + self.nsteps[1] * k)
    }

    /// Inverse of [`index`]: flat -> (i, j, k).
    #[inline]
    pub fn node_ijk(&self, flat: usize) -> (usize, usize, usize) {
        let nx = self.nsteps[0];
        let ny = self.nsteps[1];
        let i = flat % nx;
        let j = (flat / nx) % ny;
        let k = flat / (nx * ny);
        (i, j, k)
    }

    /// Signed neighbour lookup used by stencil operators. Returns `None` if
    /// the offset neighbour would fall outside the grid.
    #[inline]
    pub fn neighbour_index(
        &self,
        i: usize,
        j: usize,
        k: usize,
        d: (i32, i32, i32),
    ) -> Option<usize> {
        let ni = i as i64 + d.0 as i64;
        let nj = j as i64 + d.1 as i64;
        let nk = k as i64 + d.2 as i64;
        if ni < 0
            || nj < 0
            || nk < 0
            || ni >= self.nsteps[0] as i64
            || nj >= self.nsteps[1] as i64
            || nk >= self.nsteps[2] as i64
        {
            None
        } else {
            Some(self.index(ni as usize, nj as usize, nk as usize))
        }
    }

    /// World-space point -> the grid's own local (unrotated) frame:
    /// `local = origin + R^T (world - origin)`. Identity (returns `world`
    /// unchanged) when `rotation` is `None`. Every position-based lookup
    /// (`trilinear_coeffs`, `trilinear_gradient`) routes world-space input
    /// through this first so the rest of their axis-aligned math never has
    /// to know about rotation at all.
    pub fn to_local(&self, world: [f64; 3]) -> [f64; 3] {
        match &self.rotation {
            None => world,
            Some(r) => {
                let d = [
                    world[0] - self.origin[0],
                    world[1] - self.origin[1],
                    world[2] - self.origin[2],
                ];
                let mut local = self.origin;
                for i in 0..3 {
                    // (R^T d)[i] = sum_j R[j][i] * d[j]
                    local[i] += r[0][i] * d[0] + r[1][i] * d[1] + r[2][i] * d[2];
                }
                local
            }
        }
    }

    /// World-space *direction* -> local frame (rotation only, no
    /// translation): `local_dir = R^T world_dir`. Used to bring
    /// orientation-data vectors (gradient/tangent/norm targets, naturally
    /// specified in world space) into the same local frame that
    /// [`trilinear_gradient`](Self::trilinear_gradient) computes derivatives
    /// in, so that `v_world . grad_world(f) == direction_to_local(v_world) .
    /// grad_local(f)` for any field `f`.
    pub fn direction_to_local(&self, world_dir: [f64; 3]) -> [f64; 3] {
        match &self.rotation {
            None => world_dir,
            Some(r) => {
                let mut local = [0.0; 3];
                for i in 0..3 {
                    local[i] =
                        r[0][i] * world_dir[0] + r[1][i] * world_dir[1] + r[2][i] * world_dir[2];
                }
                local
            }
        }
    }

    /// Shared local-frame cell lookup used by both `trilinear_coeffs` and
    /// `trilinear_gradient`: given an already-local-frame point, returns the
    /// lower-corner cell index `i0` and the fractional local coordinates `t`
    /// (each clamped to `[0, 1]`), or `None` if the point falls outside the
    /// grid's bounding box.
    fn cell_and_local_coords(&self, point: [f64; 3]) -> Option<([usize; 3], [f64; 3])> {
        let mut i0 = [0usize; 3];
        let mut t = [0.0f64; 3];
        for axis in 0..3 {
            let n = self.nsteps[axis];
            if n < 2 {
                return None;
            }
            let step = self.step_vector[axis];
            let rel = (point[axis] - self.origin[axis]) / step;
            let max_index = (n - 1) as f64;
            // Small epsilon tolerance so points exactly on the upper/lower
            // boundary are accepted rather than rejected by fp rounding.
            //
            // Written as "reject unless provably inside" (rather than the
            // logically-equivalent-for-finite-values `rel < -EPS || rel >
            // max_index + EPS`) so a NaN coordinate is rejected too: every
            // comparison against NaN is false, so the `<`/`>` form would
            // silently fall through here and then have `rel.clamp(..)`
            // coerce NaN to 0.0 a few lines down, i.e. treat a garbage
            // coordinate as if it were a real point at the grid's lower
            // corner instead of dropping it like any other out-of-bounds
            // point.
            const EPS: f64 = 1e-9;
            if !(rel >= -EPS && rel <= max_index + EPS) {
                return None;
            }
            let clamped = rel.clamp(0.0, max_index);
            let mut lo = clamped.floor() as usize;
            if lo >= n - 1 {
                lo = n - 2;
            }
            i0[axis] = lo;
            t[axis] = (clamped - lo as f64).clamp(0.0, 1.0);
        }
        Some((i0, t))
    }

    /// Trilinear interpolation coefficients for a point: the 8 corner node
    /// indices of the containing cell and their weights (summing to 1).
    /// Returns `None` if the point lies outside the grid bounding box.
    /// `point` is in world space; internally transformed to the grid's local
    /// frame first (a no-op unless a `rotation` is set).
    pub fn trilinear_coeffs(&self, point: [f64; 3]) -> Option<([usize; 8], [f64; 8])> {
        let point = self.to_local(point);
        let (i0, t) = self.cell_and_local_coords(point)?;

        let mut indices = [0usize; 8];
        let mut coefs = [0.0f64; 8];
        for corner in 0..8 {
            let di = corner & 1;
            let dj = (corner >> 1) & 1;
            let dk = (corner >> 2) & 1;
            let i = i0[0] + di;
            let j = i0[1] + dj;
            let k = i0[2] + dk;
            let wx = if di == 0 { 1.0 - t[0] } else { t[0] };
            let wy = if dj == 0 { 1.0 - t[1] } else { t[1] };
            let wz = if dk == 0 { 1.0 - t[2] } else { t[2] };
            indices[corner] = self.index(i, j, k);
            coefs[corner] = wx * wy * wz;
        }
        Some((indices, coefs))
    }

    /// Trilinear *gradient* shape function: the 8 corner node indices (same
    /// order/indices [`trilinear_coeffs`] would return for the same point)
    /// and, for each corner, `d(coef)/d(x_local), d(coef)/d(y_local),
    /// d(coef)/d(z_local)` -- i.e. the analytic derivative of the same
    /// trilinear basis `trilinear_coeffs` computes, divided by
    /// `step_vector`. This is the Rust analogue of the Python
    /// `StructuredGrid.get_element_gradient_for_location`'s `T` matrix
    /// (shape `(3, 8)` per point; here returned per-corner as `[f64; 3]` for
    /// cache-friendliness). `point` is in world space; internally
    /// transformed to the grid's local frame first, so the returned
    /// gradient components are also in the grid's local axes -- callers
    /// combining this with a world-space direction vector (gradient/norm/
    /// tangent constraints) must first rotate that vector into the local
    /// frame via [`direction_to_local`](Self::direction_to_local).
    pub fn trilinear_gradient(&self, point: [f64; 3]) -> Option<([usize; 8], [[f64; 3]; 8])> {
        let point = self.to_local(point);
        let (i0, t) = self.cell_and_local_coords(point)?;

        let mut indices = [0usize; 8];
        let mut grad = [[0.0f64; 3]; 8];
        for corner in 0..8 {
            let di = corner & 1;
            let dj = (corner >> 1) & 1;
            let dk = (corner >> 2) & 1;
            let i = i0[0] + di;
            let j = i0[1] + dj;
            let k = i0[2] + dk;
            indices[corner] = self.index(i, j, k);

            let wx = if di == 0 { 1.0 - t[0] } else { t[0] };
            let wy = if dj == 0 { 1.0 - t[1] } else { t[1] };
            let wz = if dk == 0 { 1.0 - t[2] } else { t[2] };
            let dwx = if di == 0 { -1.0 } else { 1.0 };
            let dwy = if dj == 0 { -1.0 } else { 1.0 };
            let dwz = if dk == 0 { -1.0 } else { 1.0 };

            grad[corner][0] = dwx * wy * wz / self.step_vector[0];
            grad[corner][1] = wx * dwy * wz / self.step_vector[1];
            grad[corner][2] = wx * wy * dwz / self.step_vector[2];
        }
        Some((indices, grad))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small_grid() -> StructuredGrid3D {
        StructuredGrid3D::new([0.0, 0.0, 0.0], [3, 3, 3], [1.0, 1.0, 1.0])
    }

    #[test]
    fn index_roundtrip() {
        let grid = small_grid();
        for k in 0..3 {
            for j in 0..3 {
                for i in 0..3 {
                    let flat = grid.index(i, j, k);
                    assert_eq!(grid.node_ijk(flat), (i, j, k));
                }
            }
        }
    }

    #[test]
    fn trilinear_at_node_is_one_hot() {
        let grid = small_grid();
        let (idx, coefs) = grid.trilinear_coeffs([1.0, 1.0, 1.0]).unwrap();
        let centre = grid.index(1, 1, 1);
        for (i, c) in idx.iter().zip(coefs.iter()) {
            if *i == centre {
                assert!((c - 1.0).abs() < 1e-12);
            } else {
                assert!(c.abs() < 1e-12);
            }
        }
    }

    #[test]
    fn trilinear_outside_bounds_is_none() {
        let grid = small_grid();
        assert!(grid.trilinear_coeffs([-1.0, 0.0, 0.0]).is_none());
        assert!(grid.trilinear_coeffs([5.0, 0.0, 0.0]).is_none());
    }

    #[test]
    fn trilinear_coefficients_sum_to_one() {
        let grid = small_grid();
        let (_, coefs) = grid.trilinear_coeffs([0.3, 1.7, 0.9]).unwrap();
        let sum: f64 = coefs.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn trilinear_gradient_matches_finite_difference() {
        // Central-difference check of trilinear_gradient against
        // trilinear_coeffs for a handful of interior points.
        let grid = StructuredGrid3D::new([0.0, 0.0, 0.0], [6, 6, 6], [1.0, 1.0, 1.0]);
        let points = [[2.3, 1.1, 3.7], [0.5, 4.9, 0.2], [3.0, 3.0, 3.0]];
        let eps = 1e-6;
        for p in points {
            let (idx, grad) = grid.trilinear_gradient(p).unwrap();
            let (idx_c, _coefs_c) = grid.trilinear_coeffs(p).unwrap();
            assert_eq!(idx, idx_c);
            // Evaluate a synthetic scalar field (value(node) = node index *
            // 0.37) via finite differences on trilinear_coeffs and compare
            // to the analytic derivative from trilinear_gradient.
            for axis in 0..3 {
                let mut p_plus = p;
                let mut p_minus = p;
                p_plus[axis] += eps;
                p_minus[axis] -= eps;
                let (idx_p, coefs_p) = grid.trilinear_coeffs(p_plus).unwrap();
                let (idx_m, coefs_m) = grid.trilinear_coeffs(p_minus).unwrap();
                let field = |i: usize| i as f64 * 0.37;
                let val_p: f64 = idx_p
                    .iter()
                    .zip(coefs_p.iter())
                    .map(|(i, c)| field(*i) * c)
                    .sum();
                let val_m: f64 = idx_m
                    .iter()
                    .zip(coefs_m.iter())
                    .map(|(i, c)| field(*i) * c)
                    .sum();
                let fd = (val_p - val_m) / (2.0 * eps);
                let analytic: f64 = idx
                    .iter()
                    .zip(grad.iter())
                    .map(|(i, g)| field(*i) * g[axis])
                    .sum();
                assert!(
                    (fd - analytic).abs() < 1e-4,
                    "axis={axis} fd={fd} analytic={analytic}"
                );
            }
        }
    }

    #[test]
    fn rotated_grid_matches_unrotated_grid_at_local_frame_point() {
        // A grid rotated 37 degrees about its origin, queried at a
        // world-space point, should produce exactly the coefficients/indices
        // an *unrotated* grid (same origin/nsteps/step_vector) would produce
        // at the corresponding local-frame point.
        let origin = [1.0, 2.0, 0.0];
        let unrotated = StructuredGrid3D::new(origin, [5, 5, 5], [1.0, 1.0, 1.0]);
        let rotated = unrotated.with_rotation_xy_degrees(37.0);

        let local_point = [2.4, 3.1, 1.7];
        let angle = 37f64.to_radians();
        let (s, c) = angle.sin_cos();
        let d = [
            local_point[0] - origin[0],
            local_point[1] - origin[1],
            local_point[2] - origin[2],
        ];
        // world = origin + R * d
        let world = [
            origin[0] + c * d[0] - s * d[1],
            origin[1] + s * d[0] + c * d[1],
            origin[2] + d[2],
        ];

        let expected = unrotated.trilinear_coeffs(local_point).unwrap();
        let actual = rotated.trilinear_coeffs(world).unwrap();
        assert_eq!(expected.0, actual.0);
        for (e, a) in expected.1.iter().zip(actual.1.iter()) {
            assert!((e - a).abs() < 1e-9, "expected={e} actual={a}");
        }
    }
}
