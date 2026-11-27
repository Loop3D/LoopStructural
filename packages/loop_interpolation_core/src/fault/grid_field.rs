//! A scalar field on a regular 3D grid.

/// A scalar field on a regular 3D grid with trilinear interpolation.
#[derive(Debug, Clone)]
pub struct GridField {
    /// The world position of the first grid node.
    pub origin: [f64; 3],
    /// The number of nodes along each axis.
    pub nsteps: [usize; 3],
    /// The distance between nodes along each axis.
    pub step_vector: [f64; 3],
    /// The origin of the local frame, if the grid has one.
    pub local_origin: Option<[f64; 3]>,
    /// The rotation from the world frame to the local frame, if the grid has one.
    pub local_rotation: Option<[[f64; 3]; 3]>,
    /// The node values.
    pub values: Vec<f64>,
}

impl GridField {
    /// Return the field value at `p`. Return NaN outside the grid or for NaN input.
    pub fn value(&self, p: [f64; 3]) -> f64 {
        match self.lookup(p) {
            Some((corners, t)) => {
                let w = trilinear(t);
                (0..8).map(|c| w[c] * self.values[corners[c]]).sum()
            }
            None => f64::NAN,
        }
    }

    /// Return the field gradient at `p`. Return NaN outside the grid or for NaN input.
    /// The result is in world axes.
    pub fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let Some((corners, t)) = self.lookup(p) else {
            return [f64::NAN; 3];
        };
        let w = gradient_weights(t);
        let mut g = [0.0; 3];
        for (axis, ga) in g.iter_mut().enumerate() {
            *ga = (0..8)
                .map(|c| w[axis][c] * self.values[corners[c]])
                .sum::<f64>()
                / self.step_vector[axis];
        }
        // Rotate the gradient from the local frame to the world frame (transpose of R).
        match self.local_rotation {
            Some(r) => std::array::from_fn(|i| (0..3).map(|j| r[j][i] * g[j]).sum()),
            None => g,
        }
    }

    /// Find the 8 corner node indexes and the local cell coordinates of a world point.
    /// Return None when the point is outside the grid or has a NaN coordinate.
    fn lookup(&self, p: [f64; 3]) -> Option<([usize; 8], [f64; 3])> {
        // World to local frame: R (p - local_origin).
        let mut q = p;
        if let Some(o) = self.local_origin {
            q = std::array::from_fn(|i| p[i] - o[i]);
        }
        if let Some(r) = self.local_rotation {
            q = std::array::from_fn(|i| (0..3).map(|j| r[i][j] * q[j]).sum());
        }
        let mut cell = [0usize; 3];
        let mut t = [0.0; 3];
        for a in 0..3 {
            let ncells = self.nsteps[a] - 1;
            let max = self.origin[a] + ncells as f64 * self.step_vector[a];
            // The boundary has a small tolerance, as in the Python grid.
            let tol = f64::EPSILON * 32.0 * (max - self.origin[a]).abs().max(1.0);
            // The negated form is also false for NaN.
            if !(q[a] >= self.origin[a] - tol && q[a] <= max + tol) {
                return None;
            }
            let f = ((q[a] - self.origin[a]) / self.step_vector[a]).floor();
            // A point on the maximum boundary belongs to the last cell.
            cell[a] = (f.max(0.0) as usize).min(ncells - 1);
            let local = (q[a] - (self.origin[a] + self.step_vector[a] * cell[a] as f64))
                / self.step_vector[a];
            t[a] = local.clamp(0.0, 1.0);
        }
        let (nx, ny) = (self.nsteps[0], self.nsteps[1]);
        let g = cell[0] + nx * cell[1] + nx * ny * cell[2];
        let offsets = [
            0,
            1,
            nx,
            nx + 1,
            nx * ny,
            nx * ny + 1,
            nx * ny + nx,
            nx * ny + nx + 1,
        ];
        Some((offsets.map(|o| g + o), t))
    }
}

/// Return the trilinear weights of the 8 cell corners (x runs fastest).
fn trilinear(t: [f64; 3]) -> [f64; 8] {
    let [x, y, z] = t;
    [
        (1.0 - x) * (1.0 - y) * (1.0 - z),
        x * (1.0 - y) * (1.0 - z),
        (1.0 - x) * y * (1.0 - z),
        x * y * (1.0 - z),
        (1.0 - x) * (1.0 - y) * z,
        x * (1.0 - y) * z,
        (1.0 - x) * y * z,
        x * y * z,
    ]
}

/// Return the gradient weights of the 8 cell corners for each axis, in cell units.
fn gradient_weights(t: [f64; 3]) -> [[f64; 8]; 3] {
    let [x, y, z] = t;
    [
        [
            (1.0 - z) * (y - 1.0),
            (1.0 - y) * (1.0 - z),
            -y * (1.0 - z),
            y * (1.0 - z),
            -(1.0 - y) * z,
            (1.0 - y) * z,
            -y * z,
            y * z,
        ],
        [
            (x - 1.0) * (1.0 - z),
            -x * (1.0 - z),
            (1.0 - x) * (1.0 - z),
            x * (1.0 - z),
            -(1.0 - x) * z,
            -x * z,
            (1.0 - x) * z,
            x * z,
        ],
        [
            -(1.0 - x) * (1.0 - y),
            -x * (1.0 - y),
            -(1.0 - x) * y,
            -x * y,
            (1.0 - x) * (1.0 - y),
            x * (1.0 - y),
            (1.0 - x) * y,
            x * y,
        ],
    ]
}
