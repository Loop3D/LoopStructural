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
    pub fn value(&self, _p: [f64; 3]) -> f64 {
        unimplemented!("T010")
    }

    /// Return the field gradient at `p`. Return NaN outside the grid or for NaN input.
    pub fn gradient(&self, _p: [f64; 3]) -> [f64; 3] {
        unimplemented!("T010")
    }
}
