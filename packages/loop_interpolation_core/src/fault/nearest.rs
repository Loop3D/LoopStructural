//! Exact nearest-point search in 3D.

/// A set of 3D points for nearest-point search.
#[derive(Debug, Clone)]
pub struct KdTree {
    points: Vec<[f64; 3]>,
}

impl KdTree {
    /// Make a tree from `points`.
    pub fn new(points: &[[f64; 3]]) -> KdTree {
        KdTree {
            points: points.to_vec(),
        }
    }

    /// Return the index of the point nearest to `p`. For equal distances, return the lowest index.
    pub fn nearest(&self, _p: [f64; 3]) -> usize {
        let _ = &self.points;
        unimplemented!("T023")
    }
}
