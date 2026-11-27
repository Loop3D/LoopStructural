//! Exact nearest-point search in 3D.

/// A set of 3D points for nearest-point search.
///
/// The tree is implicit. Each range of the arrays has its median point in the
/// middle. The two halves of the range are the two subtrees.
#[derive(Debug, Clone)]
pub struct KdTree {
    /// The points, in tree order.
    pts: Vec<[f64; 3]>,
    /// The index in the input of each point of `pts`.
    ids: Vec<usize>,
    /// The split axis of each point of `pts`.
    axes: Vec<u8>,
}

fn build(pts: &mut [([f64; 3], usize)], axes: &mut [u8]) {
    if pts.len() < 2 {
        return;
    }
    // Split on the axis with the largest spread.
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for (p, _) in pts.iter() {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let mut axis = 0;
    for k in 1..3 {
        if hi[k] - lo[k] > hi[axis] - lo[axis] {
            axis = k;
        }
    }
    let mid = pts.len() / 2;
    pts.select_nth_unstable_by(mid, |a, b| a.0[axis].total_cmp(&b.0[axis]));
    axes[mid] = axis as u8;
    let (pl, pr) = pts.split_at_mut(mid);
    let (al, ar) = axes.split_at_mut(mid);
    build(pl, al);
    build(&mut pr[1..], &mut ar[1..]);
}

impl KdTree {
    /// Make a tree from `points`.
    pub fn new(points: &[[f64; 3]]) -> KdTree {
        let mut tmp: Vec<([f64; 3], usize)> = points
            .iter()
            .copied()
            .enumerate()
            .map(|(i, p)| (p, i))
            .collect();
        let mut axes = vec![0u8; tmp.len()];
        build(&mut tmp, &mut axes);
        KdTree {
            pts: tmp.iter().map(|t| t.0).collect(),
            ids: tmp.iter().map(|t| t.1).collect(),
            axes,
        }
    }

    /// Return the index of the point nearest to `p`. For equal distances, return the lowest index.
    /// Return `usize::MAX` if the tree has no points.
    pub fn nearest(&self, p: [f64; 3]) -> usize {
        let mut best = (f64::INFINITY, usize::MAX);
        self.search(0, self.pts.len(), p, &mut best);
        best.1
    }

    fn search(&self, lo: usize, hi: usize, p: [f64; 3], best: &mut (f64, usize)) {
        if lo >= hi {
            return;
        }
        let mid = lo + (hi - lo) / 2;
        let q = self.pts[mid];
        let d = (0..3).map(|k| (q[k] - p[k]) * (q[k] - p[k])).sum::<f64>();
        // A NaN query gives no `<` match: the first point is kept.
        if d < best.0 || (d == best.0 && self.ids[mid] < best.1) || best.1 == usize::MAX {
            *best = (d, self.ids[mid]);
        }
        let axis = self.axes[mid] as usize;
        let diff = p[axis] - q[axis];
        let (near, far) = if diff < 0.0 {
            ((lo, mid), (mid + 1, hi))
        } else {
            ((mid + 1, hi), (lo, mid))
        };
        self.search(near.0, near.1, p, best);
        // Strict `>`: an equal-distance point with a lower index must stay reachable.
        if diff * diff <= best.0 {
            self.search(far.0, far.1, p, best);
        }
    }
}
