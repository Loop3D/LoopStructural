//! Tests of the exact nearest-point search.

use loop_interpolation_core::fault::nearest::KdTree;

/// Small LCG with a fixed seed. Return a number in [0, 1).
fn lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*state >> 11) as f64 / (1u64 << 53) as f64
}

fn d2(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| (a[i] - b[i]) * (a[i] - b[i])).sum()
}

/// Brute-force search. Strict `<` keeps the lowest index on ties.
fn brute(points: &[[f64; 3]], p: [f64; 3]) -> usize {
    let mut best = usize::MAX;
    let mut best_d = f64::INFINITY;
    for (i, q) in points.iter().enumerate() {
        let d = d2(*q, p);
        if d < best_d {
            best_d = d;
            best = i;
        }
    }
    best
}

#[test]
fn matches_brute_force() {
    let mut s = 12345u64;
    let points: Vec<[f64; 3]> = (0..2000)
        .map(|_| [lcg(&mut s) * 10.0, lcg(&mut s) * 10.0, lcg(&mut s) * 10.0])
        .collect();
    let tree = KdTree::new(&points);
    for _ in 0..500 {
        let q = [
            lcg(&mut s) * 12.0 - 1.0,
            lcg(&mut s) * 12.0 - 1.0,
            lcg(&mut s) * 12.0 - 1.0,
        ];
        assert_eq!(tree.nearest(q), brute(&points, q));
    }
}

#[test]
fn equal_distance_gives_lowest_index() {
    // Many points on a lattice: the queries have many equal-distance points.
    let mut points = Vec::new();
    for x in 0..6 {
        for y in 0..6 {
            for z in 0..6 {
                points.push([x as f64, y as f64, z as f64]);
            }
        }
    }
    let tree = KdTree::new(&points);
    for q in [
        [0.5, 0.5, 0.5],
        [2.5, 2.5, 2.5],
        [1.0, 1.5, 2.5],
        [-1.0, 3.0, 3.0],
    ] {
        assert_eq!(tree.nearest(q), brute(&points, q));
    }
    let dup = [[1.0, 1.0, 1.0]; 5];
    assert_eq!(KdTree::new(&dup).nearest([0.0, 0.0, 0.0]), 0);
}

#[test]
fn one_point_and_empty() {
    let tree = KdTree::new(&[[1.0, 2.0, 3.0]]);
    assert_eq!(tree.nearest([-9.0, 4.0, 0.0]), 0);
    assert_eq!(KdTree::new(&[]).nearest([0.0, 0.0, 0.0]), usize::MAX);
}
