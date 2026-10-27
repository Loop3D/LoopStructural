//! Tests of the memory estimate and the budget check (story 4).
//! The tests use the counting allocator of `common`.

mod common;

use std::sync::{Mutex, MutexGuard};

use loop_interpolation_core::{CoreError, Problem, Settings};

/// The counting allocator is global for this test binary.
/// Tests run in parallel threads. Each test takes this lock,
/// so that only one test at a time measures allocation.
static SERIAL: Mutex<()> = Mutex::new(());

fn lock() -> MutexGuard<'static, ()> {
    // A failed test poisons the lock. Ignore the poison.
    SERIAL.lock().unwrap_or_else(|e| e.into_inner())
}

/// Small deterministic generator (LCG). It returns a number in [0, 1).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

/// Make `n` flat random points inside a grid of `nodes` per axis.
fn points(rng: &mut Lcg, n: usize, nodes: usize) -> Vec<f64> {
    let extent = (nodes - 1) as f64;
    (0..3 * n).map(|_| rng.next() * extent).collect()
}

/// Make a problem like R1: a cubic grid, value rows and normal rows.
fn r1_problem(nodes: usize, n_value: usize, n_normal: usize) -> Problem {
    let mut p = Problem::new([0.0; 3], [nodes; 3], [1.0; 3]).expect("grid must be valid");
    let mut rng = Lcg(12345);
    let pts = points(&mut rng, n_value, nodes);
    let values: Vec<f64> = (0..n_value)
        .map(|i| pts[3 * i] + 0.5 * pts[3 * i + 1] + 0.25 * pts[3 * i + 2])
        .collect();
    p.add_value(&pts, &values, &vec![1.0; n_value])
        .expect("value rows must be valid");
    let pts = points(&mut rng, n_normal, nodes);
    let normals: Vec<f64> = (0..n_normal).flat_map(|_| [1.0, 0.5, 0.25]).collect();
    p.add_normal(&pts, &normals, &vec![1.0; n_normal])
        .expect("normal rows must be valid");
    p
}

const MB: u64 = 1_000_000;

// The estimate for 100^3 with R1 data is between 40 MB and 100 MB.
// All parts are more than 0.
#[test]
fn estimate_100_cubed_is_in_range() {
    let _guard = lock();
    let p = r1_problem(100, 2000, 500);
    let e = p
        .estimate_memory(&Settings::default())
        .expect("estimate must succeed");
    assert!(e.vectors > 0, "vectors must be more than 0");
    assert!(e.constraints > 0, "constraints must be more than 0");
    assert!(e.working > 0, "working must be more than 0");
    assert!(
        e.total >= 40 * MB && e.total <= 100 * MB,
        "total {} is outside 40 MB to 100 MB",
        e.total
    );
}

// The estimate grows in a straight line with the node count.
// 60^3 has 8 times the nodes of 30^3.
#[test]
fn estimate_scales_in_a_straight_line() {
    let _guard = lock();
    let s = Settings::default();
    let small = r1_problem(30, 2000, 500)
        .estimate_memory(&s)
        .expect("estimate must succeed")
        .total;
    let large = r1_problem(60, 2000, 500)
        .estimate_memory(&s)
        .expect("estimate must succeed")
        .total;
    let ratio = large as f64 / small as f64;
    assert!(
        (7.5..=8.5).contains(&ratio),
        "ratio {ratio} is outside 7.5 to 8.5"
    );
}

// The estimate call allocates less than 1 MB.
#[test]
fn estimate_call_allocates_less_than_1_mb() {
    let _guard = lock();
    let p = r1_problem(100, 2000, 500);
    let s = Settings::default();
    common::reset_peak();
    let start = common::current();
    let result = p.estimate_memory(&s);
    let used = common::peak().saturating_sub(start);
    assert!(result.is_ok(), "estimate must succeed");
    assert!(used < 1_000_000, "estimate allocated {used} bytes");
}

// The constructor on a 200^3 grid allocates less than 1 MB.
// The constructor must not build the grid tables.
#[test]
fn problem_new_does_not_allocate() {
    let _guard = lock();
    common::reset_peak();
    let start = common::current();
    let p = Problem::new([0.0; 3], [200; 3], [1.0; 3]);
    let used = common::peak().saturating_sub(start);
    assert!(p.is_ok(), "constructor must succeed");
    assert!(used < 1_000_000, "constructor allocated {used} bytes");
}

// A budget of 1 MB is too small for 100^3.
// `start` returns OverBudget before it allocates the vectors.
#[test]
fn start_refuses_a_budget_that_is_too_small() {
    let _guard = lock();
    let budget = 1_000_000u64;
    let p = r1_problem(100, 2000, 500);
    let settings = Settings {
        budget_bytes: Some(budget),
        ..Settings::default()
    };
    common::reset_peak();
    let start = common::current();
    let result = p.start(settings.clone());
    let used = common::peak().saturating_sub(start);
    let max_nodes = match result {
        Err(CoreError::OverBudget {
            estimate_bytes,
            budget_bytes,
            max_nodes_per_axis,
        }) => {
            assert_eq!(budget_bytes, budget);
            assert!(estimate_bytes > budget, "estimate must be over the budget");
            max_nodes_per_axis
        }
        Err(e) => panic!("expected OverBudget, got {e:?}"),
        Ok(_) => panic!("expected OverBudget, got a session"),
    };
    assert!(
        (used as u64) < budget,
        "start allocated {used} bytes; budget is {budget}"
    );
    assert!(max_nodes >= 1, "max_nodes_per_axis must be at least 1");

    // A cubic grid of that size must fit in the budget.
    // A grid has at least 3 nodes on each axis.
    if max_nodes >= 3 {
        let small = r1_problem(max_nodes, 2000, 500);
        let e = small
            .estimate_memory(&Settings::default())
            .expect("estimate must succeed");
        assert!(
            e.total <= budget,
            "estimate {} of {max_nodes}^3 is over the budget",
            e.total
        );
    }
}

// Pair bounds are counted. 100 000 pair bounds add at least
// 100 000 * 200 bytes to `constraints`.
#[test]
fn pair_bounds_are_counted_in_constraints() {
    let _guard = lock();
    let nodes = 30;
    let n = 100_000usize;
    let s = Settings::default();
    let mut p = r1_problem(nodes, 100, 100);
    let before = p.estimate_memory(&s).expect("estimate must succeed");

    let mut rng = Lcg(777);
    let a = points(&mut rng, n, nodes);
    let b = points(&mut rng, n, nodes);
    p.add_pair_bound(&a, &b, &vec![-1.0; n], &vec![1.0; n], &vec![1.0; n])
        .expect("pair bounds must be valid");
    let after = p.estimate_memory(&s).expect("estimate must succeed");

    let added = after.constraints.saturating_sub(before.constraints);
    assert!(
        added >= 100_000 * 200,
        "constraints grew by {added} bytes; at least 20000000 are needed"
    );
}
