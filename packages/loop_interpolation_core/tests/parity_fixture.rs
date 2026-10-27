//! Parity test (T026): the Rust solver must agree with the matrix solver.
//!
//! The file `fixtures/fixture_20.bin` has a 20^3 problem and the expected
//! field of the Python matrix solver. The format is in plan.md, section D8.
//! This test reads the file with the standard library only.

use loop_interpolation_core::{
    BorderRegularisationWeights, Problem, RegularisationWeights, Settings,
};

/// The largest allowed difference, as a part of the range of the field.
const MAX_REL_DIFF: f64 = 1e-6;

/// The ridge factor that the ported interpolator uses (`ridge_factor` default).
/// `Problem` cannot set the ridge, so the file must have this value.
const CORE_RIDGE: f64 = 1e-8;

/// The data of the fixture file, in the order of the file.
struct Fixture {
    origin: Vec<f64>,
    counts: Vec<f64>,
    step: Vec<f64>,
    reg: Vec<f64>,
    border: Vec<f64>,
    ridge: Vec<f64>,
    value_points: Vec<f64>,
    value_values: Vec<f64>,
    value_weights: Vec<f64>,
    gradient_points: Vec<f64>,
    gradient_vectors: Vec<f64>,
    gradient_weights: Vec<f64>,
    normal_points: Vec<f64>,
    normal_vectors: Vec<f64>,
    normal_weights: Vec<f64>,
    tangent_points: Vec<f64>,
    tangent_vectors: Vec<f64>,
    tangent_weights: Vec<f64>,
    interface_points: Vec<f64>,
    interface_ids: Vec<f64>,
    interface_weight: Vec<f64>,
    expected_field: Vec<f64>,
}

/// Read the next block: a `u32` length and then that number of `f64` values.
fn read_block(bytes: &[u8], pos: &mut usize) -> Vec<f64> {
    let len_end = *pos + 4;
    assert!(
        len_end <= bytes.len(),
        "fixture is too short (block length)"
    );
    let n = u32::from_le_bytes(bytes[*pos..len_end].try_into().unwrap()) as usize;
    *pos = len_end;
    let data_end = *pos + 8 * n;
    assert!(data_end <= bytes.len(), "fixture is too short (block data)");
    let values = bytes[*pos..data_end]
        .chunks_exact(8)
        .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
        .collect();
    *pos = data_end;
    values
}

fn load_fixture() -> Fixture {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/fixtures/fixture_20.bin");
    let bytes = std::fs::read(path).expect("cannot read fixture_20.bin");
    assert!(bytes.len() >= 4, "fixture is too short (magic)");
    assert_eq!(&bytes[0..4], b"LIF1", "wrong magic in fixture");
    let mut pos = 4usize;
    let mut next = || read_block(&bytes, &mut pos);
    let fx = Fixture {
        origin: next(),
        counts: next(),
        step: next(),
        reg: next(),
        border: next(),
        ridge: next(),
        value_points: next(),
        value_values: next(),
        value_weights: next(),
        gradient_points: next(),
        gradient_vectors: next(),
        gradient_weights: next(),
        normal_points: next(),
        normal_vectors: next(),
        normal_weights: next(),
        tangent_points: next(),
        tangent_vectors: next(),
        tangent_weights: next(),
        interface_points: next(),
        interface_ids: next(),
        interface_weight: next(),
        expected_field: next(),
    };
    assert_eq!(pos, bytes.len(), "fixture has extra bytes at the end");
    fx
}

/// Make a weight vector with one value for each of `n` points.
/// If the file has one value only, repeat it.
fn per_point(weights: &[f64], n: usize) -> Vec<f64> {
    if weights.len() == n {
        weights.to_vec()
    } else if weights.len() == 1 {
        vec![weights[0]; n]
    } else {
        panic!("weights have {} values; {} are needed", weights.len(), n);
    }
}

#[test]
fn rust_solver_agrees_with_matrix_solver_fixture() {
    let fx = load_fixture();

    assert_eq!(fx.origin.len(), 3);
    assert_eq!(fx.counts.len(), 3);
    assert_eq!(fx.step.len(), 3);
    assert_eq!(fx.reg.len(), 6);
    assert_eq!(fx.border.len(), 3);
    assert_eq!(fx.ridge.len(), 1);
    assert_eq!(fx.interface_weight.len(), 1);
    // The ridge of the file must be the ridge of the ported interpolator.
    assert_eq!(
        fx.ridge[0], CORE_RIDGE,
        "the ridge of the file is not the ridge of the core"
    );

    let counts = [
        fx.counts[0] as usize,
        fx.counts[1] as usize,
        fx.counts[2] as usize,
    ];
    let n_nodes = counts[0] * counts[1] * counts[2];
    assert_eq!(fx.expected_field.len(), n_nodes);

    let mut problem = Problem::new(
        [fx.origin[0], fx.origin[1], fx.origin[2]],
        counts,
        [fx.step[0], fx.step[1], fx.step[2]],
    )
    .expect("cannot make the problem");

    // The weights of the file are in core meaning. Pass them as they are.
    problem.set_regularisation(RegularisationWeights {
        dxx: fx.reg[0],
        dyy: fx.reg[1],
        dzz: fx.reg[2],
        dxy: fx.reg[3],
        dxz: fx.reg[4],
        dyz: fx.reg[5],
    });
    problem.set_border_regularisation(BorderRegularisationWeights {
        dx: fx.border[0],
        dy: fx.border[1],
        dz: fx.border[2],
    });

    let n_value = fx.value_points.len() / 3;
    problem
        .add_value(
            &fx.value_points,
            &fx.value_values,
            &per_point(&fx.value_weights, n_value),
        )
        .expect("value constraints");

    let n_gradient = fx.gradient_points.len() / 3;
    problem
        .add_gradient(
            &fx.gradient_points,
            &fx.gradient_vectors,
            &per_point(&fx.gradient_weights, n_gradient),
        )
        .expect("gradient constraints");

    let n_normal = fx.normal_points.len() / 3;
    problem
        .add_normal(
            &fx.normal_points,
            &fx.normal_vectors,
            &per_point(&fx.normal_weights, n_normal),
        )
        .expect("normal constraints");

    let n_tangent = fx.tangent_points.len() / 3;
    problem
        .add_tangent(
            &fx.tangent_points,
            &fx.tangent_vectors,
            &per_point(&fx.tangent_weights, n_tangent),
        )
        .expect("tangent constraints");

    problem
        .add_interface(
            &fx.interface_points,
            &fx.interface_ids,
            fx.interface_weight[0],
        )
        .expect("interface constraints");

    let settings = Settings {
        tol: 1e-10,
        max_iter: 50_000,
        ..Settings::default()
    };
    let (field, report) = problem.solve(settings).expect("solve failed");
    assert!(report.converged, "the solver did not converge: {report:?}");

    let got = field.values_f64();
    assert_eq!(got.len(), n_nodes);

    let lo = fx
        .expected_field
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let hi = fx
        .expected_field
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let range = hi - lo;
    assert!(range > 0.0, "the expected field has no range");

    let (worst_index, max_diff) = got
        .iter()
        .zip(fx.expected_field.iter())
        .map(|(a, b)| (a - b).abs())
        .enumerate()
        .fold(
            (0usize, 0.0f64),
            |acc, (i, d)| if d > acc.1 { (i, d) } else { acc },
        );
    let rel = max_diff / range;
    assert!(
        rel <= MAX_REL_DIFF,
        "parity failed: max difference {max_diff:e} at node {worst_index} is {rel:e} of the range {range:e} (limit {MAX_REL_DIFF:e}); check the mapping"
    );
}
