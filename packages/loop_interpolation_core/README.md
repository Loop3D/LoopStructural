# loop-interpolation-core

## Purpose

This package is the native core of the low-memory interpolator. It is a
Python extension module, built with maturin. It is named
`loop_interpolation_core` on import. The package solves the finite
difference interpolation problem without building the sparse constraint
rows in Python.

The package is optional. The `loop-interpolation` package works without it.

## Specification

The requirements and the design are in
[specs/001-low-memory-interpolator/spec.md](../../specs/001-low-memory-interpolator/spec.md).

## Forms

The crate builds in three forms:

- Rust library (rlib). Use it from other Rust crates.
- Python module (PyO3). It uses the abi3 stable ABI, so one wheel works for
  Python 3.9 and later.
- Browser module (wasm-bindgen). It runs in a Web Worker. The worker code is
  in `web/`.

## Feature flags

| Feature | Effect |
|---|---|
| `parallel` | Uses rayon threads. Off for cargo by default. On for the Python wheel. |
| `python` | Builds the PyO3 bindings. |
| `wasm` | Builds the wasm-bindgen bindings. |

With no feature, the code runs on one thread and has no dependency.

## Build commands

### Python module

```text
pip install maturin
maturin develop --release -m packages/loop_interpolation_core/Cargo.toml
```

The `pyproject.toml` sets the features `python` and `parallel`. If you pass
`--features`, you must list both: `--features python,parallel`. The option
`--features parallel` alone replaces the list. The build then fails with an
"install cffi" message.

### Browser module

```text
bash packages/loop_interpolation_core/scripts/build_wasm.sh
```

You need these tools:

- The `wasm32-unknown-unknown` Rust target.
- A `wasm-bindgen` CLI with the same version as the `wasm-bindgen` crate in
  `Cargo.lock`. The script stops if the versions differ.
- `wasm-opt` is optional. If it is missing, the script skips this step.

The script writes the output to `packages/loop_interpolation_core/pkg/`.

### Tests and lint

Run these commands from `packages/loop_interpolation_core`:

```text
cargo test --release
cargo test --release --features parallel

cargo clippy --all-targets -- -D warnings
cargo clippy --all-targets --features python -- -D warnings
cargo clippy --features wasm --target wasm32-unknown-unknown -- -D warnings
```

To run the browser tests:

```text
npm --prefix packages/loop_interpolation_core/tests_web ci
npm --prefix packages/loop_interpolation_core/tests_web test
```

## Contract

- Inputs are flat arrays. Points are stored as 3n values, one x, y, z set
  for each point.
- The grid uses grid coordinates only. The grid has no rotation. The caller
  transforms the coordinates.
- Node order is x fastest.
- Errors name the constraint that failed. For example, `value[3]`.
- A solve is a session. You step the session to advance the solve. To cancel
  it, drop the session.
- The memory estimate is available before the solve.
- Set `budget_bytes` to limit the memory of a solve.

## Differences from loop-rs

- Grid rotation is removed. The functions `rotation`, `to_local` and
  `direction_to_local` are not in this crate. The `serde` feature is also
  removed.
- Solve sessions replace the closure-style solvers. A session has step and
  cancel.
- Memory tracking is new. The crate uses `try_reserve`, a byte budget and a
  memory estimate.
- Curvature rows are added only at nodes that are interior on all three
  axes. The border rows are centred second derivatives `[1, -2, 1]` on
  layers 1 and n-2. The gradient strike and dip basis is the Python
  strike and dip basis. So the field equals the field of the Python matrix
  solver.
- Only the ADMM solver is in this crate. There is no active-set solver and
  no hybrid solver. A feasibility check runs before the solve. The solver
  stops when the violation is at most 1e-4 of the field range.
- Single precision (`f32`) is available. The residual check runs in double
  precision.

## Measured values

- Memory: 63 bytes per node for R1, a 100 x 100 x 100 grid with 2000 value
  rows and 500 gradient rows. The value is the peak heap size divided by
  the number of nodes. It does not include the copy of the field that goes
  to JavaScript, which takes 8 bytes per node.
- Comparison with the Python matrix solver on a 100 x 100 x 100 grid: the
  native solver uses 84 MB where the matrix solver uses 2445 MB (29 times
  less). On one thread it is 10 to 30 percent slower. With the `parallel`
  feature (the default of the Python wheel) it is about 3 times faster. The
  numbers are in `specs/001-low-memory-interpolator/measurements.md`.
- Browser size: the gzip size of the `.wasm` file is about 92 KB. The
  limit is 2 MB.

## Source

Source: loop-rs crates/interpolation-fd, commit bff51342da8ebff7ca3c87e07c1d4fbd31279c8b, branch lachlan/drillhole
