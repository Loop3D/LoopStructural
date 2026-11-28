# loop_interpolation_core: browser package

This folder contains the Web Worker that runs the solver. The solver is
WebAssembly. The worker file is `worker.js`. The type declarations are in
`index.d.ts`.

## Files

| File | Content |
|---|---|
| `web/worker.js` | The Worker. It loads the WebAssembly module and handles the messages. |
| `web/index.d.ts` | The TypeScript types of the messages and of the WebAssembly classes. |
| `pkg/` | The WebAssembly module and its JavaScript glue. The build makes it. |

The worker loads `../pkg/`. Serve `web/` and `pkg/` from the same origin.

## Build

1. Install the Rust target: `rustup target add wasm32-unknown-unknown`.
2. Install the wasm-bindgen CLI. Its version must equal the `wasm-bindgen`
   version in `Cargo.lock`.
3. From the package folder, run:

   ```sh
   bash scripts/build_wasm.sh
   ```

The script makes `pkg/`. You do not need `wasm-opt`. If `wasm-opt` is
installed, the script uses it. If not, the script skips it.

## Load the worker

```js
const worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
worker.postMessage({ type: 'init' });
```

The main thread does not load the WebAssembly module. Only the worker loads it.

## No cross-origin isolation

The package uses no threads and no `SharedArrayBuffer`. Your page does not
need the headers `Cross-Origin-Opener-Policy` and `Cross-Origin-Embedder-Policy`.

## Coordinates: transform them first

The solver takes grid coordinates only. It does not rotate or move the data
(spec FR-020). Your page must change all points and vectors into grid
coordinates before it sends them. Do this for any rotation of the model.

## Messages

Messages to the worker:

| `type` | Purpose | Reply |
|---|---|---|
| `init` | Start the module load. | None. |
| `estimate` | Get the memory estimate. | `result` with `estimate`. |
| `solve` | Start a solve. | `progress` messages, then `result` or `error`. |
| `cancel` | Stop a solve. `id` is optional. | `cancelled`. |

Messages from the worker:

| `type` | Content |
|---|---|
| `progress` | `id`, `stage`, `iteration`, `residual`, `maxViolation`. |
| `result` | Solve: `field`, `memoryBytes`, `done`, `converged`, `peakBytes`, `iterations`, `residual`, `maxViolation`, `warnings`. Estimate: `estimate`. |
| `cancelled` | `id`. |
| `error` | `id` and `kind`. Other fields depend on the kind. |

Only one solve runs at a time. A second `solve` or `estimate` gets the error
kind `SolveRunning`.

### Example: solve

```js
worker.postMessage(
  {
    type: 'solve',
    id: 1,
    grid: {
      origin: new Float64Array([0, 0, 0]),
      counts: new Uint32Array([20, 20, 20]),
      step: new Float64Array([1, 1, 1]),
    },
    regularisation: { dxx: 1, dyy: 1, dzz: 1 },
    constraints: {
      value: { points, values, weights },
    },
    settings: { precision: 'f64', tol: 1e-8, maxIter: 50000 },
  },
  [points.buffer, values.buffer, weights.buffer],
);
```

`points`, `values` and `weights` are `Float64Array` objects in grid
coordinates. The transfer list moves the buffers to the worker without a copy.
After the transfer, the page cannot use these arrays.

### Example: cancel

```js
worker.postMessage({ type: 'cancel', id: 1 });
```

The worker stops the solve at the next step. It then sends `cancelled` with
the same `id`.

### Error kinds

| `kind` | Meaning |
|---|---|
| `OverBudget` | The estimate is above `budgetBytes`. The message has `estimate`, `budget` and `maxNodesPerAxis`. |
| `SolveRunning` | A solve already runs. |
| `InvalidInput` | A value is wrong. The message has `constraint` and `reason`. |
| `GridTooSmall` | An axis has too few nodes. |
| `GridTooLarge` | The grid has too many nodes. |
| `NoConstraints` | The message has no constraints. |
| `OutOfMemory` | The solver could not get memory. |
| `Infeasible` | The constraints cannot all be met. The message has `worst`. |
| `PointOutsideGrid` | A point is outside the grid. |
| `Unsupported` | A setting is not supported. |
| `Panic` | The WebAssembly module stopped with a trap. The worker loads it again. |
| `WorkerLoadError` | The module did not load. |
| `InternalError` | Another error. The message has `message`. |

## Types

Import the types in TypeScript:

```ts
import type { WorkerInMessage, WorkerOutMessage } from './index.js';
```

The file `index.d.ts` also declares the WebAssembly classes `Problem`,
`Session`, `Estimate` and `StepResult`, and the function `memoryBytes`.
