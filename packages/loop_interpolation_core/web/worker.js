// Worker wrapper for loop_interpolation_core (task T042, plan D6).
//
// This file is a module Worker. It loads the WebAssembly package from
// ../pkg/ (made by wasm-bindgen --target web). The main thread never loads
// the module.
//
// Messages in:
//   { type:'init' }
//   { type:'estimate', id, grid, constraints, settings }
//   { type:'solve', id, grid, regularisation, constraints, settings }
//   { type:'cancel', id? }
// Messages out:
//   { type:'progress', id, stage, iteration, residual, maxViolation }
//   { type:'result', id, field, memoryBytes, done, converged, peakBytes,
//     iterations, residual, warnings }
//   { type:'result', id, estimate:{ total, vectors, constraints, working } }
//   { type:'cancelled', id }
//   { type:'error', id, kind, ... }
//
// Only one solve runs at a time. A second solve or estimate gets the error
// kind 'SolveRunning'.

// Time for one step() call, in milliseconds.
const TARGET_STEP_MS = 100;
// Smallest time between two progress messages, in milliseconds.
const PROGRESS_MIN_MS = 100;
// Limits of the iteration count for one step() call.
const MIN_ITERS = 1;
const MAX_ITERS = 5000;
const START_ITERS = 5;

// State of the module loader.
let wasmApi = null; // the imported module, or null when it is not loaded
let loadPromise = null; // Promise while the module loads
let needsReinit = false; // true after a WebAssembly trap
let loadCount = 0; // makes a new module URL for each re-instantiation

// State of the solve.
let solving = false;
let cancelFlag = false;
let currentId = undefined;

// Messages that arrive while the module loads.
const queue = [];
let draining = false;

// Load (or load again) the WebAssembly module.
function ensureLoaded() {
  if (wasmApi && !needsReinit) return Promise.resolve(wasmApi);
  if (loadPromise) return loadPromise;
  const suffix = loadCount === 0 ? '' : `?reload=${loadCount}`;
  loadCount += 1;
  loadPromise = (async () => {
    // A new query string gives a new module instance, with new memory.
    const mod = await import(`../pkg/loop_interpolation_core.js${suffix}`);
    await mod.default();
    wasmApi = mod;
    needsReinit = false;
    return mod;
  })();
  loadPromise = loadPromise.then(
    (mod) => {
      loadPromise = null;
      return mod;
    },
    (err) => {
      loadPromise = null;
      wasmApi = null;
      throw err;
    },
  );
  return loadPromise;
}

// Turn a thrown value into the fields of an error message.
function describeError(err) {
  const text = err && err.message !== undefined ? String(err.message) : String(err);
  // A WebAssembly trap (panic=abort) shows as 'unreachable'.
  const isTrap =
    (typeof WebAssembly !== 'undefined' &&
      WebAssembly.RuntimeError &&
      err instanceof WebAssembly.RuntimeError) ||
    /unreachable/i.test(text);
  if (isTrap) {
    return { trap: true, fields: { kind: 'Panic', message: text } };
  }
  try {
    const parsed = JSON.parse(text);
    if (parsed && typeof parsed === 'object' && typeof parsed.kind === 'string') {
      return { trap: false, fields: parsed };
    }
  } catch (_e) {
    // The text is not JSON. Use the fallback below.
  }
  return { trap: false, fields: { kind: 'InternalError', message: text } };
}

// Post an error message for a thrown value.
// "extra" can replace fields (for example the estimate of OverBudget).
function postError(id, err, extra) {
  const info = describeError(err);
  if (info.trap) {
    // The module state is not safe now. Load it again for the next message.
    needsReinit = true;
    wasmApi = null;
  }
  const fields = Object.assign({}, info.fields, extra || {});
  fields.type = 'error';
  fields.id = id;
  self.postMessage(fields);
}

function toF64(a) {
  if (a instanceof Float64Array) return a;
  return Float64Array.from(a || []);
}

// Add all constraints of a message to the problem.
function addConstraints(problem, c) {
  if (!c) return;
  const has = (fam, key) => fam && fam[key] && fam[key].length > 0;
  if (has(c.value, 'points')) {
    problem.addValue(toF64(c.value.points), toF64(c.value.values), toF64(c.value.weights));
  }
  if (has(c.gradient, 'points')) {
    problem.addGradient(toF64(c.gradient.points), toF64(c.gradient.vectors), toF64(c.gradient.weights));
  }
  if (has(c.normal, 'points')) {
    problem.addNormal(toF64(c.normal.points), toF64(c.normal.vectors), toF64(c.normal.weights));
  }
  if (has(c.tangent, 'points')) {
    problem.addTangent(toF64(c.tangent.points), toF64(c.tangent.vectors), toF64(c.tangent.weights));
  }
  if (has(c.interface, 'points')) {
    problem.addInterface(toF64(c.interface.points), toF64(c.interface.ids), Number(c.interface.weight));
  }
  if (has(c.valueBound, 'points')) {
    const b = c.valueBound;
    problem.addValueBound(toF64(b.points), toF64(b.lower), toF64(b.upper), toF64(b.weights));
  }
  if (has(c.pairBound, 'pointsA')) {
    const b = c.pairBound;
    problem.addPairBound(
      toF64(b.pointsA), toF64(b.pointsB), toF64(b.lower), toF64(b.upper), toF64(b.weights),
    );
  }
}

// Build the Problem from the buffers of a message.
function buildProblem(api, msg, withRegularisation) {
  const g = msg.grid;
  const origin = toF64(g.origin);
  const counts = g.counts instanceof Uint32Array ? g.counts : Uint32Array.from(g.counts);
  const step = toF64(g.step);
  const problem = new api.Problem(origin, counts, step);
  try {
    if (withRegularisation && msg.regularisation) {
      problem.setRegularisation(msg.regularisation);
    }
    addConstraints(problem, msg.constraints);
  } catch (err) {
    problem.free();
    throw err;
  }
  return problem;
}

// Copy the Estimate object of wasm-bindgen to a plain object.
function plainEstimate(est) {
  const out = {
    total: Number(est.total),
    vectors: Number(est.vectors),
    constraints: Number(est.constraints),
    working: Number(est.working),
  };
  if (typeof est.free === 'function') est.free();
  return out;
}

function delay0() {
  return new Promise((r) => setTimeout(r, 0));
}

// Handle an estimate message.
function handleEstimate(api, msg) {
  let problem = null;
  try {
    problem = buildProblem(api, msg, true);
    const estimate = plainEstimate(problem.estimate(msg.settings));
    self.postMessage({ type: 'result', id: msg.id, estimate });
  } catch (err) {
    postError(msg.id, err);
  } finally {
    if (problem && !needsReinit) {
      try {
        problem.free();
      } catch (_e) {
        // Ignore. The problem is already freed or the module is not valid.
      }
    }
  }
}

// Handle a solve message. The function runs until the solve ends.
async function handleSolve(api, msg) {
  const id = msg.id;
  const settings = msg.settings || {};
  solving = true;
  cancelFlag = false;
  currentId = id;
  let session = null;
  let problem = null;
  try {
    problem = buildProblem(api, msg, true);
    // Estimate first. start() consumes the problem, but the OverBudget
    // error needs the estimate as an object.
    const estimate = plainEstimate(problem.estimate(settings));
    try {
      session = problem.start(settings);
    } catch (err) {
      const info = describeError(err);
      if (info.fields.kind === 'OverBudget') {
        const extra = { estimate };
        if (typeof info.fields.budget !== 'number') {
          extra.budget = Number(settings.budgetBytes);
        }
        if (typeof info.fields.maxNodesPerAxis !== 'number') {
          extra.maxNodesPerAxis = 1;
        }
        postError(id, err, extra);
        return;
      }
      throw err;
    } finally {
      // start() consumed the problem. Do not use it again.
      problem = null;
    }

    let maxIters = START_ITERS;
    let lastProgress = 0;
    let last = null;
    let first = true;
    for (;;) {
      if (cancelFlag) {
        session.free();
        session = null;
        self.postMessage({ type: 'cancelled', id });
        return;
      }
      const t0 = performance.now();
      last = session.step(maxIters);
      const t1 = performance.now();

      // Scale maxIters so that one step takes about TARGET_STEP_MS.
      const dt = Math.max(t1 - t0, 1);
      const scale = Math.min(Math.max(TARGET_STEP_MS / dt, 0.25), 4);
      maxIters = Math.min(MAX_ITERS, Math.max(MIN_ITERS, Math.round(maxIters * scale)));

      if (last.done) break;

      // Always post after the first step. Later: not more often than
      // every PROGRESS_MIN_MS.
      if (first || t1 - lastProgress >= PROGRESS_MIN_MS) {
        self.postMessage({
          type: 'progress',
          id,
          stage: last.stage,
          iteration: last.iteration,
          residual: last.residual,
          maxViolation: last.violation,
        });
        lastProgress = t1;
        first = false;
      }
      // Let the worker handle a cancel message.
      await delay0();
    }

    const field = session.field();
    const memory = api.memoryBytes();
    const warnings = Array.from(last.warnings || []);
    const out = {
      type: 'result',
      id,
      field,
      memoryBytes: memory,
      done: true,
      converged: last.converged,
      peakBytes: last.peakBytes,
      iterations: last.iteration,
      residual: last.residual,
      maxViolation: last.violation,
      warnings,
    };
    session.free();
    session = null;
    self.postMessage(out, [field.buffer]);
  } catch (err) {
    postError(id, err);
  } finally {
    if (session && !needsReinit) {
      try {
        session.free();
      } catch (_e) {
        // Ignore. The session is already freed.
      }
    }
    if (problem && !needsReinit) {
      try {
        problem.free();
      } catch (_e) {
        // Ignore.
      }
    }
    solving = false;
    cancelFlag = false;
    currentId = undefined;
  }
}

// Handle one message when the module is loaded.
function dispatch(api, msg) {
  switch (msg.type) {
    case 'init':
      // The module is loaded. There is no reply.
      return;
    case 'cancel':
      if (solving && (msg.id === undefined || msg.id === currentId)) {
        cancelFlag = true;
      } else if (!solving) {
        self.postMessage({ type: 'cancelled', id: msg.id });
      }
      return;
    case 'estimate':
      if (solving) {
        self.postMessage({ type: 'error', id: msg.id, kind: 'SolveRunning' });
        return;
      }
      handleEstimate(api, msg);
      return;
    case 'solve':
      if (solving) {
        self.postMessage({ type: 'error', id: msg.id, kind: 'SolveRunning' });
        return;
      }
      // Set the flag now, so that the next message sees the running solve.
      solving = true;
      handleSolve(api, msg);
      return;
    default:
      self.postMessage({
        type: 'error',
        id: msg.id,
        kind: 'InvalidInput',
        message: `unknown message type: ${String(msg.type)}`,
      });
  }
}

// Take messages from the queue in order. The module loads first.
async function drain() {
  if (draining) return;
  draining = true;
  try {
    while (queue.length > 0) {
      // A cancel or a refusal does not need a new module. All others do.
      const msg = queue.shift();
      let api;
      try {
        api = await ensureLoaded();
      } catch (err) {
        if (msg.type !== 'init') {
          postError(msg.id, err, { kind: 'WorkerLoadError' });
        } else {
          postError(undefined, err, { kind: 'WorkerLoadError' });
        }
        continue;
      }
      try {
        dispatch(api, msg);
      } catch (err) {
        postError(msg.id, err);
      }
    }
  } finally {
    draining = false;
  }
}

self.onmessage = (event) => {
  const msg = event.data;
  if (!msg || typeof msg.type !== 'string') return;
  // Handle cancel at once when a solve runs. The module is already loaded.
  if (msg.type === 'cancel' && solving) {
    if (msg.id === undefined || msg.id === currentId) cancelFlag = true;
    return;
  }
  queue.push(msg);
  drain();
};

// Load the module at start. A load error shows at the first message.
ensureLoaded().catch(() => {});
