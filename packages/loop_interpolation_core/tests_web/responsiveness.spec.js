// Story 3, story 4, SC-006, FR-010: the page stays responsive while the
// Worker solves a 60^3 problem. The tests also check cancel, the budget
// error, the refusal of a second solve, memoryBytes and progress messages.
//
// Message shapes that these tests expect (T041 and T042 must follow them):
//   in : { type:'init' }
//   in : { type:'solve', id, grid:{origin,counts,step}, regularisation:{...},
//          constraints:{value,gradient,normal,tangent,interface},
//          settings:{precision,tol,maxIter,budgetBytes?} }
//   in : { type:'estimate', id, grid, constraints, settings }
//        D6 and T042 do not fix the reply. This file chooses the reply
//        { type:'result', id, estimate:{total,vectors,constraints,working} }.
//   in : { type:'cancel', id? }
//   out: { type:'progress', id, ... }
//   out: { type:'result', id, field, memoryBytes, ... }
//   out: { type:'cancelled', id }
//   out: { type:'error', id, kind, ... }
//        kind 'OverBudget' also has estimate, budget, maxNodesPerAxis.
//        kind 'SolveRunning' answers a second solve while one runs.
// @ts-check
const { test, expect } = require('@playwright/test');

// This code runs in the page. It makes a flat random problem like
// tests/reference_problems.rs: an LCG and a plane field f = a . x.
// It sets window.makeSolve(id, n, rows, settings) -> { msg, transfer }.
function installGenerator() {
  window.makeSolve = (id, n, rows, settings) => {
    let state = 12345;
    // LCG: returns a number in [0, 1).
    const rnd = () => {
      state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
      return state / 4294967296;
    };
    const a = [0.7, -0.4, 0.5];
    const points = new Float64Array(3 * rows);
    const values = new Float64Array(rows);
    const gpoints = new Float64Array(3 * rows);
    const gvec = new Float64Array(3 * rows);
    for (let i = 0; i < rows; i++) {
      for (let k = 0; k < 3; k++) {
        points[3 * i + k] = 1 + rnd() * (n - 3);
        gpoints[3 * i + k] = 1 + rnd() * (n - 3);
        gvec[3 * i + k] = a[k];
      }
      values[i] = a[0] * points[3 * i] + a[1] * points[3 * i + 1] + a[2] * points[3 * i + 2];
    }
    const empty = () => new Float64Array(0);
    const msg = {
      type: 'solve',
      id,
      grid: {
        origin: new Float64Array([0, 0, 0]),
        counts: new Uint32Array([n, n, n]),
        step: new Float64Array([1, 1, 1]),
      },
      regularisation: { dxx: 1, dyy: 1, dzz: 1, dxy: 1, dxz: 1, dyz: 1, dx: 0.1, dy: 0.1, dz: 0.1 },
      constraints: {
        value: { points, values, weights: new Float64Array(rows).fill(1) },
        gradient: { points: gpoints, vectors: gvec, weights: new Float64Array(rows).fill(1) },
        normal: { points: empty(), vectors: empty(), weights: empty() },
        tangent: { points: empty(), vectors: empty(), weights: empty() },
        interface: { points: empty(), ids: empty(), weight: 1 },
      },
      settings: Object.assign({ precision: 'f64', tol: 1e-8, maxIter: 50000 }, settings || {}),
    };
    // The buffers go to the worker without a copy.
    const transfer = [];
    transfer.push(msg.grid.origin.buffer, msg.grid.counts.buffer, msg.grid.step.buffer);
    for (const fam of Object.values(msg.constraints)) {
      for (const k of Object.keys(fam)) {
        if (fam[k] instanceof Float64Array) transfer.push(fam[k].buffer);
      }
    }
    return { msg, transfer };
  };
  // Record the arrival time of each worker message.
  window.stamps = [];
  window.loopHarness.worker.addEventListener('message', (e) => {
    window.stamps.push({ t: performance.now(), type: e.data.type, id: e.data.id, kind: e.data.kind });
  });
  // Main-thread timer: the largest gap between two ticks, every 10 ms.
  window.gapMax = 0;
  window.gapLast = 0;
  window.startTimer = () => {
    window.gapMax = 0;
    window.gapLast = performance.now();
    window.timerHandle = setInterval(() => {
      const now = performance.now();
      window.gapMax = Math.max(window.gapMax, now - window.gapLast);
      window.gapLast = now;
    }, 10);
  };
  window.stopTimer = () => clearInterval(window.timerHandle);
}

// Post a generated solve message from the page.
async function postSolve(page, id, n, rows, settings) {
  await page.evaluate(([i, nn, r, s]) => {
    const { msg, transfer } = window.makeSolve(i, nn, r, s);
    window.loopHarness.post(msg, transfer);
  }, [id, n, rows, settings || null]);
}

// Wait for a worker message with a given id and one of the given types.
// A load error of the worker also ends the wait, so a missing worker fails fast.
function waitType(page, id, types, ms) {
  return page.evaluate(([i, t, m]) => window.loopHarness.waitFor(
    (x) => (x.id === i && t.includes(x.type)) || x.kind === 'WorkerLoadError', m).then((x) => ({
    type: x.type, id: x.id, kind: x.kind, estimate: x.estimate, budget: x.budget,
    maxNodesPerAxis: x.maxNodesPerAxis, memoryBytes: x.memoryBytes,
    hasField: !!x.field, message: x.message,
  })), [id, types, ms]);
}

test.beforeEach(async ({ page }) => {
  await page.goto('/tests_web/harness.html');
  await page.waitForFunction(() => window.harnessReady === true);
  await page.evaluate(installGenerator);
  await page.evaluate(() => window.loopHarness.post({ type: 'init' }));
});

const N = 60;
const ROWS = 2000;

test('(a) the main thread has no gap above 100 ms during a 60^3 solve', async ({ page }) => {
  await page.evaluate(() => window.startTimer());
  await postSolve(page, 1, N, ROWS);
  const done = await waitType(page, 1, ['result', 'error'], 110000);
  await page.evaluate(() => window.stopTimer());
  expect(done.type).toBe('result');
  const gap = await page.evaluate(() => window.gapMax);
  expect(gap).toBeLessThan(100);
});

test('(b) cancel after 1 s is answered within 1 s and a new solve works', async ({ page }) => {
  await postSolve(page, 1, N, ROWS);
  await page.waitForTimeout(1000);
  const t0 = await page.evaluate(() => {
    window.loopHarness.post({ type: 'cancel', id: 1 });
    return performance.now();
  });
  const cancelled = await waitType(page, 1, ['cancelled', 'result', 'error'], 5000);
  const t1 = await page.evaluate(() => {
    const s = window.stamps.filter((x) => x.id === 1 && x.type === 'cancelled');
    return s.length ? s[s.length - 1].t : -1;
  });
  expect(cancelled.type).toBe('cancelled');
  expect(t1 - t0).toBeLessThan(1000);

  // A new solve (small grid) must work after the cancel.
  await postSolve(page, 2, 12, 100);
  const done = await waitType(page, 2, ['result', 'error'], 60000);
  expect(done.type).toBe('result');
});

test('(c) a budget of 1 MB gives an OverBudget error', async ({ page }) => {
  await postSolve(page, 1, N, ROWS, { budgetBytes: 1e6 });
  const err = await waitType(page, 1, ['error', 'result'], 30000);
  expect(err.type).toBe('error');
  expect(err.kind).toBe('OverBudget');
  expect(err.estimate).toBeTruthy();
  expect(typeof err.estimate.total).toBe('number');
  expect(err.budget).toBe(1e6);
  expect(typeof err.maxNodesPerAxis).toBe('number');
  expect(err.maxNodesPerAxis).toBeGreaterThan(0);
});

test('estimate message gives the estimate in a result message', async ({ page }) => {
  await page.evaluate(([n, r]) => {
    const { msg, transfer } = window.makeSolve(1, n, r);
    msg.type = 'estimate';
    window.loopHarness.post(msg, transfer);
  }, [N, ROWS]);
  const res = await page.evaluate(() => window.loopHarness.waitFor(
    (x) => x.id === 1 && (x.type === 'result' || x.type === 'error'), 30000)
    .then((x) => ({ type: x.type, estimate: x.estimate })));
  expect(res.type).toBe('result');
  for (const key of ['total', 'vectors', 'constraints', 'working']) {
    expect(typeof res.estimate[key]).toBe('number');
  }
  expect(res.estimate.total).toBeGreaterThan(0);
});

test('(d) a second solve while one runs gives SolveRunning; the first ends normally', async ({ page }) => {
  await postSolve(page, 1, N, ROWS);
  // Wait until the first solve has started (first progress message).
  await waitType(page, 1, ['progress', 'result', 'error'], 30000);
  await postSolve(page, 2, 12, 100);
  const refused = await waitType(page, 2, ['error', 'result'], 30000);
  expect(refused.type).toBe('error');
  expect(refused.kind).toBe('SolveRunning');
  const first = await waitType(page, 1, ['result', 'error', 'cancelled'], 110000);
  expect(first.type).toBe('result');
});

test('(e) the result carries memoryBytes', async ({ page }) => {
  await postSolve(page, 1, 20, 300);
  const done = await waitType(page, 1, ['result', 'error'], 60000);
  expect(done.type).toBe('result');
  expect(typeof done.memoryBytes).toBe('number');
  expect(done.memoryBytes).toBeGreaterThan(0);
});

test('(f) progress messages arrive at least once each second', async ({ page }) => {
  const t0 = await page.evaluate(() => performance.now());
  await postSolve(page, 1, N, ROWS);
  const done = await waitType(page, 1, ['result', 'error'], 110000);
  expect(done.type).toBe('result');
  const times = await page.evaluate(() => window.stamps
    .filter((x) => x.id === 1 && (x.type === 'progress' || x.type === 'result'))
    .map((x) => x.t));
  expect(times.length).toBeGreaterThan(1);
  // Gaps: start to first progress, between progress messages, last to result.
  let prev = t0;
  let worst = 0;
  for (const t of times) {
    worst = Math.max(worst, t - prev);
    prev = t;
  }
  expect(worst).toBeLessThan(1000);
});
