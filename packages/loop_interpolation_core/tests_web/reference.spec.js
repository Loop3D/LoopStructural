// Browser reference run for R1 (SC-001, SC-003, SC-006, SC-007).
// The spec is skipped unless the environment variable RUN_REFERENCE=1 is set.
//
// Keep in step with tests/reference_problems.rs: the function `build` and the
// LCG there are the source of the formulas below. If one file changes,
// change the other file in the same way.
//
// R1: grid 100^3, 2000 value rows and 500 gradient rows, plane field
// f = x + 0.5 y + 0.25 z, weights 1, tol 1e-8, max_iter 5000, precision f64.
// The Rust test `r1_write_field` writes the reference field to
// target/r1_field.f64 (little-endian f64, x fastest). The static server
// serves it as /target/r1_field.f64.
//
// R2 (200 active pair bounds) is NOT in this file. The pair-bound path needs
// the free field first. It is not needed for T052.
//
// SC-007 (60 s) is report only. The release check uses the reference machine.
// @ts-check
const { test, expect } = require('@playwright/test');
const fs = require('fs');
const path = require('path');

const FIELD_FILE = path.join(__dirname, '..', 'target', 'r1_field.f64');
const N = 100;
const N_VALUE = 2000;
const N_NORMAL = 500;

// This code runs in the page. It makes the R1 solve message.
// The LCG uses BigInt, to get the same 64-bit wrap as Rust.
function installR1Generator() {
  window.makeR1 = (id, n, nValue, nNormal) => {
    let state = 12345n;
    // Same as Lcg::next in reference_problems.rs.
    const rnd = () => {
      state = BigInt.asUintN(64, state * 6364136223846793005n + 1442695040888963407n);
      return Number(state >> 11n) / 9007199254740992; // 2^53
    };
    const extent = n - 1;
    const points = new Float64Array(3 * nValue);
    const values = new Float64Array(nValue);
    for (let i = 0; i < nValue; i++) {
      const q = [rnd() * extent, rnd() * extent, rnd() * extent];
      values[i] = q[0] + 0.5 * q[1] + 0.25 * q[2];
      points.set(q, 3 * i);
    }
    const gpoints = new Float64Array(3 * nNormal);
    const gvec = new Float64Array(3 * nNormal);
    for (let i = 0; i < nNormal; i++) {
      gpoints.set([rnd() * extent, rnd() * extent, rnd() * extent], 3 * i);
      gvec.set([1.0, 0.5, 0.25], 3 * i);
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
      // The Rust Problem uses weight 1.0 for all regularisation terms.
      regularisation: { dxx: 1, dyy: 1, dzz: 1, dxy: 1, dxz: 1, dyz: 1, dx: 1, dy: 1, dz: 1 },
      constraints: {
        value: { points, values, weights: new Float64Array(nValue).fill(1) },
        gradient: { points: gpoints, vectors: gvec, weights: new Float64Array(nNormal).fill(1) },
        normal: { points: empty(), vectors: empty(), weights: empty() },
        tangent: { points: empty(), vectors: empty(), weights: empty() },
        interface: { points: empty(), ids: empty(), weight: 1 },
      },
      settings: { precision: 'f64', tol: 1e-8, maxIter: 5000 },
    };
    const transfer = [msg.grid.origin.buffer, msg.grid.counts.buffer, msg.grid.step.buffer];
    for (const fam of Object.values(msg.constraints)) {
      for (const k of Object.keys(fam)) {
        if (fam[k] instanceof Float64Array) transfer.push(fam[k].buffer);
      }
    }
    // The first value-row point is a check of the LCG (it must be in [0, 99]^3).
    return { msg, transfer, firstPoint: Array.from(points.slice(0, 3)) };
  };
  window.stamps = [];
  window.loopHarness.worker.addEventListener('message', (e) => {
    window.stamps.push({ t: performance.now(), type: e.data.type, id: e.data.id });
  });
}

test.describe('R1 reference run', () => {
  test.skip(!process.env.RUN_REFERENCE, 'set RUN_REFERENCE=1 to run the reference problems');

  test('R1 100^3 solve: memory, field difference, progress', async ({ page }) => {
    test.setTimeout(15 * 60 * 1000);

    if (!fs.existsSync(FIELD_FILE)) {
      throw new Error(
        `Reference field file is missing: ${FIELD_FILE}. Run: cargo test --manifest-path ` +
        'packages/loop_interpolation_core/Cargo.toml --release --test reference_problems ' +
        'r1_write_field -- --ignored --nocapture',
      );
    }
    expect(fs.statSync(FIELD_FILE).size).toBe(8 * N * N * N);

    await page.goto('/tests_web/harness.html');
    await page.waitForFunction(() => window.harnessReady === true);
    await page.evaluate(installR1Generator);
    await page.evaluate(() => window.loopHarness.post({ type: 'init' }));

    const start = await page.evaluate(([n, nv, ng]) => {
      const { msg, transfer, firstPoint } = window.makeR1(1, n, nv, ng);
      window.firstPoint = firstPoint;
      const t = performance.now();
      window.loopHarness.post(msg, transfer);
      return { t, firstPoint };
    }, [N, N_VALUE, N_NORMAL]);
    for (const c of start.firstPoint) {
      expect(c).toBeGreaterThanOrEqual(0);
      expect(c).toBeLessThanOrEqual(N - 1);
    }

    const done = await page.evaluate((ms) => window.loopHarness.waitFor(
      (x) => x.id === 1 && (x.type === 'result' || x.type === 'error'), ms).then((x) => {
      window.fieldResult = x.field;
      return {
        type: x.type, kind: x.kind, message: x.message, memoryBytes: x.memoryBytes,
        iterations: x.iterations, converged: x.converged, residual: x.residual,
        length: x.field ? x.field.length : 0, end: performance.now(),
      };
    }), 14 * 60 * 1000);
    expect(done.type, JSON.stringify(done)).toBe('result');
    const seconds = (done.end - start.t) / 1000;
    expect(done.length).toBe(N * N * N);

    // Compare with the reference field. The page reads it from the static server.
    const cmp = await page.evaluate(async () => {
      const buf = await (await fetch('/target/r1_field.f64')).arrayBuffer();
      const r = new Float64Array(buf);
      const f = window.fieldResult;
      let maxDiff = 0;
      let lo = Infinity;
      let hi = -Infinity;
      for (let i = 0; i < r.length; i++) {
        maxDiff = Math.max(maxDiff, Math.abs(f[i] - r[i]));
        lo = Math.min(lo, r[i]);
        hi = Math.max(hi, r[i]);
      }
      return { maxDiff, range: hi - lo };
    });
    const ratio = cmp.maxDiff / cmp.range;

    const times = await page.evaluate(() => window.stamps
      .filter((x) => x.id === 1 && (x.type === 'progress' || x.type === 'result'))
      .map((x) => x.t));
    let prev = start.t;
    let worst = 0;
    for (const t of times) {
      worst = Math.max(worst, t - prev);
      prev = t;
    }

    console.log(`R1 memoryBytes=${done.memoryBytes}`);
    console.log(`R1 time_s=${seconds.toFixed(1)} (SC-007 target 60 s, report only)`);
    console.log(`R1 iterations=${done.iterations} converged=${done.converged}`);
    console.log(`R1 max_diff=${cmp.maxDiff} range=${cmp.range} ratio=${ratio}`);
    console.log(`R1 progress messages=${times.length} worst_gap_ms=${worst.toFixed(0)}`);

    expect(done.converged).toBe(true);
    expect(done.memoryBytes).toBeLessThanOrEqual(100_000_000); // SC-001
    expect(ratio).toBeLessThanOrEqual(1e-6); // SC-003
    expect(times.length).toBeGreaterThan(1);
    expect(worst).toBeLessThan(1000); // SC-006
  });
});
