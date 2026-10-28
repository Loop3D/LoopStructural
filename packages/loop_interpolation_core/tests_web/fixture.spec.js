// Story 1, scenarios 1 to 4 and FR-017: the Worker solves the 20^3 fixture
// in Chromium and the field agrees with the reference field of the Python
// matrix solver. The page has NO cross-origin isolation (FR-006).
// @ts-check
const { test, expect } = require('@playwright/test');

test('worker solves the 20^3 fixture and agrees with the reference', async ({ page }) => {
  const errors = [];
  page.on('pageerror', (e) => errors.push(String(e)));

  const response = await page.goto('/tests_web/harness.html');
  expect(response && response.ok()).toBeTruthy();

  // FR-006: no COOP/COEP headers, so the page is not cross-origin isolated.
  const headers = response ? response.headers() : {};
  expect(headers['cross-origin-opener-policy']).toBeUndefined();
  expect(headers['cross-origin-embedder-policy']).toBeUndefined();
  expect(await page.evaluate(() => crossOriginIsolated)).toBe(false);

  await page.waitForFunction(() => window.harnessReady === true);

  // Fetch the fixture in the page and give the bytes to runFixture.
  const outcome = await page.evaluate(async () => {
    const res = await fetch('/fixtures/fixture_20.bin');
    if (!res.ok) throw new Error('cannot fetch the fixture: ' + res.status);
    const bytes = await res.arrayBuffer();
    return await window.runFixture(bytes);
  });

  expect(outcome.length).toBe(20 * 20 * 20);
  expect(outcome.range).toBeGreaterThan(0);
  // Maximum difference: at most 1e-6 of the range of the reference field.
  expect(outcome.maxDiff).toBeLessThanOrEqual(1e-6 * outcome.range);
  expect(errors).toEqual([]);
});
