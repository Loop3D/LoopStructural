// Playwright configuration for the browser tests of loop_interpolation_core.
// Chromium only. The server is a plain static server (Python standard
// library). It sends NO cross-origin isolation headers (FR-006): the tests
// must pass without SharedArrayBuffer and without threads.
// @ts-check
const { defineConfig } = require('@playwright/test');

const PORT = 4173;

module.exports = defineConfig({
  testDir: '.',
  testMatch: '*.spec.js',
  timeout: 120000,
  fullyParallel: false,
  workers: 1,
  reporter: 'list',
  use: {
    baseURL: `http://127.0.0.1:${PORT}`,
  },
  projects: [{ name: 'chromium', use: { browserName: 'chromium' } }],
  webServer: {
    // The served root is the package root (the folder above tests_web).
    // The command runs in the folder of this file.
    command: `python -m http.server ${PORT} --bind 127.0.0.1 --directory ..`,
    url: `http://127.0.0.1:${PORT}/tests_web/harness.html`,
    reuseExistingServer: false,
    timeout: 30000,
  },
});
