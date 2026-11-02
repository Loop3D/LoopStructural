# Legacy golden replay checkpoint

The frozen, pre-performance baseline contains 21 model archives: 16 synthetic
fault/unconformity combinations, four relationship cases, and the supplied
plugin model. T005 committed and audited the archives and manifest in
`2697c1d1`. The manifest records `baseline_revision=ddc333ea`,
`prerequisite_fix_revision=e9b72293`, and
`prerequisite_parent_revision=b22517f9`; it also records the SHA-256 digest for
each reference archive. `test_golden.py` checks those digests before replay,
then rebuilds each model from saved fields, captures outputs, and compares
them with the immutable reference. It does not invoke the explicit writer or a
solver.

The replay covers feature values and gradients, stratigraphic ids, ordered
fault restorations, and surface summaries. Each model capture is an independent
pytest case with its own elapsed-time property; a fast catalogue test checks
that all 21 model names appear in the manifest. The separate edge tests cover
empty/outside/NaN value queries, region masking, reverse and disabled fault
application, zero displacement, and one lazy first-use build through
`GeologicalModel.evaluate_model`. The NaN-coordinate value query returns the
legacy value `0.0`, alongside NaN for the outside query. The single-fault
reverse control uses two fixed points and frozen values from `ddc333ea`. The
sensitivity check perturbs fault movement by 1% only inside a scoped monkeypatch
and requires the frozen
`f1_u0` full capture to fail with a model/output and max-difference/point
diagnostic. The comparator may report an affected feature output before the
restored-point output.

The initial replay exposed an incorrect all-NaN expectation; it was corrected
to the frozen legacy result. Final verification passed all 108 tests in 108.43
seconds (112.345966 seconds wall time), including the four D0 gradient
regressions. All 21 archive comparisons, the 1% sensitivity check, and both new
empty-input regressions passed. The 120-second acceptance gate passed. Pytest
reported one existing plugin zero-direction `RuntimeWarning`.

The test harness batches surface queries and reuses grid lookups. No production
evaluation behavior changed, and the reference archives remain unchanged from
baseline `ddc333ea`. Scoped sensitivity recording restores the patched method,
and the surface summary checks preserve the original mask, finite bounds, and
fallback behavior.

Faulted-gradient NaN/exact-fault-plane inputs are intentionally excluded
because the legacy refinement loop can continue refining without bound.
Finite faulted-gradient regression coverage remains in
`tests/unit/modelling/test_faulted_gradient_output.py`.
