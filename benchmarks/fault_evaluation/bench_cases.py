"""pytest-benchmark suite for fault evaluation.

Only the runner selects this file, by explicit path. The runner passes the
pass kind, the case list and the round count in environment variables.
"""

# ruff: noqa: ANN401  (the wrapped objects are untyped model and fixture objects)

from __future__ import annotations

import json
import os
import time
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pytest

from benchmarks.fault_evaluation.instrumentation import EvaluationInstrumentation, RssSampler
from LoopStructural.modelling.features._geological_feature import GeologicalFeature
from tests.fault_evaluation.models import build_model

PASS_ENV = "FAULT_EVAL_PASS"
CASES_ENV = "FAULT_EVAL_CASES"
ROUNDS_ENV = "FAULT_EVAL_ROUNDS"
CHECKPOINT_ENV = "FAULT_EVAL_CHECKPOINT"
POINT_SEED = 1


def _groups(model: Any) -> list[Any]:
    """Return the foliation groups in creation order. The last one is the oldest."""
    groups = [feature for feature in model.features if type(feature) is GeologicalFeature]
    if not groups:
        raise ValueError("The model has no foliation group")
    return groups


_OPERATIONS: dict[str, Callable[[Any, np.ndarray], Any]] = {
    "oldest_value": lambda model, points: _groups(model)[-1].evaluate_value(points),
    "all_group_values": lambda model, points: [
        group.evaluate_value(points) for group in _groups(model)
    ],
    "model_ids": lambda model, points: model.evaluate_model(points),
    "oldest_gradient": lambda model, points: _groups(model)[-1].evaluate_gradient(points),
    "model_gradient": lambda model, points: model.evaluate_model_gradient(points),
}


def run_case(
    benchmark: Any,
    model: Any,
    points: np.ndarray,
    case: dict[str, Any],
    *,
    rounds: int,
    after_setup: Callable[[], None] | None = None,
) -> None:
    """Time one case. Cache clearing and warm priming stay outside the timed call.

    Every round clears the cache first. A warm case then evaluates once to fill
    the cache. A model with no cache is not primed. ``after_setup`` runs last, just before the timed call.
    """
    operation = _OPERATIONS[case["operation"]]
    cache = getattr(model, "evaluation_cache", None)
    warm = case["state"] == "warm"

    def call() -> Any:
        return operation(model, points)

    def setup() -> None:
        if cache is not None:
            cache.clear()
        # A model with no cache has nothing to fill, so a priming call does no work.
        if warm and cache is not None:
            call()
        if after_setup is not None:
            after_setup()

    benchmark.extra_info.update(case)
    # Warmup rounds would turn a cold case into a cache hit.
    benchmark.pedantic(call, setup=setup, rounds=rounds, warmup_rounds=0)


def _load_cases() -> list[dict[str, Any]]:
    path = os.environ.get(CASES_ENV)
    if not path:
        return []
    return json.loads(Path(path).read_text(encoding="utf8"))


@lru_cache(maxsize=1)
def _model(name: str) -> Any:
    """Build a model once. The cases of one model are neighbours in the list."""
    return build_model(name)


def _points(model: Any, count: int) -> np.ndarray:
    # Same recipe as tests.fault_evaluation.models.sample_points, without a rebuild.
    box = model.bounding_box
    rng = np.random.default_rng(POINT_SEED)
    return rng.random((count, 3)) * (box.maximum - box.origin) + box.origin


def write_checkpoint(benchmark: Any, case: dict[str, Any], elapsed_seconds: float) -> None:
    """Append one JSON line for a finished case to the checkpoint file.

    Write nothing if the checkpoint variable is absent. Flush and fsync the line,
    so a stop of the process cannot lose a finished case.
    """
    path = os.environ.get(CHECKPOINT_ENV)
    if not path:
        return
    # benchmark.stats is a Metadata object. Its .stats attribute holds the numbers.
    stats = benchmark.stats.stats
    line = {
        "case_id": case["case_id"],
        "stats": {
            "rounds": int(stats.rounds),
            "data": [float(value) for value in stats.data],
            "median": float(stats.median),
        },
        "extra_info": dict(benchmark.extra_info),
        "elapsed_seconds": float(elapsed_seconds),
    }
    with open(path, "a", encoding="utf8") as handle:
        handle.write(json.dumps(line, default=str) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


@pytest.mark.parametrize("case", _load_cases(), ids=lambda case: case["case_id"])
def test_fault_evaluation(benchmark: Any, case: dict[str, Any]) -> None:
    """Measure one case in the pass that the runner selected."""
    started = time.perf_counter()
    pass_kind = os.environ[PASS_ENV]
    rounds = int(os.environ[ROUNDS_ENV])
    model = _model(case["model"])
    points = _points(model, case["points"])
    if pass_kind == "timing":
        run_case(benchmark, model, points, case, rounds=rounds)
    elif pass_kind == "counters":
        with EvaluationInstrumentation() as counters:
            run_case(benchmark, model, points, case, rounds=1, after_setup=counters.reset)
        benchmark.extra_info["counters"] = counters.snapshot()
    elif pass_kind == "rss":
        sampler = RssSampler()
        try:
            run_case(benchmark, model, points, case, rounds=1, after_setup=sampler.start)
        finally:
            benchmark.extra_info["peak_extra_rss_bytes"] = sampler.stop()
    else:
        raise ValueError(f"Unknown pass kind {pass_kind!r}")
    write_checkpoint(benchmark, case, time.perf_counter() - started)
