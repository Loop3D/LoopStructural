"""The 94-fault model of spec 007 (SC-001, SC-004): build it, time it, check its results.

Usage::

    python -m benchmarks.fault_evaluation.many_faults --output RESULT.json \
        [--save-reference REF.npz | --compare REF.npz] [--points 500 2000 5000]

The model is deterministic (seed 100): 54 isolated faults in the east of a
10 x 10 x 3 km box, 40 long faults of two crossing strike families in the west
(each cut by at most its 3 youngest crossing faults; chains up to 20 deep), and
two stratigraphic groups with one unconformity. The build solves every
interpolator and takes about 7 minutes, so this is a local benchmark only: it
is not part of the standard suite and no test builds the model.

The result holds, for each point count, the ``evaluate_model`` time and fault
calls, and the mean time of one isolated and one cutting fault (outermost
``apply_to_points`` calls only). ``--save-reference`` stores the model ids and
the values of each group; ``--compare`` checks them bit for bit (SC-001).
"""

# ruff: noqa: ANN001, ANN002, ANN003, ANN201, ANN202, ANN204, D103  (a measurement script)

from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd

from LoopStructural import GeologicalModel
from LoopStructural.modelling.core.stratigraphic_column import StratigraphicColumn
from LoopStructural.modelling.features.fault._fault_segment import FaultSegment

SEED = 100
ORIGIN = np.array([0.0, 0.0, -3000.0])
MAXIMUM = np.array([10000.0, 10000.0, 0.0])
MAX_CUTTERS = 3  # each fault is cut by at most its 3 youngest crossing faults
GROUPS = ("group_1", "group_0")


def _unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def _plane(dip_direction, dip):
    dd, d = np.radians(dip_direction), np.radians(dip)
    normal = [np.sin(d) * np.sin(dd), np.sin(d) * np.cos(dd), np.cos(d)]
    down_dip = [np.cos(d) * np.sin(dd), np.cos(d) * np.cos(dd), -np.sin(d)]
    strike = [np.sin(dd - np.pi / 2), np.cos(dd - np.pi / 2), 0.0]
    return _unit(normal), _unit(down_dip), _unit(strike)


def _traces_cross(a, b, buffer):
    """Two traces (centre, strike, half length) cross or come within ``buffer``."""
    pa, da, la = a
    pb, db, lb = b
    samples_a = pa + np.linspace(-la, la, 60)[:, None] * da
    samples_b = pb + np.linspace(-lb, lb, 60)[:, None] * db
    dist = np.linalg.norm(samples_a[:, None, :] - samples_b[None, :, :], axis=2)
    return bool(dist.min() < buffer)


def fault_specs(n_isolated=60, n_long=40):
    """Return the faults, youngest first. ``cut_by`` holds the indices of the cutting faults."""
    rng = np.random.default_rng(SEED)
    specs = []
    for i in range(n_long):
        dip_direction = (120.0 if i % 2 == 0 else 30.0) + rng.normal(0, 8)
        dip = rng.uniform(60, 75)
        centre = np.array([rng.uniform(1800, 3200), rng.uniform(2000, 8000), -1500.0])
        specs.append(
            dict(
                kind="cutting",
                dip_direction=dip_direction,
                dip=dip,
                centre=centre,
                half_length=rng.uniform(1200, 1900),
                displacement=rng.uniform(40, 120),
            )
        )
    # Isolated faults on a jittered 6 x 13 grid of 700 m cells from x = 5600.
    cells = [(i, j) for i in range(6) for j in range(13)]
    rng.shuffle(cells)
    for i, j in cells:
        if sum(s["kind"] == "isolated" for s in specs) == n_isolated:
            break
        centre = np.array(
            [5600 + i * 700 + rng.uniform(-80, 80), 350 + j * 700 + rng.uniform(-80, 80), -1500.0]
        )
        specs.append(
            dict(
                kind="isolated",
                dip_direction=rng.uniform(0, 360),
                dip=rng.uniform(55, 85),
                centre=centre,
                half_length=rng.uniform(150, 250),
                displacement=rng.uniform(10, 30),
            )
        )
    for s in specs:
        s["normal"], s["slip"], s["strike"] = _plane(s["dip_direction"], s["dip"])
        s["major_axis"] = s["half_length"]
        s["intermediate_axis"] = 1500.0 if s["kind"] == "cutting" else 600.0
        s["minor_axis"] = s["displacement"] * 6
    traces = [(s["centre"][:2], s["strike"][:2], s["half_length"]) for s in specs]
    # Drop the isolated faults that touch any other fault.
    keep = [
        k
        for k, s in enumerate(specs)
        if s["kind"] != "isolated"
        or not any(
            _traces_cross(traces[k], traces[m], s["minor_axis"] + specs[m]["minor_axis"])
            for m in range(len(specs))
            if m != k
        )
    ]
    specs = [specs[k] for k in keep]
    traces = [traces[k] for k in keep]
    # Index 0 is the youngest. A younger cutting fault that crosses an older one cuts it.
    for k, s in enumerate(specs):
        s["name"] = f"fault_{k:03d}"
        s["cut_by"] = []
        if s["kind"] == "cutting":
            s["cut_by"] = [
                m
                for m in range(k)
                if specs[m]["kind"] == "cutting" and _traces_cross(traces[k], traces[m], 200.0)
            ][-MAX_CUTTERS:]
    return specs


def _strat_data():
    rng = np.random.default_rng(SEED + 1)
    rows = []
    for name, base, tilt in (("group_0", -2200.0, 0.02), ("group_1", -900.0, -0.03)):
        for value in (0.0, 200.0, 400.0):
            xy = rng.uniform(0, 10000, (80, 2))
            z = base + value + tilt * (xy[:, 0] - 5000)
            rows += [
                {"X": x, "Y": y, "Z": zz, "val": value, "feature_name": name}
                for (x, y), zz in zip(xy, z)
            ]
        rows.append(
            {"X": 5000.0, "Y": 5000.0, "Z": base, "gx": -tilt, "gy": 0.0, "gz": 1.0, "feature_name": name}
        )
    return pd.DataFrame(rows)


def _fault_data(specs):
    """Three trace points (value 0) and the normal at the centre of each fault."""
    rows = []
    for s in specs:
        for t in (-0.5, 0.0, 0.5):
            x, y, z = s["centre"] + t * s["half_length"] * s["strike"]
            rows.append({"X": x, "Y": y, "Z": z, "val": 0.0, "feature_name": s["name"]})
        x, y, z = s["centre"]
        gx, gy, gz = s["normal"]
        rows.append({"X": x, "Y": y, "Z": z, "gx": gx, "gy": gy, "gz": gz, "feature_name": s["name"]})
    return pd.DataFrame(rows)


def _set_column(model):
    column = StratigraphicColumn()
    column.clear(basement=False)
    intervals = []
    for k in (0, 1):  # oldest first
        for name, low, high in (
            (f"unit_{k}_high", 300.0, 1e5),
            (f"unit_{k}_mid", 100.0, 300.0),
            (f"unit_{k}_low", -1e5, 100.0),
        ):
            intervals.append((column.add_unit(name, thickness=high - low), low, high))
        if k == 0:
            column.add_unconformity("unconformity_1")
    for unit, low, high in intervals:
        unit.min_value, unit.max_value = low, high
    column.group_mapping.update({"Group_0": "group_1", "Group_1": "group_0"})
    model.stratigraphic_column = column


def build(nelements=3000):
    """Build and solve the model. Return the model and the fault specs."""
    specs = fault_specs()
    model = GeologicalModel(ORIGIN, MAXIMUM)
    model.data = pd.concat((_strat_data(), _fault_data(specs)), ignore_index=True)
    created = {}
    for k, s in enumerate(specs):
        kwargs = dict(
            fault_center=s["centre"],
            fault_normal_vector=s["normal"],
            fault_slip_vector=s["slip"],
            major_axis=s["major_axis"],
            minor_axis=s["minor_axis"],
            intermediate_axis=s["intermediate_axis"],
            nelements=nelements,
        )
        if s["cut_by"]:
            kwargs["faults"] = [created[m] for m in s["cut_by"]]
        created[k] = model.create_and_add_fault(s["name"], s["displacement"], **kwargs)
    young = model.create_and_add_foliation("group_1", nelements=nelements * 3)
    model.add_unconformity(young, -900.0)
    model.create_and_add_foliation("group_0", nelements=nelements * 3)
    _set_column(model)
    model.update()
    return model, specs


def points(count):
    return np.random.default_rng(1).uniform(ORIGIN, MAXIMUM, (count, 3))


class _FaultTimer:
    """Count every fault call and time the outermost ``apply_to_points`` call of each fault."""

    def __init__(self):
        self.calls, self.seconds, self._depth = {}, {}, 0
        self._original = FaultSegment.apply_to_points

    def __enter__(self):
        timer, original = self, self._original

        def timed(fault, *args, **kwargs):
            timer._depth += 1
            start = time.perf_counter()
            try:
                return original(fault, *args, **kwargs)
            finally:
                timer._depth -= 1
                timer.calls[fault.name] = timer.calls.get(fault.name, 0) + 1
                if timer._depth == 0:
                    elapsed = time.perf_counter() - start
                    timer.seconds[fault.name] = timer.seconds.get(fault.name, 0.0) + elapsed

        FaultSegment.apply_to_points = timed
        return self

    def __exit__(self, *exc):
        FaultSegment.apply_to_points = self._original


def measure(model, specs, counts):
    kind = {s["name"]: s["kind"] for s in specs}
    result, arrays = {}, {}
    for n in counts:
        xyz = points(n)
        with _FaultTimer() as timer:
            start = time.perf_counter()
            ids = model.evaluate_model(xyz)
            seconds = time.perf_counter() - start
        per = {"isolated": [], "cutting": []}
        for name, value in timer.seconds.items():
            per[kind[name]].append(1000 * value)
        result[str(n)] = dict(
            evaluate_model_seconds=seconds,
            fault_calls=sum(timer.calls.values()),
            isolated_mean_ms=float(np.mean(per["isolated"])),
            cutting_mean_ms=float(np.mean(per["cutting"])),
            cutting_max_ms=float(np.max(per["cutting"])),
        )
        arrays[f"model_ids_{n}"] = np.asarray(ids)
        for group in GROUPS:
            arrays[f"{group}_{n}"] = model.evaluate_feature_value(group, xyz)
        print(f"{n:6d} points: {seconds:6.2f} s, {result[str(n)]}", flush=True)
    return result, arrays


def compare(arrays, reference_path):
    """Return the names of the arrays that differ from the reference (NaN equals NaN)."""
    with np.load(reference_path) as reference:
        names = sorted(set(reference.files) | set(arrays))
        return [
            name
            for name in names
            if name not in reference.files
            or name not in arrays
            or not np.array_equal(reference[name], arrays[name], equal_nan=True)
        ]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--points", type=int, nargs="+", default=[500, 2000, 5000])
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--save-reference", type=Path)
    group.add_argument("--compare", type=Path)
    options = parser.parse_args(argv)
    logging.disable(logging.WARNING)
    start = time.perf_counter()
    model, specs = build()
    result = dict(
        faults=len(specs),
        isolated=sum(s["kind"] == "isolated" for s in specs),
        cutting=sum(s["kind"] == "cutting" for s in specs),
        build_seconds=time.perf_counter() - start,
    )
    print(f"built {result['faults']} faults in {result['build_seconds']:.0f} s", flush=True)
    result["cases"], arrays = measure(model, specs, options.points)
    if options.save_reference:
        np.savez_compressed(options.save_reference, **arrays)
    if options.compare:
        result["differences"] = compare(arrays, options.compare)
        print("bit-identical" if not result["differences"] else f"differ: {result['differences']}")
    options.output.write_text(json.dumps(result, indent=1), encoding="utf8")
    return 1 if result.get("differences") else 0


if __name__ == "__main__":
    raise SystemExit(main())
