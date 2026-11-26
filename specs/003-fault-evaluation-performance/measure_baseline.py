"""Baseline for spec 003: count fault applications and time evaluate_value
for models with F faults and U unconformities (see research.md)."""

import logging
import time

import numpy as np
import pandas as pd
from LoopStructural import GeologicalModel
from LoopStructural.modelling.features.fault import FaultSegment

logging.disable(logging.CRITICAL)

COUNTS = {'apply_to_points': 0, 'points': 0}
_original = FaultSegment.apply_to_points


def _counted(self, points, *args, **kwargs):
    COUNTS['apply_to_points'] += 1
    COUNTS['points'] += len(points)
    return _original(self, points, *args, **kwargs)


FaultSegment.apply_to_points = _counted

ORIGIN, MAXIMUM = np.array([0.0, 0.0, -1000.0]), np.array([1000.0, 1000.0, 0.0])


def build(n_faults, n_unconformities, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_faults):
        x = 150 + 700 * (i + 0.5) / max(n_faults, 1)
        name = f'fault_{i}'
        for y in (100.0, 500.0, 900.0):
            rows.append({'X': x, 'Y': y, 'Z': -500.0, 'val': 0.0, 'feature_name': name})
        rows.append({'X': x, 'Y': 500.0, 'Z': -500.0, 'gx': 1.0, 'gy': 0.2 * i, 'gz': 0.0,
                     'feature_name': name})
    groups = [f'group_{k}' for k in range(n_unconformities + 1)]
    for k, name in enumerate(groups):
        tilt = 0.05 * k
        for value in (0.0, 50.0, 100.0):
            xy = rng.random((30, 2)) * 1000
            z = -900 + 150 * k + value + tilt * xy[:, 0]
            for (x, y), zz in zip(xy, z):
                rows.append({'X': x, 'Y': y, 'Z': zz, 'val': value, 'feature_name': name})
        rows.append({'X': 500.0, 'Y': 500.0, 'Z': -900 + 150 * k, 'gx': -tilt, 'gy': 0.0, 'gz': 1.0,
                     'feature_name': name})
    model = GeologicalModel(ORIGIN, MAXIMUM)
    model.data = pd.DataFrame(rows)
    for i in range(n_faults):
        model.create_and_add_fault(f'fault_{i}', 40, nelements=3000)
    # youngest group last, with an unconformity at the top of each older group
    for k, name in enumerate(groups):
        feature = model.create_and_add_foliation(name, nelements=3000)
        if k < len(groups) - 1:
            model.add_unconformity(feature, 100.0)
    model.update()
    return model, groups


def measure(n_faults, n_unconformities, n_points=5000):
    model, groups = build(n_faults, n_unconformities)
    points = np.random.default_rng(1).random((n_points, 3)) * (MAXIMUM - ORIGIN) + ORIGIN
    youngest = model.get_feature_by_name(groups[0])
    oldest = model.get_feature_by_name(groups[-1])
    result = {}
    for label, call in (
        ('oldest group', lambda: oldest.evaluate_value(points)),
        ('all groups', lambda: [model.get_feature_by_name(g).evaluate_value(points) for g in groups]),
    ):
        COUNTS['apply_to_points'] = COUNTS['points'] = 0
        start = time.perf_counter()
        call()
        result[label] = (time.perf_counter() - start, COUNTS['apply_to_points'], COUNTS['points'])
    del youngest
    return result


print('F  U | oldest group: s, fault applications | all groups: s, fault applications')
for n_faults in (0, 1, 2, 4):
    for n_unconformities in (0, 1, 2, 3):
        r = measure(n_faults, n_unconformities)
        a, b = r['oldest group'], r['all groups']
        print(f'{n_faults}  {n_unconformities} | {a[0]:6.2f} s {a[1]:4d} | {b[0]:6.2f} s {b[1]:4d}', flush=True)
