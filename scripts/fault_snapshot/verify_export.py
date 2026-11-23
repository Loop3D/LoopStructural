"""Check an export with a small numpy reader: grid nodes + geometry must reproduce the reference
frame values, and the manifest must be complete. Usage: python verify_export.py OUT_DIR/model"""
import json, sys
from pathlib import Path
import numpy as np

def load(path, shape=None):
    a = np.fromfile(path, dtype="<f8")
    return a.reshape(shape) if shape else a

def trilinear(grid, root, pts):
    nodes = load(root / grid["nodes"])
    a = np.array(grid["local_from_world"]["A"]).reshape(3, 3); b = np.array(grid["local_from_world"]["b"])
    o, h, n = np.array(grid["origin"]), np.array(grid["step"]), np.array(grid["ncells"])
    l = pts @ a.T + b
    out = np.full(len(pts), np.nan)
    ok = np.all((l >= o - 1e-12) & (l <= o + n * h + 1e-12), axis=1)
    c = np.clip(np.floor((l[ok] - o) / h).astype(int), 0, n - 1)
    t = np.clip((l[ok] - (o + h * c)) / h, 0, 1)
    nx, ny = n[0] + 1, n[1] + 1
    v = 0
    for k in range(8):
        bx, by, bz = k & 1, (k >> 1) & 1, k >> 2
        w = np.where(bx, t[:, 0], 1 - t[:, 0]) * np.where(by, t[:, 1], 1 - t[:, 1]) * np.where(bz, t[:, 2], 1 - t[:, 2])
        v = v + w * nodes[(c[:, 0] + bx) + nx * ((c[:, 1] + by) + ny * (c[:, 2] + bz))]
    out[ok] = v
    return out

root = Path(sys.argv[1])
m = json.loads((root / "manifest.json").read_text())
pts = load(root / "reference/points.bin").reshape(-1, 3)
worst = 0.0
for f in m["faults"]:
    for c in range(3):
        ref = load(root / f"reference/fault_{f['name']}_frame_g{'xyz'[c]}.bin")
        got = trilinear(f["frame"][c], root, pts)
        assert np.array_equal(np.isnan(ref), np.isnan(got)), (f["name"], c, "nan pattern")
        worst = max(worst, np.nanmax(np.abs(ref - got)))
for s in m["series"]:
    ref = load(root / f"reference/series_{s['name']}_value.bin")
    # no faults: direct. With faults the restored points are needed: use the reference restore chain
    k = len(s["faults"])
    rp = pts if k == 0 else load(root / f"reference/series_{s['name']}_restored_after_{k-1}.bin").reshape(-1, 3)
    got = trilinear(s["grid"], root, rp)
    nanm = np.isnan(ref)
    # discrete interpolator: NaN input -> 0 value; outside grid -> NaN
    worst = max(worst, np.nanmax(np.abs(ref[~nanm] - np.nan_to_num(got)[~nanm])))
print(f"{root.name}: reproduces reference, max |diff| = {worst:.2e}; gradient scale = {m['gradient']['scale']:.6f}")
assert worst < 1e-10
