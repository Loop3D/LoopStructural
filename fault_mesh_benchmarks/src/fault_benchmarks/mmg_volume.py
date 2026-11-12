"""Labelled tetrahedral meshes of the benchmark cases with MMG3D level-set cuts.

The volume is cut by a sequence of MMG3D ``-ls`` calls. Each call inserts the
zero set of ONE smooth field, and only in the materials (tetrahedron
references) that are listed as active in an ``LSReferences`` table. All other
materials are ``nosplit``: MMG subdivides them only where it must keep the
mesh conforming, and it keeps their reference.

A Python :class:`MaterialTable` records, for each MMG reference, the
geological meaning of that material:

* the side of each fault that is already inserted (hanging wall +1,
  footwall -1, or ``None`` where the fault does not exist);
* the series region (the material is above a younger unconformity, or not
  assigned yet);
* the range of stratigraphic intervals that the material can still contain.

Every cut splits each active reference into two fresh references (negative
and positive side of the field), so after the last cut every reference has
exactly one fault block and one unit.

Cut order (youngest structure first, so that each structure is inserted
only in the region where it exists):

1. for each series (youngest first): its faults that are not yet inserted,
   in ``fault_build_order`` (younger and host faults first), then its
   unconformity surface, then all remaining material is given to the series
   if it has no unconformity;
2. the horizons of each series, one MMG call per isovalue and fault block.

A fault field is evaluated with the faults that displace it forced to one
side (``fields.forced_sides``), one call per side pattern. An abutting fault
is inserted only in the materials on the host side where it exists. A fault
that an unconformity truncates is inserted only below that unconformity.

MMG cannot end an interface inside the volume. A fault with a tip is
therefore inserted as (a) its full extension ``gx = 0`` and (b) the tip
surface ``t = 0`` on the footwall side. The geological footwall block is
``gx < 0 and t > 0``. The extension and the tip cap stay in the mesh as
interfaces between materials with the same geology; they are not exported
as fault surface.

Surfaces are the tetrahedron faces between two materials: a fault face
separates two materials with different sides of that fault; a horizon face
separates two materials of the same fault block whose intervals are on
either side of the isovalue. The OBJ files use the volume mesh vertices, so
they conform with the volume and with each other.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .builder import BuiltCase
from .fields import force_map, forced_sides, raw_frame, tip_function
from .spec import profile_tip_axes

CUBE_FACE_REFS = {(0, 0.0): 1, (0, 1.0): 2, (1, 0.0): 3, (1, 1.0): 4, (2, 0.0): 5, (2, 1.0): 6}
ISO_REF = 10  # MMG CLI reference of a new level-set interface

# --------------------------------------------------------------------------
# Mesh container and Medit IO
# --------------------------------------------------------------------------


@dataclass
class TetMesh:
    points: np.ndarray
    tets: np.ndarray
    refs: np.ndarray
    triangles: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), dtype=np.int64))
    tri_refs: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))

    def volumes(self) -> np.ndarray:
        p = self.points[self.tets]
        return np.einsum(
            "ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0]
        ) / 6.0

    def centroids(self) -> np.ndarray:
        return self.points[self.tets].mean(axis=1)

    def orient(self) -> TetMesh:
        neg = self.volumes() < 0
        self.tets[neg, :2] = self.tets[neg, 1::-1]
        return self


def kuhn_cube(cells: int) -> TetMesh:
    """Unit cube split into ``cells**3`` cubes of 6 tetrahedra (Kuhn/Freudenthal)."""
    from itertools import permutations

    n = cells + 1
    ax = np.linspace(0.0, 1.0, n)
    pts = np.stack(np.meshgrid(ax, ax, ax, indexing="ij"), axis=-1).reshape(-1, 3)
    base = np.stack(np.meshgrid(*(np.arange(cells),) * 3, indexing="ij"), axis=-1).reshape(-1, 3)
    tets = []
    for order in permutations(range(3)):
        q = base.copy()
        chain = [np.ravel_multi_index(q.T, (n, n, n))]
        for axis in order:
            q = q.copy()
            q[:, axis] += 1
            chain.append(np.ravel_multi_index(q.T, (n, n, n)))
        tets.append(np.stack(chain, axis=1))
    tets = np.concatenate(tets)
    mesh = TetMesh(pts, tets, np.ones(len(tets), dtype=np.int64)).orient()
    faces = boundary_faces(mesh)
    mesh.triangles = faces
    mesh.tri_refs = cube_face_refs(pts[faces])
    return mesh


def cube_face_refs(tri_pts: np.ndarray) -> np.ndarray:
    out = np.zeros(len(tri_pts), dtype=np.int64)
    for (axis, value), ref in CUBE_FACE_REFS.items():
        on = np.all(np.abs(tri_pts[:, :, axis] - value) < 1e-12, axis=1)
        out[on] = ref
    return out


def write_medit(path: Path, mesh: TetMesh) -> None:
    with open(path, "w", newline="\n") as f:
        f.write("MeshVersionFormatted 2\nDimension 3\n")
        f.write(f"Vertices\n{len(mesh.points)}\n")
        np.savetxt(f, np.c_[mesh.points, np.zeros(len(mesh.points))], fmt="%.17g %.17g %.17g %d")
        f.write(f"Tetrahedra\n{len(mesh.tets)}\n")
        np.savetxt(f, np.c_[mesh.tets + 1, mesh.refs], fmt="%d")
        if len(mesh.triangles):
            f.write(f"Triangles\n{len(mesh.triangles)}\n")
            np.savetxt(f, np.c_[mesh.triangles + 1, mesh.tri_refs], fmt="%d")
        f.write("End\n")


_WIDTHS = {
    "Vertices": 4, "Tetrahedra": 5, "Triangles": 4, "Edges": 3, "Corners": 1,
    "RequiredVertices": 1, "RequiredTriangles": 1, "RequiredTetrahedra": 1,
    "RequiredEdges": 1, "Ridges": 1, "Normals": 3, "NormalAtVertices": 2,
    "Tangents": 3, "TangentAtVertices": 2, "Quadrilaterals": 5, "Prisms": 7,
}


def read_medit(path: Path) -> TetMesh:
    words = Path(path).read_text().split()
    data = {}
    i = 0
    while i < len(words):
        key = words[i]
        i += 1
        if key == "End":
            break
        if key in ("MeshVersionFormatted", "Dimension"):
            i += 1
            continue
        if key not in _WIDTHS:
            raise ValueError(f"unsupported Medit section {key!r} in {path}")
        n = int(words[i])
        i += 1
        w = _WIDTHS[key]
        dtype = float if key in ("Vertices", "Normals", "Tangents") else np.int64
        data[key] = np.array(words[i : i + n * w], dtype=dtype).reshape(n, w)
        i += n * w
    if "Prisms" in data and len(data["Prisms"]):
        raise ValueError("MMG wrote prisms; only tetrahedra are supported")
    v, t = data["Vertices"], data["Tetrahedra"]
    tri = data.get("Triangles", np.zeros((0, 4), dtype=np.int64))
    return TetMesh(v[:, :3].copy(), t[:, :4] - 1, t[:, 4].copy(), tri[:, :3] - 1,
                   tri[:, 3].copy()).orient()


# --------------------------------------------------------------------------
# Face topology
# --------------------------------------------------------------------------

_TET_FACES = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])


def face_table(mesh: TetMesh):
    """Unique faces, and for each face its one or two incident tetrahedra.

    Returns (faces (F,3) sorted vertex ids, tet_a (F,), tet_b (F,) with -1 for
    a boundary face, counts (F,)).
    """
    nt = len(mesh.tets)
    f = mesh.tets[:, _TET_FACES].reshape(-1, 3)
    owner = np.repeat(np.arange(nt), 4)
    fs = np.sort(f, axis=1)
    order = np.lexsort((fs[:, 2], fs[:, 1], fs[:, 0]))
    fs, owner = fs[order], owner[order]
    new = np.ones(len(fs), dtype=bool)
    new[1:] = np.any(fs[1:] != fs[:-1], axis=1)
    start = np.flatnonzero(new)
    counts = np.diff(np.append(start, len(fs)))
    faces = fs[start]
    tet_a = owner[start]
    tet_b = np.full(len(faces), -1, dtype=np.int64)
    two = counts >= 2
    tet_b[two] = owner[start[two] + 1]
    return faces, tet_a, tet_b, counts


def boundary_faces(mesh: TetMesh) -> np.ndarray:
    faces, _, _, counts = face_table(mesh)
    return faces[counts == 1]


# --------------------------------------------------------------------------
# MMG runner
# --------------------------------------------------------------------------


def find_mmg(explicit: str | None = None) -> str:
    if explicit:
        return str(Path(explicit).resolve())
    if os.environ.get("MMG3D"):
        return os.environ["MMG3D"]
    here = Path(__file__).resolve()
    for parent in here.parents:
        cand = parent / "fault_mmg_experiment" / ".mmg-runtime" / "bin" / "mmg3d_O3.exe"
        if cand.exists():
            return str(cand)
    for name in ("mmg3d_O3", "mmg3d"):
        if shutil.which(name):
            return shutil.which(name)
    raise FileNotFoundError("MMG3D not found: install conda-forge::mmgsuite or set MMG3D")


@dataclass
class MmgOptions:
    executable: str | None = None
    adapt: bool = True
    """True: MMG remeshes after each cut (-hsiz/-hausd). False: pure cut."""
    hsiz: float = 0.05
    hausd: float = 0.002
    hgrad: float = 1.4
    hmin: float | None = None
    extra: tuple = ()
    timeout: float = 900.0
    keep_files: bool = False


def run_mmg_ls(mesh: TetMesh, values: np.ndarray, mapping: dict, workdir: Path,
               opts: MmgOptions) -> tuple[TetMesh, dict]:
    """One MMG3D level-set call. ``mapping``: ref -> (neg, pos) or None (nosplit)."""
    workdir.mkdir(parents=True, exist_ok=True)
    write_medit(workdir / "in.mesh", mesh)
    with open(workdir / "in.sol", "w", newline="\n") as f:
        f.write(f"MeshVersionFormatted 2\nDimension 3\nSolAtVertices\n{len(values)}\n1 1\n")
        np.savetxt(f, values, fmt="%.17g")
        f.write("End\n")
    lines = ["LSReferences", str(len(mapping))]
    for ref, pair in sorted(mapping.items()):
        lines.append(f"{ref} nosplit" if pair is None else f"{ref} {pair[0]} {pair[1]}")
    (workdir / "in.mmg3d").write_text("\n".join(lines) + "\n")
    cmd = [find_mmg(opts.executable), "-in", "in.mesh", "-sol", "in.sol", "-out", "out.mesh",
           "-ls", "0", "-f", "in.mmg3d", "-v", "1"]
    if opts.adapt:
        cmd += ["-hsiz", str(opts.hsiz), "-hausd", str(opts.hausd), "-hgrad", str(opts.hgrad)]
        if opts.hmin is not None:
            cmd += ["-hmin", str(opts.hmin)]
    else:
        cmd += ["-noinsert", "-nomove", "-noswap", "-hausd", "1e-8"]
    cmd += list(opts.extra)
    t0 = time.perf_counter()
    res = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, errors="replace",
                         timeout=opts.timeout, check=False)
    log = res.stdout + res.stderr
    (workdir / "mmg.log").write_text(log)
    rec = {
        "command": " ".join(Path(c).name if i == 0 else c for i, c in enumerate(cmd)),
        "returncode": res.returncode,
        "seconds": round(time.perf_counter() - t0, 3),
        "in_vertices": int(len(mesh.points)),
        "in_tets": int(len(mesh.tets)),
        "warnings": sorted({ln.strip() for ln in log.splitlines()
                            if "Warning" in ln or "Error" in ln or "WARNING" in ln}),
    }
    if res.returncode != 0 or not (workdir / "out.mesh").exists():
        raise RuntimeError(f"MMG failed ({res.returncode}): see {workdir / 'mmg.log'}")
    out = read_medit(workdir / "out.mesh")
    rec["out_vertices"] = int(len(out.points))
    rec["out_tets"] = int(len(out.tets))
    if not opts.keep_files:
        for name in ("in.mesh", "in.sol", "out.mesh", "out.sol"):
            (workdir / name).unlink(missing_ok=True)
    return out, rec


def run_mmg_remesh(mesh: TetMesh, workdir: Path, opts: MmgOptions,
                   optim: bool = False) -> tuple[TetMesh, dict]:
    """Plain MMG adaptation (no level set) that keeps material interfaces."""
    workdir.mkdir(parents=True, exist_ok=True)
    write_medit(workdir / "in.mesh", mesh)
    cmd = [find_mmg(opts.executable), "-in", "in.mesh", "-out", "out.mesh", "-v", "1",
           "-hausd", str(opts.hausd), "-hgrad", str(opts.hgrad)]
    cmd += ["-optim"] if optim else ["-hsiz", str(opts.hsiz)]
    t0 = time.perf_counter()
    res = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, errors="replace",
                         timeout=opts.timeout, check=False)
    (workdir / "mmg.log").write_text(res.stdout + res.stderr)
    if res.returncode != 0 or not (workdir / "out.mesh").exists():
        raise RuntimeError(f"MMG remesh failed ({res.returncode}): see {workdir / 'mmg.log'}")
    out = read_medit(workdir / "out.mesh")
    rec = {"command": " ".join(cmd[1:]), "returncode": res.returncode,
           "seconds": round(time.perf_counter() - t0, 3), "out_tets": int(len(out.tets))}
    if not opts.keep_files:
        for name in ("in.mesh", "out.mesh", "out.sol"):
            (workdir / name).unlink(missing_ok=True)
    return out, rec


# --------------------------------------------------------------------------
# Material table
# --------------------------------------------------------------------------


@dataclass
class Material:
    ref: int
    sides: dict = field(default_factory=dict)
    """fault -> +1 / -1 geological side (None: fault absent here)."""
    gx: dict = field(default_factory=dict)
    """fault -> sign of the raw gx cut (differs from ``sides`` beyond a tip)."""
    inside_tip: dict = field(default_factory=dict)
    region: str | None = None
    """Series that owns this material; None while it is below every unconformity so far."""
    interval: tuple[int, int] | None = None
    """(lo, hi) range of interval indices still possible in the region series."""
    parent: int | None = None
    cut: str | None = None
    """Name of the cut that made this material."""

    def child(self, ref: int, cut: str) -> Material:
        return Material(ref, dict(self.sides), dict(self.gx), dict(self.inside_tip),
                        self.region, self.interval, self.ref, cut)

    def to_json(self) -> dict:
        return {"ref": self.ref, "sides": self.sides, "gx": self.gx,
                "inside_tip": self.inside_tip, "region": self.region,
                "interval": list(self.interval) if self.interval else None,
                "parent": self.parent, "cut": self.cut}


class MaterialTable:
    def __init__(self):
        self.mats: dict[int, Material] = {1: Material(1)}
        self.history: dict[int, Material] = {1: Material(1)}
        self._next = 1000

    def new_ref(self) -> int:
        self._next += 1
        return self._next

    def split(self, active: list[int], cut: str, neg_update, pos_update) -> dict:
        """Mapping for MMG: each active ref -> two fresh refs; all others nosplit."""
        mapping: dict[int, tuple[int, int] | None] = {r: None for r in self.mats}
        for r in active:
            parent = self.mats[r]
            n, p = self.new_ref(), self.new_ref()
            mapping[r] = (n, p)
            for ref, upd in ((n, neg_update), (p, pos_update)):
                m = parent.child(ref, cut)
                upd(m)
                self.mats[ref] = m
                self.history[ref] = m
        return mapping

    def prune(self, present: np.ndarray) -> None:
        keep = set(np.unique(present).tolist())
        self.mats = {r: m for r, m in self.mats.items() if r in keep}


# --------------------------------------------------------------------------
# Cut planning and execution
# --------------------------------------------------------------------------


@dataclass
class VolumeOptions:
    cells: int = 20
    mmg: MmgOptions = field(default_factory=MmgOptions)
    final_optim: bool = True
    """Run one MMG -optim pass on the final mesh (no level set)."""


class VolumeMesher:
    def __init__(self, built: BuiltCase, workdir: Path, opts: VolumeOptions):
        self.built = built
        self.spec = built.spec
        self.workdir = Path(workdir)
        self.opts = opts
        self.mesh = kuhn_cube(opts.cells)
        self.table = MaterialTable()
        self.records: list[dict] = []
        self.inserted_faults: list[str] = []
        self.series_order = [s.name for s in self.spec.series]
        self.n_calls = 0

    # -- helpers ---------------------------------------------------------
    def _active_vertices(self, active_refs) -> np.ndarray:
        sel = np.isin(self.mesh.refs, list(active_refs))
        return np.unique(self.mesh.tets[sel])

    def _cut(self, name: str, active: list[int], fn, neg_update, pos_update) -> None:
        """Insert fn = 0 in the active materials (fn evaluated at their vertices)."""
        if not active:
            return
        vids = self._active_vertices(active)
        values = np.ones(len(self.mesh.points))
        vals = np.asarray(fn(self.mesh.points[vids]), dtype=float)
        bad = ~np.isfinite(vals)
        if bad.any():
            raise ValueError(f"{name}: {int(bad.sum())} non-finite field values")
        values[vids] = vals
        rec = {"cut": name, "active_refs": len(active), "active_vertices": int(len(vids)),
               "min": float(vals.min()), "max": float(vals.max())}
        if vals.min() >= 0 or vals.max() <= 0:
            # No sign change: no MMG call, relabel every active material.
            positive = not (vals.max() <= 0 and vals.min() < 0)
            for r in active:
                m = self.table.mats[r]
                (pos_update if positive else neg_update)(m)
                m.cut = name
            rec["mmg"] = "skipped (no sign change)"
            rec["side"] = "positive" if positive else "negative"
            self.records.append(rec)
            return
        mapping = self.table.split(active, name, neg_update, pos_update)
        self.n_calls += 1
        wd = self.workdir / f"{self.n_calls:03d}_{_slug(name)}"
        out, mrec = run_mmg_ls(self.mesh, values, mapping, wd, self.opts.mmg)
        unknown = set(np.unique(out.refs).tolist()) - set(self.table.mats)
        if unknown:
            raise RuntimeError(f"{name}: MMG returned unknown references {sorted(unknown)}")
        out.tri_refs[out.tri_refs == ISO_REF] = 100 + self.n_calls
        self.mesh = out
        self.table.prune(out.refs)
        rec.update(mrec)
        rec["materials_after"] = len(self.table.mats)
        self.records.append(rec)

    def _refs(self, pred) -> list[int]:
        present = set(np.unique(self.mesh.refs).tolist())
        return [r for r, m in self.table.mats.items() if r in present and pred(m)]

    def _patterns(self, names: list[str], pred) -> dict[tuple, list[int]]:
        """Group refs that satisfy ``pred`` by their side pattern over ``names``."""
        groups: dict[tuple, list[int]] = {}
        for r in self._refs(pred):
            m = self.table.mats[r]
            key = tuple(m.sides.get(n) for n in names)
            if any(k is None for k in key):
                raise RuntimeError(f"material {r}: side undefined for {names} ({m.sides})")
            groups.setdefault(key, []).append(r)
        return groups

    # -- faults ----------------------------------------------------------
    def insert_fault(self, name: str) -> None:
        fs = self.spec.fault(name)
        fault = self.built.faults[name]
        parents = list(fs.faulted_by)
        region_ok = (lambda m: m.region is None
                     or name in self.spec.series_by_name(m.region).faults)

        def exists(m):
            if not region_ok(m):
                return False
            if fs.abuts is not None:
                host, side = fs.abuts
                return m.sides.get(host) == side
            return True

        # Where the fault cannot exist in its region: canonical hanging-wall side.
        for r in self._refs(lambda m: region_ok(m) and not exists(m)):
            m = self.table.mats[r]
            m.sides[name] = 1
            m.gx[name] = None

        tip_y, tip_z = profile_tip_axes(fs.profile)
        has_tip = tip_y or tip_z
        for key, refs in self._patterns(parents, exists).items():
            pattern = dict(zip(parents, key))

            def gx_fn(p, pattern=pattern):
                with forced_sides(force_map(self.built, pattern)):
                    return raw_frame(fault, p)

            def neg(m):
                m.sides[name] = -1
                m.gx[name] = -1

            def pos(m):
                m.sides[name] = 1
                m.gx[name] = 1

            self._cut(f"fault {name} {_pat(pattern)}", refs, gx_fn, neg, pos)
        if has_tip:
            # Footwall side only: beyond the tip (t < 0) the footwall is geologically
            # the hanging-wall block (the fault displaces nothing there).
            for key, refs in self._patterns(
                    parents, lambda m: exists(m) and m.gx.get(name) == -1).items():
                pattern = dict(zip(parents, key))

                def t_fn(p, pattern=pattern):
                    with forced_sides(force_map(self.built, pattern)):
                        return tip_function(self.built, name, p)

                def outside(m):
                    m.sides[name] = 1
                    m.inside_tip[name] = False

                def inside(m):
                    m.sides[name] = -1
                    m.inside_tip[name] = True

                self._cut(f"tip {name} {_pat(pattern)}", refs, t_fn, outside, inside)
        self.inserted_faults.append(name)

    # -- series ----------------------------------------------------------
    def _isos(self, sname):
        ss = self.spec.series_by_name(sname)
        order = np.argsort(ss.isovalues)
        return ss, np.asarray(ss.isovalues, dtype=float)[order], order

    def _series_fn(self, sname, pattern, iso):
        feat = self.built.series[sname]

        def fn(p):
            with forced_sides(force_map(self.built, pattern)):
                return np.asarray(feat.evaluate_value(p, ignore_regions=True), dtype=float) - iso
        return fn

    def insert_unconformity(self, sname: str) -> None:
        ss, isos, _ = self._isos(sname)
        u = float(ss.unconformity)
        lo = int(np.searchsorted(isos, u, side="right"))
        n = len(isos)
        for key, refs in self._patterns(list(ss.faults), lambda m: m.region is None).items():
            pattern = dict(zip(ss.faults, key))

            def above(m, lo=lo):
                m.region = sname
                m.interval = (lo, n)

            self._cut(f"unconformity {sname} {_pat(pattern)}", refs,
                      self._series_fn(sname, pattern, u), lambda m: None, above)

    def assign_rest(self, sname: str) -> None:
        n = len(self.spec.series_by_name(sname).isovalues)
        for r in self._refs(lambda m: m.region is None):
            m = self.table.mats[r]
            m.region = sname
            m.interval = (0, n)

    def insert_horizons(self, sname: str) -> None:
        ss, isos, _ = self._isos(sname)
        for k, iso in enumerate(isos):
            if ss.unconformity is not None and iso <= ss.unconformity:
                continue  # below or at the unconformity: not a cut inside this series

            def can_contain(m, k=k):
                return m.region == sname and m.interval[0] <= k < m.interval[1]

            for key, refs in self._patterns(list(ss.faults), can_contain).items():
                pattern = dict(zip(ss.faults, key))

                def below(m, k=k):
                    m.interval = (m.interval[0], k)

                def above(m, k=k):
                    m.interval = (k + 1, m.interval[1])

                self._cut(f"horizon {sname} {iso:g} {_pat(pattern)}", refs,
                          self._series_fn(sname, pattern, float(iso)), below, above)

    # -- driver ----------------------------------------------------------
    def run(self) -> TetMesh:
        order = [f.name for f in self.spec.fault_build_order()]
        in_series = {f for s in self.spec.series for f in s.faults}
        for name in order:
            if name not in in_series:
                self.insert_fault(name)
        for ss in self.spec.series:
            for name in order:
                if name in ss.faults and name not in self.inserted_faults:
                    self.insert_fault(name)
            if ss.unconformity is not None:
                self.insert_unconformity(ss.name)
            else:
                self.assign_rest(ss.name)
        for ss in self.spec.series:
            self.insert_horizons(ss.name)
        if self.opts.final_optim and self.opts.mmg.adapt:
            self.n_calls += 1
            out, rec = run_mmg_remesh(self.mesh, self.workdir / f"{self.n_calls:03d}_optim",
                                      self.opts.mmg, optim=True)
            if set(np.unique(out.refs).tolist()) - set(self.table.mats):
                raise RuntimeError("optim pass returned unknown references")
            self.mesh = out
            rec["cut"] = "final optim"
            self.records.append(rec)
        return self.mesh


def _pat(pattern: dict) -> str:
    return "[" + ",".join(f"{k}{'+' if v > 0 else '-'}" for k, v in pattern.items()) + "]" \
        if pattern else "[all]"


def _slug(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-." else "_" for c in s)[:60]


# --------------------------------------------------------------------------
# Labels, surfaces and export
# --------------------------------------------------------------------------


def unit_offsets(spec) -> dict[str, int]:
    off, out = 0, {}
    for ss in spec.series:
        out[ss.name] = off
        off += len(ss.isovalues) + 1
    return out


def cell_labels(mesher: VolumeMesher) -> dict[str, np.ndarray]:
    spec, mats, refs = mesher.spec, mesher.table.mats, mesher.mesh.refs
    offs = unit_offsets(spec)
    uniq, inv = np.unique(refs, return_inverse=True)
    sname = [s.name for s in spec.series]

    def per_ref(fn, dtype):
        return np.asarray([fn(mats[int(r)]) for r in uniq], dtype=dtype)[inv]

    def unit(m):
        if m.region is None or m.interval is None or m.interval[0] != m.interval[1]:
            return -1
        return offs[m.region] + m.interval[0]

    out = {
        "material": refs.astype(np.int32),
        "unit": per_ref(unit, np.int32),
        "series": per_ref(lambda m: sname.index(m.region) if m.region else -1, np.int32),
    }
    for fs in spec.faults:
        out[f"side:{fs.name}"] = per_ref(lambda m, n=fs.name: m.sides.get(n) or 0, np.int8)
    for ss in spec.series:
        # Fault block of this series: index of the side pattern, -1 outside the series.
        pats = sorted({tuple(m.sides.get(f) for f in ss.faults)
                       for m in mats.values() if m.region == ss.name})

        def blk(m, ss=ss, pats=pats):
            if m.region != ss.name:
                return -1
            return pats.index(tuple(m.sides.get(f) for f in ss.faults))
        out[f"block:{ss.name}"] = per_ref(blk, np.int32)
    return out


@dataclass
class Surface:
    name: str
    kind: str
    faces: np.ndarray  # (F,3) vertex ids, oriented
    blocks: np.ndarray  # (F,) index into ``patterns``
    patterns: list  # side pattern (fault -> +1/-1) of each block
    meta: dict = field(default_factory=dict)


def extract_surfaces(mesher: VolumeMesher) -> list[Surface]:
    spec, mesh, mats = mesher.spec, mesher.mesh, mesher.table.mats
    faces, ta, tb, counts = face_table(mesh)
    if counts.max() > 2:
        raise RuntimeError("non-manifold volume: a face has more than two tetrahedra")
    inner = tb >= 0
    faces, ta, tb = faces[inner], ta[inner], tb[inner]
    ra, rb = mesh.refs[ta], mesh.refs[tb]
    diff = ra != rb
    faces, ta, tb, ra, rb = faces[diff], ta[diff], tb[diff], ra[diff], rb[diff]
    cent = mesh.centroids()
    pairs = np.unique(np.stack([ra, rb], axis=1), axis=0)
    pair_index = {(int(a), int(b)): i for i, (a, b) in enumerate(pairs)}
    pid = np.array([pair_index[(int(a), int(b))] for a, b in zip(ra, rb)])

    out: list[Surface] = []
    order = [f.name for f in spec.fault_build_order()]
    claimed = np.zeros(len(faces), dtype=bool)

    def orient(sel, flip_mask):
        """Faces oriented so the normal points from tet b to tet a (flip where asked)."""
        fc = faces[sel].copy()
        p = mesh.points
        n = np.cross(p[fc[:, 1]] - p[fc[:, 0]], p[fc[:, 2]] - p[fc[:, 0]])
        d = cent[ta[sel]] - cent[tb[sel]]
        wrong = np.einsum("ij,ij->i", n, d) < 0
        wrong ^= flip_mask
        fc[wrong] = fc[wrong][:, [0, 2, 1]]
        return fc

    for name in order:
        hit = np.zeros(len(pairs), dtype=bool)
        a_is_hw = np.zeros(len(pairs), dtype=bool)
        labels = []
        blk = np.zeros(len(pairs), dtype=np.int64)
        for i, (a, b) in enumerate(pairs):
            ma, mb = mats[int(a)], mats[int(b)]
            sa, sb = ma.sides.get(name), mb.sides.get(name)
            ga, gb = ma.gx.get(name, sa), mb.gx.get(name, sb)
            if sa is None or sb is None or sa == sb or ga is None or gb is None or ga == gb:
                continue
            hit[i] = True
            a_is_hw[i] = sa > 0
            lab = {f: ma.sides.get(f) for f in spec.fault(name).faulted_by
                   if ma.sides.get(f) is not None}
            if lab not in labels:
                labels.append(lab)
            blk[i] = labels.index(lab)
        sel = hit[pid] & ~claimed
        if not sel.any():
            continue
        claimed |= sel
        # normal toward the hanging wall: from b to a when a is the hanging wall
        fc = orient(sel, ~a_is_hw[pid[sel]])
        out.append(Surface(f"fault_{name}", "fault", fc, blk[pid[sel]], labels,
                           {"fault": name}))

    for ss in spec.series:
        isos = np.asarray(ss.isovalues, dtype=float)
        srt = np.argsort(isos)
        older = set(mesher.series_order[mesher.series_order.index(ss.name) + 1:]) | {None}
        for j, k_orig in enumerate(srt):
            iso = isos[k_orig]
            hit = np.zeros(len(pairs), dtype=bool)
            a_up = np.zeros(len(pairs), dtype=bool)
            blk = np.zeros(len(pairs), dtype=np.int64)
            labels = []
            for i, (a, b) in enumerate(pairs):
                ma, mb = mats[int(a)], mats[int(b)]
                for up, lo, flag in ((ma, mb, True), (mb, ma, False)):
                    if up.region != ss.name or up.interval is None or up.interval[0] != j + 1:
                        continue
                    if any(up.sides.get(f) != lo.sides.get(f) for f in ss.faults):
                        continue
                    below_same = (lo.region == ss.name and lo.interval is not None
                                  and lo.interval[1] == j)
                    below_older = (ss.unconformity is not None
                                   and abs(iso - ss.unconformity) < 1e-12
                                   and lo.region in older)
                    if below_same or below_older:
                        hit[i] = True
                        a_up[i] = flag
                        lab = {f: up.sides.get(f) for f in ss.faults}
                        if lab not in labels:
                            labels.append(lab)
                        blk[i] = labels.index(lab)
            sel = hit[pid] & ~claimed
            if not sel.any():
                continue
            fc = orient(sel, ~a_up[pid[sel]])
            out.append(Surface(f"horizon_{ss.name}_{int(k_orig)}", "horizon", fc,
                               blk[pid[sel]], labels,
                               {"series": ss.name, "isovalue": float(iso),
                                "unconformity": ss.unconformity is not None
                                and abs(iso - ss.unconformity) < 1e-12}))
    return out


def write_obj(path: Path, points: np.ndarray, faces: np.ndarray, header: str = "") -> None:
    used, inv = np.unique(faces, return_inverse=True)
    with open(path, "w", newline="\n") as f:
        if header:
            f.write("".join(f"# {ln}\n" for ln in header.splitlines()))
        np.savetxt(f, points[used], fmt="v %.10f %.10f %.10f")
        np.savetxt(f, inv.reshape(-1, 3) + 1, fmt="f %d %d %d")


def tet_quality(mesh: TetMesh) -> np.ndarray:
    """Normalised mean-ratio-like quality used by MMG (1 = regular tetrahedron)."""
    p = mesh.points[mesh.tets]
    d = p[:, [0, 0, 0, 1, 1, 2]] - p[:, [1, 2, 3, 2, 3, 3]]
    s = np.einsum("ijk,ijk->i", d, d)
    return 6.0 * mesh.volumes() / np.maximum(s, 1e-300) ** 1.5 * 12.0 * np.sqrt(3.0)


def min_dihedral_deg(mesh: TetMesh) -> np.ndarray:
    p = mesh.points[mesh.tets]
    fn = []
    for f in _TET_FACES:
        a, b, c = p[:, f[0]], p[:, f[1]], p[:, f[2]]
        n = np.cross(b - a, c - a)
        fn.append(n / np.maximum(np.linalg.norm(n, axis=1), 1e-300)[:, None])
    ang = []
    for i in range(4):
        for j in range(i + 1, 4):
            c = np.clip(-np.einsum("ij,ij->i", fn[i], fn[j]), -1, 1)
            ang.append(np.degrees(np.arccos(c)))
    return np.min(np.stack(ang, axis=1), axis=1)


def save_vtk(path: Path, mesh: TetMesh, cell_data: dict) -> None:
    import pyvista as pv

    cells = np.hstack([np.full((len(mesh.tets), 1), 4), mesh.tets]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(len(mesh.tets), pv.CellType.TETRA), mesh.points)
    for k, v in cell_data.items():
        grid.cell_data[k] = v
    grid.save(str(path))


def save_vtp(path: Path, points: np.ndarray, faces: np.ndarray, cell_data: dict) -> None:
    import pyvista as pv

    used, inv = np.unique(faces, return_inverse=True)
    tri = inv.reshape(-1, 3)
    poly = pv.PolyData(points[used], np.hstack([np.full((len(tri), 1), 3), tri]).ravel())
    for k, v in cell_data.items():
        poly.cell_data[k] = v
    poly.save(str(path))


# --------------------------------------------------------------------------
# QA and case output
# --------------------------------------------------------------------------


def _tri_area(points, faces) -> np.ndarray:
    p = points[faces]
    return 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)


def _edge_counts(faces):
    e = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    _, c = np.unique(e, axis=0, return_counts=True)
    return c


def _read_obj_area(path: Path) -> float | None:
    if not path.exists():
        return None
    v, f = [], []
    for ln in path.read_text().splitlines():
        if ln.startswith("v "):
            v.append([float(x) for x in ln.split()[1:4]])
        elif ln.startswith("f "):
            f.append([int(x.split("/")[0]) - 1 for x in ln.split()[1:4]])
    if not f:
        return 0.0
    return float(_tri_area(np.asarray(v), np.asarray(f)).sum())


def _distance_error(fn, pts, h=1e-5) -> np.ndarray:
    """|f| / |grad f| (first-order distance to the zero set), gradient by central differences."""
    v = np.asarray(fn(pts), dtype=float)
    g = np.zeros((len(pts), 3))
    for a in range(3):
        d = np.zeros(3)
        d[a] = h
        g[:, a] = (np.asarray(fn(pts + d), dtype=float)
                   - np.asarray(fn(pts - d), dtype=float)) / (2 * h)
    return np.abs(v) / np.maximum(np.linalg.norm(g, axis=1), 1e-12)


def surface_qa(mesher: VolumeMesher, s: Surface, reference_dir: Path | None) -> dict:
    pts = mesher.mesh.points
    area = _tri_area(pts, s.faces)
    ec = _edge_counts(s.faces)
    rec = {
        "kind": s.kind, **s.meta, "n_triangles": int(len(s.faces)),
        "n_vertices": int(len(np.unique(s.faces))), "area": float(area.sum()),
        "boundary_edges": int((ec == 1).sum()), "nonmanifold_edges": int((ec > 2).sum()),
        "degenerate_triangles": int((area < 1e-12).sum()),
        "blocks": [_pat(p) for p in s.patterns],
    }
    errs = []
    for b, pattern in enumerate(s.patterns):
        vids = np.unique(s.faces[s.blocks == b])
        if s.kind == "fault":
            fault = mesher.built.faults[s.meta["fault"]]

            def fn(p, pattern=pattern, fault=fault):
                with forced_sides(force_map(mesher.built, pattern)):
                    return raw_frame(fault, p)
        else:
            fn = mesher._series_fn(s.meta["series"], pattern, s.meta["isovalue"])
        errs.append(_distance_error(fn, pts[vids]))
    e = np.concatenate(errs) if errs else np.zeros(0)
    rec["vertex_distance_to_implicit_surface"] = {
        "max": float(e.max()) if len(e) else 0.0, "mean": float(e.mean()) if len(e) else 0.0,
        "p99": float(np.percentile(e, 99)) if len(e) else 0.0}
    if reference_dir is not None:
        ref = _read_obj_area(reference_dir / f"{s.name}.obj")
        rec["reference_area"] = ref
        if ref:
            rec["area_relative_difference"] = float(area.sum() / ref - 1.0)
    return rec


def volume_qa(mesher: VolumeMesher, surfaces: list[Surface], labels: dict) -> dict:
    mesh, mats = mesher.mesh, mesher.table.mats
    vol = mesh.volumes()
    q = tet_quality(mesh)
    dih = min_dihedral_deg(mesh)
    faces, ta, tb, counts = face_table(mesh)
    bnd = faces[counts == 1]
    on_cube = np.zeros(len(bnd), dtype=bool)
    bp = mesh.points[bnd]
    for axis in range(3):
        for val in (0.0, 1.0):
            on_cube |= np.all(np.abs(bp[:, :, axis] - val) < 1e-9, axis=1)
    inner = counts == 2
    ra, rb = mesh.refs[ta[inner]], mesh.refs[tb[inner]]
    fi = faces[inner][ra != rb]
    ra, rb = ra[ra != rb], rb[ra != rb]

    def geo(m):
        return (m.region, m.interval, tuple(sorted((k, v) for k, v in m.sides.items())))

    classified = set()
    for s in surfaces:
        classified.update(map(tuple, np.sort(s.faces, axis=1).tolist()))
    same_geo = np.array([geo(mats[int(a)]) == geo(mats[int(b)]) for a, b in zip(ra, rb)],
                        dtype=bool)
    is_cls = np.array([tuple(f) in classified for f in np.sort(fi, axis=1).tolist()],
                      dtype=bool)
    unexplained = ~same_geo & ~is_cls
    units = labels["unit"]
    return {
        "n_vertices": int(len(mesh.points)), "n_tets": int(len(mesh.tets)),
        "volume": float(vol.sum()), "min_tet_volume": float(vol.min()),
        "nonpositive_tets": int((vol <= 0).sum()),
        "faces_with_more_than_two_tets": int((counts > 2).sum()),
        "boundary_faces_not_on_cube": int((~on_cube).sum()),
        "unlabelled_tets": int((units < 0).sum()),
        "materials": int(len(np.unique(mesh.refs))),
        "units_present": sorted(int(u) for u in np.unique(units)),
        "unit_volumes": {int(u): float(vol[units == u].sum()) for u in np.unique(units)},
        "quality": {"min": float(q.min()), "median": float(np.median(q)),
                    "below_0.1": int((q < 0.1).sum()), "below_0.01": int((q < 0.01).sum())},
        "min_dihedral_deg": {"min": float(dih.min()), "below_5": int((dih < 5).sum()),
                             "below_1": int((dih < 1).sum())},
        "interface_faces": {
            "geological": int(is_cls.sum()),
            "non_geological_same_geology": int(same_geo.sum()),
            "non_geological_area": float(_tri_area(mesh.points, fi[same_geo]).sum()),
            "unexplained": int(unexplained.sum()),
            "unexplained_area": float(_tri_area(mesh.points, fi[unexplained]).sum()),
        },
    }


def write_case(mesher: VolumeMesher, out: Path, reference_dir: Path | None = None) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    (out / "surfaces").mkdir(exist_ok=True)
    mesh = mesher.mesh
    labels = cell_labels(mesher)
    labels["quality"] = tet_quality(mesh)
    labels["min_dihedral_deg"] = min_dihedral_deg(mesh)
    surfaces = extract_surfaces(mesher)
    save_vtk(out / "volume.vtu", mesh, labels)
    write_medit(out / "volume.mesh", mesh)
    surf_rep = {}
    for s in surfaces:
        header = (f"{s.name}: {s.kind} from LoopStructural benchmark {mesher.spec.name}\n"
                  f"MMG3D level-set volume mesh; vertices shared with volume.mesh")
        write_obj(out / "surfaces" / f"{s.name}.obj", mesh.points, s.faces, header)
        save_vtp(out / "surfaces" / f"{s.name}.vtp", mesh.points, s.faces, {"block": s.blocks})
        surf_rep[s.name] = surface_qa(mesher, s, reference_dir)
    offs = unit_offsets(mesher.spec)
    unit_table = []
    for ss in mesher.spec.series:
        isos = sorted(ss.isovalues)
        for i in range(len(isos) + 1):
            unit_table.append({"unit": offs[ss.name] + i, "series": ss.name,
                               "lower": isos[i - 1] if i > 0 else None,
                               "upper": isos[i] if i < len(isos) else None})
    mats = {}
    for r, m in sorted(mesher.table.mats.items()):
        d = m.to_json()
        d["unit"] = ((offs[m.region] + m.interval[0]) if m.region and m.interval
                     and m.interval[0] == m.interval[1] else -1)
        d["n_tets"] = int((mesh.refs == r).sum())
        mats[str(r)] = d
    (out / "materials.json").write_text(json.dumps(
        {"units": unit_table, "materials": mats,
         "lineage": {str(r): m.to_json() for r, m in sorted(mesher.table.history.items())}},
        indent=1, default=_json_default))
    report = {
        "case": mesher.spec.name,
        "options": {"cells": mesher.opts.cells,
                    **{k: v for k, v in vars(mesher.opts.mmg).items() if k != "extra"}},
        "mmg_calls": mesher.n_calls,
        "volume": volume_qa(mesher, surfaces, labels),
        "surfaces": surf_rep,
        "missing_surfaces": _missing(surfaces, reference_dir),
        "cuts": mesher.records,
    }
    (out / "report.json").write_text(json.dumps(report, indent=1, default=_json_default))
    return report


def _missing(surfaces, reference_dir) -> list[str]:
    if reference_dir is None or not reference_dir.exists():
        return []
    have = {s.name for s in surfaces}
    return sorted(p.stem for p in reference_dir.glob("*.obj") if p.stem not in have)


def _json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)
