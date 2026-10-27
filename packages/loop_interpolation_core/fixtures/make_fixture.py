"""Make and check the cross-build fixture ``fixture_20.bin``.

Use:
    python -I make_fixture.py            write fixture_20.bin
    python -I make_fixture.py --check    read fixture_20.bin, solve again with the
                                         Python matrix solver, compare

The problem is a 20 x 20 x 20 ``FiniteDifferenceInterpolator`` problem with
value, gradient, normal, tangent and interface data (fixed seed 20240607).
The expected field is the reference solution: the matrix solver with
``solver="cg"`` and a tight tolerance (rtol = atol = 1e-10).

File format (little endian). Bytes 0 to 3 are the ASCII text ``LIF1``. Then a
list of blocks. Each block is a ``u32`` length N and N values of ``f64``.
Block order:

    origin(3) counts(3) step(3) reg(6: dxx dyy dzz dxy dxz dyz)
    border(3: dx dy dz) ridge(1)
    value_points value_values value_weights
    gradient_points gradient_vectors gradient_weights
    normal_points normal_vectors normal_weights
    tangent_points tangent_vectors tangent_weights
    interface_points interface_ids interface_weight(1)
    expected_field(N = product of counts)

Points and vectors are flat arrays (3 values for each point). The node order
of ``expected_field`` is the node order of ``interpolator.c`` of the matrix
solver (x is the fastest index).

MAPPING from the Python data columns to the core meaning
--------------------------------------------------------
All weights in the file are in core meaning: the core applies the weight
verbatim to a row that it has already scaled to unit L2 norm. The Python
matrix solver (``loop_interpolation.FiniteDifferenceInterpolator``) gives
these equivalents. Python column layouts are x y z | data | weight.

* Grid. ``counts`` is the number of NODES for each axis. The Python grid is
  ``StructuredGrid(origin, nsteps_cells=counts - 1, step_vector=step)``. The
  node count is ``counts``. No rotation.
* Regularisation. ``dxx dyy dzz`` are the same in Python and core. The
  cross-derivative weights are divided by 4 inside
  ``StructuredGrid.get_operators``. The file holds the core meaning, so
  ``dxy_file = dxy_python / 4`` (and the same for ``dxz``, ``dyz``). Python
  setup uses ``dxy = 4 * dxy_file``.
* Border. ``dx dy dz`` are the same in Python and core. The Python border rows
  are the second-derivative rows on the node layer 1 and layer n - 2 of each
  axis (``assemble_borders``).
* Ridge. ``ridge`` is ``ridge_factor`` (Python default 1e-8). It adds
  ``ridge**2 * I`` to the normal equations.
* Value. Python row weight = ``cpw * column_weight`` (column 4). File weight =
  ``cpw * column_weight``. The generator uses ``cpw = 1``, so the file weight
  is the column weight.
* Gradient (direction only, unit vector, two rows for each point: strike and
  dip). Python row weight = ``gpw``. The weight column is NOT used by Python.
  File weight = ``gpw``. All gradient weights in the file are equal.
* Normal (vector with length, three rows for each point, one for each axis).
  Python row weight = ``npw / 3 * column_weight`` (column 6). File weight =
  ``npw / 3 * column_weight``. The generator uses ``npw = 3``, so the file
  weight is the column weight. The vector is NOT made a unit vector.
* Tangent (``grad . t = 0``, vector made unit length, one row for each point).
  Python row weight = ``tpw``. The weight column is NOT used by Python. File
  weight = ``tpw``. All tangent weights in the file are equal.
* Interface (equal value for all points with one id). Python row weight =
  ``ipw``. Python must use ``interface_pair_mode = "star"``: each point is
  paired with the first point of its id group. This is the pairing of the core.
  File weight = ``ipw`` (one value).
* Rows outside the grid are dropped. The generator puts all data inside the
  grid.
* Solver. ``solve_system("cg", solver_kwargs={"rtol": 1e-10, "atol": 1e-10,
  "maxiter": 50000})``. The field is ``interpolator.c``.
"""

from __future__ import annotations

import struct
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FIXTURE = HERE / "fixture_20.bin"
MAGIC = b"LIF1"
SEED = 20240607

# Put the packages of this repository first on the import path. This also works
# with ``python -I``, which ignores PYTHONPATH.
for _src in (
    HERE.parents[1] / "loop_common" / "src",
    HERE.parents[1] / "loop_interpolation" / "src",
):
    if _src.is_dir():
        sys.path.insert(0, str(_src))

CG_KWARGS = {"rtol": 1e-10, "atol": 1e-10, "maxiter": 50000}

ORIGIN = np.array([10.0, -5.0, 3.0])
COUNTS = np.array([20, 20, 20])
STEP = np.array([1.0, 1.0, 1.0])
# dxx dyy dzz dxy dxz dyz in core meaning (cross terms already divided by 4)
REG = np.array([1.0, 1.0, 1.0, 0.25, 0.25, 0.25])
BORDER = np.array([1.0, 1.0, 1.0])
RIDGE = 1e-8
GRADIENT_WEIGHT = 1.0
TANGENT_WEIGHT = 1.0
INTERFACE_WEIGHT = 1.0

BLOCK_NAMES = (
    "origin counts step reg border ridge "
    "value_points value_values value_weights "
    "gradient_points gradient_vectors gradient_weights "
    "normal_points normal_vectors normal_weights "
    "tangent_points tangent_vectors tangent_weights "
    "interface_points interface_ids interface_weight "
    "expected_field"
).split()


def _field(p):
    """Smooth scalar field that makes the test data."""
    x = (p[:, 0] - ORIGIN[0]) / 19.0
    y = (p[:, 1] - ORIGIN[1]) / 19.0
    z = (p[:, 2] - ORIGIN[2]) / 19.0
    return 3.0 * x + 1.5 * np.sin(2.0 * y) + 2.0 * z * z + 0.5 * x * y


def _field_gradient(p):
    x = (p[:, 0] - ORIGIN[0]) / 19.0
    y = (p[:, 1] - ORIGIN[1]) / 19.0
    z = (p[:, 2] - ORIGIN[2]) / 19.0
    gx = 3.0 + 0.5 * y
    gy = 3.0 * np.cos(2.0 * y) + 0.5 * x
    gz = 4.0 * z
    return np.column_stack([gx, gy, gz]) / 19.0


def make_inputs() -> dict:
    """Make the input arrays (no expected field). All data is inside the grid."""
    rng = np.random.default_rng(SEED)
    lo = ORIGIN + 1.5 * STEP
    hi = ORIGIN + (COUNTS - 1 - 1.5) * STEP

    def pts(n):
        return rng.uniform(lo, hi, size=(n, 3))

    vp = pts(60)
    gp = pts(20)
    npts = pts(15)
    tp = pts(15)
    ip = pts(12)
    g_dir = _field_gradient(gp)
    g_dir = g_dir / np.linalg.norm(g_dir, axis=1)[:, None]
    n_vec = _field_gradient(npts)
    t_grad = _field_gradient(tp)
    # tangent vector: a unit vector that is orthogonal to the field gradient
    helper = np.tile([0.0, 0.0, 1.0], (len(tp), 1))
    t_vec = np.cross(t_grad, helper)
    t_vec /= np.linalg.norm(t_vec, axis=1)[:, None]
    ids = np.repeat(np.arange(4, dtype=float), 3)
    return {
        "origin": ORIGIN,
        "counts": COUNTS.astype(float),
        "step": STEP,
        "reg": REG,
        "border": BORDER,
        "ridge": np.array([RIDGE]),
        "value_points": vp.ravel(),
        "value_values": _field(vp),
        "value_weights": rng.uniform(0.5, 2.0, len(vp)),
        "gradient_points": gp.ravel(),
        "gradient_vectors": g_dir.ravel(),
        "gradient_weights": np.full(len(gp), GRADIENT_WEIGHT),
        "normal_points": npts.ravel(),
        "normal_vectors": n_vec.ravel(),
        "normal_weights": rng.uniform(0.5, 2.0, len(npts)),
        "tangent_points": tp.ravel(),
        "tangent_vectors": t_vec.ravel(),
        "tangent_weights": np.full(len(tp), TANGENT_WEIGHT),
        "interface_points": ip.ravel(),
        "interface_ids": ids,
        "interface_weight": np.array([INTERFACE_WEIGHT]),
    }


def solve_reference(d: dict) -> np.ndarray:
    """Solve the problem of ``d`` with the Python matrix solver (cg)."""
    from loop_common.supports import StructuredGrid  # noqa: PLC0415
    from loop_interpolation import FiniteDifferenceInterpolator  # noqa: PLC0415

    counts = d["counts"].astype(int)
    grid = StructuredGrid(origin=d["origin"], nsteps_cells=counts - 1, step_vector=d["step"])
    assert tuple(np.asarray(grid.nsteps).astype(int)) == tuple(counts), grid.nsteps
    interp = FiniteDifferenceInterpolator(grid)

    def rows(name, k):
        return np.asarray(d[name], dtype=float).reshape(-1, k)

    interp.set_value_constraints(
        np.column_stack([rows("value_points", 3), d["value_values"], d["value_weights"]])
    )
    gw = d["gradient_weights"]
    assert np.all(gw == gw[0]), "Python ignores the gradient weight column"
    interp.set_gradient_constraints(
        np.column_stack([rows("gradient_points", 3), rows("gradient_vectors", 3), gw])
    )
    interp.set_normal_constraints(
        np.column_stack([rows("normal_points", 3), rows("normal_vectors", 3), d["normal_weights"]])
    )
    tw = d["tangent_weights"]
    assert np.all(tw == tw[0]), "Python ignores the tangent weight column"
    interp.set_tangent_constraints(
        np.column_stack([rows("tangent_points", 3), rows("tangent_vectors", 3), tw])
    )
    interp.set_interface_constraints(
        np.column_stack([rows("interface_points", 3), d["interface_ids"]])
    )
    interp.ridge_factor = float(d["ridge"][0])
    r = d["reg"]
    b = d["border"]
    interp.setup_interpolator(
        dxx=r[0],
        dyy=r[1],
        dzz=r[2],
        dxy=4.0 * r[3],
        dxz=4.0 * r[4],
        dyz=4.0 * r[5],
        dx=b[0],
        dy=b[1],
        dz=b[2],
        cpw=1.0,
        gpw=float(gw[0]),
        npw=3.0,
        tpw=float(tw[0]),
        ipw=float(d["interface_weight"][0]),
        interface_pair_mode="star",
    )
    ok = interp.solve_system("cg", solver_kwargs=dict(CG_KWARGS))
    if not ok:
        raise RuntimeError("the Python cg solve did not converge")
    return np.asarray(interp.c, dtype=float).copy()


def write(path: Path, d: dict, expected: np.ndarray) -> None:
    with open(path, "wb") as f:
        f.write(MAGIC)
        for name in BLOCK_NAMES:
            a = expected if name == "expected_field" else d[name]
            a = np.ascontiguousarray(a, dtype="<f8").ravel()
            f.write(struct.pack("<I", a.size))
            f.write(a.tobytes())


def read(path: Path) -> dict:
    raw = path.read_bytes()
    if raw[:4] != MAGIC:
        raise ValueError("bad magic")
    pos = 4
    out = {}
    for name in BLOCK_NAMES:
        (n,) = struct.unpack_from("<I", raw, pos)
        pos += 4
        out[name] = np.frombuffer(raw, dtype="<f8", count=n, offset=pos).copy()
        pos += 8 * n
    if pos != len(raw):
        raise ValueError("extra bytes at the end of the file")
    return out


def main(argv: list[str]) -> int:
    if "--check" in argv:
        d = read(FIXTURE)
        n = int(np.prod(d["counts"]))
        if d["expected_field"].size != n:
            print("expected_field has the wrong size")
            return 1
        field = solve_reference(d)
        span = float(np.ptp(d["expected_field"]))
        diff = float(np.max(np.abs(field - d["expected_field"])))
        print(f"max difference {diff:.3e}, field range {span:.3e}, ratio {diff / span:.3e}")
        return 0 if diff <= 1e-6 * span else 1
    d = make_inputs()
    expected = solve_reference(d)
    write(FIXTURE, d, expected)
    print(f"wrote {FIXTURE} ({FIXTURE.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
