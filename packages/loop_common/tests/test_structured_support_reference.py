"""Pin the behaviour of the structured interpolation supports.

``StructuredGrid``, ``TetMesh`` and ``P2TetMesh`` share ``BaseStructuredSupport``. This test
records a wide set of their outputs (geometry, node/element arrays, point location, evaluation,
regularisation operators, setters, ``onGeometryChange`` calls, constructor validation) on a few
grids and compares them with ``data/structured_support_reference.npz``.

The reference was generated from the code *before* ``BaseStructuredSupport`` was changed to extend
``StructuredGrid3DGeometry``, so it checks that refactor (and later ones) did not change behaviour.

Regenerate after an intended change, and review the diff:

    python packages/loop_common/tests/test_structured_support_reference.py --update
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
from loop_common.supports import P2TetMesh, StructuredGrid, TetMesh

REFERENCE = Path(__file__).parent / "data" / "structured_support_reference.npz"

# (origin, cells per axis, step vector). Step vectors are floats on purpose.
GRIDS = {
    "unit": ([0.0, 0.0, 0.0], [4, 4, 4], [1.0, 1.0, 1.0]),
    "aniso": ([10.0, -3.0, 2.0], [6, 4, 3], [0.5, 2.0, 1.5]),
    "small": ([0.0, 0.0, 0.0], [2, 2, 2], [1.0, 1.0, 1.0]),
}
CLASSES = {"StructuredGrid": StructuredGrid, "TetMesh": TetMesh, "P2TetMesh": P2TetMesh}
GEOMETRY_ATTRIBUTES = (
    "nsteps",
    "nsteps_cells",
    "maximum",
    "origin",
    "step_vector",
    "n_nodes",
    "n_elements",
    "volume",
    "nodes",
    "elements",
    "rotation_xy",
    "element_size",
    "element_scale",
)
REGULARISATION = {"dxx": 1, "dyy": 1, "dzz": 1, "dxy": 1, "dxz": 1, "dyz": 1}


def _build(cls, origin, cells, step):
    if cls is P2TetMesh:
        return cls(origin=origin, nsteps=cells, step_vector=step)
    return cls(origin=origin, nsteps_cells=cells, step_vector=step)


def _flatten(value, prefix, out):
    """Store ``value`` under ``prefix`` (tuples and dicts become ``prefix.<index or key>``)."""
    if isinstance(value, dict):
        for key in sorted(value):
            _flatten(value[key], f"{prefix}.{key}", out)
    elif isinstance(value, (tuple, list)) and not all(np.isscalar(v) for v in value):
        for i, item in enumerate(value):
            _flatten(item, f"{prefix}.{i}", out)
    else:
        out[prefix] = np.asarray(value)


def _record(out, key, call):
    try:
        _flatten(call(), key, out)
    except Exception as exc:  # noqa: BLE001 - the exception type is part of the behaviour
        out[key] = np.asarray(f"EXC:{type(exc).__name__}")


def collect() -> dict[str, np.ndarray]:
    """Run every probe and return ``{key: array}``."""
    out: dict[str, np.ndarray] = {}
    rng = np.random.default_rng(3)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for cls_name, cls in CLASSES.items():
            for grid_name, (origin, cells, step) in GRIDS.items():
                prefix = f"{cls_name}/{grid_name}"
                grid = _build(cls, origin, cells, step)
                calls = []
                original = type(grid).onGeometryChange
                type(grid).onGeometryChange = lambda self, _o=original, _c=calls: (
                    _c.append(1),
                    _o(self),
                )[1]
                try:
                    _probe(out, prefix, grid, origin, cells, step, rng, calls)
                finally:
                    type(grid).onGeometryChange = original
        for cells in ([0, 3, 3], [-1, 3, 3]):
            _record(
                out,
                f"validation/{cells}",
                lambda cells=cells: StructuredGrid(
                    origin=[0.0, 0.0, 0.0], nsteps_cells=cells, step_vector=[1.0, 1.0, 1.0]
                ),
            )
    return out


def _probe(out, prefix, grid, origin, cells, step, rng, calls):
    for name in GEOMETRY_ATTRIBUTES:
        _record(out, f"{prefix}/{name}", lambda name=name: getattr(grid, name))
    _record(out, f"{prefix}/to_dict", grid.to_dict)

    low = np.array(origin) - 1.0
    high = np.array(origin) + np.array(cells) * np.array(step) + 1.0
    anywhere = rng.uniform(low, high, (120, 3))
    # points exactly on the boundary: every corner and the centre of every face
    corners = np.array(np.meshgrid(*zip(origin, grid.maximum), indexing="ij")).reshape(3, -1).T
    centre = 0.5 * (np.array(origin) + grid.maximum)
    faces = []
    for axis in range(3):
        for bound in (origin[axis], grid.maximum[axis]):
            point = centre.copy()
            point[axis] = bound
            faces.append(point)
    anywhere = np.vstack([anywhere, corners, np.array(faces)])
    interior = rng.uniform(np.array(origin) + 1e-3, grid.maximum - 1e-3, (80, 3))
    _record(out, f"{prefix}/inside", lambda: grid.inside(anywhere))
    _record(out, f"{prefix}/position_to_cell_index", lambda: grid.position_to_cell_index(anywhere))
    _record(out, f"{prefix}/position_to_cell_corners", lambda: grid.position_to_cell_corners(interior))
    _record(out, f"{prefix}/position_to_cell_vertices", lambda: grid.position_to_cell_vertices(interior))
    _record(out, f"{prefix}/check_position_1d", lambda: grid.check_position(np.array([1.0, 2.0, 3.0])))
    _record(out, f"{prefix}/check_position_list", lambda: grid.check_position([[1.0, 2.0, 3.0]]))
    _record(out, f"{prefix}/rotate", lambda: grid.rotate(anywhere))

    cell_index = np.arange(grid.n_elements)
    node_index = np.arange(grid.n_nodes)
    _record(out, f"{prefix}/global_index_to_cell_index", lambda: grid.global_index_to_cell_index(cell_index))
    node_ijk = grid.global_index_to_node_index(node_index)
    _record(out, f"{prefix}/global_index_to_node_index", lambda: node_ijk)
    _record(out, f"{prefix}/node_indexes_to_position", lambda: grid.node_indexes_to_position(node_ijk[:20]))
    if cls_has_ijk(node_ijk):
        _record(out, f"{prefix}/global_node_indices", lambda: grid.global_node_indices(node_ijk))
        _record(
            out,
            f"{prefix}/global_cell_indices",
            lambda: grid.global_cell_indices(grid.global_index_to_cell_index(cell_index)),
        )

    values = rng.normal(size=grid.n_nodes)
    _record(out, f"{prefix}/evaluate_value", lambda: grid.evaluate_value(interior, values))
    _record(out, f"{prefix}/evaluate_gradient", lambda: grid.evaluate_gradient(interior, values))
    _record(out, f"{prefix}/get_element_for_location", lambda: grid.get_element_for_location(interior))
    _record(
        out,
        f"{prefix}/get_element_gradient_for_location",
        lambda: grid.get_element_gradient_for_location(interior),
    )
    _record(out, f"{prefix}/neighbour_global_indexes", lambda: grid.neighbour_global_indexes())
    if hasattr(grid, "get_operators"):
        _record(
            out,
            f"{prefix}/get_operators",
            lambda: {k: v[:2] for k, v in grid.get_operators(REGULARISATION).items()},
        )
    _record(out, f"{prefix}/vtk_n_cells", lambda: np.array([grid.vtk().n_cells]))

    # setters, applied in sequence; each records the state it leaves behind
    before = len(calls)

    def state():
        return {"nsteps": grid.nsteps, "step_vector": grid.step_vector, "origin": grid.origin}

    grid.set_nelements(200)
    _record(out, f"{prefix}/after_set_nelements", lambda: {**state(), "n_elements": grid.n_elements})
    grid.nsteps = np.array(grid.nsteps) + 1
    _record(out, f"{prefix}/after_nsteps", state)
    grid.step_vector = np.array(grid.step_vector) * 2
    _record(out, f"{prefix}/after_step_vector", state)
    grid.origin = np.array(grid.origin) - 1
    _record(out, f"{prefix}/after_origin", state)
    grid.maximum = np.array(grid.maximum) + 3
    _record(out, f"{prefix}/after_maximum", lambda: {**state(), "maximum": grid.maximum})
    out[f"{prefix}/onGeometryChange_calls"] = np.array([len(calls) - before])
    grid.rotation_xy = 30
    _record(out, f"{prefix}/after_rotation", lambda: grid.rotation_xy)


def cls_has_ijk(node_ijk) -> bool:
    """The tet meshes return something other than an (N, 3) array here; the ijk helpers do not apply."""
    return isinstance(node_ijk, np.ndarray) and node_ijk.ndim == 2 and node_ijk.shape[-1] == 3


def _assert_same(key, got, want):
    if want.dtype.kind in "US":
        assert got.dtype.kind in "US" and str(got) == str(want), f"{key}: {got} != {want}"
        return
    assert got.dtype.kind not in "US", f"{key}: now raises {got}"
    assert got.shape == want.shape, f"{key}: shape {got.shape} != {want.shape}"
    assert got.dtype.kind == want.dtype.kind, f"{key}: dtype {got.dtype} != {want.dtype}"
    if want.dtype.kind == "f":
        np.testing.assert_allclose(got, want, rtol=1e-9, atol=1e-12, err_msg=key)
    else:
        np.testing.assert_array_equal(got, want, err_msg=key)


@pytest.fixture(scope="module")
def current():
    return collect()


@pytest.fixture(scope="module")
def reference():
    with np.load(REFERENCE) as data:
        return {key: data[key] for key in data.files}


def test_same_probes_are_recorded(current, reference):
    assert sorted(current) == sorted(reference)


@pytest.mark.parametrize("cls_name", list(CLASSES))
def test_structured_support_behaviour_is_unchanged(cls_name, current, reference):
    keys = [k for k in reference if k.startswith(f"{cls_name}/")]
    assert keys
    for key in keys:
        _assert_same(key, current[key], reference[key])


def test_constructor_validation_is_unchanged(current, reference):
    for key in [k for k in reference if k.startswith("validation/")]:
        _assert_same(key, current[key], reference[key])


if __name__ == "__main__":
    if sys.argv[1:] != ["--update"]:
        raise SystemExit(__doc__)
    REFERENCE.parent.mkdir(exist_ok=True)
    result = collect()
    np.savez_compressed(REFERENCE, **result)
    print(f"wrote {len(result)} arrays to {REFERENCE}")
