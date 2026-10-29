"""Pin the behaviour of the 2D structured interpolation support.

Records a wide set of ``StructuredGrid2D`` outputs (geometry, node/element arrays, point location
including exact boundary points, bilinear weights, value/gradient evaluation, element lookup,
regularisation operators, index conversions, ``set_nelements``) on a few grids and compares them
with ``data/structured_grid_2d_reference.npz``.

The reference was generated from the code *before* ``StructuredGrid2D`` was changed to extend
``StructuredGrid2DGeometry``. Regenerate after an intended change, and review the diff:

    python packages/loop_common/tests/test_structured_grid_2d_reference.py --update
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
from loop_common.supports import StructuredGrid2D

REFERENCE = Path(__file__).parent / "data" / "structured_grid_2d_reference.npz"

# (origin, nodes per axis, step vector). StructuredGrid2D takes node counts.
GRIDS = {
    "unit": ([0.0, 0.0], [5, 5], [1.0, 1.0]),
    "aniso": ([10.0, -3.0], [7, 4], [0.5, 2.0]),
    "minimal": ([0.0, 0.0], [2, 2], [3.0, 1.0]),
}
ATTRIBUTES = (
    "nsteps",
    "nsteps_cells",
    "maximum",
    "origin",
    "step_vector",
    "n_nodes",
    "n_elements",
    "element_size",
    "n_cell_x",
    "n_cell_y",
    "dim",
    "dimension",
    "nodes",
    "elements",
    "barycentre",
)
OPERATOR_WEIGHTS = {"dxx": 1.0, "dyy": 1.0, "dxy": 1.0}


def _flatten(value, prefix, out):
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
    out: dict[str, np.ndarray] = {}
    rng = np.random.default_rng(5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name, (origin, nodes, step) in GRIDS.items():
            grid = StructuredGrid2D(origin=origin, nsteps=nodes, step_vector=step)
            _probe(out, name, grid, origin, nodes, step, rng)
    return out


def _probe(out, prefix, grid, origin, nodes, step, rng):
    for attribute in ATTRIBUTES:
        _record(out, f"{prefix}/{attribute}", lambda attribute=attribute: getattr(grid, attribute))
    _record(out, f"{prefix}/str_type", lambda: np.asarray(int(grid.type)))

    low = np.array(origin) - 1.0
    high = np.array(origin) + np.array(nodes) * np.array(step) + 1.0
    anywhere = rng.uniform(low, high, (120, 2))
    # points exactly on the boundary of the nodes' extent and of `maximum`
    node_max = np.array(origin) + (np.array(nodes) - 1) * np.array(step)
    corners = np.array(np.meshgrid(*zip(origin, node_max), indexing="ij")).reshape(2, -1).T
    faces = []
    centre = 0.5 * (np.array(origin) + node_max)
    for axis in range(2):
        for bound in (origin[axis], node_max[axis]):
            point = centre.copy()
            point[axis] = bound
            faces.append(point)
    anywhere = np.vstack([anywhere, corners, np.array(faces)])
    interior = rng.uniform(np.array(origin) + 1e-3, node_max - 1e-3, (80, 2))

    _record(out, f"{prefix}/inside", lambda: grid.inside(anywhere))
    _record(out, f"{prefix}/position_to_cell_index", lambda: grid.position_to_cell_index(anywhere))
    _record(out, f"{prefix}/check_position_1d", lambda: grid.check_position(np.array([1.0, 2.0])))
    _record(out, f"{prefix}/check_position_3d", lambda: grid.check_position(np.zeros((2, 2, 2))))
    _record(out, f"{prefix}/position_to_cell_corners", lambda: grid.position_to_cell_corners(anywhere))
    _record(out, f"{prefix}/position_to_cell_vertices", lambda: grid.position_to_cell_vertices(interior))
    _record(out, f"{prefix}/position_to_local_coordinates", lambda: grid.position_to_local_coordinates(interior))
    _record(out, f"{prefix}/position_to_dof_coefs", lambda: grid.position_to_dof_coefs(interior))
    _record(out, f"{prefix}/bilinear", lambda: grid.bilinear(rng.uniform(0, 1, (10, 2))))

    cell_index = np.arange(grid.n_elements)
    node_index = np.arange(grid.n_nodes)
    _record(out, f"{prefix}/cell_centres", lambda: grid.cell_centres(cell_index))
    _record(out, f"{prefix}/global_index_to_cell_index", lambda: grid.global_index_to_cell_index(cell_index))
    node_ijk = grid.global_index_to_node_index(node_index)
    _record(out, f"{prefix}/global_index_to_node_index", lambda: node_ijk)
    cell_ijk = grid.global_index_to_cell_index(cell_index)
    _record(out, f"{prefix}/cell_corner_indexes", lambda: grid.cell_corner_indexes(cell_ijk))
    _record(out, f"{prefix}/global_node_indices", lambda: grid.global_node_indices(node_ijk))
    _record(out, f"{prefix}/global_cell_indices", lambda: grid.global_cell_indices(cell_ijk))
    _record(out, f"{prefix}/node_indexes_to_position", lambda: grid.node_indexes_to_position(node_ijk[:10]))
    _record(out, f"{prefix}/neighbour_global_indexes", lambda: grid.neighbour_global_indexes())

    values = rng.normal(size=grid.n_nodes)
    _record(out, f"{prefix}/evaluate_value", lambda: grid.evaluate_value(anywhere, values))
    _record(out, f"{prefix}/evaluate_gradient", lambda: grid.evaluate_gradient(interior, values))
    _record(out, f"{prefix}/get_element_for_location", lambda: grid.get_element_for_location(interior))
    _record(
        out,
        f"{prefix}/get_element_gradient_for_location",
        lambda: grid.get_element_gradient_for_location(interior),
    )
    _record(
        out,
        f"{prefix}/get_operators",
        lambda: dict(grid.get_operators(OPERATOR_WEIGHTS)),
    )
    _record(out, f"{prefix}/vtk_n_cells", lambda: np.array([grid.vtk(0.0).n_cells]))
    _record(out, f"{prefix}/set_nelements", lambda: grid.set_nelements(100))


def _assert_same(key, got, want):
    if want.dtype.kind in "US":
        assert got.dtype.kind in "US" and str(got) == str(want), f"{key}: {got} != {want}"
        return
    assert got.dtype.kind not in "US", f"{key}: now raises {got}"
    assert got.shape == want.shape, f"{key}: shape {got.shape} != {want.shape}"
    assert got.dtype.kind == want.dtype.kind, f"{key}: dtype {got.dtype} != {want.dtype}"
    if want.dtype.kind == "f":
        np.testing.assert_allclose(got, want, rtol=1e-9, atol=1e-12, err_msg=key, equal_nan=True)
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


@pytest.mark.parametrize("grid_name", list(GRIDS))
def test_structured_grid_2d_behaviour_is_unchanged(grid_name, current, reference):
    keys = [k for k in reference if k.startswith(f"{grid_name}/")]
    assert keys
    for key in keys:
        _assert_same(key, current[key], reference[key])


if __name__ == "__main__":
    if sys.argv[1:] != ["--update"]:
        raise SystemExit(__doc__)
    REFERENCE.parent.mkdir(exist_ok=True)
    result = collect()
    np.savez_compressed(REFERENCE, **result)
    print(f"wrote {len(result)} arrays to {REFERENCE}")
