import numpy as np
from loop_common.supports import UnStructuredTetMesh
from loop_common.math import rng
from os.path import dirname

file_path = dirname(__file__)


def _brute_force_tetra(nodes, elements, points):
    """Reference point-in-tetra lookup that tests every element directly,
    used to check the aabb-accelerated lookup for correctness."""
    vertices = nodes[elements, :]
    pos = points[:, :]
    vap = pos[:, None, :] - vertices[None, :, 0, :]
    vbp = pos[:, None, :] - vertices[None, :, 1, :]
    vab = vertices[None, :, 1, :] - vertices[None, :, 0, :]
    vac = vertices[None, :, 2, :] - vertices[None, :, 0, :]
    vad = vertices[None, :, 3, :] - vertices[None, :, 0, :]
    vbc = vertices[None, :, 2, :] - vertices[None, :, 1, :]
    vbd = vertices[None, :, 3, :] - vertices[None, :, 1, :]

    va = np.einsum("ikj, ikj->ik", vbp, np.cross(vbd, vbc, axisa=2, axisb=2)) / 6.0
    vb = np.einsum("ikj, ikj->ik", vap, np.cross(vac, vad, axisa=2, axisb=2)) / 6.0
    vc = np.einsum("ikj, ikj->ik", vap, np.cross(vad, vab, axisa=2, axisb=2)) / 6.0
    vd = np.einsum("ikj, ikj->ik", vap, np.cross(vab, vac, axisa=2, axisb=2)) / 6.0
    v = np.einsum("ikj, ikj->ik", vab, np.cross(vac, vad, axisa=2, axisb=2)) / 6.0
    c = np.zeros((pos.shape[0], va.shape[1], 4))
    c[:, :, 0] = va / v
    c[:, :, 1] = vb / v
    c[:, :, 2] = vc / v
    c[:, :, 3] = vd / v
    found = np.all(c >= 0, axis=2)
    inside = np.any(found, axis=1)
    tetra_idx = np.argmax(found, axis=1)
    return inside, tetra_idx


def _load_mesh():
    nodes = np.loadtxt("{}/nodes.txt".format(file_path))
    elements = np.loadtxt("{}/elements.txt".format(file_path))
    elements = np.array(elements, dtype="int64")
    neighbours = np.loadtxt("{}/neighbours.txt".format(file_path))
    return nodes, elements, neighbours


def test_get_elements():
    nodes, elements, neighbours = _load_mesh()
    mesh = UnStructuredTetMesh(nodes, elements, neighbours)
    points = rng.random((100, 3))
    verts, c, tetra, inside = mesh.get_element_for_location(points)

    _, tetra_idx = _brute_force_tetra(nodes, elements, points)

    # check if the calculated tetra from the mesh method using aabb
    # is the same as using the barycentric coordinates on all elelemts for
    # all points
    assert np.all(elements[tetra_idx[inside]] - elements[tetra[inside]] == 0)


def test_get_elements_outside_bounds():
    # points partly outside the mesh's bounding box exercise the aabb grid's
    # "outside" branch, which previously misaligned the compacted candidate
    # indices back onto the query points and silently corrupted results.
    nodes, elements, neighbours = _load_mesh()
    mesh = UnStructuredTetMesh(nodes, elements, neighbours)

    local_rng = np.random.default_rng(42)
    n = 600
    inside_pts = mesh.minimum + local_rng.random((n // 2, 3)) * (mesh.maximum - mesh.minimum)
    outside_pts = mesh.maximum + 10 + local_rng.random((n // 2, 3)) * 5
    points = np.vstack([inside_pts, outside_pts])
    local_rng.shuffle(points)

    verts, c, tetra, inside = mesh.get_element_for_location(points)
    brute_inside, brute_tetra = _brute_force_tetra(nodes, elements, points)

    assert np.array_equal(inside, brute_inside)
    assert np.all(elements[tetra[inside]] - elements[brute_tetra[inside]] == 0)
    # barycentric coordinates of located points should sum to 1
    assert np.allclose(c[inside].sum(axis=1), 1.0)


def test_get_elements_chunk_boundary():
    # get_element_for_location processes points in blocks of 1e4; use more
    # than one block, with some points outside the mesh, to make sure the
    # per-block results are stitched back together at the right offsets.
    nodes, elements, neighbours = _load_mesh()
    mesh = UnStructuredTetMesh(nodes, elements, neighbours)

    local_rng = np.random.default_rng(7)
    n = 25_000
    inside_pts = mesh.minimum + local_rng.random((n * 3 // 4, 3)) * (mesh.maximum - mesh.minimum)
    outside_pts = mesh.maximum + 10 + local_rng.random((n // 4, 3)) * 5
    points = np.vstack([inside_pts, outside_pts])
    local_rng.shuffle(points)

    verts, c, tetra, inside = mesh.get_element_for_location(points)
    brute_inside, brute_tetra = _brute_force_tetra(nodes, elements, points)

    assert np.array_equal(inside, brute_inside)
    assert np.all(elements[tetra[inside]] - elements[brute_tetra[inside]] == 0)


def test_get_elements_small_mesh_real_world_scale():
    # meshes with fewer than 2000 elements used to fall back to a hardcoded
    # 1-unit aabb grid regardless of the mesh's actual coordinate scale,
    # silently breaking lookups for any mesh not roughly 1 unit across.
    nodes, elements, neighbours = _load_mesh()
    keep = np.unique(elements[:500].ravel())
    remap = -np.ones(nodes.shape[0], dtype=np.int64)
    remap[keep] = np.arange(keep.shape[0])
    small_elements = remap[elements[:500]]
    small_neighbours = np.where(
        np.isin(neighbours[:500], np.arange(500)), neighbours[:500], -1
    )
    small_nodes = nodes[keep] * 1000.0  # push coordinates to a "real world" scale

    mesh = UnStructuredTetMesh(small_nodes, small_elements, small_neighbours)
    assert mesh.n_elements < 2000

    points = small_nodes[small_elements[:, :4]].mean(axis=1)  # element barycentres
    verts, c, tetra, inside = mesh.get_element_for_location(points)

    assert np.all(inside)
    assert np.allclose(c.sum(axis=1), 1.0)
