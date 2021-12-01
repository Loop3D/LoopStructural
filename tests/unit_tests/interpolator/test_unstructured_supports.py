import numpy as np
from LoopStructural.interpolators.unstructured_tetra import UnStructuredTetMesh
from os.path import dirname
file_path = dirname(__file__)

def test_get_elements():
    nodes = np.loadtxt('{}/nodes.txt'.format(file_path))
    elements = np.loadtxt('{}/elements.txt'.format(file_path))
    elements=np.array(elements,dtype='int64')
    neighbours = np.loadtxt('{}/neighbours.txt'.format(file_path))
    
    mesh = UnStructuredTetMesh(nodes,elements,neighbours)
    points = np.random.random((100,3))
    verts, c, tetra, inside = mesh.get_element_for_location(points)

    vertices = nodes[elements,:]
    pos = points[:,:]
    vap = pos[:,None, :] - vertices[None, :, 0, :]
    vbp = pos[:,None, :] - vertices[None, :, 1, :]
    #         # vcp = p - points[:, 2, :]
    #         # vdp = p - points[:, 3, :]
    vab = vertices[None, :, 1, :] - vertices[None, :,  0, :]
    vac = vertices[None, :, 2, :] - vertices[None, :,  0, :]
    vad = vertices[None, :, 3, :] - vertices[None, :,  0, :]
    vbc = vertices[None, :, 2, :] - vertices[None, :,  1, :]
    vbd = vertices[None, :, 3, :] - vertices[None, :, 1, :]

    va = np.einsum('ikj, ikj->ik', vbp, np.cross(vbd, vbc, axisa=2, axisb=2)) / 6.
    vb = np.einsum('ikj, ikj->ik', vap, np.cross(vac, vad, axisa=2, axisb=2)) / 6.
    vc = np.einsum('ikj, ikj->ik', vap, np.cross(vad, vab, axisa=2, axisb=2)) / 6.
    vd = np.einsum('ikj, ikj->ik', vap, np.cross(vab, vac, axisa=2, axisb=2)) / 6.
    v = np.einsum('ikj, ikj->ik', vab, np.cross(vac, vad, axisa=2, axisb=2)) / 6.
    # c = np.zeros((va.shape[0], 4))
    # c[:,  0] = va / v
    # c[:,  1] = vb / v
    # c[:,  2] = vc / v
    # c[:,  3] = vd / v
    c = np.zeros((pos.shape[0],va.shape[1],  4))
    c[:,:,  0] = va / v
    c[:,:,  1] = vb / v
    c[:,:,  2] = vc / v

    c[:,:,  3] = vd / v
    row,col = np.where(np.all(c>=0,axis=2))

    tetra_idx = col
    
    # check if the calculated tetra from the mesh method using aabb
    # is the same as using the barycentric coordinates on all elelemts for 
    # all points
    assert np.all(elements[tetra_idx]-tetra==0)