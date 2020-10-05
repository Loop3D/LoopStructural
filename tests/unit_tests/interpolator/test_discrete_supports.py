from LoopStructural.interpolators import TetMesh
from LoopStructural.interpolators import StructuredGrid
import numpy as np

## structured grid tests
def test_create_structured_grid():
    grid = StructuredGrid()
    
def test_create_structured_grid_origin_nsteps():
    grid = StructuredGrid(origin=np.zeros(3),nsteps=np.array([5,5,5]))
    assert grid.n_nodes == 5*5*5
    assert np.sum(grid.maximum- np.ones(3)*5) == 0  

def test_create_structured_grid_origin_nsteps():
    grid = StructuredGrid(origin=np.zeros(3),
                        nsteps=np.array([10,10,10]),
                        step_vector=np.array([0.1,0.1,0.1]))
    assert np.sum(grid.step_vector - np.array([0.1,0.1,0.1])) == 0
    assert np.sum(grid.maximum - np.ones(3)) == 0
  
def test_evaluate_value():
    grid = StructuredGrid()
    grid.update_property('X',grid.nodes[:,0])
    assert np.sum(grid.barycentre()[:,0] - 
    grid.evaluate_value(grid.barycentre(),'X')) ==0

def test_evaluate_gradient():
    grid = StructuredGrid()
    grid.update_property('Y',grid.nodes[:,1])
    vector = np.mean(grid.evaluate_gradient(grid.barycentre(),'Y'),axis=0)
    # vector/=np.linalg.norm(vector)
    assert np.sum(vector-np.array([0,grid.step_vector[1],0])) == 0
    
def test_get_element():
    grid = StructuredGrid()
    point = grid.barycentre()[[0],:]
    idc, inside = grid.position_to_cell_corners(point)
    bary = np.mean(grid.nodes[idc,:],axis=0)
    assert np.sum(point-bary) == 0  

def test_global_to_local_coordinates():
    grid = StructuredGrid()
    point = np.array([[1.2,1.5,1.7]])
    lx,ly, lz = grid.position_to_local_coordinates(point)
    assert(np.isclose(lx[0],.2))
    assert(np.isclose(ly[0], .5))
    assert(np.isclose(lz[0], .7))

def test_get_element_outside():
    grid = StructuredGrid()
    point = np.array([grid.origin - np.ones(3)])
    idc, inside = grid.position_to_cell_corners(point)
    assert inside[0] == False
## structured tetra tests
def test_create_testmesh():
    grid = TetMesh()
    
def test_create_testmesh_origin_nsteps():
    grid = TetMesh(origin=np.zeros(3),nsteps=np.array([5,5,5]))
    assert grid.n_nodes == 5*5*5
    assert np.sum(grid.maximum- np.ones(3)*5) == 0  

def test_create_tetmesh_origin_nsteps():
    grid = TetMesh(origin=np.zeros(3),
                        nsteps=np.array([10,10,10]),
                        step_vector=np.array([0.1,0.1,0.1]))
    assert np.sum(grid.step_vector - np.array([0.1,0.1,0.1])) == 0
    assert np.sum(grid.maximum - np.ones(3)) == 0


# def test_trilinear():

# def test_position_to_local():

# def test_neighbours():

def test_evaluate_value_tetmesh():
    grid = TetMesh()
    grid.update_property('X',grid.nodes[:,0])
    assert np.sum(grid.barycentre()[:,0] - 
    grid.evaluate_value(grid.barycentre(),'X')) ==0

def test_evaluate_gradient_tetmesh():
    grid = TetMesh()
    grid.update_property('Y',grid.nodes[:,1])
    vector = np.mean(grid.evaluate_gradient(grid.barycentre(),'Y'),axis=0)
    # vector/=np.linalg.norm(vector)
    assert np.sum(vector-np.array([0,grid.step_vector[1],0])) == 0