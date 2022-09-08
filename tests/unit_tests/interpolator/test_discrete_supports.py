from LoopStructural.interpolators import TetMesh
from LoopStructural.interpolators import StructuredGrid, StructuredGrid2D
import numpy as np

## structured grid tests
def test_create_structured_grid():
    grid = StructuredGrid()


def test_create_structured_grid_origin_nsteps():
    grid = StructuredGrid(origin=np.zeros(3), nsteps=np.array([5, 5, 5]))
    assert grid.n_nodes == 5 * 5 * 5
    assert np.sum(grid.maximum - np.ones(3) * 5) == 0


def test_create_structured_grid_origin_nsteps():
    grid = StructuredGrid(
        origin=np.zeros(3),
        nsteps=np.array([10, 10, 10]),
        step_vector=np.array([0.1, 0.1, 0.1]),
    )
    assert np.sum(grid.step_vector - np.array([0.1, 0.1, 0.1])) == 0
    assert np.sum(grid.maximum - np.ones(3)) == 0


def test_evaluate_value():
    grid = StructuredGrid()
    assert (
        np.sum(
            grid.barycentre[:, 0]
            - grid.evaluate_value(grid.barycentre, grid.nodes[:, 0])
        )
        == 0
    )


def test_evaluate_gradient():
    grid = StructuredGrid()
    # test by setting the scalar field to the y coordinate
    vector = grid.evaluate_gradient(grid.barycentre, grid.nodes[:, 1])
    assert np.sum(vector - np.array([0, 1, 0])) == 0

    # same test but for a bigger grid, making sure scaling for cell is ok
    grid = StructuredGrid(step_vector=np.array([100, 100, 100]))
    vector = grid.evaluate_gradient(grid.barycentre, grid.nodes[:, 1])
    assert np.sum(vector - np.array([0, 1, 0])) == 0

def test_outside_box():
    grid = StructuredGrid()
    # test by setting the scalar field to the y coordinate
    inside = grid.inside(grid.barycentre+5)
    assert np.all(~inside == np.any((grid.barycentre+5) > grid.maximum,axis=1))
    inside = grid.inside(grid.barycentre-5)
    assert np.all(~inside == np.any((grid.barycentre-5) < grid.origin,axis=1))

    cix, ciy, ciz = grid.position_to_cell_index(grid.barycentre-5)
    assert np.all(cix[inside] < grid.nsteps_cells[0])
    assert np.all(ciy[inside] < grid.nsteps_cells[1])
    assert np.all(ciz[inside] < grid.nsteps_cells[2])
    cornersx, cornersy, cornersz = grid.cell_corner_indexes(cix, ciy, ciz)
    assert np.all(cornersx[inside] < grid.nsteps[0])
    assert np.all(cornersy[inside] < grid.nsteps[1])
    assert np.all(cornersz[inside] < grid.nsteps[2])
    globalidx = grid.global_indicies(np.dstack([cornersx, cornersy, cornersz]).T)
    # print(globalidx[inside],grid.n_nodes,inside)
    assert np.all(globalidx[inside] < grid.n_nodes)
    inside = grid.inside(grid.barycentre-5)
    inside, grid.position_to_cell_corners(grid.barycentre-5)
    vector = grid.evaluate_gradient(grid.barycentre-5, grid.nodes[:, 1])
    assert np.sum(np.mean(vector[inside,:],axis=0) - np.array([0, 1, 0])) == 0
    vector = grid.evaluate_gradient( grid.nodes, grid.nodes[:, 1])

def test_evaluate_gradient2():
    # this test is the same as above but we will use a random vector
    np.random.seed(0)
    for i in range(10):
        step = np.random.uniform(0, 100)
        grid = StructuredGrid(step_vector=np.array([step, step, step]))

        # define random vector
        n = np.random.random(3)
        n /= np.linalg.norm(n)
        distance = (
            n[0] * grid.nodes[:, 0] + n[1] * grid.nodes[:, 1] + n[2] * grid.nodes[:, 2]
        )
        vector = grid.evaluate_gradient(
            np.random.uniform(1, 8, size=(100, 3)), distance
        )
        assert np.all(np.isclose(np.sum(vector - n[None, :], axis=1), 0)) == True
    assert i == 9


def test_get_element():
    grid = StructuredGrid()
    point = grid.barycentre[[0], :]
    idc, inside = grid.position_to_cell_corners(point)
    bary = np.mean(grid.nodes[idc, :], axis=0)
    assert np.sum(point - bary) == 0


def test_global_to_local_coordinates():
    grid = StructuredGrid()
    point = np.array([[1.2, 1.5, 1.7]])
    lx, ly, lz = grid.position_to_local_coordinates(point)
    assert np.isclose(lx[0], 0.2)
    assert np.isclose(ly[0], 0.5)
    assert np.isclose(lz[0], 0.7)


def test_get_element_outside():
    grid = StructuredGrid()
    point = np.array([grid.origin - np.ones(3)])
    idc, inside = grid.position_to_cell_corners(point)
    assert inside[0] == False

## structured grid 2d tests
def test_create_structured_grid2d():
    grid = StructuredGrid2D()
    
def test_create_structured_grid2d_origin_nsteps():
    grid = StructuredGrid2D(origin=np.zeros(2),nsteps=np.array([5,5]))
    assert grid.n_nodes == 5*5
    assert np.sum(grid.maximum- np.ones(2)*5) == 0  

def test_create_structured_grid2d_origin_nsteps():
    grid = StructuredGrid2D(origin=np.zeros(2),
                        nsteps=np.array([10,10]),
                        step_vector=np.array([0.1,0.1]))
    assert np.sum(grid.step_vector - np.array([0.1,0.1])) == 0
    assert np.sum(grid.maximum - np.ones(2)) == 0
  
def test_evaluate_value_2d():
    grid = StructuredGrid2D()
    grid.update_property('X',grid.nodes[:,0])
    assert np.sum(grid.barycentre[:,0] - 
    grid.evaluate_value(grid.barycentre,'X')) ==0

def test_evaluate_gradient_2d():
    grid = StructuredGrid2D()
    grid.update_property('Y',grid.nodes[:,1])
    vector = np.mean(grid.evaluate_gradient(grid.barycentre,'Y'),axis=0)
    # vector/=np.linalg.norm(vector)
    assert np.sum(vector-np.array([0,grid.step_vector[1]])) == 0
    
def test_get_element_2d():
    grid = StructuredGrid2D()
    point = grid.barycentre[[0],:]
    idc, inside = grid.position_to_cell_corners(point)
    bary = np.mean(grid.nodes[idc,:],axis=0)
    assert np.sum(point-bary) == 0  

def test_global_to_local_coordinates2d():
    grid = StructuredGrid2D()
    point = np.array([[1.2,1.5,1.7]])
    lx,ly = grid.position_to_local_coordinates(point)
    assert(np.isclose(lx[0],.2))
    assert(np.isclose(ly[0], .5))

def test_get_element_outside2d():
    grid = StructuredGrid2D()
    point = np.array([grid.origin - np.ones(2)])
    idc, inside = grid.position_to_cell_corners(point)
    assert inside[0] == False
## structured tetra tests
def test_create_testmesh():
    grid = TetMesh()


def test_create_testmesh_origin_nsteps():
    grid = TetMesh(origin=np.zeros(3), nsteps=np.array([5, 5, 5]))
    assert grid.n_nodes == 6 * 6 * 6
    assert np.sum(grid.maximum - np.ones(3) * 5) == 0


def test_create_tetmesh_origin_nsteps():
    grid = TetMesh(
        origin=np.zeros(3),
        nsteps=np.array([10, 10, 10]),
        step_vector=np.array([0.1, 0.1, 0.1]),
    )
    assert np.sum(grid.step_vector - np.array([0.1, 0.1, 0.1])) == 0
    assert np.sum(grid.maximum - np.ones(3)) == 0


# def test_trilinear():

# def test_position_to_local():

# def test_neighbours():


def test_evaluate_value_tetmesh():
    grid = TetMesh()
    assert (
        np.sum(
            grid.barycentre[:, 0]
            - grid.evaluate_value(grid.barycentre, grid.nodes[:, 0])
        )
        == 0
    )


def test_evaluate_gradient_tetmesh():
    grid = TetMesh()
    vector = np.mean(
        grid.evaluate_gradient(grid.barycentre, grid.nodes[:, 1]), axis=0
    )
    # vector/=np.linalg.norm(vector)
    assert np.sum(vector - np.array([0, grid.step_vector[1], 0])) == 0


def test_change_origin():
    grid = StructuredGrid(origin=np.zeros(3), nsteps=np.array([5, 5, 5]))
    grid.origin = np.array([-1, -1, -1])
    assert np.all(grid.origin == np.array([-1, -1, -1]))
    assert np.all(grid.nsteps == np.array([6, 6, 6]))
    assert np.all(grid.step_vector == np.ones(3))


def test_change_maximum():
    grid = StructuredGrid(origin=np.zeros(3), nsteps=np.array([5, 5, 5]))
    grid.maximum = np.array([7, 7, 7])
    assert np.all(np.isclose(grid.nsteps,np.array([8, 8, 8])))
    assert np.all(np.isclose(grid.maximum, np.array([7, 7, 7])))
    assert np.all(np.isclose(grid.step_vector,np.ones(3)))


def test_change_maximum_and_origin():
    grid = StructuredGrid(origin=np.zeros(3), nsteps=np.array([5, 5, 5]))
    grid.origin = np.array([-1.0, -1.0, -1.0])
    assert np.all(np.isclose(grid.origin,np.array([-1, -1, -1])))
    assert np.all(np.isclose(grid.nsteps, np.array([6, 6, 6])))
    assert np.all(np.isclose(grid.step_vector, np.ones(3)))
    grid.maximum = np.array([7.0, 7.0, 7.0])
    assert np.all(np.isclose(grid.nsteps, np.array([9, 9, 9])))
    assert np.all(np.isclose(grid.maximum, np.array([7.0, 7.0, 7.0])))
    assert np.all(np.isclose(grid.step_vector,np.ones(3)))

if __name__ == "__main__":
    test_create_structured_grid()
    test_create_structured_grid2d()
    test_create_structured_grid2d_origin_nsteps()
    test_change_maximum()
    test_change_maximum_and_origin()
    test_change_origin()
    test_create_structured_grid_origin_nsteps()
    test_create_testmesh()
    test_evaluate_gradient()
    test_evaluate_gradient2()
    test_outside_box()
    test_create_testmesh_origin_nsteps()
    