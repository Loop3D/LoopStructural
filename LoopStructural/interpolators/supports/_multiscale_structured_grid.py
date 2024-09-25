from ._3d_structured_grid import StructuredGrid
import numpy as np
from ...utils import getLogger

logger = getLogger(__name__)


class MultiScaleStructuredGrid(StructuredGrid):
    def __init__(
        self,
        origin=np.zeros(3),
        nsteps=np.array([10, 10, 10]),
        step_vector=np.ones(3),
        rotation_xy=None,
        levels=1,
    ):
        """_summary_

        Parameters
        ----------
        origin : _type_, optional
            _description_, by default np.zeros(3)
        nsteps : _type_, optional
            _description_, by default np.array([10, 10, 10])
        step_vector : _type_, optional
            _description_, by default np.ones(3)
        rotation_xy : _type_, optional
            _description_, by default None
        """
        # we need to check that the nsteps is divisible by 2**levels
        # otherwise we should optimise the nsteps to be divisible by 2**levels
        if levels == 0:
            raise ValueError('levels must be greater than 0, or use a StructuredGrid')
        if levels > 10:
            logger.warning("levels is greater than 10, this may result in a very large grid")
        if not np.all(nsteps % 2**levels == 0):
            # jsut increase the number of steps so that it is divisible by 2**levels
            # it will make the mesh bigger, but this can be considered later.
            nsteps = (np.ceil(nsteps / 2**levels) * 2**levels).astype(int)

        super().__init__(origin, nsteps, step_vector, rotation_xy)
        self.levels = {}
        for l in range(1, levels + 1):
            self.levels[l] = {
                'nsteps': (nsteps / 2 ** (levels - l)).astype(int),
                'step_vector': step_vector * 2 ** (levels - l),
            }
        self.current_level = 1
        self.max_level = levels
        self.set_level(self.current_level)

    def next_level(self):
        if self.current_level < self.max_level:
            self.set_level(self.current_level + 1)
            return True
        else:
            return False

    def set_level(self, level):
        self.nsteps = self.levels[level]['nsteps']
        # self.step_vector = self.levels[level]['step_vector']
        self.current_level = level

    def propagate_property(self, property, to_level):
        from_level = self.current_level
        # set the grid to be the shape of the new level
        # and get the node locations
        self.set_level(to_level)
        new_nodes = self.nodes
        self.set_level(from_level)  # elf.levels[from_level]['nsteps']
        # calculate the global indices of the shared nodes (every 2**level nodes)
        # we dont' want to interpolate these just reuse the values
        xi = np.arange(0, self.levels[to_level]['nsteps'][0], 2 ** (to_level - from_level))
        yi = np.arange(0, self.levels[to_level]['nsteps'][1], 2 ** (to_level - from_level))
        zi = np.arange(0, self.levels[to_level]['nsteps'][2], 2 ** (to_level - from_level))
        ii, jj, kk = np.meshgrid(xi, yi, zi, indexing='ij')
        new_nodes = np.vstack(
            [ii.flatten(order='f'), jj.flatten(order='f'), kk.flatten(order='f')]
        ).T
        self.nsteps = self.levels[to_level]['nsteps']
        gi = self.global_node_indices(new_nodes).astype(int)
        self.set_level(from_level)  # reset to the original level
        # get property value on all nodes
        new_property = self.evaluate_value(
            new_nodes,
            property,
        )
        # set all of the shared nodes to the previous value
        new_property[gi] = property
        return new_property

    @classmethod
    def from_structured_grid(cls, grid, levels=1):
        new_grid = cls(
            origin=grid.origin,
            nsteps=grid.nsteps,
            step_vector=grid.step_vector,
            rotation_xy=grid.rotation_xy,
            levels=levels,
        )
        return new_grid
