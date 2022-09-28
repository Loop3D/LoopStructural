from LoopStructural.interpolators import StructuredGrid
import numpy as np


class ScalarField:
    """A scalar field defined by a regular grid and values"""

    def __init__(
        self,
        values,
        step_vector,
        nsteps,
        origin=np.zeros(3),
        name="scalar field",
        mode="node",
    ):
        """[summary]

        Parameters
        ----------
        values : numpy array
            nodes values of regular grid
        step_vector : numpy array
            step vector for x,y,z step
        nsteps : numpy array
            number of steps
        name : string, optional
            name of the feature for the visualisation
        mode : string, optional
            property stored on nodes or as cell values, 'node' or 'cell'
        """
        self.values = values
        self.grid = StructuredGrid(origin, nsteps, step_vector)
        self.name = name
        self.mode = mode

    @property
    def nodes(self):
        return self.grid.nodes

    def evaluate_value(self, xyz):
        """Evaluate the scalar field at locations

        Parameters
        ----------
        xyz : numpy array
            locations in real coordinates

        Returns
        -------
        numpy array
            interpolated values, same shape as xyz
        """
        xyz = self.grid.check_position(xyz)
        if self.mode == "node":
            return self.grid.evaluate_value(xyz, self.values)
        if self.mode == "grid":
            indexes = np.array(self.grid.position_to_cell_index(xyz))
            inside = np.all(indexes >= 0, axis=0)
            inside = np.logical_and(
                np.all(indexes < self.grid.nsteps_cells[:, None], axis=0), inside
            )
            v = np.zeros(xyz.shape[0], dtype=float)
            v[:] = np.nan
            v[inside] = self.values[self.grid.global_cell_indicies(indexes[:, inside])]

            return v

    def __call__(self, xyz):
        return self.evaluate_value(xyz)

    def min(self):
        return np.min(self.values)

    def max(self):
        return np.max(self.values)  # 14
