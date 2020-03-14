import logging

import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)


class MapView:
    def __init__(self, model = None, bounding_box=None, nsteps=None, **kwargs):
        """

        Parameters
        ----------
        origin - lower left
        maximum - upper right
        nsteps - number of cells
        kwargs
        """
        self.bounding_box = bounding_box
        self.nsteps = nsteps
        if model is not None:
            self.bounding_box = model.bounding_box
            self.nsteps = model.nsteps

        x = np.linspace(self.bounding_box[0,0], self.bounding_box[1,0], self.nsteps[0])
        y = np.linspace(self.bounding_box[0,1], self.bounding_box[1,1], self.nsteps[1])
        self.xx, self.yy = np.meshgrid(x, y, indexing='ij')
        self.xx = self.xx.flatten()
        self.yy = self.yy.flatten()
        self.fig, self.ax = plt.subplots(1, figsize=(10, 10))
        self.ax.set_aspect('equal', adjustable='box')

    def add_data(self, feature, val=True, grad=True, **kwargs):
        """
        Adds the data associated to the feature to the plot
        Parameters
        ----------
        feature geological feature
        kwargs are passed to matplotlib functions and draw strike
        Returns
        -------

        """

        ori_data = []
        gradient_data = feature.support.interpolator.get_gradient_constraints()
        if gradient_data.shape[0] > 0:
            ori_data.append(gradient_data)
        norm_data = feature.support.interpolator.get_norm_constraints()
        if norm_data.shape[0] > 0:
            ori_data.append(norm_data)
        cmap = kwargs.pop('cmap','rainbow')
        # if single colour then specify kwarg, otherwise use point value
        if val:
            value_data = feature.support.interpolator.get_value_constraints()
            point_colour = kwargs.pop('point_colour',None)
            if point_colour is None:
                self.ax.scatter(value_data[:, 0], value_data[:, 1], c=value_data[:,3],
                            vmin=feature.min(), vmax=feature.max(),cmap=cmap)
            if point_colour is not None:
                self.ax.scatter(value_data[:, 0], value_data[:, 1], c=point_colour)
        # points = strati.support.interpolator.get_gradient_control()
        if grad:
            symb_colour = kwargs.pop('symb_colour','black')
            gradient_data = np.hstack(ori_data)
            gradient_data[:, 3:5] /= np.linalg.norm(gradient_data[:, 3:5], axis=1)[:, None]
            t = gradient_data[:, [4, 3]] * np.array([1, -1]).T
            n = gradient_data[:, 3:5]
            t *= 0.01
            n *= 0.005
            p1 = gradient_data[:, [0, 1]] - t
            p2 = gradient_data[:, [0, 1]] + t
            # plt.scatter(val[:,0],val[:,1],c='black')
            self.ax.plot([p1[:, 0], p2[:, 0]], [p1[:, 1], p2[:, 1]], symb_colour)
            p1 = gradient_data[:, [0, 1]]
            p2 = gradient_data[:, [0, 1]] + n
            self.ax.plot([p1[:, 0], p2[:, 0]], [p1[:, 1], p2[:, 1]], symb_colour)
            # plt.gca().
            # strike = normal_vector_to_strike_and_dip(gradient_data[:, 3:6])
            # for i in range(len(strike)):
            #     self.draw_strike(gradient_data[i, 0], gradient_data[i, 1],
            #                      -strike[i, 0], **kwargs)

    def add_scalar_field(self, feature, z=0, **kwargs):
        """
        Draw the
        Parameters
        ----------
        feature
        z
        kwargs

        Returns
        -------

        """
        zz = np.zeros(self.xx.shape)
        zz[:] = z
        v = feature.evaluate_value(np.array([self.xx, self.yy, zz]).T)
        self.ax.imshow(v.reshape(self.nsteps).T,
                       extent=[self.bounding_box[0,0], self.bounding_box[1,0], self.bounding_box[0,1],
                               self.bounding_box[1,1]],
                       vmin=feature.min(), vmax=feature.max(),
                       **kwargs)

    def add_contour(self, feature, values, z=0):
        zz = np.zeros(self.xx.shape)
        zz[:] = z
        v = feature.evaluate_value(np.array([self.xx, self.yy, zz]).T)
        self.ax.contour(np.rot90(v.reshape(self.nsteps),1), levels=values,
                       extent=[self.bounding_box[0, 0], self.bounding_box[1, 0], self.bounding_box[0, 1],
                               self.bounding_box[1, 1]],
                       )

        pass
