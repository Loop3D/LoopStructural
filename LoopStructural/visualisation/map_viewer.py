import matplotlib.pyplot as plt
import numpy as np
from ..utils.helper import normal_vector_to_strike_and_dip
from ..utils.utils import strike_symbol

class MapView:
    def __init__(self, origin, maximum, nsteps,**kwargs):
        """

        Parameters
        ----------
        origin - lower left
        maximum - upper right
        nsteps - number of cells
        kwargs
        """
        self.origin = origin
        self.maximum = maximum
        self.nsteps = nsteps
        x = np.linspace(self.origin[0],self.maximum[0],self.nsteps[0])
        y = np.linspace(self.origin[1],self.maximum[1],self.nsteps[1])
        self.xx, self.yy = np.meshgrid(x,y,indexing='ij')
        self.xx=self.xx.flatten()
        self.yy = self.yy.flatten()
        self.fig, self.ax = plt.subplots(1, figsize=(10, 10))

    def draw_strike(self, x, y, strike, scale=100, colour='black'):
        """
        Draw a strike symbol on the map
        Parameters
        ----------
        x
        y
        strike
        scale
        colour

        Returns
        -------

        """
        rotated, r2 = strike_symbol(-strike)
        rotated *= scale
        r2 *= scale
        self.ax.plot([x, x + rotated[0]], [y, y + rotated[1]], colour)
        self.ax.plot([x - rotated[0], x], [y - rotated[1], y], colour)
        self.ax.plot([x, x + r2[0]], [y, y + r2[1]], colour)

    def add_data(self, feature, **kwargs):
        """
        Adds the data associated to the feature to the plot
        Parameters
        ----------
        feature geological feature
        kwargs are passed to matplotlib functions and draw strike
        Returns
        -------

        """
        gradient_data = feature.support.interpolator.get_gradient_control()
        value_data = feature.support.interpolator.get_control_points()
        self.ax.scatter(value_data[:,0],value_data[:,1],c=value_data[:,3],
                        vmin=feature.min(),vmax=feature.max(),**kwargs)
        #points = strati.support.interpolator.get_gradient_control()
        strike = normal_vector_to_strike_and_dip(gradient_data[:, 3:])
        for i in range(len(strike)):
            self.draw_strike(gradient_data[i, 0], gradient_data[i, 1], -strike[i, 0],**kwargs)

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
        v = feature.evaluate_value(np.array([self.xx, self.yy,zz]).T)
        self.ax.imshow(v.reshape(self.nsteps),
                   extent=[self.origin[0],self.maximum[0],self.origin[1],self.maximum[1]],
                       vmin=feature.min(),vmax=feature.max(),**kwargs)

    def add_contour(self, feature):
        pass