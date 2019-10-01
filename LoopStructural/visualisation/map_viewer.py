import matplotlib.pyplot as plt
import numpy as np
from ..utils.helper import normal_vector_to_strike_and_dip


class MapView:
    def __init__(self,origin,maximum,nsteps,**kwargs):
        self.origin = origin
        self.maximum = maximum
        x = np.linspace(self.origin[0],self.maximum[0],self.nsteps[0])
        y = np.linspace(self.origin[1],self.maximum[1],self.nsteps[1])
        self.xx, self.yy = np.meshgrid(x,y,indexing='ij')
        self.xx=self.xx.flatten()
        self.yy = self.yy.flatten()
    def add_data_to_plot(self, feature, **kwargs):
        """
        Adds the data associated to the feature to the plot
        Parameters
        ----------
        feature geological feature

        Returns
        -------

        """
        gradient_data = feature.interpolator.get_gradient_control()
        value_data = feature.interpolator.get_control_points()
        plt.scatter(value_data[:,0],value_data[:,1],c=value_data[:,3])

        def strike_symbol(strike):
            R = np.zeros((2, 2))
            R[0, 0] = np.cos(np.deg2rad(-strike))
            R[0, 1] = -np.sin(np.deg2rad(-strike))
            R[1, 0] = np.sin(np.deg2rad(-strike))
            R[1, 1] = np.cos(np.deg2rad(-strike))
            R = np.zeros((2, 2))
            R[0, 0] = np.cos(np.deg2rad(-strike))
            R[0, 1] = -np.sin(np.deg2rad(-strike))
            R[1, 0] = np.sin(np.deg2rad(-strike))
            R[1, 1] = np.cos(np.deg2rad(-strike))

            vec = np.array([0, 1])
            rotated = R @ vec
            vec2 = np.array([-0.5, 0])
            r2 = R @ vec2
            return rotated, r2

        def plot_strike(x, y, strike, scale=100, colour='black'):
            rotated, r2 = strike_symbol(-strike)
            rotated *= scale
            r2 *= scale
            plt.plot([x, x + rotated[0]], [y, y + rotated[1]], colour)
            plt.plot([x - rotated[0], x], [y - rotated[1], y], colour)
            plt.plot([x, x + r2[0]], [y, y + r2[1]], colour)

        fig, ax = plt.subplots(1, figsize=(10, 10))

        #points = strati.support.interpolator.get_gradient_control()
        strike = normal_vector_to_strike_and_dip(points[:, 3:])
        for i in range(len(strike)):
            plot_strike(points[i, 0], points[i, 1], -strike[i, 0])
        # ax.scatter(points[:,0],points[:,1],c=points[:,3])
        for i in range(1, 3):
            points = fault_frame.features[0].support.interpolator.get_gradient_control()
            strike = normal_vector_to_strike_and_dip(points[:, 3:])
            for i in range(len(strike)):
                plot_strike(points[i, 0], points[i, 1], strike[i, 0] + 180, colour='red')
    def plot_feature(self,feature):
        v = feature.evaluate_value(np.array(self.xx,self.yy).T)
        plt.imshow(v.reshape(self.nsteps))

    def contour_feature(self, feature):
        pass