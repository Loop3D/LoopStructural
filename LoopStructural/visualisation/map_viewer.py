import matplotlib.pyplot as plt
import numpy as np

class MapView:
    def __init__(self,origin,maximum,nsteps,**kwargs):
        self.origin = origin
        self.maximum = maximum
        x = np.linspace(self.origin[0],self.maximum[0],self.nsteps[0])
        y = np.linspace(self.origin[1],self.maximum[1],self.nsteps[1])
        self.xx, self.yy = np.meshgrid(x,y,indexing='ij')
        self.xx=self.xx.flatten()
        self.yy = self.yy.flatten()
    def add_data_to_plot(self,feature):
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
    def plot_feature(self,feature):
        v = feature.evaluate_value(np.array(self.xx,self.yy).T)
        plt.imshow(v.reshape(self.nsteps))

    def contour_feature(self,feature):
        pass