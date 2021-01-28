import logging

import matplotlib.pyplot as plt
import numpy as np

from LoopStructural.utils import getLogger
logger = getLogger(__name__)


class RotationAnglePlotter:
    def __init__(self, feature, axis=True):
        """

        """
        self.fig, self.ax = plt.subplots(2, 2, figsize=(30, 15))
        self.ax[0][0].set_ylim(-90, 90)
        self.ax[1][0].set_ylim(-90, 90)
        self.feature = feature

    def plot(self, x, y, ix, iy, symb,**kwargs):
        """

        Parameters
        ----------
        x : np.array
            vector of x
        y
        ix
        iy
        symb

        Returns
        -------

        """
        return self.ax[iy][ix].plot(x, y, symb,**kwargs)
    def default_titles(self):
        self.ax[0][0].set_title('Fold Axis S-Plot')
        self.ax[0][1].set_title('Fold Axis S-Variogram')
        self.ax[1][0].set_title('Fold Limb S-Plot')
        self.ax[1][1].set_title('Fold Limb S-Variogram')

        self.ax[1][1].set_xlabel('Variogram Steps')
        self.ax[1][1].set_ylabel('Fold Limb S-Variogram')
        self.ax[1][0].set_ylabel('Fold Limb Rotation Angle')
        self.ax[1][0].set_xlabel('Fold Frame Axial Surface Field')

        self.ax[0][1].set_xlabel('Variogram Steps')
        self.ax[0][1].set_ylabel('Fold Axis S-Variogram')
        self.ax[0][0].set_ylabel('Fold Axis Rotation Angle')
        self.ax[0][0].set_xlabel('Fold Frame Axis Direction Field')
    def add_fold_limb_data(self, symb="bo",**kwargs):
        fold_frame = self.feature.fold.fold_limb_rotation.fold_frame_coordinate
        rotation = self.feature.fold.fold_limb_rotation.rotation_angle
        return self.plot(fold_frame, rotation, 0, 1, symb,**kwargs)
        

    def add_fold_limb_curve(self, symb='r-',**kwargs):
        x = np.linspace(self.feature.fold.foldframe[0].min(),self.feature.fold.foldframe[0].max(),100)
        return self.plot(x,self.feature.fold.fold_limb_rotation(x), 0, 1, symb,**kwargs)

    def add_axis_svariogram(self, symb='bo',**kwargs):
        svariogram = self.feature.fold.fold_axis_rotation.svario
        if svariogram:
            svariogram.calc_semivariogram()
            return self.plot(svariogram.lags, svariogram.variogram, 1, 0, symb,**kwargs)

    def add_limb_svariogram(self, symb='bo', **kwargs):
        svariogram = self.feature.fold.fold_limb_rotation.svario
        if svariogram:
            svariogram.calc_semivariogram()
            return self.plot(svariogram.lags, svariogram.variogram, 1, 1, symb,**kwargs)

    def add_fold_axis_data(self, symb='bo',**kwargs):
        fold_frame = self.feature.fold.fold_axis_rotation.fold_frame_coordinate
        rotation = self.feature.fold.fold_axis_rotation.rotation_angle
        return self.plot(fold_frame, rotation, 0, 0, symb,**kwargs)

    def add_fold_axis_curve(self, symb='r-',**kwargs):
        x = np.linspace(self.feature.fold.foldframe[1].min(),self.feature.fold.foldframe[1].max(),100)
        return self.plot(x,self.feature.fold.fold_axis_rotation(x), 0, 0, symb,**kwargs)


