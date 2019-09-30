import matplotlib.pyplot as plt


class RotationAnglePlotter:
    def __init__(self):
        self.fig, self.ax = plt.subplots(2, 2, figsize=(15, 15))
        self.ax[0][0].set_ylim(-90, 90)
        self.ax[1][0].set_ylim(-90, 90)

    def plot(self,x,y,ix,iy,symb):
        self.ax[iy][ix].plot(x,y,symb)

    def add_fold_limb_data(self,rotation,fold_frame):
        self.plot(fold_frame,rotation,0,1,"bo")

    def add_fold_limb_curve(self,rotation,fold_frame):
        self.plot(fold_frame,rotation,0,1, 'r-')

    def add_axis_svariogram(self,svariogram):
        svariogram.calc_semivariogram()
        self.plot(svariogram.lags,svariogram.variogram, 1, 0,'bo')

    def add_limb_svariogram(self,svariogram):
        svariogram.calc_semivariogram()

        self.plot(svariogram.lags,svariogram.variogram, 1, 1,'bo')

    def add_fold_axis_data(self,rotation,fold_frame):
        self.plot(fold_frame,rotation,0,0,"bo")
    def add_fold_axis_curve(self,rotation,fold_frame):
        self.plot(fold_frame,rotation,0,0,"r-")

# ax[0][0].set_ylim(-90,90)
# svario = SVariogram(s1gy,far)
# # guess = [20000]
# guess = svario.find_wavelengths()
# guess/=2.
# ax[0][1].plot(svario.lags, svario.variogram, 'bo')
# ax[0][1].axvline(guess[0])
# far_tan = np.tan(np.deg2rad(far))
# rbf_fold_axis = Rbf(s1gy,np.zeros(s1gy.shape),np.zeros(s1gy.shape),far_tan,
#                     function='gaussian',
#                     epsilon=guess[0],
#                     smooth=0.05)
# xi = np.linspace(f1_frame.features[1].min(),
#                  f1_frame.features[1].max(),1000)
# ax[0][0].plot(xi, np.rad2deg(np.arctan(rbf_fold_axis(xi,np.zeros(1000),np.zeros(1000)))))
# ax[0][0].plot(s1gy,far,'bo')
#
# def fold_axis_rotation(x):
#     return np.rad2deg(np.arctan(rbf_fold_axis(x,np.zeros(x.shape),np.zeros(x.shape))))
# fold.fold_axis_rotation = fold_axis_rotation
#
# axis = fold.get_fold_axis_orientation(xyz)
# axis/=np.linalg.norm(axis,axis=1)[:,None]
# flr = f1_frame.calculate_fold_limb_rotation(np.hstack([xyz,s0g]),axis=axis)
# ##quick figure
# ax[1][0].plot(s1, flr, 'bo')
# ax[1][0].set_ylim(-90, 90)
# # plt.show()
# ax[1][0].set_ylim(-90, 90)
# svario = SVariogram(s1,flr)
# guess = svario.find_wavelengths()
# guess[0] = 5000.
# guess/=2.
# ax[1][1].plot(svario.lags,svario.variogram,'bo')
# ax[1][1].axvline(guess[0])
#
#
# flr_tan = np.tan(np.deg2rad(flr))
# rbf_fold_limb = Rbf(s1,np.zeros(s1.shape),np.zeros(s1.shape),flr_tan,
#                     function='gaussian',
#                     epsilon=guess[0],
#                     smooth=0.05)
# xi = np.linspace(f1_frame.features[0].min(),f1_frame.features[0].max(),1000)
# ax[1][0].plot(xi,np.rad2deg(np.arctan(rbf_fold_limb(xi,np.zeros(1000),np.zeros(1000)))))
# # plt.plot(s1,flr,'bo')
# # plt.ylim(-90,90)