"""
2a. Modelling folds
====================
This tutorial shows how LoopStructural improves the modelling of folds by
using an accurate parameterisation of fold geometry, by:

1. modelling a folded surface without structural geology - i.e. using only
   data points and letting the interpolator's regularisation shape the
   surface between them, and
2. modelling the same surface using structural geology, which involves
   describing a local fold frame, calculating fold rotation angles, and
   constructing folded foliations using fold geostatistics within the
   fold frame coordinate system.
"""

######################################################################
# Imports
# -------

import pandas as pd

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_noddy_single_fold
from LoopStructural.visualisation import Loop3DView, RotationAnglePlotter

######################################################################
# Structural geology of folds
# ----------------------------
# Folds are one of the most common features found in deformed rocks and
# are defined by the location of higher curvature. The geometry of the
# folded surface can be characterised by three geometrical elements:
#
# 1. the fold hinge is the point of maximum curvature along folded surface
# 2. the axial surface is a surfaces that passes through all curvature
#    points in all folded foliations
# 3. the fold axis is the intersection of the folded foliation and the
#    axial surface
#
# Modelling folded surfaces using standard implicit algorithms is
# challenging because the implicit modelling methods are generally trying
# to minimise the resulting curvature of the surface. To model folded
# surfaces the geologist will need to characterise the geometry of the
# folded surface in high detail.


######################################################################
# Modelling folded surfaces without structural geology
# ----------------------------------------------------
#
# In the following section we will attempt to model a synthetic fold shape
# that is defined by a sinusoidal folded surface. For simplicity we will
# consider the fold as cylindrical and therefore only consider the fold in
# a 2D plane. The data set has been sampled from a model generated using
# Noddy and represents a fold with a wavelength of ~10km and amplitude of
# ~2km.
#
# The orientation of the structure has been sampled within the model
# volume (10km,7km,5km) at 500m intervals.
#
# **The aim of this exercise is to investigate how standard implicit
# modelling techniques are fundamentally limited when trying to model
# folded surfaces.**
#
# 1. Load data from sample datasets
# 2. Visualise data
# 3. Look at varying degrees of sampling e.g. 200 points, 100 points, 10
#    points.
# 4. Look at using data points ONLY from a map surface
#


######################################################################
# Modelling folded surfaces using loop structural
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# load the sample data
data, boundary_points = load_noddy_single_fold()
data.head()


######################################################################
# The input dataset was generated using Noddy by sampling the orientation
# of a structure on a regular grid. We have loaded it into a pandas
# DataFrame, this is basically an excel spreadsheet for python. Above are
# the first 5 rows of the dataset and as we can see it is regularly
# sampled with data points being sampled regularly along the :math:`z`,
# :math:`y` and :math:`x` axes. In order to avoid artefacts due to the
# sampling errors we will shuffle the data. We can do this using the
# ``random`` column in the DataFrame (ensuring everyone has the same
# data).
#

data = data.sort_values(
    "random"
)  # sort the data by a random int then we can select N random points
data.head()


######################################################################
# The data points are now randomly ordered and can now be subsampled by
# choosing the first N samples from the dataframe
#
# .. code:: python
#
#    data[:100]
#
# returns the first 100 data points from the array
#


######################################################################
# Testing data density
# ~~~~~~~~~~~~~~~~~~~~
#
# The number of points used to build the model is controlled by
# ``npoints`` below - try changing it and re-running to see how the shape
# of the interpolated fold degrades as fewer points are used, since
# without a fold frame the interpolator only has the regularisation term
# to constrain the surface between observations.
#
# **The black arrows are the normal vector to the folded surface**
#
npoints = 20
model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
model.data = data[:npoints]
stratigraphy = model.create_and_add_foliation(
    "s0", interpolatortype="PLI", nelements=5000, buffer=0.3, cgw=0.1
)
viewer = Loop3DView(model, background="white")
viewer.plot_data(stratigraphy)
viewer.plot_surface(stratigraphy, value=10)
viewer.show()


######################################################################
# Modelling folds using structural geology
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# The following section will describe how the fold axis, fold axial
# surface and fold vergence can be used to help constrain the shape of the
# folded surface. To do this we need to build a fold frame which is
# curvilinear coordinate system based around the fold axis and the fold
# axial surface.
#
# There are three coordinates to the fold frame:
#
# * coordinate 0 is the axial surface of the fold and is parallel to the
#   axial foliation
# * coordinate 1 is the fold axis direction field and is orthogonal to the
#   axial foliation
# * coordinate 2 is orthogonal to both the fold axis direction field and
#   axial foliation and is roughly parallel to the extension direction of
#   the fold
#
# Three direction vectors are defined by the normalised gradient of these
# fields: :math:`e_0` (red), :math:`e_1` (green), :math:`e_2` (blue).
#
# The orientation of the folded foliation can be defined by rotating
# :math:`e_1` around :math:`e_0` by the fold axis rotation angle
# :math:`\alpha_P` to give the orientation of the fold axis. The
# orientation of the folded foliation can then be defined by rotating the
# plane defined by the fold axis and :math:`e_0` around the fold axis by
# the fold limb rotation angle :math:`\alpha_L`.
#
# Calculating the fold rotation angles
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# The rotation angles can be calculated for observations of the folded
# foliation and associated lineations. For example, the fold axis rotation
# angle is found by calculating the angle between the gradient of the fold
# axis direction field and the intersection lineations. The fold limb
# rotation angle is found by finding the angle needed to rotate the
# folded foliation to be parallel to the plane of the axial foliation.
# The wavelength can be specified by the user or, in some cases, estimated
# from the s-variogram of the fold frame coordinate system.
#
mdata = pd.concat([data[:npoints], data[data["feature_name"] == "s1"]])
model = GeologicalModel(boundary_points[0, :], boundary_points[1, :])
model.data = mdata
fold_frame = model.create_and_add_fold_frame(
    "s1",
    interpolatortype="PLI",
    nelements=10000,
    buffer=0.5,
)
stratigraphy = model.create_and_add_folded_foliation(
    "s0",
    fold_frame=fold_frame,
    nelements=10000,
    fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
    limb_wl=12000,
    buffer=0.5,
)
viewer = Loop3DView(model, background="white")
viewer.plot_surface(
    fold_frame[0],
    value=10,
    colour="blue",
    opacity=0.5,
)
viewer.plot_data(stratigraphy)
viewer.plot_surface(stratigraphy, value=10)
viewer.show()

###########################################
# Plotting the fold rotation angles
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The fold limb rotation angle can be plotted against the fold frame
# coordinate to show the calculated data (points), the fitted rotation
# curve, and the S-variogram used to estimate the fold wavelength.
rotation_plots = RotationAnglePlotter(stratigraphy)
rotation_plots.add_fold_limb_data()
rotation_plots.add_fold_limb_curve()
rotation_plots.add_limb_svariogram()
rotation_plots.fig.show()
