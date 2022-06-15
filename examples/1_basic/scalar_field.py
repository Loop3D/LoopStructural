"""
1f. Implicit modelling: signed distance field
=============================================

This example explores the relationship between the value of the scalar
field and the geometry that is modelled.

"""

from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


######################################################################
# Changing magnitude of the gradient norm
# ---------------------------------------
#
# The first example demonstrates how the value of the scalar field is
# related to the gradient norm, using a single value constraint and a
# gradient norm constraint where the magnitude of the norm is increased. A
# larger gradient norm means that the distance between the same isovalues
# is decreased, as the function changes value at a faster rate.
#

images = {}
for v in [1, 5, 1 / 5]:
    v_data = np.zeros((1, 4))
    v_data[0, :3] += 0.1
    data = pd.DataFrame(v_data, columns=["X", "Y", "Z", "val"])
    data["feature_name"] = "test"
    data["feature_name"] = "test"
    data["nx"] = np.nan
    data["ny"] = np.nan
    data["nz"] = np.nan
    data.loc[3, :] = [0, 0, 0, np.nan, "test", v, 0, 0]
    # data.loc[3,['nx','ny','nz']]/=np.linalg.norm(data.loc[3,['nx','ny','nz']])
    # data.loc[4,:] = [0,0,1,np.nan,'test',1,0,0]
    model = GeologicalModel(np.zeros(3), np.ones(3) * 10)
    model.data = data
    model.create_and_add_foliation("test", nelements=1e4, interpolatortype="FDI")
    view = LavaVuModelViewer(model)
    view.add_isosurface(model["test"], slices=[0, 1], name="test")
    view.add_data(model["test"])
    view.rotate([-92.68915557861328, 2.879497528076172, 1.5840799808502197])
    view.xmin = 0
    view.ymin = 0
    view.zmin = 0
    view.xmax = 10
    view.ymax = 10
    view.zmax = 10
    images[v] = view.image_array()
fig, ax = plt.subplots(1, 3, figsize=(30, 10))
ax[0].imshow(images[1])
ax[1].imshow(images[5])
ax[2].imshow(images[1 / 5])
ax[0].axis("off")
ax[1].axis("off")
ax[2].axis("off")
ax[0].set_title(f"A. Magnitude of gradient norm {1.0} ")
ax[1].set_title(f"B. Magnitude of gradient norm 5.0 ")
ax[2].set_title(f"C. Magnitude of gradient norm {1/5.} ")


######################################################################
# Changing magnitude of gradient norm with constant value constraints
# -------------------------------------------------------------------
#
# Changing the magnitude of the gradient norm changes the distances
# between the isosurfaces 0 and 1. What happens when we add value
# constraints? With a gradient norm magnitude of 1, adding a point with
# the value of 1 at the coordinates 1,1,1
#

images = {}
for v in [1, 5, 1 / 5]:
    v_data = np.zeros((1, 4))
    v_data[0, :3] += 0.1
    data = pd.DataFrame(v_data, columns=["X", "Y", "Z", "val"])
    data["feature_name"] = "test"
    data["feature_name"] = "test"
    data["nx"] = np.nan
    data["ny"] = np.nan
    data["nz"] = np.nan
    data.loc[1, :] = [1.1, 1.1, 1.1, 1.0, "test", np.nan, np.nan, np.nan]
    data.loc[3, :] = [5, 5, 5, np.nan, "test", v, 0, 0]

    # data.loc[3,['nx','ny','nz']]/=np.linalg.norm(data.loc[3,['nx','ny','nz']])
    # data.loc[4,:] = [0,0,1,np.nan,'test',1,0,0]
    model = GeologicalModel(np.zeros(3), np.ones(3) * 10)
    model.data = data
    model.create_and_add_foliation("test", nelements=1e4, interpolatortype="FDI")
    view = LavaVuModelViewer(model)
    view.add_isosurface(model["test"], slices=[0, 1], name="test")
    # view.add_isosurface(model['test'],nslices=5,name='test2',colour='blue')

    view.add_data(model["test"])
    view.rotate([-92.68915557861328, 2.879497528076172, 1.5840799808502197])
    view.xmin = 0
    view.ymin = 0
    view.zmin = 0
    view.xmax = 10
    view.ymax = 10
    view.zmax = 10
    images[v] = view.image_array()
fig, ax = plt.subplots(1, 3, figsize=(30, 10))
ax[0].imshow(images[1])
ax[1].imshow(images[5])
ax[2].imshow(images[1 / 5])
ax[0].axis("off")
ax[1].axis("off")
ax[2].axis("off")
ax[0].set_title(f"A. Magnitude of gradient norm {1.0} ")
ax[1].set_title(f"B. Magnitude of gradient norm 5.0 ")
ax[2].set_title(f"C. Magnitude of gradient norm {1/5.} ")


######################################################################
# -  The model in A. has a norm constraint where the magnitude norm is 1.0
#    and the distance between the value points 0 and 1 is 1.0m. In this
#    case the gradient norm is consistent with the data and the resulting
#    geometry is as we would expect.
# -  The model in B. has a norm constraint with the magnitude norm of 5.0
#    meaning the scalar field will grow quickly with the distance between
#    the isovalues 0 and 1 being shorter.
# -  The model in C. has a norm constraint with the magnitude of the
#    vector of 0.2, this means that the scalar field is going to grow
#    slower with the distance between 0 and 1 being larger.
#
# The resulting surfaces in B and C are oblique to the Z,Y plane as this
# makes the projected distance between the value points closer to the
# constraint values. In this example the interpolation is minimised by not
# fitting the direction of the gradient norm but minimising the magnitude
# of the gradient norm. This is an important take away for building a 3D
# model as if the gradient norm of the norm constraints are not consistent
# with the value constraints the model may have unexpected geometries.
#


######################################################################
# Using gradient constraints
# --------------------------
#
# An alternative approach if the norm of the scalar field is not known is
# to use gradient constraints. This type of constraint finds two vectors
# perpendicular to the norm (strike and dip vectors) and adds two
# constraints for the gradient of the implicit function to be
# perpendicular to these constraints. This type of constraint only
# constrains the orientation of the scalar field and not the direction or
# magnitude of the norm. In the following examples it does not matter what
# the magnitude of the norm is for these constraints the models are the
# same.
#

images = {}
for v in [1, 5, 1 / 5]:
    v_data = np.zeros((1, 4))
    v_data[0, :3] += 0.1
    data = pd.DataFrame(v_data, columns=["X", "Y", "Z", "val"])
    data["feature_name"] = "test"
    data["feature_name"] = "test"
    data["gx"] = np.nan
    data["gy"] = np.nan
    data["gz"] = np.nan
    data.loc[1, :] = [1.1, 1.1, 1.1, 1.0, "test", np.nan, np.nan, np.nan]
    data.loc[3, :] = [5, 5, 5, np.nan, "test", v, 0, 0]

    # data.loc[3,['nx','ny','nz']]/=np.linalg.norm(data.loc[3,['nx','ny','nz']])
    # data.loc[4,:] = [0,0,1,np.nan,'test',1,0,0]
    model = GeologicalModel(np.zeros(3), np.ones(3) * 10)
    model.data = data
    model.create_and_add_foliation("test", nelements=1e4, interpolatortype="FDI")
    view = LavaVuModelViewer(model)
    view.add_isosurface(model["test"], slices=[0, 1], name="test")
    # view.add_isosurface(model['test'],nslices=5,name='test2',colour='blue')

    view.add_data(model["test"])
    view.rotate([-92.68915557861328, 2.879497528076172, 1.5840799808502197])
    view.xmin = 0
    view.ymin = 0
    view.zmin = 0
    view.xmax = 10
    view.ymax = 10
    view.zmax = 10
    images[v] = view.image_array()
fig, ax = plt.subplots(1, 3, figsize=(30, 10))
ax[0].imshow(images[1])
ax[1].imshow(images[5])
ax[2].imshow(images[1 / 5])
ax[0].axis("off")
ax[1].axis("off")
ax[2].axis("off")
ax[0].set_title(f"A. Magnitude of vector {1.0} ")
ax[1].set_title(f"B. Magnitude of vector 5.0 ")
ax[2].set_title(f"C. Magnitude of vector {1/5.} ")


######################################################################
# Dome structures
# ---------------
#
# The geometry of dome structures are very sensitive to the difference in
# the value of the scalar field. In this example there are two rings of
# value constraints where the centre ring has a scalar field value of 0
# and the outer ring has a scalar field value ranging of 1, 5 or 9. There
# is also a gradient norm constraint where the with different scalar field
# values
#

images = {}
interpolator = "FDI"
# for interpolator in ['PLI']:
images[interpolator] = {}
for val in [1, 5, 9]:  # range(1,10,1):
    a = np.linspace(0, 360, 30)
    radius_mult = 10
    x = 0 + radius_mult * 3 * np.cos(np.deg2rad(a))
    y = 0 + radius_mult * 3 * np.sin(np.deg2rad(a))
    z = np.zeros_like(a)
    x2 = 0 + radius_mult * 5.5 * np.cos(np.deg2rad(a))
    y2 = 0 + radius_mult * 5.5 * np.sin(np.deg2rad(a))
    z2 = np.zeros_like(a)

    data = np.hstack([np.vstack([x, y, z, z]), np.vstack([x2, y2, z2, z2 + val])])

    dataframe = pd.DataFrame(data.T, columns=["X", "Y", "Z", "val"])
    dataframe["feature_name"] = "test"
    dataframe["nx"] = np.nan
    dataframe["ny"] = np.nan
    dataframe["nz"] = np.nan
    dataframe.loc[len(dataframe), :] = [0, 0, 0, np.nan, "test", 0, 0, 1 / val]
    model = GeologicalModel(
        -radius_mult * np.ones(3) * np.array([10, 10, 1]),
        radius_mult * np.ones(3) * np.array([10, 10, 1]),
        rescale=False,
    )
    model.data = dataframe
    model.create_and_add_foliation(
        "test", nelements=5e4, interpolatortype=interpolator
    )  # ,gpw=1,cpw=1,cgw=0.05)

    view.clear()
    view.model = model
    view.add_isosurface(model.features[0], nslices=5, name="test")
    view.add_data(model.features[0])
    view.rotation = [-64.91520690917969, -46.954345703125, -13.14844036102295]
    images[interpolator][val] = view.image_array()
fig, ax = plt.subplots(1, 3, figsize=(30, 10))
ax[0].imshow(images["FDI"][1])
ax[1].imshow(images["FDI"][5])
ax[2].imshow(images["FDI"][9])
ax[0].axis("off")
ax[1].axis("off")
ax[2].axis("off")
ax[0].set_title(f"A. Value range {1.0} ")
ax[1].set_title(f"B. Value range 5.0 ")
ax[2].set_title(f"C. Value range 9.0 ")


######################################################################
# Setting magnitude of norm to be 1
# ---------------------------------
#

images = {}
interpolator = "FDI"
# for interpolator in ['PLI']:
images[interpolator] = {}
for val in [1, 5, 9]:  # range(1,10,1):
    a = np.linspace(0, 360, 30)
    radius_mult = 10
    x = 0 + radius_mult * 3 * np.cos(np.deg2rad(a))
    y = 0 + radius_mult * 3 * np.sin(np.deg2rad(a))
    z = np.zeros_like(a)
    x2 = 0 + radius_mult * 5.5 * np.cos(np.deg2rad(a))
    y2 = 0 + radius_mult * 5.5 * np.sin(np.deg2rad(a))
    z2 = np.zeros_like(a)

    data = np.hstack([np.vstack([x, y, z, z]), np.vstack([x2, y2, z2, z2 + val])])

    dataframe = pd.DataFrame(data.T, columns=["X", "Y", "Z", "val"])
    dataframe["feature_name"] = "test"
    dataframe["nx"] = np.nan
    dataframe["ny"] = np.nan
    dataframe["nz"] = np.nan
    dataframe.loc[len(dataframe), :] = [0, 0, 0, np.nan, "test", 0, 0, 1]
    model = GeologicalModel(
        -radius_mult * np.ones(3) * np.array([10, 10, 1]),
        radius_mult * np.ones(3) * np.array([10, 10, 1]),
        rescale=False,
    )
    model.data = dataframe
    model.create_and_add_foliation(
        "test", nelements=5e4, interpolatortype=interpolator
    )  # ,gpw=1,cpw=1,cgw=0.05)

    view.clear()
    view.model = model
    view.add_isosurface(model.features[0], nslices=5, name="test")
    view.add_data(model.features[0])
    view.rotation = [-64.91520690917969, -46.954345703125, -13.14844036102295]
    images[interpolator][val] = view.image_array()
fig, ax = plt.subplots(1, 3, figsize=(30, 10))
ax[0].imshow(images["FDI"][1])
ax[1].imshow(images["FDI"][5])
ax[2].imshow(images["FDI"][9])
ax[0].axis("off")
ax[1].axis("off")
ax[2].axis("off")
ax[0].set_title(f"A. Value range {1.0} ")
ax[1].set_title(f"B. Value range 5.0 ")
ax[2].set_title(f"C. Value range 9.0 ")


######################################################################
# Using gradient constraints
# --------------------------
#
# If the norm direction is not known and gradient constraints are used the
# solution is the same for all values.
#

images = {}
interpolator = "FDI"
# for interpolator in ['PLI']:
images[interpolator] = {}
for val in [1, 5, 9]:  # range(1,10,1):
    a = np.linspace(0, 360, 30)
    radius_mult = 10
    x = 0 + radius_mult * 3 * np.cos(np.deg2rad(a))
    y = 0 + radius_mult * 3 * np.sin(np.deg2rad(a))
    z = np.zeros_like(a)
    x2 = 0 + radius_mult * 5.5 * np.cos(np.deg2rad(a))
    y2 = 0 + radius_mult * 5.5 * np.sin(np.deg2rad(a))
    z2 = np.zeros_like(a)

    data = np.hstack([np.vstack([x, y, z, z]), np.vstack([x2, y2, z2, z2 + val])])

    dataframe = pd.DataFrame(data.T, columns=["X", "Y", "Z", "val"])
    dataframe["feature_name"] = "test"
    dataframe["gx"] = np.nan
    dataframe["gy"] = np.nan
    dataframe["gz"] = np.nan
    dataframe.loc[len(dataframe), :] = [0, 0, 0, np.nan, "test", 0, 0, 1]
    model = GeologicalModel(
        -radius_mult * np.ones(3) * np.array([10, 10, 1]),
        radius_mult * np.ones(3) * np.array([10, 10, 1]),
        rescale=False,
    )
    model.data = dataframe
    model.create_and_add_foliation(
        "test", nelements=5e4, interpolatortype=interpolator
    )  # ,gpw=1,cpw=1,cgw=0.05)

    view.clear()
    view.model = model
    view.add_isosurface(model.features[0], nslices=5, name="test")
    view.add_data(model.features[0])
    view.rotation = [-64.91520690917969, -46.954345703125, -13.14844036102295]
    images[interpolator][val] = view.image_array()
fig, ax = plt.subplots(1, 3, figsize=(30, 10))
ax[0].imshow(images["FDI"][1])
ax[1].imshow(images["FDI"][5])
ax[2].imshow(images["FDI"][9])
ax[0].axis("off")
ax[1].axis("off")
ax[2].axis("off")
ax[0].set_title(f"A. Value range {1.0} ")
ax[1].set_title(f"B. Value range 5.0 ")
ax[2].set_title(f"C. Value range 9.0 ")
