# ### Imports

from LoopStructural.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from LoopStructural.interpolators.discrete_fold_interpolator import DiscreteFoldInterpolator as DFI
from LoopStructural.modelling.features.geological_feature import GeologicalFeatureInterpolator
from LoopStructural.modelling.features.faulted_geological_feature import FaultedGeologicalFeature
from LoopStructural.modelling.structural_frame import StructuralFrameBuilder, StructuralFrame
from LoopStructural.modelling.fold.foldframe import FoldFrame
from LoopStructural.modelling.fold.fold import FoldEvent
from LoopStructural.modelling.fold.svariogram import SVariogram
from LoopStructural.modelling.fault.fault_segment import FaultSegment
from LoopStructural.supports.tet_mesh import TetMesh
from LoopStructural.visualisation.model_visualisation import LavaVuModelViewer
from LoopStructural.visualisation.rotation_angle_plotter import RotationAnglePlotter
# import other libraries
import geopandas
import numpy as np
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt


# ### Load data from shapefile
# We use geopandas to load the objects from a shapefile. In this case we use a geopackage with three
# shapefile layers:
# * points - value data
# * orientations - orientation data
# * bounding box - simple defines the map view of the model area
# The bounding box
viewer = LavaVuModelViewer(background="white")

points = geopandas.read_file(
    '/home/lgrose/lachlan.grose@monash.edu/Loop/case_study_data/faulted_fold/data2.gpkg',layer='points')
orientations = geopandas.read_file(
    '/home/lgrose/lachlan.grose@monash.edu/Loop/case_study_data/faulted_fold/data2.gpkg',layer='orientations')
print(orientations)
boundary_points = np.zeros((2,3))
boundary_points[0,0] = 0#-10 #np.min(coords[:,0])-10
boundary_points[0,1] = 0#-10 #np.min(coords[:,1])
boundary_points[0,2] = -5000#0#-40#minz
boundary_points[1,0] = 7000 #np.max(coords[:,0])
boundary_points[1,1] = 10000 #np.max(coords[:,1])
boundary_points[1,2] = 500#-50000#-minz*0.1
mesh = TetMesh()
mesh.setup_mesh(boundary_points, n_tetra=20000,)
fault_frame_interpolator = PLI(mesh)
fault_frame_builder = StructuralFrameBuilder(
    interpolator=fault_frame_interpolator,
    mesh=mesh,
    name='fault_frame')
# Interfacing with dataframes should be done using a convenience wrapper function
for i, r in orientations.iterrows():
    if r['label'] == 'fault':
        xy = r['geometry'].xy

        z = 0#5000
        if 'z' in r:
            z = r['z']
        fault_frame_builder.add_strike_and_dip([xy[0][0],xy[1][0],z],
                                                     r['strike']+180,
                                                     -r['dip'],
                                                     itype='gx')
    if r['label'] == 'fault_slip':
        xy = r['geometry'].xy
        z = 0#5000
        if 'z' in r:
            z = r['z']
        print(r)
        fault_frame_builder.add_strike_and_dip([xy[0][0],xy[1][0],z],r['strike']-90,90.,itype='gy')

for i, r in points.iterrows():
    if r['label'] == 'fault':
        xy = r['geometry'].xy
        z = 0#5000
        if 'z' in r:
            z = r['z']
        print(r)
        fault_frame_builder.add_point([xy[0][0],xy[1][0],z],r['value'],itype='gx')
# #         fault_frame_builder.add_point([xy[0][0],xy[1][0],z],r['value'],itype='gy')
#
# # # We define weights for the orthogonal constraint and the regularisation constraints. The solver to
# # # use to solve the least squares system. Possible solvers include
# # # * **chol** cholesky decompsition
# # # * **lu** - lower upper decomposition
# # # * **cg** - conjugate gradient
# # # * **bicg** - biconjugate gradient
# #
ogw = 3000
ogw /= mesh.n_elements
cgw = 6000
solver='lu'
fault_frame = fault_frame_builder.build(
    # frame=FoldFrame,
    solver=solver,
    gxxgy=2 * ogw,
    gxxgz=2 * ogw,
    gyxgz=ogw,
    gxcg=cgw,
    gycg=cgw,
    gzcg=cgw,
    shape='rectangular',)
#
fold_frame_interpolator  = PLI(mesh)
fold_frame_builder = StructuralFrameBuilder(
    interpolator=fold_frame_interpolator,
    mesh=mesh,
    name='fold_frame'
)
for i, r in orientations.iterrows():
    if r['label'] == 'fold':
        xy = r['geometry'].xy
        z = 0
        if 'z' in r:
            z = r['z']
        fold_frame_builder.add_strike_and_dip([xy[0][0],xy[1][0],z],r['strike'],r['dip'],itype='gx')
for i, r in orientations.iterrows():
    if r['label'] == 'fold_axis':
        xy = r['geometry'].xy
        z = 0
        if 'z' in r:
            z = r['z']
        fold_frame_builder.add_strike_and_dip([xy[0][0],xy[1][0],z],r['strike'],r['dip'],itype='gy')
for i, r in points.iterrows():
    if r['label'] == 'fold':
        xy = r['geometry'].xy
        z = 0
        if 'z' in r:
            z = r['z']
        fault_frame_builder.add_point([xy[0][0],xy[1][0],z],r['value'],itype='gx')
#
#
ogw = 3000
ogw /= mesh.n_elements
cgw = 6000
solver='lu'
fold_frame = fold_frame_builder.build(
    frame=FoldFrame,
    solver=solver,
    gxxgy=2 * ogw,
    gxxgz=2 * ogw,
    gyxgz=ogw,
    gxcg=cgw,
    gycg=cgw,
    gzcg=cgw,
    shape='rectangular',)
#
fault = FaultSegment(fault_frame,
                     displacement=1000)
fault.apply_to_data(fold_frame.features[0].data)
fold_frame_features = []
for f in fold_frame.features:
    fold_frame_features.append(FaultedGeologicalFeature(f, fault))
faulted_fold_frame = FoldFrame("faulted_frame", fold_frame_features)
for i, r in points.iterrows():
    if r['label'] == 's1':
        xy = r['geometry'].xy
        z = 0
        if 'z' in r:
            z = r['z']
#
# # ### Create a fold event linked to the fold frame
# # We need to create an empty fold event that links our fold frame to the fold event so that it can
# #  given to the fold interpolator.
fold = FoldEvent(faulted_fold_frame,None,None)
# #
# # # ### Create a DiscreteFoldInterpolator object
# # # The DiscreteFoldInterpolator is a daughter class of the PiecewiseLinearInterpolator that
# # # uses a fold event and a mesh to define additional constraints in the least squares system.
stratigraphy_interpolator = DFI(mesh, fold)
# #
# # # ### Build the stratigraphy geological feature
# # # We can build the stratigraphy geological feature using the fold interpolator object and
# # # then linking the observations from the shapefile to the interpolator. .
# #
stratigraphy_builder = GeologicalFeatureInterpolator(stratigraphy_interpolator, name="folded_stratigraphy")
for i, r in orientations.iterrows():
    if r['label'] == 's0':
        xy = r['geometry'].xy
        z = 0
        if 'z' in r:
            z = r['z']
        stratigraphy_builder.add_strike_and_dip([xy[0][0],xy[1][0],z],r['strike'],r['dip'])
strati = stratigraphy_builder.build()
# viewer.plot_vector_data(strati.support.interpolator.get_gradient_control()[:,:3],
#                         strati.support.interpolator.get_gradient_control()[:,3:],
#                         "strati_grad2",
#                         colour='green')
strati_g = strati.support.interpolator.get_gradient_control()
faulted_strati = FaultedGeologicalFeature(strati, fault)
xyz = stratigraphy_builder.interpolator.get_gradient_control()[:,:3]
s0g = stratigraphy_builder.interpolator.get_gradient_control()[:,3:]
l1 = faulted_fold_frame.calculate_intersection_lineation(np.hstack([xyz,s0g]))
far = faulted_fold_frame.calculate_fold_axis_rotation(np.hstack([xyz,l1]))
s1 = faulted_fold_frame.features[0].evaluate_value(xyz)
s1gy = faulted_fold_frame.features[1].evaluate_value(xyz)
axis_svariogram = SVariogram(s1gy,far)
guess = axis_svariogram.find_wavelengths()
# guess/=2
rotation_plots = RotationAnglePlotter()
rotation_plots.add_fold_axis_data(far,s1gy)
rotation_plots.add_axis_svariogram(axis_svariogram)
def fold_axis_rotation(x):
    v = np.zeros(x.shape)
    v[:] = 0
    return v
fold.fold_axis_rotation = fold_axis_rotation
#
# # ### Evaluate the fold limb rotation angle
# # The fold axis can be queried for any location in the model. The axis is a unit vector and
# # does not need to be normalised here. THe fold limb rotation angle is calculated by finding
# # the angle between the folded foliation and the axial foliation in the fold frame. The calculated
# # axis can be used by passing an N,3 array of the fold axis for every N locations. Or if no
# # axis is specified the local intersection lineation is used.
axis = fold.get_fold_axis_orientation(xyz)
flr = faulted_fold_frame.calculate_fold_limb_rotation(np.hstack([xyz,s0g]),axis=axis)
limb_svariogram = SVariogram(s1,flr)
guess = np.array(limb_svariogram.find_wavelengths())
# guess[0] = 5000.
guess/=2.
flr_tan = np.tan(np.deg2rad(flr))
rbf_fold_limb = Rbf(s1,np.zeros(s1.shape),np.zeros(s1.shape),flr_tan,
                    function='gaussian',
                    epsilon=guess[1],
                    smooth=0.05)
xi = np.linspace(faulted_fold_frame.features[0].min(),faulted_fold_frame.features[0].max(),1000)
def fold_limb_rotation(x):
    return np.rad2deg(np.arctan(rbf_fold_limb(x,np.zeros(x.shape),np.zeros(x.shape))))

# The fold rotation angle function is added to the **FoldEvent**
# # fold.fold_limb_rotation = fold_limb_rotation
rotation_plots.add_fold_limb_data(flr,s1)
rotation_plots.add_limb_svariogram(limb_svariogram)
# # rotation_plots.add_fold_limb_curve(np.rad2deg(np.arctan(rbf_fold_limb(xi,np.zeros(1000),np.zeros(1000)))), xi)
# #
# # plt.show()
# # rotation_plots.add_fold_axis_curve(np.rad2deg(np.arctan(rbf_fold_axis(xi,np.zeros(1000),np.zeros(1000)))), xi)
rotation_plots.add_fold_limb_data(flr,s1)
rotation_plots.add_limb_svariogram(limb_svariogram)
# # rotation_plots.add_fold_limb_curve(np.rad2deg(np.arctan(rbf_fold_limb(xi,np.zeros(1000),np.zeros(1000)))), xi)
plt.show()
# # # ### Set up the fold geometry
# # # **get_gradient_control()** returns an array N,6 array with the xyz and normal vector
# # # components of the gradient control points. We can use these arrays to calculate the
# # # fold axis using the intersection lineation between the axial foliation and the normal
# # # to the folded foliation. The fold axis rotation angle can then be calculated by finding
# # # the angle between the fold axis direction field and the intersection lineation. A
# # # SVariogram can then be used to automatically pick the wavelength of the fold
# #
# # xyz = stratigraphy_builder.interpolator.get_gradient_control()[:,:3]
# # s0g = stratigraphy_builder.interpolator.get_gradient_control()[:,3:]
# # l1 = f1_frame.calculate_intersection_lineation(np.hstack([xyz,s0g]))
# # far = f1_frame.calculate_fold_axis_rotation(np.hstack([xyz,l1]))
# # s1 = f1_frame.features[0].evaluate_value(xyz)
# # s1gy = f1_frame.features[1].evaluate_value(xyz)
# # axis_svariogram = SVariogram(s1gy,far)
# # guess = axis_svariogram.find_wavelengths()
# # guess/=2
# #
# #
# # # Rather than interpolate the angle directly we interpolate the gradient
# # # of the curve. Which can be calculated by finding the tangent of the
# # # angle. This means that the interpolated value will not be > 90 or < -90.
# # # We use SciPy's RBF to interpolate the fold rotation angles. A smoothing parameter
# # # prevents overfitting of the interpolated angles.
# #
# # far_tan = np.tan(np.deg2rad(far))
# # rbf_fold_axis = Rbf(s1gy,np.zeros(s1gy.shape),np.zeros(s1gy.shape),far_tan,
# #                     function='gaussian',
# #                     epsilon=guess[0],
# #                     smooth=0.05)
# # xi = np.linspace(f1_frame.features[1].min(),
# #                  f1_frame.features[1].max(), 1000)
# #
# # # We create a wrapper function for the fold axis rotation angle that can be given to the **FoldEvent**
# # # object so that it can calculate the fold rotation angle for any value of the fold frame.
# #
# # def fold_axis_rotation(x):
# #     return np.rad2deg(np.arctan(rbf_fold_axis(x,np.zeros(x.shape),np.zeros(x.shape))))
# # fold.fold_axis_rotation = fold_axis_rotation
# #
# # # ### Evaluate the fold limb rotation angle
# # # The fold axis can be queried for any location in the model. The axis is a unit vector and
# # # does not need to be normalised here. THe fold limb rotation angle is calculated by finding
# # # the angle between the folded foliation and the axial foliation in the fold frame. The calculated
# # # axis can be used by passing an N,3 array of the fold axis for every N locations. Or if no
# # # axis is specified the local intersection lineation is used.
# # axis = fold.get_fold_axis_orientation(xyz)
# # flr = f1_frame.calculate_fold_limb_rotation(np.hstack([xyz,s0g]),axis=axis)
# # limb_svariogram = SVariogram(s1,flr)
# # guess = limb_svariogram.find_wavelengths()
# # guess[0] = 5000.
# # guess/=2.
# # flr_tan = np.tan(np.deg2rad(flr))
# # rbf_fold_limb = Rbf(s1,np.zeros(s1.shape),np.zeros(s1.shape),flr_tan,
# #                     function='gaussian',
# #                     epsilon=guess[0],
# #                     smooth=0.05)
# # xi = np.linspace(f1_frame.features[0].min(),f1_frame.features[0].max(),1000)
# # def fold_limb_rotation(x):
# #     return np.rad2deg(np.arctan(rbf_fold_limb(x,np.zeros(x.shape),np.zeros(x.shape))))
# #
# # # The fold rotation angle function is added to the **FoldEvent**
# # fold.fold_limb_rotation = fold_limb_rotation
# #
# # # ### Plot rotation angles
# # rotation_plots = RotationAnglePlotter()
# # rotation_plots.add_fold_axis_data(far,s1gy)
# # rotation_plots.add_axis_svariogram(axis_svariogram)
# # rotation_plots.add_fold_axis_curve(np.rad2deg(np.arctan(rbf_fold_axis(xi,np.zeros(1000),np.zeros(1000)))), xi)
# # rotation_plots.add_fold_limb_data(flr,s1)
# # rotation_plots.add_limb_svariogram(limb_svariogram)
# # rotation_plots.add_fold_limb_curve(np.rad2deg(np.arctan(rbf_fold_limb(xi,np.zeros(1000),np.zeros(1000)))), xi)
# #
# # # ### Define interpolation weights
# # # The DiscreteFoldInterpolator has the ability to weigh the different geometrical properties of the
# # # fold differently in the least squares system. Much the same as data points and regularisation are
# # # handled differently for standard interpolation.
# # # * **fold_orientation** is the weight of the dot product between a
# # # * **fold_axis** is the weight of the dot product between the fold axis and the interpolated foliation
# # # * **fold_normalisation**
# # # * **fold_regularisation** is the fold regularisation term that ensures the fold is smooth along strike
# # # These weights need to be added to a dictionary and then passed to the feature builder along with the fold
# # fold_weights = {}
# # fold_weights['fold_orientation'] = 50.
# # fold_weights['fold_axis'] = 3.
# # fold_weights['fold_normalisation'] = 1.
# # fold_weights['fold_regularisation'] = 10.
# # folded_stratigraphy = stratigraphy_builder.build(solver=solver,
# #                                                  fold_weights=fold_weights,
# #                                                  fold=fold)
# #
# # # ### Visualising Results
# # # We can visualise the surfaces for isovalues of the fold frame or interpolated surface.
# # # Normal vectors to the scalar fields can be visualised for evaluation points (if too many
# # # points are chosen it can be difficult to see the locations. We suggest slicing numpy
# # # arrays to create regularly sampled locations.
# #
viewer = LavaVuModelViewer(background="white")
viewer.plot_isosurface(fault_frame.features[0],  colour='green', isovalue=0)
viewer.plot_isosurface(fault_frame.features[1],  colour='blue')
viewer.plot_vector_data(fault_frame.features[0].support.interpolator.get_gradient_control()[:,:3],
                        fault_frame.features[0].support.interpolator.get_gradient_control()[:,3:],
                        "gx_grad",colour='green')

viewer.plot_vector_data(fault_frame.features[1].support.interpolator.get_gradient_control()[:,:3],
                        fault_frame.features[1].support.interpolator.get_gradient_control()[:,3:],
                        "gy_grad",
                        colour='blue')

viewer.plot_vector_data(strati.support.interpolator.get_gradient_control()[:,:3],
                        strati.support.interpolator.get_gradient_control()[:,3:],
                        "strati_grad",
                        colour='red')
# viewer.plot_vector_data(strati_g[:,:3],
#                         strati_g[:,3:],
#                         "strati_grad1",
#                         colour='black')
viewer.plot_isosurface(fold_frame.features[0],  colour='black',isovalue=0)
viewer.plot_isosurface(faulted_fold_frame.features[0],  colour='purple',isovalue=0)
viewer.plot_value_data(fault_frame.features[0].support.interpolator.get_control_points()[:,:3],
                       fault_frame.features[0].support.interpolator.get_control_points()[:,3],
                       'fault_data')
# viewer.plot_value_data(strati.support.interpolator.get_control_points()[:,:3],
#                        strati.support.interpolator.get_control_points()[:,3],
#                        'strati_data')
# # viewer.plot_isosurface(fold_frame.features[1],  colour='red')
# # viewer.plot_vector_data(fold_frame.features[0].support.interpolator.get_gradient_control()[:,:3],
# #                         fold_frame.features[0].support.interpolator.get_gradient_control()[:,3:],
# #                         "gx_grad")
# # viewer.plot_vector_data(fold_frame.features[1].support.interpolator.get_gradient_control()[:,:3],
# #                         fold_frame.features[1].support.interpolator.get_gradient_control()[:,3:],
# #                         "gy_grad",
# #                         colour='red')
# # viewer.plot_vector_data(stratigraphy_builder.interpolator.get_gradient_control()[:,:3],
# #                         stratigraphy_builder.interpolator.get_gradient_control()[:,3:],
# #                         "s0",
# #                         colour='pink',
# #                         size=4)
# # viewer.plot_isosurface(folded_stratigraphy,
# #                        colour='purple',
# #                        nslices=10,
# #                        # paint_with=f1_frame.features[0]
# #                        )
locations = mesh.barycentre[::20,:]
# # viewer.plot_vector_field(f1_frame.features[2], locations=locations, colour='red')
# viewer.plot_vector_field(fault_frame.features[1], locations=locations, colour='blue')
# # viewer.plot_vector_field(f1_frame.features[0], locations=locations, colour='blue')
viewer.interactive()
# #
