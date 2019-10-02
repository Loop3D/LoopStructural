# ### Imports

from LoopStructural.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from LoopStructural.interpolators.discrete_fold_interpolator import DiscreteFoldInterpolator as DFI
from LoopStructural.modelling.features.geological_feature import GeologicalFeatureInterpolator
from LoopStructural.modelling.features.structural_frame import StructuralFrameBuilder
from LoopStructural.modelling.fold.foldframe import FoldFrame
from LoopStructural.modelling.fold.fold import FoldEvent
from LoopStructural.modelling.fold.svariogram import SVariogram
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

points = geopandas.read_file('data.gpkg',layer='points')
orientations = geopandas.read_file('data.gpkg',layer='orientations')
model_area = geopandas.read_file('data.gpkg',layer='bounding_box')

# If we have a look at the orientation and points data we can see that the points have a geometry
# column as well as type and itype.


geom = model_area['geometry']#.shapes()
coords = np.array(geom[0].exterior.coords)#[0]
minz = -(np.max(coords[:,0])-np.min(coords[:,0]))/2.

# ### Build mesh
boundary_points = np.zeros((2,3))
boundary_points[0,0] = np.min(coords[:,0])-10
boundary_points[0,1] = np.min(coords[:,1])
boundary_points[0,2] = minz
boundary_points[1,0] = np.max(coords[:,0])
boundary_points[1,1] = np.max(coords[:,1])
boundary_points[1,2] = -minz*0.1
mesh = TetMesh()
mesh.setup_mesh(boundary_points, n_tetra=20000,)


# ### Build a fold frame
# Fold frame is built from the observations of foliations and fold axis and consists of 2 interpolated
# scalar fields that are orthoognal to each other and one computed vector field that is orthogonal
# to both interpolated fields. The StructuralFrameBuilder needs to be given an interpolator object to use as
# a template to interpolate the curvilinear coordinate systems. The StructuralFrameBuilder creates copies
# of the interpolator for the different fields that need to be interpolated and assigns property names to
# the different fields being interpolated.
fold_frame_interpolator = PLI(mesh)
fold_frame_builder = StructuralFrameBuilder(
    interpolator=fold_frame_interpolator,
    mesh=mesh,
    name='F1_fold_frame')
# Interfacing with dataframes should be done using a convenience wrapper function
for i, r in orientations.iterrows():
    if r['type'] == 's1':
        xy = r['geometry'].xy
        z = 0
        if 'z' in r:
            z = r['z']
        fold_frame_builder.add_strike_and_dip([xy[0][0],xy[1][0],z],r['strike'],r['dip'],itype=r['itype'])
for i, r in points.iterrows():
    if r['type'] == 's1':
        xy = r['geometry'].xy
        z = 0
        if 'z' in r:
            z = r['z']
        fold_frame_builder.add_point([xy[0][0],xy[1][0],z],r['value'],itype=r['itype'])

# We define weights for the orthogonal constraint and the regularisation constraints. The solver to
# use to solve the least squares system. Possible solvers include
# * **chol** cholesky decompsition
# * **lu** - lower upper decomposition
# * **cg** - conjugate gradient
# * **bicg** - biconjugate gradient

ogw = 3000
ogw /= mesh.n_elements
cgw = 6000
solver='lu'
f1_frame = fold_frame_builder.build(
    frame=FoldFrame,
    solver=solver,
    gxxgy=2 * ogw,
    gxxgz=2 * ogw,
    gyxgz=ogw,
    gxcg=cgw,
    gycg=cgw,
    gzcg=cgw,
    shape='rectangular',)

# ### Create a fold event linked to the fold frame
# We need to create an empty fold event that links our fold frame to the fold event so that it can
#  given to the fold interpolator.
fold = FoldEvent(f1_frame,None,None)

# ### Create a DiscreteFoldInterpolator object
# The DiscreteFoldInterpolator is a daughter class of the PiecewiseLinearInterpolator that
# uses a fold event and a mesh to define additional constraints in the least squares system.
fold_interpolator = DFI(mesh,fold)

# ### Build the stratigraphy geological feature
# We can build the stratigraphy geological feature using the fold interpolator object and
# then linking the observations from the shapefile to the interpolator. .

stratigraphy_builder = GeologicalFeatureInterpolator(fold_interpolator, name="folded_stratigraphy")
for i, r in orientations.iterrows():
    if r['type'] == 's0':
        xy = r['geometry'].xy
        z = 0
        if 'z' in r:
            z = r['z']
        stratigraphy_builder.add_strike_and_dip([xy[0][0],xy[1][0],z],r['strike'],r['dip'])
for i, r in points.iterrows():
    if r['type'] == 's0':
        xy = r['geometry'].xy
        z = 0
        if 'z' in r:
            z = r['z']
        stratigraphy_builder.add_point([xy[0][0],xy[1][0],z],r['value'])

# ### Set up the fold geometry
# **get_gradient_control()** returns an array N,6 array with the xyz and normal vector
# components of the gradient control points. We can use these arrays to calculate the
# fold axis using the intersection lineation between the axial foliation and the normal
# to the folded foliation. The fold axis rotation angle can then be calculated by finding
# the angle between the fold axis direction field and the intersection lineation. A
# SVariogram can then be used to automatically pick the wavelength of the fold

xyz = stratigraphy_builder.interpolator.get_gradient_control()[:,:3]
s0g = stratigraphy_builder.interpolator.get_gradient_control()[:,3:]
l1 = f1_frame.calculate_intersection_lineation(np.hstack([xyz,s0g]))
far = f1_frame.calculate_fold_axis_rotation(np.hstack([xyz,l1]))
s1 = f1_frame.features[0].evaluate_value(xyz)
s1gy = f1_frame.features[1].evaluate_value(xyz)
axis_svariogram = SVariogram(s1gy,far)
guess = axis_svariogram.find_wavelengths()
guess/=2


# Rather than interpolate the angle directly we interpolate the gradient
# of the curve. Which can be calculated by finding the tangent of the
# angle. This means that the interpolated value will not be > 90 or < -90.
# We use SciPy's RBF to interpolate the fold rotation angles. A smoothing parameter
# prevents overfitting of the interpolated angles.

far_tan = np.tan(np.deg2rad(far))
rbf_fold_axis = Rbf(s1gy,np.zeros(s1gy.shape),np.zeros(s1gy.shape),far_tan,
                    function='gaussian',
                    epsilon=guess[0],
                    smooth=0.05)
xi = np.linspace(f1_frame.features[1].min(),
                 f1_frame.features[1].max(), 1000)

# We create a wrapper function for the fold axis rotation angle that can be given to the **FoldEvent**
# object so that it can calculate the fold rotation angle for any value of the fold frame.

def fold_axis_rotation(x):
    return np.rad2deg(np.arctan(rbf_fold_axis(x,np.zeros(x.shape),np.zeros(x.shape))))
fold.fold_axis_rotation = fold_axis_rotation

# ### Evaluate the fold limb rotation angle
# The fold axis can be queried for any location in the model. The axis is a unit vector and
# does not need to be normalised here. THe fold limb rotation angle is calculated by finding
# the angle between the folded foliation and the axial foliation in the fold frame. The calculated
# axis can be used by passing an N,3 array of the fold axis for every N locations. Or if no
# axis is specified the local intersection lineation is used.
axis = fold.get_fold_axis_orientation(xyz)
flr = f1_frame.calculate_fold_limb_rotation(np.hstack([xyz,s0g]),axis=axis)
limb_svariogram = SVariogram(s1,flr)
guess = limb_svariogram.find_wavelengths()
guess[0] = 5000.
guess/=2.
flr_tan = np.tan(np.deg2rad(flr))
rbf_fold_limb = Rbf(s1,np.zeros(s1.shape),np.zeros(s1.shape),flr_tan,
                    function='gaussian',
                    epsilon=guess[0],
                    smooth=0.05)
xi = np.linspace(f1_frame.features[0].min(),f1_frame.features[0].max(),1000)
def fold_limb_rotation(x):
    return np.rad2deg(np.arctan(rbf_fold_limb(x,np.zeros(x.shape),np.zeros(x.shape))))

# The fold rotation angle function is added to the **FoldEvent**
fold.fold_limb_rotation = fold_limb_rotation

# ### Plot rotation angles
rotation_plots = RotationAnglePlotter()
rotation_plots.add_fold_axis_data(far,s1gy)
rotation_plots.add_axis_svariogram(axis_svariogram)
rotation_plots.add_fold_axis_curve(np.rad2deg(np.arctan(rbf_fold_axis(xi,np.zeros(1000),np.zeros(1000)))), xi)
rotation_plots.add_fold_limb_data(flr,s1)
rotation_plots.add_limb_svariogram(limb_svariogram)
rotation_plots.add_fold_limb_curve(np.rad2deg(np.arctan(rbf_fold_limb(xi,np.zeros(1000),np.zeros(1000)))), xi)
plt.show()
# ### Define interpolation weights
# The DiscreteFoldInterpolator has the ability to weigh the different geometrical properties of the
# fold differently in the least squares system. Much the same as data points and regularisation are
# handled differently for standard interpolation.
# * **fold_orientation** is the weight of the dot product between a
# * **fold_axis** is the weight of the dot product between the fold axis and the interpolated foliation
# * **fold_normalisation**
# * **fold_regularisation** is the fold regularisation term that ensures the fold is smooth along strike
# These weights need to be added to a dictionary and then passed to the feature builder along with the fold
fold_weights = {}
fold_weights['fold_orientation'] = 50.
fold_weights['fold_axis'] = 3.
fold_weights['fold_normalisation'] = 1.
fold_weights['fold_regularisation'] = 10.

# create the geological feature make sure that the constant gradient is 0
folded_stratigraphy = stratigraphy_builder.build(solver=solver,
                                                 fold_weights=fold_weights,
                                                 fold=fold,
                                                 cgw=0)

# ### Visualising Results
# We can visualise the surfaces for isovalues of the fold frame or interpolated surface.
# Normal vectors to the scalar fields can be visualised for evaluation points (if too many
# points are chosen it can be difficult to see the locations. We suggest slicing numpy
# arrays to create regularly sampled locations.

viewer = LavaVuModelViewer(background="white")
viewer.add_isosurface(f1_frame.features[0], colour='green')
viewer.add_isosurface(f1_frame.features[1], colour='blue')
viewer.add_isosurface(folded_stratigraphy,
                      colour='purple',
                      nslices=10,
                      # paint_with=f1_frame.features[0]
                      )
locations = mesh.barycentre[::20,:]
viewer.add_vector_field(f1_frame.features[2], locations=locations, colour='red')
viewer.add_vector_field(f1_frame.features[1], locations=locations, colour='green')
viewer.add_vector_field(f1_frame.features[0], locations=locations, colour='blue')
viewer.interactive()

