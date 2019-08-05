import geopandas
import numpy as np
import matplotlib.pyplot as plt
from FME.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from FME.interpolators.discrete_fold_interpolator import DiscreteFoldInterpolator as DFI
from FME.modelling.features.geological_feature import GeologicalFeatureInterpolator
from FME.supports.tet_mesh import TetMesh
from FME.visualisation.model_visualisation import LavaVuModelViewer
from FME.modelling.structural_frame import StructuralFrameBuilder
from FME.modelling.fold.foldframe import FoldFrame
from FME.modelling.fold.fold import FoldEvent
from FME.utils.helper import strike_dip_vector
points = geopandas.read_file('data.gpkg',layer='points')
orientations = geopandas.read_file('data.gpkg',layer='orientations')
model_area = geopandas.read_file('data.gpkg',layer='bounding_box')

geom = model_area['geometry']#.shapes()
coords = np.array(geom[0].exterior.coords)#[0]
minz = -(np.max(coords[:,0])-np.min(coords[:,0]))/2.

boundary_points = np.zeros((2,3))
boundary_points[0,0] = np.min(coords[:,0])-10
boundary_points[0,1] = np.min(coords[:,1])
boundary_points[0,2] = minz
boundary_points[1,0] = np.max(coords[:,0])
boundary_points[1,1] = np.max(coords[:,1])
boundary_points[1,2] = -minz*0.1
mesh = TetMesh()
mesh.setup_mesh(boundary_points,n_tetra=50000,)

fold_frame_interpolator = PLI(mesh)
fold_frame_builder = StructuralFrameBuilder(
    interpolator=fold_frame_interpolator,
    mesh=mesh,
    name='F1_fold_frame')

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
ogw = 3000
ogw /= mesh.n_elements
cgw = 6000
# cgw = cgw / mesh.n_elements
solver='cgp'

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

##processing geopandas data so it works
nans = np.zeros(len(orientations))
nans[:] = np.nan
orientations['x'] = nans
orientations['y'] = nans
orientations['z'] = nans
for i, r in orientations.iterrows():
    orientations.at[i, 'x'] = r['geometry'].xy[0][0]
    orientations.at[i, 'y'] = r['geometry'].xy[1][0]
    orientations.at[i, 'z'] = 0.0
s0 = orientations[orientations['type'] == 's0']
## end geopandas
## get the xyz coordinates of the data and find the fold frame orientation
xyz = s0.loc[:, ['x', 'y', 'z']].as_matrix()
#convert from strike and dip to
s0g = strike_dip_vector(s0.loc[:, 'strike'], s0.loc[:, 'dip'])
s0g /= np.linalg.norm(s0g,axis=1)[:,None]
l1 = f1_frame.calculate_intersection_lineation(np.hstack([xyz,s0g]))
far = f1_frame.calculate_fold_axis_rotation(np.hstack([xyz,l1]))
s1 = f1_frame.features[0].evaluate_value(xyz)
s1gy = f1_frame.features[1].evaluate_value(xyz)

##quick figure
# fig, ax = plt.subplots(1,2,figsize=(15,5))
# ax[0].plot(s1gy,far,'bo')
# ax[0].set_ylim(-90,90)
# plt.show()
# plt.savefig("test.png")
#ax[0].ylim(-90,90)
# svario = s_variogram(s1gy,far)
# svario.setup()
guess = [20000]
# guess = svario.find_wavelengths()
# ax[1].plot(svario.h,svario.var,'bo')
# #ax[1].axvline(guess[0])
from scipy.interpolate import Rbf
far_tan = np.tan(np.deg2rad(far))
rbf_fold_axis = Rbf(s1gy,np.zeros(s1gy.shape),np.zeros(s1gy.shape),far_tan,function='gaussian',epsilon=guess[0],smooth=.1)
xi = np.linspace(f1_frame.features[1].min(),
                 f1_frame.features[1].max(),1000)
plt.plot(xi,np.rad2deg(np.arctan(rbf_fold_axis(xi,np.zeros(1000),np.zeros(1000)))))
plt.plot(s1gy,far,'bo')
plt.ylim(-90,90)
def fold_axis_rotation(x):
    return np.rad2deg(np.arctan(rbf_fold_axis(x,np.zeros(x.shape),np.zeros(x.shape))))
# plt.show()


fold = FoldEvent(f1_frame,fold_axis_rotation,None)
axis = fold.get_fold_axis_orientation(xyz)
axis/=np.linalg.norm(axis,axis=1)[:,None]
flr = f1_frame.calculate_fold_limb_rotation(np.hstack([xyz,s0g]),axis=axis)
##quick figure
fig, ax = plt.subplots(2,2,figsize=(15,10))
ax[0][0].plot(s1,flr,'bo')
ax[0][0].set_ylim(-90,90)
# plt.show()
#ax[0].ylim(-90,90)
# svario = s_variogram(s1,flr)
# svario.setup()
# guess = svario.find_wavelengths()
# ax[1].plot(svario.h,svario.var,'bo')
#ax[1].axvline(guess[1])
guess = np.array(guess)
guess[0] = 2500

flr_tan = np.tan(np.deg2rad(flr))
rbf_fold_limb = Rbf(s1,np.zeros(s1.shape),np.zeros(s1.shape),flr_tan,function='gaussian',epsilon=guess[0],smooth=.05)
xi = np.linspace(f1_frame.features[0].min(),f1_frame.features[0].max(),1000)
plt.plot(xi,np.rad2deg(np.arctan(rbf_fold_limb(xi,np.zeros(1000),np.zeros(1000)))))
plt.plot(s1,flr,'bo')
plt.ylim(-90,90)
def fold_limb_rotation(x):
    return np.rad2deg(np.arctan(rbf_fold_limb(x,np.zeros(x.shape),np.zeros(x.shape))))
# plt.show()
fold.fold_limb_rotation = fold_limb_rotation
fold_interpolator = DFI(mesh,fold)
stratigraphy_builder = GeologicalFeatureInterpolator(fold_interpolator, name="folded_stratigraphy")
# for i, r in orientations.iterrows():
#     if r['type'] == 's0':
#         xy = r['geometry'].xy
#         z = 0
#         if 'z' in r:
#             z = r['z']
#         stratigraphy_builder.add_strike_and_dip([xy[0][0],xy[1][0],z],r['strike'],r['dip'])
# for i, r in points.iterrows():
#     if r['type'] == 's0':
#         xy = r['geometry'].xy
#         z = 0
#         if 'z' in r:
#             z = r['z']
#         stratigraphy_builder.add_point([xy[0][0],xy[1][0],z],r['value'])
stratigraphy_builder.add_point(mesh.pca.inverse_transform([0,0,0]),0.)
fold_weights = {}
fold_weights['fold_orientation'] = 50.
fold_weights['fold_axis'] = 3.
fold_weights['fold_normalisation'] = 1.
fold_weights['fold_regularisation'] = 10.
folded_stratigraphy = stratigraphy_builder.build(solver=solver,fold_weights=fold_weights,fold=fold)

viewer = LavaVuModelViewer(background="white")
# viewer.plot_isosurface(f1_frame.features[0],  colour='green')
# viewer.plot_isosurface(f1_frame.features[1],  colour='blue')
viewer.plot_isosurface(folded_stratigraphy, colour='purple')
locations = mesh.barycentre[::20,:]
# viewer.plot_vector_field(f1_frame.features[2], locations=locations, colour='red')

viewer.lv.interactive()