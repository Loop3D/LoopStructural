import geopandas
import numpy as np

from FME.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
from FME.supports.tet_mesh import TetMesh
from FME.modelling.geological_feature import GeologicalFeature, FaultedGeologicalFeature, GeologicalFeatureBuilder
from FME.visualisation.model_visualisation import LavaVuModelViewer
from FME.modelling.structural_frame import StructuralFrameBuilder, StructuralFrame

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
mesh.setup_mesh(boundary_points,n_tetra=20000,)

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
    solver=solver,
    gxxgy=2 * ogw,
    gxxgz=2 * ogw,
    gyxgz=ogw,
    gxcg=cgw,
    gycg=cgw,
    gzcg=cgw,
    shape='rectangular'
                         )

viewer = LavaVuModelViewer(background="white")
viewer.plot_isosurface(f1_frame.features[0],  colour='green')
viewer.plot_isosurface(f1_frame.features[1],  colour='blue')
locations = mesh.barycentre[::20,:]
viewer.plot_vector_field(f1_frame.features[2], locations=locations, colour='red')

viewer.lv.interactive()