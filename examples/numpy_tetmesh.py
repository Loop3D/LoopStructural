from LoopStructural.interpolators.structured_tetra import TetMesh
from LoopStructural.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator as PLI
# from LoopStructural.supports.tet_mesh import TetMesh
from LoopStructural.core.features.geological_feature_builder import GeologicalFeatureInterpolator

import numpy as np
import pandas as pd
import glob

boundary_points = np.zeros((2,3))
boundary_points[0,0] = 548800
boundary_points[0,1] = 7816600
boundary_points[0,2] = -11010
boundary_points[1,0] = 552500
boundary_points[1,1] = 7822000
boundary_points[1,2] = -8400
steps = np.ones(3).astype(int)

steps[:] = 10

step_vector = boundary_points[1,:] - boundary_points[0,:]
# build the mesh
mesh = TetMesh(boundary_points[0,:],steps,step_vector)
# mesh.setup_mesh(boundary_points, n_tetra=50000,)

# link mesh to the interpolator
interpolator = PLI(mesh)

dips = pd.read_csv('../notebooks/data/Dips.csv',delimiter=';')


# import all of the csv into the same dataframe use glob to find all files matching pattern
dfs = []
for f in glob.glob('../notebooks/data/*Points.csv'):
    dfs.append(pd.read_csv(f,delimiter=';'))
points = pd.concat(dfs,axis=0,ignore_index=True)

dfs = []
for f in glob.glob('../notebooks/data/*Section.csv'):
    dfs.append(pd.read_csv(f,delimiter=';'))
sections = pd.concat(dfs,axis=0,ignore_index=True)
stratigraphy_builder = GeologicalFeatureInterpolator(
    interpolator=interpolator,
    name='stratigraphy')
solver = 'cg'
for i, r in points.iterrows():
    stratigraphy_builder.add_point([r['X'], r['Y'], r['Z']],
                                   r['Strati'])  # xy[0][0],xy[1][0],z],r['value'],itype=r['itype'])
for i, r in sections.iterrows():
    stratigraphy_builder.add_point([r['X'], r['Y'], r['Z']],
                                   r['Strati'])  # xy[0][0],xy[1][0],z],r['value'],itype=r['itype'])
for i, r in dips.iterrows():
    stratigraphy_builder.add_planar_constraint([r['X'], r['Y'], r['Z']], [r['OrientX'], r['OrientY'], r['OrientZ']])

stratigraphy = stratigraphy_builder.build(solver='cg',cpw=1,cgw=0.1)