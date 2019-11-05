import pandas as pd
import numpy as np
## one fold data
boundary_points = np.zeros((2,3))
boundary_points[0,0] = 0
boundary_points[0,1] = 0
boundary_points[0,2] = 5000
boundary_points[1,0] = 10000
boundary_points[1,1] = 7000
boundary_points[1,2] = 10000

data = pd.read_csv('../notebooks/data/onefold.csv',delimiter='\t')
data = data[data['event'] == 0]
data['type'] = 's0'

data = data[data['x'] > 200]
data = data[data['x'] < 9800]
data = data[data['y'] > 200]
data = data[data['y'] < 6800]
data = data[data['z'] > 5200]
data = data[data['z'] < 9800]
data = data.rename(columns={"x":"X","y":"Y","z":"Z"})
fold_frame_data = pd.DataFrame(np.array([[1000,1000,5200,0,90,0,'s1'],
                                         [1000,1000,5200,90,90,1,'s1']]),columns=
                               ['X','Y','Z','strike','dip','coord','type'])
data = data.drop(columns=['Unnamed: 6','event'])
data = pd.concat([data,fold_frame_data],sort=False)
data['random'] = np.random.rand(len(data))
data[["X","Y","Z","strike","dip","coord"]] = data[["X","Y","Z","strike","dip","coord"]].apply(pd.to_numeric)
data.to_pickle('../LoopStructural/datasets/data/onefolddata.pkl')
np.savetxt('../LoopStructural/datasets/data/onefoldbb.txt', boundary_points)