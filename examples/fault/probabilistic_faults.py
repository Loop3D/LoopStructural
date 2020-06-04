"""
Probabilistic faults
====================
Imports
~~~~~~~

"""

from LoopStructural import GeologicalModel
from LoopStructural.visualisation import LavaVuModelViewer
from LoopStructural.datasets import load_intrusion
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
# %matplotlib inline
# %load_ext snakeviz

data, bb = load_intrusion()

fault_data = data[data['type']=='fault']

fault_data

model = GeologicalModel(bb[0,:],bb[1,:])
model.set_model_data(fault_data)
fault = model.create_and_add_fault('fault',
                                   -600,
                                   nelements=2000,
                                   steps=4,
                                   interpolatortype='PLI',
                                  buffer=0.3
                                  )


bedding_val = np.random.random((40,4))
bedding_val[:,0]*=bb[1,0]
bedding_val[:,1]*=bb[1,1]
bedding_val[:,2]=-600
bedding_val[:,3]=0
bedding_val = np.vstack([bedding_val,bedding_val])
bedding_val[40:,2]-=-500
bedding_val[40:,3]= -1
# print(bedding_val)
# print(fault['feature'].evaluate(model.scale(bedding_val)))
bedding_val[:,:3] = model.rescale(fault['feature'].apply_to_points(model.scale(bedding_val[:,:3])))

# print(bedding_val)

new_data = pd.DataFrame(bedding_val,columns=['X','Y','Z','val'])
new_data['type'] = 'strati'
# new_data['val'] = 0

# normal_vec = pd.DataFrame([[9000,10,10,0,0,1]],columns=['X','Y','Z','nx','ny','nz'])
# normal_vec['type'] = 'strati'

data = pd.concat([fault_data,new_data],sort=False)
data

bb

# %%snakeviz
# for i in range(10):
model = GeologicalModel(bb[0,:],bb[1,:])
model.set_model_data(data)
fault = model.create_and_add_fault('fault',
                                   600,
                                   nelements=2000,
                                   steps=4,
                                   interpolatortype='PLI',
                                  buffer=0.3,
                                   solver='pyamg'
                                  )
strati = model.create_and_add_foliation('strati',
                                        nelements=10000,
                                        interpolatortype='PLI',
                                        cgw=0.1,
                                        solver='pyamg',
                                        buffer=0.5
                                       )

viewer = LavaVuModelViewer(model)
viewer.add_isosurface(strati['feature'],isovalue=0)
viewer.add_isosurface(fault['feature'],isovalue=0)
# viewer.add_scalar_field(fault['feature'].displacementfeature)
viewer.add_vector_field(fault['feature'][1],locations=model.regular_grid()[::100,:])
viewer.add_data(strati['feature'])
# viewer.add_points(strati['feature'].get_interpolator()[:,:3],pts)
viewer.interactive()

viewer = LavaVuModelViewer(model)
viewer.add_isosurface(strati['feature'],isovalue=0)
viewer.add_isosurface(fault['feature'],isovalue=0)
# viewer.add_scalar_field(fault['feature'].displacementfeature)
viewer.add_vector_field(fault['feature'][1],locations=model.regular_grid()[::100,:])
viewer.add_data(strati['feature'])
viewer.add_points(strati['feature'].get_interpolator().get_value_constraints()[:,:3],name='pts',pointsize=10)
viewer.interactive()

model = GeologicalModel(bb[0,:],bb[1,:],reuse_supports=True)
model.set_model_data(data)
displacement = -700
viewer = LavaVuModelViewer(model)
viewer.rotate([0.019632680341601372, 88.20027923583984, -89.94925689697266])#lv.rotatex(90)

fault = model.create_and_add_fault('fault',
                                   displacement,
                                   nelements=2000,
                                   steps=4,
                                   interpolatortype='PLI',
                                  buffer=0.3,
                                   solver='lu'
                                  )
dist = []
for displacement in [-1000,-500,0,500,1000]:
#     print("displacement: {}".format(displacement))
    fault['feature'].set_displacement(displacement)
#     print('one')
    strati = model.create_and_add_foliation('strati',
                                            nelements=10000,
                                            interpolatortype='PLI',
                                            cgw=0.1,
                                            solver='lu',
                                            buffer=0.5
                                           )
    
#     print('two')
#     plt.figure()
#     print(np.sum(np.abs(strati['feature'].evaluate_value_misfit())))

# #     strati['feature'].faults_enabled = False#toggle_faults()
#     viewer.add_isosurface(strati['feature'],slices=[-1,0])
#     viewer.add_isosurface(fault['feature'][0],isovalue=0)
#     viewer.add_data(strati['feature'])
    viewer.add_points(strati['feature'].get_interpolator().get_value_constraints()[:,:3],name='points_{}'.format(displacement)
                      ,pointsize=3)
#     #viewer.lv.rotatey(90)
    viewer.display()

# viewer = LavaVuModelViewer(model)
# viewer.add_isosurface(strati['feature'],isovalue=0)
# viewer.add_isosurface(fault['feature'][0],isovalue=0)
# viewer.add_data(strati['feature'])
# viewer.rotate([0.019632680341601372, 88.20027923583984, -89.94925689697266])#lv.rotatex(90)
# #viewer.lv.rotatey(90)
# viewer.display()

plt.plot([-1000,-500,0,500,1000],dist,'bo')

# viewer.lv['xyzrotate']

def planeFit(points):
    """
    p, n = planeFit(points)

    Given an array, points, of shape (d,...)
    representing points in d-dimensional space,
    fit an d-dimensional plane to the points.
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.
    """
    import numpy as np
    from numpy.linalg import svd
    #points = points.T
    #print('p',points.shape)
#     points = np.reshape(points, (np.shape(points)[0], -1)) # Collapse trialing dimensions
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1], points.shape[0])
    ctr = points.mean(axis=1)
    x = points - ctr[:,np.newaxis]
    M = np.dot(x, x.T) # Could also use np.cov(x) here.
    U,S,V = svd(M)
    normal = V[-1]
    d = -np.sum(normal*ctr)
    return np.hstack([normal,[d]])

def planeDistance(points):
    params = planeFit(points)
    a, b, c, d = params
    x, y, z = points
    length = np.sqrt(a**2 + b**2 + c**2)
    return (np.abs(a * x + b * y + c * z + d) / length).mean()



import emcee

def log_prior(theta):
    displacement, sigma, mu = theta
#     mu = 0
    if sigma <= 0:
        return -np.inf
    if mu <= 0:
        return -np.inf
    prob = -np.log(1.0/(np.sqrt(2*np.pi)*sigma))-0.5*(displacement-mu)**2/sigma**2 - np.log(sigma) - np.log(mu)
#     print()
    return -np.log(1.0/(np.sqrt(2*np.pi)*sigma))-0.5*(displacement-mu)**2/sigma**2 - np.log(sigma) - np.log(mu)

import dill as pickle

model = GeologicalModel(bb[0,:],bb[1,:],reuse_supports=True)
model.set_model_data(data)
fault = model.create_and_add_fault('fault',
                                   0,
                                   nelements=2000,
                                   steps=1,
                                   interpolatortype='PLI',
                                  buffer=0.8,
                                   solver='pyamg'
                                  )
def log_likelihood(theta):
    displacement, sigma, mu = theta
#     print("displacement: {}".format(displacement))
    fault['feature'].set_displacement(displacement)

    #strati['feature'].get_interpolator().data_added = False
    strati = model.create_and_add_foliation('strati',
                                            nelements=2000,
                                            interpolatortype='PLI',
                                            cgw=0.1,
                                            solver='fake',
                                            buffer=1
                                           )
    strati['feature'].builder.add_data_to_interpolator()
    points = strati['feature'].get_interpolator().get_value_constraints()[:,:4]
    unique_values = np.unique(points[:,3])
    distance = np.zeros_like(unique_values).astype(float)
    for i, u in enumerate(unique_values):
        distance[i] = planeDistance(points[points[:,3] == u,:3].T)
    
        
#     plt.hist(strati['feature'].evaluate_value_misfit())
    n = len(distance)#strati['feature'].interpolator.get_value_constraints()[:,:3].shape[0]
    print(distance)
    log_like = -0.5 * np.sum(np.log(2 * np.pi * sigma ** 2) + (0 - distance) ** 2 / sigma ** 2)
    #data_added = False

#     sigma2 = 3
#     log_like = -(n/2)*np.log(2*np.pi) - (n/2)*np.log(sigma2)
#     log_like-= (1/(2*sigma2))*np.sum(np.abs(strati['feature'].evaluate_value_misfit())**2)
    
#     sigma2 = strati['feature'].evaluate_value(strati['feature'].interpolator.get_value_constraints()[:,:3]) ** 2 
#     log_like = -0.5 * np.sum((strati['feature'].evaluate_value_misfit()) ** 2 / sigma2 + np.log(sigma2))
#     print("log likelihood {}".format(log_like))
#     print("missfit {}".format(np.sum(strati['feature'].evaluate_value_misfit())))
    if ~np.isfinite(log_like):
        return -np.inf
#     pickle.dump(model,open("models/model2_sigma_{}_mu_{}_displacement_{}.pkl".format(sigma,mu,displacement),"wb"))
    return log_like                         


for d in [-600,0,600]:
    log_likelihood((d,1,1))

def log_probability(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(theta)
    return lp + log_likelihood(theta)

import emcee
start = np.array([600,0,0])
pos = start + np.array([1e3,1e2,1e2]) * np.random.randn(50, 3)
pos[:,1:] = np.abs(pos[:,1:])
nwalkers, ndim = pos.shape

print(pos)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability)
sampler.run_mcmc(pos, 100, progress=True,tune=True)

flat_samples = sampler.get_chain(discard=20, flat=True)

plt.hist(flat_samples[:,0])

import corner
flat_samples.shape
fig = corner.corner(
    flat_samples
);

plt.plot(sampler.get_chain(flat=True)[:,0])

chain = sampler.get_chain()

chain.shape

for i in range(20):
    plt.plot(chain[:,i,0])

for i in range(20):
    plt.plot(chain[:,i,1])

