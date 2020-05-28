"""
Probabilistic folds
===================
Imports
-------

"""

#import the Forward Modelling Engine modules - LoopStructural
from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_noddy_single_fold
from LoopStructural.visualisation.model_visualisation import LavaVuModelViewer
from LoopStructural.utils.helper import strike_dip_vector, plunge_and_plunge_dir_to_vector
# from LoopStructural.visualisation.rotation_angle_plotter import RotationAnglePlotter
# import other libraries
import pandas as pd
import numpy as np
# from scipy.interpolate import Rbf
# import matplotlib.pyplot as plt
import logging
logging.getLogger().setLevel(logging.DEBUG)
# %load_ext snakeviz

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

data = data.sort_values('random') # sort the data by a random int then we can select N random points 
data.head()


######################################################################
# The data can be visualised using the lavavu 3d viewer - by first
# converting from strike and dip to normal vectors. Note that there are a
# lot of data points to display as the model volume was regularly sampled
# on a grid.
# 

npoints = 20


######################################################################
# Build model using maximum likelihood
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# val = /
mdata = pd.concat([data[:npoints],data[data['type']=='s1']])
model = GeologicalModel(boundary_points[0,:],boundary_points[1,:])
model.set_model_data(mdata)
fold_frame = model.create_and_add_fold_frame('s1',
                                             interpolatortype='PLI',
                                             nelements=10000,buffer=0.5,
                                             solver='pyamg',
                                            damp=True
                                            )
stratigraphy = model.create_and_add_folded_foliation('s0',
                                               fold_frame['feature'],
                                                nelements=10000,
                                               fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
#                                                    limb_wl=1
                                                     buffer=0.5
                                                    )
viewer = LavaVuModelViewer(model,background="white")
# viewer.add_scalar_field(model.bounding_box,(38,55,30),
#                       'box',
#                      paint_with=stratigraphy,
#                      cmap='prism')
viewer.add_isosurface(fold_frame['feature'][0],
                      colour='blue',
#                       isovalue=0.4,
                      alpha=0.5)
viewer.add_data(fold_frame['feature'][0])
viewer.add_isosurface(fold_frame['feature'][1],colour='green',alpha=0.5,isovalue=0)
# viewer.add_vector_field(fold_frame['feature'][0],locations=fold_frame['feature'][0].get_interpolator().support.barycentre())
viewer.add_data(fold_frame['feature'][1])

# viewer.add_data(stratigraphy['feature'])
viewer.add_isosurface(stratigraphy['feature'])
viewer.interactive()
logging.getLogger().critical("aaaa")
# plt.plot(stratigraphy['foliation'],stratigraphy['limb_rotation'],'bo')
# x = np.linspace(fold_frame['feature'][0].min(),fold_frame['feature'][1].max(),100)
# plt.plot(x,stratigraphy['fold'].fold_limb_rotation(x),'r--')


######################################################################
# Useful pdfs
# ~~~~~~~~~~~
# 

def normal(value,mu, sigma):
    prob = -np.log(1.0/(np.sqrt(2*np.pi)*sigma))-0.5*(value-mu)**2/sigma**2
#     print(prob)
    return prob#-np.log(1.0/(np.sqrt(2*np.pi)*sigma))-0.5*(value-mu)**2/sigma**2

def jefferyprior(sigma):
    if sigma <= 0:
        return -np.inf
    return np.log(sigma)


######################################################################
# Define the fixed parts of the model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

from LoopStructural.modelling.fold import fourier_series
model = GeologicalModel(boundary_points[0,:],boundary_points[1,:])
model.set_model_data(mdata)
fold_frame = model.create_and_add_fold_frame('s1',
                                             interpolatortype='PLI',
                                             nelements=10000,buffer=0.5,
                                             solver='lu',
                                            damp=True
                                            )


######################################################################
# Define log_likelihood
# ~~~~~~~~~~~~~~~~~~~~~
# 
# We can move the expensive parts of the computation out of the likelihood
# function to reduce computational time
# 

wl = 10
c0 = 0
c1 = 0
c2 = 0 
fold_limb_function = lambda x: np.rad2deg(
            np.arctan(
                fourier_series(x, c0, c1, c2, wl)))

strati = model.create_and_add_folded_foliation('s0',
                                           fold_frame['feature'],
                                            nelements=10000,
                                           fold_axis=[-6.51626577e-06, -5.00013645e-01, -8.66017526e-01],
#                                                    limb_wl=1
                                                 buffer=0.5,
                                               limb_function=fold_limb_function,
                                               solver = 'fake'
                                        )
def log_likelihood(theta):
    # unpack parameters
    wl, c0, c1, c2, sigma = theta
    fold_limb_function = lambda x: np.rad2deg(
            np.arctan(
                fourier_series(x, c0, c1, c2, wl)))
    strati['fold'].fold_limb_rotation.set_function(fold_limb_function)
#     sigma = .1

    misfit = strati['fold'].fold_limb_rotation.calculate_misfit()
    
    log_like = -0.5 * np.sum(np.log(2 * np.pi * sigma ** 2) + (0 - misfit) ** 2 / sigma ** 2)

    if ~np.isfinite(log_like):
        return -np.inf
    
    return log_like


######################################################################
# Define prior
# ~~~~~~~~~~~~
# 
# Assign each of the model parameters to a PDF using the helper functions
# before.
# 

def log_prior(theta):
    wl, c0, c1, c2, sigma = theta
    lp = jefferyprior(sigma)
    lp+= normal(c0,0,2)
    lp+= normal(c1,0,2)
    lp+= normal(c2,0,2)
    lp+= normal(wl,10,3)
    return lp

def log_probability(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(theta)
    return lp + log_likelihood(theta)

from scipy.optimize import minimize
nll = lambda *args : -log_likelihood(*args)
initial = [10,1,1,1,0.0001]
soln = minimize(nll, initial)
print(soln)

x = np.linspace(-10,10,100)
plt.plot(stratigraphy['fold'].fold_limb_rotation.fold_frame_coordinate,stratigraphy['fold'].fold_limb_rotation.rotation_angle,'bo')
wl, c0, c1, c2, sig = soln['x']
plt.plot(x,np.rad2deg(
            np.arctan(
                fourier_series(x, c0, c1, c2, wl))),alpha=0.3)

import emcee
# ndim = 5
start = soln.x#stratigraphy['fold'].fold_limb_rotation.fitted_params#np.zeros(ndim)#np.array([600,0,0,0])
pos = soln.x + 1e-3 * np.random.randn(10, ndim)
nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability)
sampler.run_mcmc(pos, 200, progress=True)

flat_samples = sampler.get_chain(discard=100, flat=True)#,thin=10)

plt.hist(flat_samples[:,0])

import corner
flat_samples.shape
fig = corner.corner(
    flat_samples
);

import matplotlib.pyplot as plt

x = np.linspace(-10,10,100)
plt.plot(stratigraphy['fold'].fold_limb_rotation.fold_frame_coordinate,stratigraphy['fold'].fold_limb_rotation.rotation_angle,'bo')
for i in range(len(flat_samples)):
    wl, c0, c1, c2, sig = flat_samples[i,:]
    plt.plot(x,np.rad2deg(
                np.arctan(
                    fourier_series(x, c0, c1, c2, wl))),alpha=0.1,color='black')
    

plt.plot(stratigraphy['fold'].fold_limb_rotation.fold_frame_coordinate,stratigraphy['fold'].fold_limb_rotation.rotation_angle,'bo')

