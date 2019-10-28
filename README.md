# LoopStructural
Loop3D Geological Forward Modeling Engine.

* Python/cython implementation of a piecewise linear interpolation (Discrete Smooth Interpolator in Gocad) 
* Fold interpolation using constraints outlined in Laurent 2016 with fold geostatistical tools shown in Grose et al., 2017
* Fault interpolation 

If you want to use LoopStructural the easiest way to get started is to use a docker container and a jupyter notebook environment.  

## Using docker
Follow the installation instructions for docker [here](https://docs.docker.com/install/).

Clone this repository to your local drive and change directory to the LoopStructural folder `cd ${LOOPSTRUCTURALDIRECTORY}`
The docker container can be built by running the following command `docker build -t=loop .`.
LoopStructural can be used by running  `docker run  -i -t -p 8888:8888 -v $WORKDIR:/notebooks LoopStructural ` where **`$WORKDIR`** is the path to your jupyter notebooks. You must change this to the absolute path! This will start a notebook server running on localhost:8888 without  password or certificate. Be aware any changes made to the notebooks within the notebooks directory will be saved on your local versions.

## Installing LoopStructural locally
Follow the installation instructions [here](https://github.com/Loop3D/LoopStructural/blob/master/docs/source/installation.rst) 

## Problems
Any bugs/feature requests/comments send to lachlan.grose@monash.edu
