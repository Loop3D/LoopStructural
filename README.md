# LoopStructural
Loop3D Geological Forward Modeling Engine.

* Python/cython implementation of a piecewise linear interpolation (Discrete Smooth Interpolator in Gocad) 
* Fold interpolation using constraints outlined in Laurent 2016 with fold geostatistical tools shown in Grose et al., 2017
* Fault interpolation 

If you want to use LoopStructural the easiest way to get started is to use a docker container and a jupyter notebook environment.  

## Using docker
Follow the installation instructions for docker [here](https://docs.docker.com/install/).

Using a github client (e.g. for windows Github Desktop) Clone this repository to your local drive and change directory to the location where you cloned it using `cd THE LOCATION YOU CLONED LOOPSTRUCTURAL`
The docker container can be built by running the following command `docker build -t=loop .`.
LoopStructural can be used by running  `docker run  -i -t -p 8888:8888 loop ` This will start a jupyter notebook server running on localhost:8888 without password or certificate required. Be aware any changes made to the notebooks within the notebooks directory will **NOT** be saved on your local versions.

If you want to use your own data with the docker container you will need to link your local directory (this can be anywhere) with the docker container. To do this add `-v LOCALDIRPATH:/home/joyvan/shared_volume` to the docker command so it becomes `docker run  -i -t -p 8888:8888 -v LOCALDIRPATH:/home/joyvan/shared_volume`. **LOCALDIRPATH** is the full path to the directory you want to share.

## Installing LoopStructural locally
Follow the installation instructions [here](https://github.com/Loop3D/LoopStructural/blob/master/docs/source/installation.rst) 

## Problems
Any bugs/feature requests/comments send to lachlan.grose@monash.edu
