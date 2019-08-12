# FME
Loop3D Geological Forward Modeling Engine.

* Python/cython implementation of a piecewise linear interpolation (Discrete Smooth Interpolator in Gocad) 
* Fold interpolation using constraints outlined in Laurent 2016 with fold geostatistical tools shown in Grose et al., 2017
* Fault interpolation 

If you want to use FME the best way to use it is within a jupyter notebook environment. 

## Using docker
FME can then be used by first building the docker image `docker build -t=fme .`.
FME can be used by run  `docker run  -i -t -p 8888:8888 -v $WORKDIR:/notebooks fme ` where **`$WORKDIR`** is the path to your jupyter notebooks. You must change this to the absolute path! This will start a notebook server running on localhost:8888 without  password or certificate. Be aware any changes made to the notebooks within the notebooks directory will be saved on your local versions.

## On Linux
FME requires the following libraries to be install on linux. 

Required dependencies
* pybind11-dev
* mesa-common-dev
* mesa-utils
* libgl1mesa-dev
`sudo apt-get install pybind11-dev mesa-common-dev mesa-utils libl1mesa-dev`

To build FME first install the required modules for python `pip3 install -r requirements.txt` and then install `python3 setup.py install build_ext --inplace`

It is recommended that you use a python environment manager. 

## On windows
You can install FME on windows natively if you have a working C compiler. Otherwise the easiest way is to use the linux subsystem for linux or docker. 
In order to use lavavu you will have to use some method of fowarding X to the windows client. 

## Problems
Any bugs/feature requests/comments send to lachlan.grose@monash.edu
