# FME
Loop3D Geological Forward Modeling Engine.

* Python/cython implementation of DSI 
* Fold interpolation using constraints outlined in Laurent 2016
* Fault interpolation 

If you want to use FME the best way to use it is within a jupyter notebook environment. 
## Using docker
FME can then be used by first building the docker image `docker build -t=fme .`.
FME can be used by run  `docker run  -i -t -p 8888:8888 -v $WORKDIR:/notebooks fme ` where **`$WORKDIR`** is the path to your jupyter notebooks. You must change this to the absolute path! This will start a notebook server running on localhost:8888 without  password or certificate. Be aware any changes made to the notebooks within the notebooks directory will be saved on your local versions.

## On Linux
If you are running Linux you should be able to run python3 setup.py install build_ext --inplace - you may need to install some additional dependencies for meshypy and lavavu. It is recommended that you use a python environment manager. 

## On windows
Use the linux subsystem on windows. cd to the working directory and run python3 setup.py install build_ext --inplace. In order to use lavavu you will have to use TINI. 

## Problems
Any bugs/feature requests/comments send to lachlan.grose@monash.edu
