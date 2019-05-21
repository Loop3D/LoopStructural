# FME
Loop3D Geological Forward Modeling Engine.

* Python/cython implementation of DSI 
* Fold interpolation using constraints outlined in Laurent 2016
* Fault interpolation 

**Getting started with FME**
The easiest way to get started with **FME** is to clone this repository and build the docker image using the Dockerfile. If you do not have docker installed, find the appropriate version of docker for your system and install it.

FME can then be used by first building the docker image `docker build -t=fme .`.
FME can be used by run  `docker run  -i -t -p 8888:8888 -v $WORKDIR:/notebooks fme %DATADIR:/data` where `$WORKDIR` is the path to your jupyter notebooks and `$DATA` is the path to your data folder. 




