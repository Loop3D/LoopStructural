# FME
Loop3D Geological Forward Modeling Engine.

Currently contains a discrete tetrahedral mesh interpolation algorithm based on gocads DSI. 

* Fold interpolation using constraints outlined in Laurent 2016
* Fault interpolation 

**Getting started with FME**
The easiest way to get started with **FME** is to clone this repository. If you do not have docker installed, find the appropriate version of docker for your system and install.

FME can then be run by first calling the following command `docker build -t=fme .` which builds the docker image for FME.
FME can be used by running the docker container  `docker run  -i -t -p 8888:8888 -v $WORKDIR:/notebooks fme`


